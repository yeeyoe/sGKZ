#include "gkz/gkz.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void require_close(double actual, double expected, double tolerance,
                   const std::string& message) {
  if (std::abs(actual - expected) > tolerance) {
    throw std::runtime_error(message + ": actual=" +
                             std::to_string(actual) +
                             ", expected=" + std::to_string(expected));
  }
}

gkz::PointConfiguration square_configuration() {
  return gkz::PointConfiguration::from_points(
      {{0, 0}, {1, 0}, {0, 1}, {1, 1}});
}

using PolygonTriangulation =
    std::vector<std::array<std::size_t, 3>>;

std::vector<PolygonTriangulation> triangulate_interval(std::size_t first,
                                                       std::size_t last) {
  if (last <= first + 1) {
    return {PolygonTriangulation{}};
  }
  std::vector<PolygonTriangulation> result;
  for (std::size_t middle = first + 1; middle < last; ++middle) {
    const auto left = triangulate_interval(first, middle);
    const auto right = triangulate_interval(middle, last);
    for (const auto& left_faces : left) {
      for (const auto& right_faces : right) {
        PolygonTriangulation faces = left_faces;
        faces.insert(faces.end(), right_faces.begin(), right_faces.end());
        faces.push_back({first, middle, last});
        result.push_back(std::move(faces));
      }
    }
  }
  return result;
}

gkz::GkzVector gkz_from_polygon_triangulation(
    const gkz::PointConfiguration& configuration,
    const PolygonTriangulation& faces) {
  gkz::GkzVector vector;
  vector.area_numerators.assign(configuration.size(), 0);
  for (const auto& face : faces) {
    const auto& a = configuration.points()[face[0]];
    const auto& b = configuration.points()[face[1]];
    const auto& c = configuration.points()[face[2]];
    const gkz::WideInt signed_area =
        static_cast<gkz::WideInt>(b.x - a.x) * (c.y - a.y) -
        static_cast<gkz::WideInt>(b.y - a.y) * (c.x - a.x);
    const gkz::WideInt area = signed_area < 0 ? -signed_area : signed_area;
    require(area > 0, "Hexagon triangulation contains a degenerate face.");
    for (const std::size_t vertex : face) {
      vector.area_numerators[vertex] += area;
    }
  }
  const gkz::WideInt denominator = 3 * configuration.twice_area();
  vector.values = Eigen::VectorXd::Zero(
      static_cast<Eigen::Index>(configuration.size()));
  for (std::size_t i = 0; i < configuration.size(); ++i) {
    vector.values[static_cast<Eigen::Index>(i)] =
        static_cast<double>(vector.area_numerators[i]) /
        static_cast<double>(denominator);
  }
  vector.visible_vertices = configuration.size();
  vector.triangles = faces.size();
  vector.triangulation =
      std::make_shared<const PolygonTriangulation>(faces);
  return vector;
}

void test_gkz_duplicate_detection_uses_exact_area_numerators() {
  gkz::GkzVector first;
  first.area_numerators = {100, 200, 300};
  first.values = Eigen::Vector3d(0.1, 0.2, 0.3);

  gkz::GkzVector exact_duplicate = first;
  exact_duplicate.values[0] += 1e-14;
  require(first.has_same_area_numerators(exact_duplicate),
          "Equal area numerators were not recognized as an exact duplicate.");

  gkz::GkzVector numerically_close_but_distinct = first;
  numerically_close_but_distinct.area_numerators[0] += 1;
  numerically_close_but_distinct.values[0] += 1e-14;
  require(!first.has_same_area_numerators(numerically_close_but_distinct),
          "Distinct exact GKZ data were collapsed by floating-point "
          "closeness.");
}

void test_lattice_enumeration() {
  const auto configuration =
      gkz::PointConfiguration::from_lattice_polygon(
          {{0, 0}, {1, 0}, {1, 1}, {0, 1}}, 2);
  require(configuration.size() == 9,
          "The level-2 unit square must have 9 lattice points.");
  require(configuration.twice_area() == 8,
          "The scaled unit square must have doubled area 8.");
  require_close(configuration.centroid_x(), 1.0, 1e-15,
                "Incorrect polygon centroid x");
  require_close(configuration.centroid_y(), 1.0, 1e-15,
                "Incorrect polygon centroid y");
}

void test_square_oracle_against_both_triangulations() {
  const auto configuration = square_configuration();
  const gkz::RegularTriangulationOracle oracle(configuration);

  Eigen::Vector4d heights;
  heights << 0.3, -0.2, 0.1, 0.4;
  const auto answer = oracle.minimize(heights);
  require(answer.triangulation && answer.triangulation->size() == 2,
          "The square oracle must retain its two finite faces.");

  Eigen::Vector4d diagonal_00_11;
  diagonal_00_11 << 1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 3.0;
  Eigen::Vector4d diagonal_10_01;
  diagonal_10_01 << 1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0;
  const double expected =
      std::min(heights.dot(diagonal_00_11),
               heights.dot(diagonal_10_01));
  require_close(heights.dot(answer.values), expected, 1e-14,
                "The regular-triangulation oracle chose the wrong diagonal");

  require_close(answer.values.sum(), 1.0, 1e-14,
                "GKZ total mass is not 1");
  double moment_x = 0.0;
  double moment_y = 0.0;
  for (std::size_t i = 0; i < configuration.size(); ++i) {
    moment_x += answer.values[static_cast<Eigen::Index>(i)] *
                configuration.points()[i].x;
    moment_y += answer.values[static_cast<Eigen::Index>(i)] *
                configuration.points()[i].y;
  }
  require_close(moment_x, 0.5, 1e-14, "Incorrect GKZ x moment");
  require_close(moment_y, 0.5, 1e-14, "Incorrect GKZ y moment");

  const std::vector<CGAL::Gmpq> exact_heights = {
      CGAL::Gmpq(3), CGAL::Gmpq(-2), CGAL::Gmpq(1), CGAL::Gmpq(4)};
  const std::vector<CGAL::Gmpz> integer_heights = {CGAL::Gmpz(3),
                                                    CGAL::Gmpz(-2),
                                                    CGAL::Gmpz(1),
                                                    CGAL::Gmpz(4)};
  const auto exact_answer = oracle.minimize_exact(exact_heights);
  const auto integer_answer = oracle.minimize_exact_integer(integer_heights);
  require(exact_answer.has_same_area_numerators(integer_answer),
          "Integer exact oracle disagrees with rational exact oracle.");
}

void test_ell_pentagon_example() {
  const auto configuration = gkz::PointConfiguration::from_points(
      {{0, 0}, {1, 0}, {2, 0}, {1, 1}, {0, 1}});
  const auto ell = gkz::compute_ell(configuration);
  require(ell.exact_coefficients[0] == CGAL::Gmpq(5, 27),
          "Pentagon ell_A constant coefficient is wrong.");
  require(ell.exact_coefficients[1] == CGAL::Gmpq(0),
          "Pentagon ell_A x coefficient is wrong.");
  require(ell.exact_coefficients[2] == CGAL::Gmpq(1, 27),
          "Pentagon ell_A y coefficient is wrong.");
  const std::vector<CGAL::Gmpq> expected = {
      CGAL::Gmpq(5, 27), CGAL::Gmpq(5, 27), CGAL::Gmpq(5, 27),
      CGAL::Gmpq(2, 9), CGAL::Gmpq(2, 9)};
  for (std::size_t i = 0; i < expected.size(); ++i) {
    require_close(ell.values[static_cast<Eigen::Index>(i)],
                  expected[i].to_double(), 1e-14,
                  "Pentagon ell_A restriction is wrong");
  }
}

void test_ell_polygon_coordinates_are_normalized() {
  const auto polygon_configuration =
      gkz::PointConfiguration::from_lattice_polygon(
          {{0, 0}, {2, 0}, {1, 1}, {0, 1}}, 2);
  const auto raw_configuration = gkz::PointConfiguration::from_points(
      polygon_configuration.points());
  const auto normalized_ell = gkz::compute_ell(polygon_configuration);
  const auto raw_ell = gkz::compute_ell(raw_configuration);
  require(normalized_ell.exact_coefficients[0] ==
              raw_ell.exact_coefficients[0],
          "Normalizing A_k changed the constant coefficient of ell_A.");
  require(normalized_ell.exact_coefficients[1] ==
              2 * raw_ell.exact_coefficients[1] &&
              normalized_ell.exact_coefficients[2] ==
                  2 * raw_ell.exact_coefficients[2],
          "A_k ell_A slopes were not converted to coordinates on P.");
  require(normalized_ell.exact_values == raw_ell.exact_values,
          "Coordinate normalization changed ell_A restricted to A_k.");
}

void test_plot_data_uses_lower_envelope_at_hidden_points() {
  const auto configuration =
      gkz::PointConfiguration::from_lattice_polygon(
          {{0, 0}, {1, 0}, {1, 1}, {0, 1}}, 2);
  gkz::SolverResult result;
  result.sigma = Eigen::VectorXd::Zero(
      static_cast<Eigen::Index>(configuration.size()));
  std::size_t center = configuration.size();
  for (std::size_t i = 0; i < configuration.size(); ++i) {
    if (configuration.points()[i] == gkz::IntPoint{1, 1}) {
      center = i;
      break;
    }
  }
  require(center < configuration.size(), "Level-2 square has no center.");
  result.sigma[static_cast<Eigen::Index>(center)] = 1.0;

  const auto prefix =
      std::filesystem::temp_directory_path() / "gkz_hidden_surface_test";
  gkz::write_plot_data(prefix, configuration, result,
                       gkz::compute_ell(configuration));
  const auto surface_path = prefix.string() + "_surface.csv";
  std::ifstream stream(surface_path);
  require(static_cast<bool>(stream), "Cannot read hidden-point surface CSV.");
  std::string line;
  std::getline(stream, line);
  for (std::size_t i = 0; i <= center; ++i) {
    require(static_cast<bool>(std::getline(stream, line)),
            "Surface CSV ended before the center point.");
  }
  std::replace(line.begin(), line.end(), ',', ' ');
  std::istringstream parser(line);
  double x = 0.0;
  double y = 0.0;
  double sigma = 0.0;
  double sigma_vee = 0.0;
  double psi = 0.0;
  require(static_cast<bool>(parser >> x >> y >> sigma >> sigma_vee >> psi),
          "Cannot parse hidden-point surface row.");
  require_close(x, 0.5, 1e-15, "Hidden point normalized x is wrong");
  require_close(y, 0.5, 1e-15, "Hidden point normalized y is wrong");
  require_close(sigma, 1.0, 1e-15, "Hidden point sigma is wrong");
  require_close(sigma_vee, 0.0, 1e-15,
                "Hidden point lower-envelope value is wrong");
  require_close(psi, -4.0, 1e-15,
                "psi_k did not use sigma_k^vee at a hidden point");

  const auto subdivision_path = prefix.string() + "_subdivision.csv";
  std::ifstream subdivision(subdivision_path);
  require(static_cast<bool>(subdivision),
          "Cannot read hidden-point subdivision CSV.");
  require(static_cast<bool>(std::getline(subdivision, line)),
          "Subdivision CSV has no header.");
  std::size_t subdivision_rows = 0;
  while (std::getline(subdivision, line)) {
    ++subdivision_rows;
  }
  require(subdivision_rows == 4,
          "The coplanar square subdivision should have one quadrilateral cell.");

  std::filesystem::remove(surface_path);
  std::filesystem::remove(prefix.string() + "_triangles.csv");
  std::filesystem::remove(subdivision_path);
  std::filesystem::remove(prefix.string() + "_ell.csv");
}

void test_exact_plot_data_uses_exact_lower_envelope() {
  const auto configuration =
      gkz::PointConfiguration::from_lattice_polygon(
          {{0, 0}, {1, 0}, {1, 1}, {0, 1}}, 2);
  gkz::SolverResult result;
  result.sigma = Eigen::VectorXd::Zero(
      static_cast<Eigen::Index>(configuration.size()));
  result.exact.certified = true;
  result.exact.sigma.assign(configuration.size(), CGAL::Gmpq(1, 3));
  std::size_t center = configuration.size();
  for (std::size_t i = 0; i < configuration.size(); ++i) {
    if (configuration.points()[i] == gkz::IntPoint{1, 1}) {
      center = i;
      break;
    }
  }
  require(center < configuration.size(), "Level-2 square has no center.");
  result.exact.sigma[center] = CGAL::Gmpq(2, 3);

  const auto prefix =
      std::filesystem::temp_directory_path() / "gkz_exact_surface_test";
  gkz::write_plot_data(prefix, configuration, result,
                       gkz::compute_ell(configuration));
  const auto surface_path = prefix.string() + "_surface.csv";
  std::ifstream stream(surface_path);
  require(static_cast<bool>(stream), "Cannot read exact surface CSV.");
  std::string line;
  std::getline(stream, line);
  for (std::size_t i = 0; i <= center; ++i) {
    require(static_cast<bool>(std::getline(stream, line)),
            "Exact surface CSV ended before the center point.");
  }
  std::replace(line.begin(), line.end(), ',', ' ');
  std::istringstream parser(line);
  double x = 0.0;
  double y = 0.0;
  double sigma = 0.0;
  double sigma_vee = 0.0;
  double psi = 0.0;
  require(static_cast<bool>(parser >> x >> y >> sigma >> sigma_vee >> psi),
          "Cannot parse exact hidden-point surface row.");
  require_close(sigma, 2.0 / 3.0, 1e-15,
                "Exact surface sigma was not converted correctly.");
  require_close(sigma_vee, 1.0 / 3.0, 1e-15,
                "Exact lower envelope was not used.");
  require_close(psi, 4.0 / 3.0, 1e-15,
                "Exact lower envelope was not used for psi.");

  const auto subdivision_path = prefix.string() + "_subdivision.csv";
  std::ifstream subdivision(subdivision_path);
  require(static_cast<bool>(subdivision),
          "Cannot read exact subdivision CSV.");
  require(static_cast<bool>(std::getline(subdivision, line)),
          "Exact subdivision CSV has no header.");
  std::size_t subdivision_rows = 0;
  while (std::getline(subdivision, line)) {
    ++subdivision_rows;
  }
  require(subdivision_rows == 4,
          "Exact coplanar square subdivision should have one quadrilateral.");

  std::filesystem::remove(surface_path);
  std::filesystem::remove(prefix.string() + "_triangles.csv");
  std::filesystem::remove(subdivision_path);
  std::filesystem::remove(prefix.string() + "_ell.csv");
}

void test_triangle_shortest_vector() {
  const auto configuration =
      gkz::PointConfiguration::from_points({{0, 0}, {1, 0}, {0, 1}});
  gkz::SolverOptions options;
  options.tolerance = 1e-13;
  options.absolute_tolerance = 1e-15;
  const auto result = gkz::ShortestGkzSolver(options).solve(configuration);
  require(result.converged, "Triangle solver did not converge.");
  require(result.exact.certified, "Triangle exact certificate failed.");
  for (Eigen::Index i = 0; i < result.sigma.size(); ++i) {
    require_close(result.sigma[i], 1.0 / 3.0, 1e-14,
                  "Triangle shortest GKZ coordinate is wrong");
    require(result.exact.sigma[static_cast<std::size_t>(i)] ==
                CGAL::Gmpq(1, 3),
            "Triangle exact shortest GKZ coordinate is wrong.");
  }
}

void test_square_shortest_vector() {
  gkz::SolverOptions options;
  options.tolerance = 1e-13;
  options.absolute_tolerance = 1e-15;
  const auto result =
      gkz::ShortestGkzSolver(options).solve(square_configuration());
  require(result.converged, "Square solver did not converge.");
  require(result.exact.certified, "Square exact certificate failed.");
  require(result.active_vectors.size() == 2,
          "Square shortest point should use its two triangulations.");
  for (Eigen::Index i = 0; i < result.sigma.size(); ++i) {
    require_close(result.sigma[i], 0.25, 1e-14,
                  "Square shortest GKZ coordinate is wrong");
    require(result.exact.sigma[static_cast<std::size_t>(i)] ==
                CGAL::Gmpq(1, 4),
            "Square exact shortest GKZ coordinate is wrong.");
  }
}

void test_six_point_regression() {
  const auto configuration = gkz::PointConfiguration::from_points(
      {{0, 0}, {0, 2}, {2, 0}, {2, 2}, {1, 1}, {6, 5}});
  gkz::SolverOptions options;
  options.tolerance = 1e-12;
  options.absolute_tolerance = 1e-14;
  options.max_iterations = 200;
  const auto result = gkz::ShortestGkzSolver(options).solve(configuration);
  require(result.converged, "Six-point solver did not converge.");
  require(result.exact.certified,
          "Six-point exact certificate failed: " + result.exact.message);

  const std::vector<CGAL::Gmpq> expected = {
      CGAL::Gmpq(15, 121),   CGAL::Gmpq(5, 33),
      CGAL::Gmpq(493, 3025), CGAL::Gmpq(1511, 9075),
      CGAL::Gmpq(1318, 9075), CGAL::Gmpq(2267, 9075)};
  require(result.exact.sigma == expected,
          "Six-point exact sigma disagrees with main.tex.");
}

void test_projection_mode_preserves_six_point_solution() {
  const auto configuration = gkz::PointConfiguration::from_points(
      {{0, 0}, {0, 2}, {2, 0}, {2, 2}, {1, 1}, {6, 5}});
  gkz::SolverOptions options;
  options.tolerance = 1e-12;
  options.absolute_tolerance = 1e-14;
  options.max_iterations = 200;
  options.projection = true;
  // Force the experimental path on this small configuration. The default
  // plateau detector is intentionally much more conservative.
  options.projection_window = 1;
  options.projection_stall_ratio = 1e-4;
  options.projection_relative_gap = 1.0;
  options.projection_probe_period = 1;
  const auto result = gkz::ShortestGkzSolver(options).solve(configuration);
  require(result.converged, "Projected six-point solver did not converge.");
  require(result.exact.certified,
          "Projected six-point exact certificate failed: " +
              result.exact.message);
  require(result.projection_start_iteration >= 0,
          "Projection mode did not activate.");
  require(result.projection_probes > 0,
          "Projection mode did not issue an exact probe.");
  require(result.projection_affine_rank > 0,
          "Projection mode did not build an affine direction.");

  const std::vector<CGAL::Gmpq> expected = {
      CGAL::Gmpq(15, 121),   CGAL::Gmpq(5, 33),
      CGAL::Gmpq(493, 3025), CGAL::Gmpq(1511, 9075),
      CGAL::Gmpq(1318, 9075), CGAL::Gmpq(2267, 9075)};
  require(result.exact.sigma == expected,
          "Projection mode changed the exact six-point shortest GKZ.");
}

void test_projection_mode_batches_a_new_vertex() {
  const auto configuration = gkz::PointConfiguration::from_points(
      {{-5, 6}, {2, 1}, {4, -1}, {1, -3}, {-3, -2}, {-4, -1}});
  gkz::SolverOptions options;
  options.tolerance = 1e-13;
  options.absolute_tolerance = 1e-15;
  options.projection = true;
  options.projection_window = 1;
  options.projection_stall_ratio = 1e-4;
  options.projection_relative_gap = 1.0;
  options.projection_probe_period = 1;
  const auto result = gkz::ShortestGkzSolver(options).solve(configuration);
  require(result.converged, "Projected hexagon solver did not converge.");
  require(result.exact.certified,
          "Projected hexagon exact certificate failed: " +
              result.exact.message);
  require(result.projection_probes > 0,
          "Projected hexagon did not issue a probe.");
  require(result.projection_new_vertices > 0,
          "Projected hexagon probe did not add a new active vertex.");
}

void test_exact_qp_with_affinely_redundant_hexagon_gkz_vectors() {
  const auto configuration = gkz::PointConfiguration::from_points(
      {{0, 0}, {2, 0}, {3, 1}, {2, 3}, {0, 3}, {-1, 1}});
  const auto triangulations = triangulate_interval(0, 5);
  require(triangulations.size() == 14,
          "A convex hexagon must have 14 triangulations.");

  std::vector<gkz::GkzVector> active;
  active.reserve(triangulations.size());
  for (const auto& triangulation : triangulations) {
    active.push_back(
        gkz_from_polygon_triangulation(configuration, triangulation));
  }

  // Sigma(A) has affine dimension at most |A|-3=3. Hence these 14 GKZ
  // vectors are necessarily affinely redundant, so the old square-system
  // KKT solve was singular.
  require(active.size() > configuration.size() - 2,
          "The test active set is not forced to be affinely redundant.");
  const auto certificate =
      gkz::certify_active_set_exact(configuration, active);
  require(certificate.certified,
          "Exact QP failed on redundant hexagon GKZ vectors: " +
              certificate.message);
  require(certificate.coefficients.size() == active.size(),
          "Exact QP returned the wrong number of coefficients.");
  require(certificate.support_value == certificate.norm_squared,
          "The hexagon exact oracle equality was not verified.");
}

}  // namespace

int main() {
  try {
    test_gkz_duplicate_detection_uses_exact_area_numerators();
    test_lattice_enumeration();
    test_square_oracle_against_both_triangulations();
    test_ell_pentagon_example();
    test_ell_polygon_coordinates_are_normalized();
    test_plot_data_uses_lower_envelope_at_hidden_points();
    test_exact_plot_data_uses_exact_lower_envelope();
    test_triangle_shortest_vector();
    test_square_shortest_vector();
    test_six_point_regression();
    test_projection_mode_preserves_six_point_solution();
    test_projection_mode_batches_a_new_vertex();
    test_exact_qp_with_affinely_redundant_hexagon_gkz_vectors();
    std::cout << "All shortest-GKZ tests passed.\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "TEST FAILURE: " << error.what() << '\n';
    return 1;
  }
}
