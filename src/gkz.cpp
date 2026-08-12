#include "gkz/gkz.hpp"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/QP_functions.h>
#include <CGAL/QP_models.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Regular_triangulation_face_base_2.h>
#include <CGAL/Regular_triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <Eigen/QR>
#include <gmpxx.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace gkz {
namespace {

using ExactKernel = CGAL::Exact_predicates_exact_constructions_kernel;
using NumericKernel = CGAL::Exact_predicates_inexact_constructions_kernel;

template <typename K>
struct TriangulationTypes {
  using VertexBase0 = CGAL::Regular_triangulation_vertex_base_2<K>;
  using VertexBase =
      CGAL::Triangulation_vertex_base_with_info_2<std::size_t, K,
                                                  VertexBase0>;
  using FaceBase = CGAL::Regular_triangulation_face_base_2<K>;
  using TriangulationData =
      CGAL::Triangulation_data_structure_2<VertexBase, FaceBase>;
  using Triangulation = CGAL::Regular_triangulation_2<K, TriangulationData>;
};
using Rational = CGAL::Gmpq;

WideInt cross(const IntPoint& origin, const IntPoint& a,
              const IntPoint& b) {
  const WideInt ax = static_cast<WideInt>(a.x) - origin.x;
  const WideInt ay = static_cast<WideInt>(a.y) - origin.y;
  const WideInt bx = static_cast<WideInt>(b.x) - origin.x;
  const WideInt by = static_cast<WideInt>(b.y) - origin.y;
  return ax * by - ay * bx;
}

WideInt abs_wide(WideInt value) { return value < 0 ? -value : value; }

long double as_long_double(WideInt value) {
  return static_cast<long double>(value);
}

CGAL::Gmpz gmpz_from_wide(WideInt value) {
  using UnsignedWide = __uint128_t;
  const UnsignedWide magnitude =
      value < 0 ? UnsignedWide(0) - static_cast<UnsignedWide>(value)
                : static_cast<UnsignedWide>(value);
  const std::uint64_t limbs[2] = {static_cast<std::uint64_t>(magnitude),
                                  static_cast<std::uint64_t>(magnitude >> 64)};
  CGAL::Gmpz result;
  mpz_import(result.mpz(), 2, -1, sizeof(std::uint64_t), 0, 0, limbs);
  if (value < 0) {
    mpz_neg(result.mpz(), result.mpz());
  }
  return result;
}

Rational rational_from_wide(WideInt value) {
  return Rational(gmpz_from_wide(value));
}

double rational_to_double(const Rational& value) {
  mpq_class converted;
  mpq_set(converted.get_mpq_t(), value.mpq());
  return converted.get_d();
}

std::vector<CGAL::Gmpz> clear_height_denominators(
    const std::vector<Rational>& heights) {
  if (heights.empty()) {
    return {};
  }
  std::vector<Rational> shifted;
  shifted.reserve(heights.size());
  const Rational gauge = heights.front();
  CGAL::Gmpz common_denominator(1);
  for (const Rational& height : heights) {
    shifted.push_back(height - gauge);
    mpz_lcm(common_denominator.mpz(), common_denominator.mpz(),
            shifted.back().denominator().mpz());
  }

  std::vector<CGAL::Gmpz> integers;
  integers.reserve(shifted.size());
  CGAL::Gmpz common_divisor(0);
  for (const Rational& height : shifted) {
    const CGAL::Gmpz value =
        height.numerator() *
        (common_denominator / height.denominator());
    integers.push_back(value);
    mpz_gcd(common_divisor.mpz(), common_divisor.mpz(), value.mpz());
  }
  if (common_divisor != 0 && common_divisor != 1) {
    for (CGAL::Gmpz& value : integers) {
      value /= common_divisor;
    }
  }
  return integers;
}

Eigen::VectorXd lower_envelope_values(
    const PointConfiguration& configuration,
    const Eigen::VectorXd& heights,
    const std::vector<std::array<std::size_t, 3>>& faces) {
  if (heights.size() != static_cast<Eigen::Index>(configuration.size())) {
    throw std::invalid_argument("Lower-envelope heights have wrong size.");
  }
  Eigen::VectorXd values = Eigen::VectorXd::Constant(
      heights.size(), std::numeric_limits<double>::quiet_NaN());
  for (const auto& face : faces) {
    for (const std::size_t index : face) {
      values[static_cast<Eigen::Index>(index)] =
          heights[static_cast<Eigen::Index>(index)];
    }
  }

  const auto& points = configuration.points();
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (std::isfinite(values[static_cast<Eigen::Index>(i)])) {
      continue;
    }
    bool found = false;
    for (const auto& face : faces) {
      const auto& a = points[face[0]];
      const auto& b = points[face[1]];
      const auto& c = points[face[2]];
      const WideInt side0 = cross(a, b, points[i]);
      const WideInt side1 = cross(b, c, points[i]);
      const WideInt side2 = cross(c, a, points[i]);
      const bool inside =
          (side0 >= 0 && side1 >= 0 && side2 >= 0) ||
          (side0 <= 0 && side1 <= 0 && side2 <= 0);
      if (!inside) {
        continue;
      }

      const WideInt denominator = cross(a, b, c);
      if (denominator == 0) {
        throw std::runtime_error("A lower-envelope face is degenerate.");
      }
      const long double inverse_denominator =
          1.0L / as_long_double(denominator);
      const long double lambda0 =
          as_long_double(cross(points[i], b, c)) * inverse_denominator;
      const long double lambda1 =
          as_long_double(cross(points[i], c, a)) * inverse_denominator;
      const long double lambda2 =
          as_long_double(cross(points[i], a, b)) * inverse_denominator;
      values[static_cast<Eigen::Index>(i)] = static_cast<double>(
          lambda0 * heights[static_cast<Eigen::Index>(face[0])] +
          lambda1 * heights[static_cast<Eigen::Index>(face[1])] +
          lambda2 * heights[static_cast<Eigen::Index>(face[2])]);
      found = true;
      break;
    }
    if (!found) {
      throw std::runtime_error(
          "A point of A is not covered by the lower triangulation.");
    }
  }
  return values;
}

std::vector<Rational> lower_envelope_values_exact(
    const PointConfiguration& configuration,
    const std::vector<Rational>& heights,
    const std::vector<std::array<std::size_t, 3>>& faces) {
  if (heights.size() != configuration.size()) {
    throw std::invalid_argument(
        "Exact lower-envelope heights have wrong size.");
  }
  std::vector<Rational> values(configuration.size(), Rational(0));
  std::vector<bool> assigned(configuration.size(), false);
  for (const auto& face : faces) {
    for (const std::size_t index : face) {
      values[index] = heights[index];
      assigned[index] = true;
    }
  }

  const auto& points = configuration.points();
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (assigned[i]) {
      continue;
    }
    bool found = false;
    for (const auto& face : faces) {
      const auto& a = points[face[0]];
      const auto& b = points[face[1]];
      const auto& c = points[face[2]];
      const WideInt side0 = cross(a, b, points[i]);
      const WideInt side1 = cross(b, c, points[i]);
      const WideInt side2 = cross(c, a, points[i]);
      const bool inside =
          (side0 >= 0 && side1 >= 0 && side2 >= 0) ||
          (side0 <= 0 && side1 <= 0 && side2 <= 0);
      if (!inside) {
        continue;
      }

      const WideInt denominator = cross(a, b, c);
      if (denominator == 0) {
        throw std::runtime_error(
            "An exact lower-envelope face is degenerate.");
      }
      values[i] =
          rational_from_wide(cross(points[i], b, c)) /
              rational_from_wide(denominator) * heights[face[0]] +
          rational_from_wide(cross(points[i], c, a)) /
              rational_from_wide(denominator) * heights[face[1]] +
          rational_from_wide(cross(points[i], a, b)) /
              rational_from_wide(denominator) * heights[face[2]];
      assigned[i] = true;
      found = true;
      break;
    }
    if (!found) {
      throw std::runtime_error(
          "A point of A is not covered by the exact lower triangulation.");
    }
  }
  return values;
}

struct SubdivisionCell {
  std::vector<std::size_t> vertices;
  long double plane_a = 0.0L;
  long double plane_b = 0.0L;
  long double plane_c = 0.0L;
};

std::array<long double, 3> face_plane(
    const PointConfiguration& configuration, const Eigen::VectorXd& heights,
    const std::array<std::size_t, 3>& face) {
  const auto& p0 = configuration.points()[face[0]];
  const auto& p1 = configuration.points()[face[1]];
  const auto& p2 = configuration.points()[face[2]];
  const long double x0 = p0.x;
  const long double y0 = p0.y;
  const long double x1 = p1.x;
  const long double y1 = p1.y;
  const long double x2 = p2.x;
  const long double y2 = p2.y;
  const long double z0 = heights[static_cast<Eigen::Index>(face[0])];
  const long double z1 = heights[static_cast<Eigen::Index>(face[1])];
  const long double z2 = heights[static_cast<Eigen::Index>(face[2])];
  const long double determinant =
      (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);
  if (determinant == 0.0L) {
    throw std::runtime_error("A subdivision face is degenerate.");
  }
  const long double a =
      (z0 * (y1 - y2) + z1 * (y2 - y0) + z2 * (y0 - y1)) /
      determinant;
  const long double b =
      (z0 * (x2 - x1) + z1 * (x0 - x2) + z2 * (x1 - x0)) /
      determinant;
  const long double c = z0 - a * x0 - b * y0;
  return {a, b, c};
}

bool planes_close(const SubdivisionCell& cell,
                  const std::array<long double, 3>& plane) {
  constexpr long double tolerance = 1e-8L;
  for (std::size_t i = 0; i < 3; ++i) {
    const long double scale =
        std::max({1.0L, std::abs(cell.plane_a), std::abs(cell.plane_b),
                  std::abs(cell.plane_c), std::abs(plane[i])});
    const long double current =
        i == 0 ? cell.plane_a : (i == 1 ? cell.plane_b : cell.plane_c);
    if (std::abs(current - plane[i]) > tolerance * scale) {
      return false;
    }
  }
  return true;
}

std::vector<std::size_t> cell_boundary(
    const PointConfiguration& configuration,
    std::vector<std::size_t> indices) {
  auto& points = configuration.points();
  std::sort(indices.begin(), indices.end(), [&](std::size_t i, std::size_t j) {
    return std::tie(points[i].x, points[i].y) <
           std::tie(points[j].x, points[j].y);
  });
  indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
  if (indices.size() < 3) {
    throw std::runtime_error("A subdivision cell has fewer than 3 vertices.");
  }

  std::vector<std::size_t> hull;
  hull.reserve(indices.size() * 2);
  for (const auto index : indices) {
    while (hull.size() >= 2 &&
           cross(points[hull[hull.size() - 2]], points[hull.back()],
                 points[index]) <= 0) {
      hull.pop_back();
    }
    hull.push_back(index);
  }
  const std::size_t lower_size = hull.size();
  for (auto it = indices.rbegin() + 1; it != indices.rend(); ++it) {
    while (hull.size() > lower_size &&
           cross(points[hull[hull.size() - 2]], points[hull.back()],
                 points[*it]) <= 0) {
      hull.pop_back();
    }
    hull.push_back(*it);
  }
  hull.pop_back();
  return hull;
}

std::vector<SubdivisionCell> subdivision_cells(
    const PointConfiguration& configuration, const Eigen::VectorXd& heights,
    const std::vector<std::array<std::size_t, 3>>& faces) {
  std::vector<SubdivisionCell> cells;
  for (const auto& face : faces) {
    const auto plane = face_plane(configuration, heights, face);
    auto cell = std::find_if(cells.begin(), cells.end(),
                             [&](const SubdivisionCell& candidate) {
                               return planes_close(candidate, plane);
                             });
    if (cell == cells.end()) {
      SubdivisionCell new_cell;
      new_cell.plane_a = plane[0];
      new_cell.plane_b = plane[1];
      new_cell.plane_c = plane[2];
      new_cell.vertices.assign(face.begin(), face.end());
      cells.push_back(std::move(new_cell));
    } else {
      cell->vertices.insert(cell->vertices.end(), face.begin(), face.end());
    }
  }
  for (auto& cell : cells) {
    cell.vertices = cell_boundary(configuration, std::move(cell.vertices));
  }
  return cells;
}

std::array<Rational, 3> exact_face_plane(
    const PointConfiguration& configuration,
    const std::vector<Rational>& heights,
    const std::array<std::size_t, 3>& face) {
  const auto& p0 = configuration.points()[face[0]];
  const auto& p1 = configuration.points()[face[1]];
  const auto& p2 = configuration.points()[face[2]];
  const Rational x0 = rational_from_wide(p0.x);
  const Rational y0 = rational_from_wide(p0.y);
  const Rational x1 = rational_from_wide(p1.x);
  const Rational y1 = rational_from_wide(p1.y);
  const Rational x2 = rational_from_wide(p2.x);
  const Rational y2 = rational_from_wide(p2.y);
  const Rational determinant =
      (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);
  if (determinant == 0) {
    throw std::runtime_error("An exact subdivision face is degenerate.");
  }
  const Rational& z0 = heights[face[0]];
  const Rational& z1 = heights[face[1]];
  const Rational& z2 = heights[face[2]];
  const Rational a =
      (z0 * (y1 - y2) + z1 * (y2 - y0) + z2 * (y0 - y1)) /
      determinant;
  const Rational b =
      (z0 * (x2 - x1) + z1 * (x0 - x2) + z2 * (x1 - x0)) /
      determinant;
  return {a, b, z0 - a * x0 - b * y0};
}

std::vector<SubdivisionCell> subdivision_cells_exact(
    const PointConfiguration& configuration,
    const std::vector<Rational>& heights,
    const std::vector<std::array<std::size_t, 3>>& faces) {
  struct ExactCell {
    std::array<Rational, 3> plane;
    std::vector<std::size_t> vertices;
  };
  std::vector<ExactCell> exact_cells;
  for (const auto& face : faces) {
    const auto plane = exact_face_plane(configuration, heights, face);
    auto cell = std::find_if(
        exact_cells.begin(), exact_cells.end(),
        [&](const ExactCell& candidate) { return candidate.plane == plane; });
    if (cell == exact_cells.end()) {
      ExactCell new_cell;
      new_cell.plane = plane;
      new_cell.vertices.assign(face.begin(), face.end());
      exact_cells.push_back(std::move(new_cell));
    } else {
      cell->vertices.insert(cell->vertices.end(), face.begin(), face.end());
    }
  }

  std::vector<SubdivisionCell> cells;
  cells.reserve(exact_cells.size());
  for (auto& exact_cell : exact_cells) {
    SubdivisionCell cell;
    cell.vertices =
        cell_boundary(configuration, std::move(exact_cell.vertices));
    cells.push_back(std::move(cell));
  }
  return cells;
}

ExactKernel::FT ft_from_wide(WideInt value) {
  return ExactKernel::FT(mpq_class(to_string(value)));
}

ExactKernel::FT ft_from_rational(const Rational& value) {
  mpq_class converted;
  mpq_set(converted.get_mpq_t(), value.mpq());
  return ExactKernel::FT(std::move(converted));
}

ExactKernel::FT ft_from_gmpz(const CGAL::Gmpz& value) {
  mpq_class converted;
  mpq_set_z(converted.get_mpq_t(), value.mpz());
  return ExactKernel::FT(std::move(converted));
}

std::vector<IntPoint> convex_hull(std::vector<IntPoint> points) {
  std::sort(points.begin(), points.end(), [](const IntPoint& a,
                                             const IntPoint& b) {
    return std::tie(a.x, a.y) < std::tie(b.x, b.y);
  });
  points.erase(std::unique(points.begin(), points.end()), points.end());
  if (points.size() < 3) {
    throw std::invalid_argument(
        "A two-dimensional point configuration needs at least 3 points.");
  }

  std::vector<IntPoint> hull;
  hull.reserve(points.size() * 2);
  for (const auto& point : points) {
    while (hull.size() >= 2 &&
           cross(hull[hull.size() - 2], hull.back(), point) <= 0) {
      hull.pop_back();
    }
    hull.push_back(point);
  }
  const std::size_t lower_size = hull.size();
  for (auto it = points.rbegin() + 1; it != points.rend(); ++it) {
    while (hull.size() > lower_size &&
           cross(hull[hull.size() - 2], hull.back(), *it) <= 0) {
      hull.pop_back();
    }
    hull.push_back(*it);
  }
  hull.pop_back();
  if (hull.size() < 3) {
    throw std::invalid_argument("The point configuration is collinear.");
  }
  return hull;
}

WideInt polygon_twice_area(const std::vector<IntPoint>& polygon) {
  WideInt area = 0;
  for (std::size_t i = 0; i < polygon.size(); ++i) {
    const auto& a = polygon[i];
    const auto& b = polygon[(i + 1) % polygon.size()];
    area += static_cast<WideInt>(a.x) * b.y -
            static_cast<WideInt>(a.y) * b.x;
  }
  return abs_wide(area);
}

bool is_inside_ccw_convex_polygon(const std::vector<IntPoint>& polygon,
                                  const IntPoint& point) {
  for (std::size_t i = 0; i < polygon.size(); ++i) {
    if (cross(polygon[i], polygon[(i + 1) % polygon.size()], point) < 0) {
      return false;
    }
  }
  return true;
}

std::vector<IntPoint> read_point_file(const std::filesystem::path& path) {
  std::ifstream stream(path);
  if (!stream) {
    throw std::runtime_error("Cannot open point file: " + path.string());
  }

  std::vector<IntPoint> points;
  std::string line;
  std::size_t line_number = 0;
  while (std::getline(stream, line)) {
    ++line_number;
    if (const auto comment = line.find('#'); comment != std::string::npos) {
      line.erase(comment);
    }
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream parser(line);
    long long x = 0;
    long long y = 0;
    if (!(parser >> x)) {
      continue;
    }
    if (!(parser >> y)) {
      throw std::runtime_error("Missing y coordinate at " + path.string() +
                               ":" + std::to_string(line_number));
    }
    std::string extra;
    if (parser >> extra) {
      throw std::runtime_error("Unexpected third field at " + path.string() +
                               ":" + std::to_string(line_number));
    }
    points.push_back({static_cast<std::int64_t>(x),
                      static_cast<std::int64_t>(y)});
  }
  return points;
}

template <typename K>
GkzVector compute_oracle(
    const PointConfiguration& configuration,
    const std::vector<typename K::FT>& cached_x,
    const std::vector<typename K::FT>& cached_y,
    const std::vector<typename K::FT>& cached_squared_radius,
    const std::vector<typename K::FT>& lifted_heights,
    bool keep_faces) {
  if (lifted_heights.size() != configuration.size()) {
    throw std::invalid_argument("Oracle height vector has the wrong size.");
  }

  using Triangulation = typename TriangulationTypes<K>::Triangulation;
  Triangulation triangulation;
  const auto& points = configuration.points();
  std::vector<std::pair<typename Triangulation::Weighted_point, std::size_t>>
      weighted_points;
  weighted_points.reserve(points.size());
  for (std::size_t i = 0; i < points.size(); ++i) {
    const typename K::FT weight = cached_squared_radius[i] - lifted_heights[i];
    weighted_points.emplace_back(
        typename Triangulation::Weighted_point(
            typename K::Point_2(cached_x[i], cached_y[i]), weight),
        i);
  }
  triangulation.insert(weighted_points.begin(), weighted_points.end());

  if (triangulation.dimension() != 2 || !triangulation.is_valid()) {
    throw std::runtime_error("CGAL returned an invalid regular triangulation.");
  }

  GkzVector result;
  result.area_numerators.assign(points.size(), 0);
  std::shared_ptr<std::vector<std::array<std::size_t, 3>>> faces;
  if (keep_faces) {
    faces = std::make_shared<std::vector<std::array<std::size_t, 3>>>();
  }
  WideInt covered_area = 0;
  for (auto face = triangulation.finite_faces_begin();
       face != triangulation.finite_faces_end(); ++face) {
    const std::size_t i0 = face->vertex(0)->info();
    const std::size_t i1 = face->vertex(1)->info();
    const std::size_t i2 = face->vertex(2)->info();
    if (i0 >= points.size() || i1 >= points.size() || i2 >= points.size()) {
      throw std::runtime_error("A triangulation vertex has invalid point info.");
    }
    const WideInt area = abs_wide(cross(points[i0], points[i1], points[i2]));
    if (area == 0) {
      throw std::runtime_error("CGAL returned a degenerate finite face.");
    }
    covered_area += area;
    result.area_numerators[i0] += area;
    result.area_numerators[i1] += area;
    result.area_numerators[i2] += area;
    ++result.triangles;
    if (faces) {
      faces->push_back({i0, i1, i2});
    }
  }

  if (covered_area != configuration.twice_area()) {
    throw std::runtime_error(
        "Regular triangulation does not cover the convex hull: covered area " +
        to_string(covered_area) + ", expected " +
        to_string(configuration.twice_area()) + ".");
  }

  const WideInt denominator = 3 * configuration.twice_area();
  const WideInt numerator_sum =
      std::accumulate(result.area_numerators.begin(),
                      result.area_numerators.end(), WideInt{0});
  if (numerator_sum != denominator) {
    throw std::runtime_error("GKZ mass invariant failed inside the oracle.");
  }

  result.values = Eigen::VectorXd::Zero(
      static_cast<Eigen::Index>(configuration.size()));
  const long double denominator_ld = as_long_double(denominator);
  for (std::size_t i = 0; i < points.size(); ++i) {
    result.values[static_cast<Eigen::Index>(i)] = static_cast<double>(
        as_long_double(result.area_numerators[i]) / denominator_ld);
  }
  result.visible_vertices = triangulation.number_of_vertices();
  result.hidden_vertices = triangulation.number_of_hidden_vertices();
  result.triangulation = std::move(faces);
  return result;
}

long double dot_long_double(const Eigen::VectorXd& a,
                            const Eigen::VectorXd& b) {
  if (a.size() != b.size()) {
    throw std::invalid_argument("Dot product size mismatch.");
  }
  long double value = 0.0L;
  for (Eigen::Index i = 0; i < a.size(); ++i) {
    value += static_cast<long double>(a[i]) * b[i];
  }
  return value;
}

class AffineHullState {
 public:
  explicit AffineHullState(double rank_tolerance)
      : rank_tolerance_(rank_tolerance) {}

  bool add(const Eigen::VectorXd& point) {
    if (!has_anchor_) {
      anchor_ = point;
      basis_.resize(point.size(), 0);
      has_anchor_ = true;
      ++observations_;
      return false;
    }
    if (point.size() != anchor_.size()) {
      throw std::invalid_argument("Affine-hull point dimension mismatch.");
    }

    Eigen::VectorXd residual = point - anchor_;
    // Reorthogonalization makes rank decisions more reliable after a long
    // sequence of nearly dependent oracle vertices.
    for (int pass = 0; pass < 2; ++pass) {
      if (basis_.cols() != 0) {
        residual.noalias() -= basis_ * (basis_.transpose() * residual);
      }
    }
    ++observations_;
    const double scale = std::max(
        {1.0, anchor_.norm(), point.norm(), (point - anchor_).norm()});
    const double residual_norm = residual.norm();
    if (residual_norm <= rank_tolerance_ * scale) {
      return false;
    }

    const Eigen::Index old_rank = basis_.cols();
    basis_.conservativeResize(Eigen::NoChange, old_rank + 1);
    basis_.col(old_rank) = residual / residual_norm;
    return true;
  }

  bool has_projection() const { return has_anchor_ && rank() > 0; }

  Eigen::VectorXd projection() const {
    if (!has_projection()) {
      throw std::logic_error("Affine hull has no nonconstant direction.");
    }
    return anchor_ - basis_ * (basis_.transpose() * anchor_);
  }

  std::size_t rank() const {
    return static_cast<std::size_t>(basis_.cols());
  }

  std::size_t observations() const { return observations_; }

 private:
  double rank_tolerance_;
  bool has_anchor_ = false;
  Eigen::VectorXd anchor_;
  Eigen::MatrixXd basis_;
  std::size_t observations_ = 0;
};

Eigen::MatrixXd make_gram(const std::vector<GkzVector>& active) {
  const Eigen::Index size = static_cast<Eigen::Index>(active.size());
  Eigen::MatrixXd gram(size, size);
  for (Eigen::Index i = 0; i < size; ++i) {
    for (Eigen::Index j = 0; j <= i; ++j) {
      const double value = active[static_cast<std::size_t>(i)]
                               .values.dot(active[static_cast<std::size_t>(j)]
                                               .values);
      gram(i, j) = value;
      gram(j, i) = value;
    }
  }
  return gram;
}

Eigen::VectorXd combine(const std::vector<GkzVector>& active,
                        const Eigen::VectorXd& coefficients) {
  if (active.empty() || coefficients.size() !=
                            static_cast<Eigen::Index>(active.size())) {
    throw std::invalid_argument("Invalid active representation.");
  }
  Eigen::VectorXd result = Eigen::VectorXd::Zero(active.front().values.size());
  for (std::size_t i = 0; i < active.size(); ++i) {
    result.noalias() += coefficients[static_cast<Eigen::Index>(i)] *
                        active[i].values;
  }
  return result;
}

Eigen::VectorXd solve_active_qp(const Eigen::MatrixXd& gram,
                                Eigen::VectorXd coefficients,
                                const SolverOptions& options) {
  const Eigen::Index m = gram.rows();
  if (m == 0 || gram.cols() != m || coefficients.size() != m) {
    throw std::invalid_argument("Invalid active-set QP dimensions.");
  }
  for (Eigen::Index i = 0; i < m; ++i) {
    if (coefficients[i] < 0.0 &&
        coefficients[i] >= -options.correction_tolerance) {
      coefficients[i] = 0.0;
    }
    if (coefficients[i] < 0.0) {
      throw std::runtime_error("Active-set QP received an infeasible start.");
    }
  }
  coefficients /= coefficients.sum();

  std::vector<Eigen::Index> working;
  for (Eigen::Index i = 0; i < m; ++i) {
    if (coefficients[i] > options.prune_tolerance) {
      working.push_back(i);
    }
  }
  if (working.empty()) {
    Eigen::Index maximum = 0;
    coefficients.maxCoeff(&maximum);
    working.push_back(maximum);
  }

  for (int step = 0; step < options.max_correction_steps; ++step) {
    const Eigen::Index r = static_cast<Eigen::Index>(working.size());
    Eigen::MatrixXd kkt = Eigen::MatrixXd::Zero(r + 1, r + 1);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(r + 1);
    rhs[r] = 1.0;
    for (Eigen::Index i = 0; i < r; ++i) {
      for (Eigen::Index j = 0; j < r; ++j) {
        kkt(i, j) = gram(working[static_cast<std::size_t>(i)],
                         working[static_cast<std::size_t>(j)]);
      }
      kkt(i, r) = -1.0;
      kkt(r, i) = 1.0;
    }

    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> decomposition(kkt);
    const Eigen::VectorXd solution = decomposition.solve(rhs);
    if ((kkt * solution - rhs).norm() > 1e-8) {
      throw std::runtime_error("Active-set KKT system is numerically singular.");
    }
    const Eigen::VectorXd affine = solution.head(r);

    if (affine.minCoeff() < -options.correction_tolerance) {
      double theta = 1.0;
      for (Eigen::Index i = 0; i < r; ++i) {
        if (affine[i] < 0.0) {
          const double current =
              coefficients[working[static_cast<std::size_t>(i)]];
          theta = std::min(theta, current / (current - affine[i]));
        }
      }
      for (Eigen::Index i = 0; i < r; ++i) {
        const auto index = working[static_cast<std::size_t>(i)];
        coefficients[index] =
            (1.0 - theta) * coefficients[index] + theta * affine[i];
        if (coefficients[index] <= options.prune_tolerance) {
          coefficients[index] = 0.0;
        }
      }
      working.erase(
          std::remove_if(working.begin(), working.end(), [&](Eigen::Index i) {
            return coefficients[i] <= options.prune_tolerance;
          }),
          working.end());
      if (working.empty()) {
        throw std::runtime_error("Active-set correction lost feasibility.");
      }
      continue;
    }

    coefficients.setZero();
    for (Eigen::Index i = 0; i < r; ++i) {
      coefficients[working[static_cast<std::size_t>(i)]] =
          std::max(0.0, affine[i]);
    }
    coefficients /= coefficients.sum();

    const Eigen::VectorXd gradient = gram * coefficients;
    const double multiplier = coefficients.dot(gradient);
    Eigen::Index entering = -1;
    double smallest_gradient = std::numeric_limits<double>::infinity();
    for (Eigen::Index candidate = 0; candidate < m; ++candidate) {
      if (std::find(working.begin(), working.end(), candidate) !=
          working.end()) {
        continue;
      }
      if (gradient[candidate] < smallest_gradient) {
        smallest_gradient = gradient[candidate];
        entering = candidate;
      }
    }
    const double kkt_tolerance =
        options.correction_tolerance * std::max(1.0, std::abs(multiplier));
    if (entering < 0 || smallest_gradient >= multiplier - kkt_tolerance) {
      return coefficients;
    }
    working.push_back(entering);
  }
  throw std::runtime_error("Active-set QP exceeded its correction-step limit.");
}

void prune_active(std::vector<GkzVector>& active,
                  Eigen::VectorXd& coefficients, Eigen::MatrixXd& gram,
                  double tolerance) {
  std::vector<Eigen::Index> keep;
  for (Eigen::Index i = 0; i < coefficients.size(); ++i) {
    if (coefficients[i] > tolerance) {
      keep.push_back(i);
    }
  }
  if (keep.empty()) {
    Eigen::Index maximum = 0;
    coefficients.maxCoeff(&maximum);
    keep.push_back(maximum);
  }

  std::vector<GkzVector> reduced;
  reduced.reserve(keep.size());
  Eigen::VectorXd reduced_coefficients(
      static_cast<Eigen::Index>(keep.size()));
  Eigen::MatrixXd reduced_gram(keep.size(), keep.size());
  for (std::size_t i = 0; i < keep.size(); ++i) {
    reduced.push_back(std::move(active[static_cast<std::size_t>(keep[i])]));
    reduced_coefficients[static_cast<Eigen::Index>(i)] = coefficients[keep[i]];
    for (std::size_t j = 0; j < keep.size(); ++j) {
      reduced_gram(static_cast<Eigen::Index>(i),
                   static_cast<Eigen::Index>(j)) = gram(keep[i], keep[j]);
    }
  }
  reduced_coefficients /= reduced_coefficients.sum();
  active = std::move(reduced);
  coefficients = std::move(reduced_coefficients);
  gram = std::move(reduced_gram);
}

bool solve_exact_linear_system(std::vector<std::vector<Rational>>& matrix,
                               std::vector<Rational>& rhs,
                               std::vector<Rational>& solution) {
  const std::size_t n = matrix.size();
  if (rhs.size() != n) {
    return false;
  }
  for (const auto& row : matrix) {
    if (row.size() != n) {
      return false;
    }
  }

  for (std::size_t column = 0; column < n; ++column) {
    std::size_t pivot = column;
    while (pivot < n && matrix[pivot][column] == 0) {
      ++pivot;
    }
    if (pivot == n) {
      return false;
    }
    if (pivot != column) {
      std::swap(matrix[pivot], matrix[column]);
      std::swap(rhs[pivot], rhs[column]);
    }
    const Rational pivot_value = matrix[column][column];
    for (std::size_t j = column; j < n; ++j) {
      matrix[column][j] /= pivot_value;
    }
    rhs[column] /= pivot_value;

    for (std::size_t row = 0; row < n; ++row) {
      if (row == column || matrix[row][column] == 0) {
        continue;
      }
      const Rational factor = matrix[row][column];
      for (std::size_t j = column; j < n; ++j) {
        matrix[row][j] -= factor * matrix[column][j];
      }
      rhs[row] -= factor * rhs[column];
    }
  }
  solution = rhs;
  return true;
}

std::array<Rational, 2> exact_polygon_centroid(
    const PointConfiguration& configuration) {
  Rational signed_area2 = 0;
  Rational first_moment_x = 0;
  Rational first_moment_y = 0;
  const auto& hull = configuration.hull();
  for (std::size_t i = 0; i < hull.size(); ++i) {
    const auto& a = hull[i];
    const auto& b = hull[(i + 1) % hull.size()];
    const Rational ax = rational_from_wide(a.x);
    const Rational ay = rational_from_wide(a.y);
    const Rational bx = rational_from_wide(b.x);
    const Rational by = rational_from_wide(b.y);
    const Rational edge_cross = ax * by - ay * bx;
    signed_area2 += edge_cross;
    first_moment_x += (ax + bx) * edge_cross;
    first_moment_y += (ay + by) * edge_cross;
  }
  if (signed_area2 == 0) {
    throw std::invalid_argument("Cannot compute ell_A on a zero-area hull.");
  }
  return {first_moment_x / (3 * signed_area2),
          first_moment_y / (3 * signed_area2)};
}

ExactCertificate solve_exact_active_qp(
    const PointConfiguration& configuration,
    const std::vector<GkzVector>& active) {
  ExactCertificate certificate;
  if (active.empty()) {
    certificate.message = "The active set is empty.";
    return certificate;
  }

  const std::size_t m = active.size();
  const std::size_t n = configuration.size();
  const Rational denominator =
      rational_from_wide(3 * configuration.twice_area());

  using ExactInteger = CGAL::Gmpz;
  using ExactProgram = CGAL::Quadratic_program<ExactInteger>;
  ExactProgram program(CGAL::EQUAL, true, ExactInteger(0), false,
                       ExactInteger(0));
  program.set_b(0, ExactInteger(1));
  for (const auto& vector : active) {
    if (vector.area_numerators.size() != n) {
      certificate.message =
          "An active GKZ vector has the wrong coordinate count.";
      return certificate;
    }
  }

  // The Gram entries are sums of n products of area numerators. These fit
  // comfortably into a 128-bit accumulator for realistic inputs, so the hot
  // path avoids per-coordinate Gmpz construction; inputs large enough to
  // threaten 2^120 fall back to the exact big-integer loop.
  WideInt max_numerator = 0;
  for (const auto& vector : active) {
    for (const WideInt numerator : vector.area_numerators) {
      max_numerator = std::max(max_numerator, numerator);
    }
  }
  const long double gram_entry_bound =
      as_long_double(max_numerator) * as_long_double(max_numerator) *
      static_cast<long double>(n);
  const bool wide_gram = gram_entry_bound < 0x1p+120L;

  for (std::size_t i = 0; i < m; ++i) {
    program.set_a(static_cast<int>(i), 0, ExactInteger(1));
    for (std::size_t j = 0; j <= i; ++j) {
      ExactInteger inner;
      if (wide_gram) {
        WideInt accumulated = 0;
        const auto& numerators_i = active[i].area_numerators;
        const auto& numerators_j = active[j].area_numerators;
        for (std::size_t coordinate = 0; coordinate < n; ++coordinate) {
          accumulated += numerators_i[coordinate] * numerators_j[coordinate];
        }
        inner = gmpz_from_wide(accumulated);
      } else {
        inner = ExactInteger(0);
        for (std::size_t coordinate = 0; coordinate < n; ++coordinate) {
          inner += gmpz_from_wide(active[i].area_numerators[coordinate]) *
                   gmpz_from_wide(active[j].area_numerators[coordinate]);
        }
      }
      // CGAL minimizes c^T lambda + 1/2 lambda^T D lambda.
      program.set_d(static_cast<int>(i), static_cast<int>(j), 2 * inner);
    }
  }

  auto solution = CGAL::solve_quadratic_program(program, ExactInteger());
  if (!solution.is_optimal()) {
    certificate.message = solution.is_infeasible()
                              ? "The exact active QP is infeasible."
                              : "The exact active QP did not return an "
                                "optimal solution.";
    return certificate;
  }
  if (!solution.solves_quadratic_program(program)) {
    certificate.message =
        "CGAL could not validate the exact active QP solution.";
    return certificate;
  }

  certificate.coefficients.reserve(m);
  auto value = solution.variable_values_begin();
  for (std::size_t i = 0; i < m; ++i, ++value) {
    certificate.coefficients.emplace_back(value->numerator(),
                                          value->denominator());
  }
  Rational coefficient_sum(0);
  for (const Rational& coefficient : certificate.coefficients) {
    if (coefficient < 0) {
      certificate.message = "The exact QP returned a negative coefficient.";
      return certificate;
    }
    coefficient_sum += coefficient;
  }
  if (coefficient_sum != 1) {
    certificate.message =
        "The exact QP coefficients do not sum to one.";
    return certificate;
  }

  std::vector<Rational> unnormalized_sigma(n, Rational(0));
  certificate.sigma.assign(n, Rational(0));
  for (std::size_t coordinate = 0; coordinate < n; ++coordinate) {
    for (std::size_t i = 0; i < m; ++i) {
      unnormalized_sigma[coordinate] +=
          certificate.coefficients[i] *
          rational_from_wide(active[i].area_numerators[coordinate]);
    }
    certificate.sigma[coordinate] =
        unnormalized_sigma[coordinate] / denominator;
  }

  certificate.norm_squared = 0;
  for (const auto& value : certificate.sigma) {
    certificate.norm_squared += value * value;
  }

  const std::vector<CGAL::Gmpz> exact_heights =
      clear_height_denominators(unnormalized_sigma);
  RegularTriangulationOracle oracle(configuration);
  const GkzVector minimizing_vertex =
      oracle.minimize_exact_integer(exact_heights, /*keep_faces=*/false);
  certificate.support_value = 0;
  for (std::size_t coordinate = 0; coordinate < n; ++coordinate) {
    certificate.support_value +=
        certificate.sigma[coordinate] *
        rational_from_wide(
            minimizing_vertex.area_numerators[coordinate]) /
        denominator;
  }
  if (certificate.support_value != certificate.norm_squared) {
    certificate.message =
        "Exact oracle found a positive Frank--Wolfe gap.";
    certificate.has_witness = true;
    certificate.witness = std::move(minimizing_vertex);
    return certificate;
  }
  certificate.certified = true;
  certificate.message =
      "Exact convex QP and equality H(x) = ||x||^2 verified.";
  return certificate;
}

std::string rational_to_string(const Rational& value) {
  std::ostringstream stream;
  if (value.denominator() == 1) {
    stream << value.numerator();
  } else {
    stream << value;
  }
  return stream.str();
}

}  // namespace

struct RegularTriangulationOracle::Cache {
  std::vector<NumericKernel::FT> numeric_x;
  std::vector<NumericKernel::FT> numeric_y;
  std::vector<NumericKernel::FT> numeric_squared_radius;
  std::vector<ExactKernel::FT> exact_x;
  std::vector<ExactKernel::FT> exact_y;
  std::vector<ExactKernel::FT> exact_squared_radius;
  long double coordinate_scale = 1.0L;
};

RegularTriangulationOracle::RegularTriangulationOracle(
    const PointConfiguration& configuration)
    : configuration_(configuration), cache_(std::make_unique<Cache>()) {
  const auto& points = configuration.points();
  cache_->numeric_x.reserve(points.size());
  cache_->numeric_y.reserve(points.size());
  cache_->numeric_squared_radius.reserve(points.size());
  cache_->exact_x.reserve(points.size());
  cache_->exact_y.reserve(points.size());
  cache_->exact_squared_radius.reserve(points.size());
  for (const auto& point : points) {
    const ExactKernel::FT px = ft_from_wide(point.x);
    const ExactKernel::FT py = ft_from_wide(point.y);
    cache_->exact_x.push_back(px);
    cache_->exact_y.push_back(py);
    cache_->exact_squared_radius.push_back(px * px + py * py);
    const double nx = static_cast<double>(point.x);
    const double ny = static_cast<double>(point.y);
    cache_->numeric_x.push_back(nx);
    cache_->numeric_y.push_back(ny);
    cache_->numeric_squared_radius.push_back(nx * nx + ny * ny);
    cache_->coordinate_scale =
        std::max(cache_->coordinate_scale,
                 static_cast<long double>(point.x) * point.x +
                     static_cast<long double>(point.y) * point.y);
  }
}

RegularTriangulationOracle::~RegularTriangulationOracle() = default;

ExactCertificate certify_active_set_exact(
    const PointConfiguration& configuration,
    const std::vector<GkzVector>& active) {
  return solve_exact_active_qp(configuration, active);
}

std::string to_string(WideInt value) {
  if (value == 0) {
    return "0";
  }
  const bool negative = value < 0;
  using UnsignedWide = __uint128_t;
  UnsignedWide magnitude = negative
                               ? static_cast<UnsignedWide>(-(value + 1)) + 1
                               : static_cast<UnsignedWide>(value);
  std::string digits;
  while (magnitude > 0) {
    digits.push_back(static_cast<char>('0' + magnitude % 10));
    magnitude /= 10;
  }
  if (negative) {
    digits.push_back('-');
  }
  std::reverse(digits.begin(), digits.end());
  return digits;
}

PointConfiguration::PointConfiguration(std::vector<IntPoint> points,
                                       std::vector<IntPoint> hull,
                                       WideInt twice_area,
                                       std::int64_t level,
                                       WideInt base_twice_area)
    : points_(std::move(points)),
      hull_(std::move(hull)),
      twice_area_(twice_area),
      level_(level),
      base_twice_area_(base_twice_area) {}

PointConfiguration PointConfiguration::from_points(
    std::vector<IntPoint> points) {
  std::set<std::pair<std::int64_t, std::int64_t>> unique;
  for (const auto& point : points) {
    if (!unique.emplace(point.x, point.y).second) {
      throw std::invalid_argument("Point configuration contains duplicates.");
    }
  }
  const auto hull = convex_hull(points);
  const WideInt area = polygon_twice_area(hull);
  if (area == 0) {
    throw std::invalid_argument("Point configuration has zero area.");
  }
  return PointConfiguration(std::move(points), hull, area);
}

PointConfiguration PointConfiguration::from_points_file(
    const std::filesystem::path& path) {
  return from_points(read_point_file(path));
}

PointConfiguration PointConfiguration::from_lattice_polygon(
    const std::vector<IntPoint>& polygon_vertices, std::int64_t k) {
  if (k <= 0) {
    throw std::invalid_argument("k must be a positive integer.");
  }
  const auto base_hull = convex_hull(polygon_vertices);
  const WideInt base_area = polygon_twice_area(base_hull);
  std::vector<IntPoint> scaled_hull;
  scaled_hull.reserve(base_hull.size());
  for (const auto& point : base_hull) {
    const WideInt x = static_cast<WideInt>(point.x) * k;
    const WideInt y = static_cast<WideInt>(point.y) * k;
    if (x < std::numeric_limits<std::int64_t>::min() ||
        x > std::numeric_limits<std::int64_t>::max() ||
        y < std::numeric_limits<std::int64_t>::min() ||
        y > std::numeric_limits<std::int64_t>::max()) {
      throw std::overflow_error("Scaled polygon coordinate overflows int64.");
    }
    scaled_hull.push_back(
        {static_cast<std::int64_t>(x), static_cast<std::int64_t>(y)});
  }

  std::int64_t minimum_x = scaled_hull.front().x;
  std::int64_t maximum_x = scaled_hull.front().x;
  std::int64_t minimum_y = scaled_hull.front().y;
  std::int64_t maximum_y = scaled_hull.front().y;
  for (const auto& point : scaled_hull) {
    minimum_x = std::min(minimum_x, point.x);
    maximum_x = std::max(maximum_x, point.x);
    minimum_y = std::min(minimum_y, point.y);
    maximum_y = std::max(maximum_y, point.y);
  }

  std::vector<IntPoint> lattice_points;
  for (std::int64_t x = minimum_x;; ++x) {
    for (std::int64_t y = minimum_y;; ++y) {
      const IntPoint point{x, y};
      if (is_inside_ccw_convex_polygon(scaled_hull, point)) {
        lattice_points.push_back(point);
      }
      if (y == maximum_y) {
        break;
      }
    }
    if (x == maximum_x) {
      break;
    }
  }
  PointConfiguration configuration = from_points(std::move(lattice_points));
  configuration.level_ = k;
  configuration.base_twice_area_ = base_area;
  return configuration;
}

PointConfiguration PointConfiguration::from_polygon_file(
    const std::filesystem::path& path, std::int64_t k) {
  return from_lattice_polygon(read_point_file(path), k);
}

double PointConfiguration::centroid_x() const {
  long double numerator = 0.0L;
  long double signed_area = 0.0L;
  for (std::size_t i = 0; i < hull_.size(); ++i) {
    const auto& a = hull_[i];
    const auto& b = hull_[(i + 1) % hull_.size()];
    const long double term =
        static_cast<long double>(a.x) * b.y -
        static_cast<long double>(a.y) * b.x;
    signed_area += term;
    numerator += (static_cast<long double>(a.x) + b.x) * term;
  }
  return static_cast<double>(numerator / (3.0L * signed_area));
}

double PointConfiguration::centroid_y() const {
  long double numerator = 0.0L;
  long double signed_area = 0.0L;
  for (std::size_t i = 0; i < hull_.size(); ++i) {
    const auto& a = hull_[i];
    const auto& b = hull_[(i + 1) % hull_.size()];
    const long double term =
        static_cast<long double>(a.x) * b.y -
        static_cast<long double>(a.y) * b.x;
    signed_area += term;
    numerator += (static_cast<long double>(a.y) + b.y) * term;
  }
  return static_cast<double>(numerator / (3.0L * signed_area));
}

AffineFunction compute_ell(const PointConfiguration& configuration) {
  const std::size_t n = configuration.size();
  if (n < 3) {
    throw std::invalid_argument("ell_A needs a two-dimensional point set.");
  }

  std::vector<std::vector<Rational>> matrix(
      3, std::vector<Rational>(3, Rational(0)));
  const Rational coordinate_denominator =
      configuration.is_polygon_level()
          ? rational_from_wide(configuration.level())
          : Rational(1);
  for (const auto& point : configuration.points()) {
    const Rational x = rational_from_wide(point.x) / coordinate_denominator;
    const Rational y = rational_from_wide(point.y) / coordinate_denominator;
    matrix[0][0] += 1;
    matrix[0][1] += x;
    matrix[0][2] += y;
    matrix[1][0] += x;
    matrix[1][1] += x * x;
    matrix[1][2] += x * y;
    matrix[2][0] += y;
    matrix[2][1] += x * y;
    matrix[2][2] += y * y;
  }

  auto centroid = exact_polygon_centroid(configuration);
  centroid[0] /= coordinate_denominator;
  centroid[1] /= coordinate_denominator;
  std::vector<Rational> rhs = {Rational(1), centroid[0], centroid[1]};
  std::vector<Rational> solution;
  if (!solve_exact_linear_system(matrix, rhs, solution) || solution.size() != 3) {
    throw std::runtime_error("The affine moment system for ell_A is singular.");
  }

  AffineFunction result;
  for (std::size_t i = 0; i < 3; ++i) {
    result.exact_coefficients[i] = solution[i];
    result.coefficients[i] = rational_to_double(solution[i]);
  }
  result.values = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(n));
  result.exact_values.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    const Rational x =
        rational_from_wide(configuration.points()[i].x) /
        coordinate_denominator;
    const Rational y =
        rational_from_wide(configuration.points()[i].y) /
        coordinate_denominator;
    const Rational value = solution[0] + solution[1] * x + solution[2] * y;
    result.exact_values.push_back(value);
    result.values[static_cast<Eigen::Index>(i)] = rational_to_double(value);
  }
  return result;
}

GkzVector RegularTriangulationOracle::minimize(
    const Eigen::VectorXd& heights, bool keep_faces) const {
  if (heights.size() != static_cast<Eigen::Index>(configuration_.size())) {
    throw std::invalid_argument("Oracle height vector has the wrong size.");
  }
  if (!heights.allFinite()) {
    throw std::invalid_argument("Oracle heights must be finite.");
  }

  const double gauge = heights[0];
  double maximum_deviation = 0.0;
  for (Eigen::Index i = 0; i < heights.size(); ++i) {
    maximum_deviation =
        std::max(maximum_deviation, std::abs(heights[i] - gauge));
  }
  double scale = 1.0;
  if (maximum_deviation > 0.0) {
    const double ratio = static_cast<double>(cache_->coordinate_scale) /
                         maximum_deviation;
    if (std::isfinite(ratio) && ratio > 0.0) {
      const int exponent = std::clamp(
          static_cast<int>(std::floor(std::log2(ratio))), -500, 500);
      scale = std::ldexp(1.0, exponent);
    }
  }

  std::vector<NumericKernel::FT> numeric_heights;
  numeric_heights.reserve(configuration_.size());
  for (Eigen::Index i = 0; i < heights.size(); ++i) {
    numeric_heights.push_back((heights[i] - gauge) * scale);
  }
  return compute_oracle<NumericKernel>(
      configuration_, cache_->numeric_x, cache_->numeric_y,
      cache_->numeric_squared_radius, numeric_heights, keep_faces);
}

GkzVector RegularTriangulationOracle::minimize_exact(
    const std::vector<CGAL::Gmpq>& heights, bool keep_faces) const {
  if (heights.size() != configuration_.size()) {
    throw std::invalid_argument("Exact oracle height vector has wrong size.");
  }
  std::vector<ExactKernel::FT> exact_heights;
  exact_heights.reserve(heights.size());
  const ExactKernel::FT gauge = ft_from_rational(heights.front());
  for (const auto& height : heights) {
    exact_heights.push_back(ft_from_rational(height) - gauge);
  }
  return compute_oracle<ExactKernel>(
      configuration_, cache_->exact_x, cache_->exact_y,
      cache_->exact_squared_radius, exact_heights, keep_faces);
}

GkzVector RegularTriangulationOracle::minimize_exact_integer(
    const std::vector<CGAL::Gmpz>& heights, bool keep_faces) const {
  if (heights.size() != configuration_.size()) {
    throw std::invalid_argument(
        "Exact integer oracle heights have wrong size.");
  }
  std::vector<ExactKernel::FT> exact_heights;
  exact_heights.reserve(heights.size());
  const ExactKernel::FT gauge = ft_from_gmpz(heights.front());
  for (const auto& height : heights) {
    exact_heights.push_back(ft_from_gmpz(height) - gauge);
  }
  return compute_oracle<ExactKernel>(
      configuration_, cache_->exact_x, cache_->exact_y,
      cache_->exact_squared_radius, exact_heights, keep_faces);
}

SolverResult ShortestGkzSolver::solve(
    const PointConfiguration& configuration) const {
  if (options_.projection_window <= 0 ||
      options_.projection_rank_stall_window <= 0 ||
      options_.projection_stall_ratio <= 0.0 ||
      options_.projection_stall_ratio > 1.0 ||
      options_.projection_relative_gap <= 0.0 ||
      options_.projection_relative_gap > 1.0 ||
      options_.projection_rank_tolerance <= 0.0) {
    throw std::invalid_argument("Invalid projection solver options.");
  }
  RegularTriangulationOracle oracle(configuration);
  const Eigen::VectorXd zero = Eigen::VectorXd::Zero(
      static_cast<Eigen::Index>(configuration.size()));

  std::vector<GkzVector> active;
  active.push_back(oracle.minimize(zero, /*keep_faces=*/false));
  Eigen::VectorXd coefficients = Eigen::VectorXd::Ones(1);
  Eigen::MatrixXd gram = make_gram(active);

  SolverResult result;
  result.projection_enabled = options_.projection;
  // The numerical EPICK oracle is fast but may return a non-minimizing
  // triangulation near a secondary-fan wall, which can understate the gap.
  // Once the numerical gap first meets the stopping threshold, the loop
  // switches to the exact EPECK oracle (the "exact endgame"): every
  // subsequent gap value and certification witness is then exact, restoring
  // the same convergence guarantees as a fully exact run.
  bool exact_endgame = false;

  const auto rational_heights = [&](const Eigen::VectorXd& heights) {
    std::vector<CGAL::Gmpq> exact_heights(configuration.size());
    for (std::size_t i = 0; i < configuration.size(); ++i) {
      exact_heights[i] = CGAL::Gmpq(heights[static_cast<Eigen::Index>(i)]);
    }
    return exact_heights;
  };
  const auto active_contains = [&](const GkzVector& vertex) {
    for (const auto& vector : active) {
      if (vector.has_same_area_numerators(vertex)) {
        return true;
      }
    }
    return false;
  };
  const auto expand_active = [&](GkzVector required_vertex,
                                 std::optional<GkzVector> optional_vertex,
                                 bool prune) {
    if (active_contains(required_vertex)) {
      throw std::runtime_error(
          "Oracle returned an exact duplicate active GKZ vector while the "
          "gap exceeds the stopping threshold.");
    }
    std::vector<GkzVector> appended;
    appended.push_back(std::move(required_vertex));
    bool optional_added = false;
    if (optional_vertex.has_value() && !active_contains(*optional_vertex) &&
        !appended.front().has_same_area_numerators(*optional_vertex)) {
      appended.push_back(std::move(*optional_vertex));
      optional_added = true;
    }

    const Eigen::Index old_size = static_cast<Eigen::Index>(active.size());
    const Eigen::Index new_size =
        old_size + static_cast<Eigen::Index>(appended.size());
    active.reserve(static_cast<std::size_t>(new_size));
    for (auto& vertex : appended) {
      active.push_back(std::move(vertex));
    }
    Eigen::MatrixXd expanded =
        Eigen::MatrixXd::Zero(new_size, new_size);
    expanded.topLeftCorner(old_size, old_size) = gram;
    for (Eigen::Index i = old_size; i < new_size; ++i) {
      for (Eigen::Index j = 0; j <= i; ++j) {
        const double value = active[static_cast<std::size_t>(i)].values.dot(
            active[static_cast<std::size_t>(j)].values);
        expanded(i, j) = value;
        expanded(j, i) = value;
      }
    }
    gram = std::move(expanded);
    Eigen::VectorXd expanded_coefficients = Eigen::VectorXd::Zero(new_size);
    expanded_coefficients.head(old_size) = coefficients;
    coefficients = solve_active_qp(gram, expanded_coefficients, options_);
    if (prune) {
      prune_active(active, coefficients, gram, options_.prune_tolerance);
    }
    return optional_added;
  };

  double initial_gap = -1.0;
  std::deque<double> recent_gaps;
  bool stable_collection_mode = false;
  int rank_stall_observations = 0;
  std::optional<AffineHullState> affine_hull;
  const std::size_t stable_rank_limit =
      configuration.size() >= 4 ? configuration.size() - 4 : 0;

  for (int iteration = 0; iteration <= options_.max_iterations; ++iteration) {
    const Eigen::VectorXd candidate = combine(active, coefficients);
    GkzVector minimizing_vertex =
        exact_endgame
            ? oracle.minimize_exact(rational_heights(candidate),
                                    /*keep_faces=*/false)
            : oracle.minimize(candidate, /*keep_faces=*/false);
    const long double norm_squared = dot_long_double(candidate, candidate);
    const long double support =
        dot_long_double(candidate, minimizing_vertex.values);
    long double gap = norm_squared - support;
    const long double roundoff =
        64.0L * std::numeric_limits<double>::epsilon() *
        std::max(1.0L, std::abs(norm_squared));
    if (gap < 0.0L && gap >= -roundoff) {
      gap = 0.0L;
    }
    if (gap < 0.0L) {
      throw std::runtime_error(
          "Oracle gap is negative beyond the floating-point error budget.");
    }

    result.iterations = iteration;
    result.norm_squared = static_cast<double>(norm_squared);
    result.gap = static_cast<double>(gap);
    result.l2_error_bound = std::sqrt(2.0 * result.gap);
    if (initial_gap < 0.0) {
      initial_gap = result.gap;
    }

    const double stopping_threshold =
        options_.absolute_tolerance +
        options_.tolerance * std::max(1e-300, result.norm_squared);
    if (result.gap <= stopping_threshold && !exact_endgame) {
      // Numerical gap met the threshold; redo this gap check with the
      // exact oracle before trusting it. This happens before the verbose
      // log so each iteration number is logged at most once.
      exact_endgame = true;
      --iteration;
      continue;
    }

    if (options_.projection && !stable_collection_mode) {
      recent_gaps.push_back(result.gap);
      const std::size_t history_size =
          static_cast<std::size_t>(options_.projection_window) + 1;
      if (recent_gaps.size() > history_size) {
        recent_gaps.pop_front();
      }
      if (recent_gaps.size() == history_size &&
          result.gap <= initial_gap * options_.projection_relative_gap &&
          result.gap >= recent_gaps.front() *
                            options_.projection_stall_ratio) {
        stable_collection_mode = true;
        affine_hull.emplace(options_.projection_rank_tolerance);
        result.projection_start_iteration = iteration;
        if (options_.verbose) {
          std::cerr << "projection event=enabled iteration=" << iteration
                    << " gap=" << std::setprecision(17) << result.gap
                    << " window_gap=" << recent_gaps.front() << '\n';
        }
      }
    }
    if (options_.verbose) {
      std::cerr << "iteration=" << iteration << " active=" << active.size()
                << " norm2=" << std::setprecision(17)
                << result.norm_squared << " gap=" << result.gap << '\n';
    }
    if (result.gap <= stopping_threshold) {
      result.converged = true;
      result.sigma = candidate;
      result.active_vectors = active;
      result.coefficients = coefficients;
      if (options_.exact_certification) {
        if (options_.exact_max_active > 0 &&
            active.size() >
                static_cast<std::size_t>(options_.exact_max_active)) {
          result.exact.message =
              "Exact certification skipped because active_size exceeds "
              "exact_max_active.";
        } else {
          result.exact = certify_active_set_exact(configuration, active);
          if (!result.exact.certified && result.exact.has_witness &&
              iteration < options_.max_iterations) {
            // Certification's exact oracle found a better vertex; use it as
            // the next Frank--Wolfe direction and keep iterating. The
            // witness must not be pruned away: the numerical QP assigns it
            // a near-zero coefficient because the current candidate is
            // already numerically optimal, but the exact QP in the next
            // certification needs it to represent the true minimizer.
            expand_active(std::move(result.exact.witness), std::nullopt,
                          /*prune=*/false);
            result.exact = ExactCertificate();
            result.converged = false;
            continue;
          }
        }
      }
      return result;
    }
    if (iteration == options_.max_iterations) {
      result.sigma = candidate;
      result.active_vectors = active;
      result.coefficients = coefficients;
      if (options_.projection && !result.stable_projection_available) {
        result.stable_projection_stop_reason =
            stable_collection_mode ? "max_iterations_before_stable_rank"
                                   : "max_iterations_before_trigger";
      }
      break;
    }

    if (stable_collection_mode) {
      const bool rank_increased = affine_hull->add(minimizing_vertex.values);
      if (rank_increased) {
        rank_stall_observations = 0;
      } else {
        ++rank_stall_observations;
      }
      if (options_.verbose) {
        std::cerr << "projection event=rank iteration=" << iteration
                  << " rank=" << affine_hull->rank()
                  << " observations=" << affine_hull->observations()
                  << " stall=" << rank_stall_observations << '\n';
      }

      const bool reached_rank_limit =
          stable_rank_limit > 0 && affine_hull->rank() >= stable_rank_limit;
      const bool rank_stalled =
          affine_hull->rank() > 0 &&
          rank_stall_observations >= options_.projection_rank_stall_window;
      if (reached_rank_limit || rank_stalled) {
        const Eigen::VectorXd projection = affine_hull->projection();
        result.stable_projection_available = true;
        result.stable_projection = projection;
        result.stable_projection_norm_squared =
            static_cast<double>(dot_long_double(projection, projection));
        result.stable_projection_stop_iteration = iteration;
        result.stable_projection_rank = affine_hull->rank();
        result.stable_projection_observations = affine_hull->observations();
        result.stable_projection_stop_reason =
            reached_rank_limit ? "rank_limit" : "rank_stall";

        GkzVector projection_vertex = oracle.minimize_exact(
            rational_heights(projection), /*keep_faces=*/false);
        const bool projection_vertex_new =
            !active_contains(projection_vertex) &&
            !minimizing_vertex.has_same_area_numerators(projection_vertex);
        if (options_.verbose) {
          std::cerr << "projection event=stable iteration=" << iteration
                    << " rank=" << affine_hull->rank()
                    << " observations=" << affine_hull->observations()
                    << " reason=" << result.stable_projection_stop_reason
                    << " p_norm2=" << std::setprecision(17)
                    << result.stable_projection_norm_squared
                    << " vertex_new=" << std::boolalpha
                    << projection_vertex_new << '\n';
        }

        expand_active(std::move(minimizing_vertex),
                      std::move(projection_vertex), /*prune=*/true);
        const Eigen::VectorXd final_candidate = combine(active, coefficients);
        const GkzVector final_vertex =
            oracle.minimize(final_candidate, /*keep_faces=*/false);
        const long double final_norm_squared =
            dot_long_double(final_candidate, final_candidate);
        long double final_gap =
            final_norm_squared - dot_long_double(final_candidate,
                                                  final_vertex.values);
        const long double final_roundoff =
            64.0L * std::numeric_limits<double>::epsilon() *
            std::max(1.0L, std::abs(final_norm_squared));
        if (final_gap < 0.0L && final_gap >= -final_roundoff) {
          final_gap = 0.0L;
        }
        if (final_gap < 0.0L) {
          throw std::runtime_error(
              "Final QP oracle gap is negative beyond the floating-point "
              "error budget.");
        }
        result.sigma = final_candidate;
        result.norm_squared = static_cast<double>(final_norm_squared);
        result.gap = static_cast<double>(final_gap);
        result.l2_error_bound = std::sqrt(2.0 * result.gap);
        result.active_vectors = active;
        result.coefficients = coefficients;
        result.final_qp_performed = true;
        result.final_qp_norm_squared = result.norm_squared;
        result.final_qp_gap = result.gap;
        return result;
      }
    }

    expand_active(std::move(minimizing_vertex), std::nullopt,
                  /*prune=*/true);
  }
  return result;
}

void write_result_csv(const std::filesystem::path& path,
                      const PointConfiguration& configuration,
                      const SolverResult& result) {
  std::ofstream stream(path);
  if (!stream) {
    throw std::runtime_error("Cannot open output file: " + path.string());
  }
  const AffineFunction ell = compute_ell(configuration);
  stream << "x,y,sigma,stable_projection,sigma_exact,ell_A,ell_A_exact\n";
  stream << std::setprecision(17);
  for (std::size_t i = 0; i < configuration.size(); ++i) {
    stream << configuration.points()[i].x << ','
           << configuration.points()[i].y << ','
           << result.sigma[static_cast<Eigen::Index>(i)] << ',';
    if (result.stable_projection_available) {
      stream << result.stable_projection[static_cast<Eigen::Index>(i)];
    }
    stream << ',';
    if (result.exact.certified) {
      stream << rational_to_string(result.exact.sigma[i]);
    }
    stream << ',' << ell.values[static_cast<Eigen::Index>(i)] << ','
           << rational_to_string(ell.exact_values[i])
           << '\n';
  }
}

void write_plot_data(const std::filesystem::path& prefix,
                     const PointConfiguration& configuration,
                     const SolverResult& result,
                     const AffineFunction& ell) {
  const RegularTriangulationOracle oracle(configuration);
  Eigen::VectorXd surface_values = result.sigma;
  Eigen::VectorXd sigma_vee;
  std::vector<Rational> exact_sigma_vee;
  GkzVector surface_triangulation;
  if (result.exact.certified) {
    surface_triangulation = oracle.minimize_exact_integer(
        clear_height_denominators(result.exact.sigma),
        /*keep_faces=*/true);
    for (std::size_t i = 0; i < configuration.size(); ++i) {
      surface_values[static_cast<Eigen::Index>(i)] =
          rational_to_double(result.exact.sigma[i]);
    }
  } else {
    surface_triangulation =
        oracle.minimize(result.sigma, /*keep_faces=*/true);
  }
  if (!surface_triangulation.triangulation) {
    throw std::runtime_error("The oracle did not return triangulation faces.");
  }
  if (result.exact.certified) {
    exact_sigma_vee = lower_envelope_values_exact(
        configuration, result.exact.sigma,
        *surface_triangulation.triangulation);
    sigma_vee = Eigen::VectorXd::Zero(
        static_cast<Eigen::Index>(configuration.size()));
    for (std::size_t i = 0; i < configuration.size(); ++i) {
      sigma_vee[static_cast<Eigen::Index>(i)] =
          rational_to_double(exact_sigma_vee[i]);
    }
  } else {
    sigma_vee = lower_envelope_values(
        configuration, surface_values,
        *surface_triangulation.triangulation);
  }

  const std::filesystem::path surface_path = prefix.string() + "_surface.csv";
  const std::filesystem::path triangles_path =
      prefix.string() + "_triangles.csv";
  const std::filesystem::path subdivision_path =
      prefix.string() + "_subdivision.csv";
  const std::filesystem::path ell_path = prefix.string() + "_ell.csv";
  std::ofstream surface(surface_path);
  std::ofstream triangles(triangles_path);
  std::ofstream subdivision(subdivision_path);
  std::ofstream ell_stream(ell_path);
  if (!surface || !triangles || !subdivision || !ell_stream) {
    throw std::runtime_error(
        "Cannot open one of the plot data files for prefix: " +
        prefix.string());
  }

  const long double coordinate_scale =
      configuration.is_polygon_level()
          ? static_cast<long double>(configuration.level())
          : 1.0L;
  const bool has_psi = configuration.is_polygon_level();
  const long double psi_factor =
      has_psi
          ? as_long_double(configuration.base_twice_area()) *
                static_cast<long double>(configuration.level()) *
                static_cast<long double>(configuration.level()) *
                static_cast<long double>(configuration.level())
          : 0.0L;
  const Rational exact_psi_factor =
      has_psi && result.exact.certified
          ? rational_from_wide(configuration.base_twice_area()) *
                Rational(configuration.level()) *
                Rational(configuration.level()) *
                Rational(configuration.level())
          : Rational(0);

  surface << "x,y,sigma,sigma_vee,psi\n" << std::setprecision(17);
  ell_stream << "x,y,ell_A,ell_A_exact\n" << std::setprecision(17);
  for (std::size_t i = 0; i < configuration.size(); ++i) {
    const auto& point = configuration.points()[i];
    const double x = static_cast<double>(point.x / coordinate_scale);
    const double y = static_cast<double>(point.y / coordinate_scale);
    const double sigma = surface_values[static_cast<Eigen::Index>(i)];
    const double envelope = sigma_vee[static_cast<Eigen::Index>(i)];
    surface << x << ',' << y << ',' << sigma << ',' << envelope << ',';
    if (has_psi) {
      if (result.exact.certified) {
        const Rational exact_psi =
            exact_psi_factor * exact_sigma_vee[i] -
            Rational(2) * Rational(configuration.level());
        surface << rational_to_double(exact_psi);
      } else {
        surface << static_cast<double>(psi_factor * envelope -
                                       2.0L * configuration.level());
      }
    }
    surface << '\n';

    ell_stream << x << ',' << y << ','
               << ell.values[static_cast<Eigen::Index>(i)] << ','
               << rational_to_string(ell.exact_values[i]) << '\n';
  }

  triangles << "i,j,l\n";
  for (const auto& face : *surface_triangulation.triangulation) {
    triangles << face[0] << ',' << face[1] << ',' << face[2] << '\n';
  }

  const auto cells = result.exact.certified
                         ? subdivision_cells_exact(
                               configuration, result.exact.sigma,
                               *surface_triangulation.triangulation)
                         : subdivision_cells(
                               configuration, surface_values,
                               *surface_triangulation.triangulation);
  subdivision << "cell,vertex,x,y\n" << std::setprecision(17);
  for (std::size_t cell_index = 0; cell_index < cells.size(); ++cell_index) {
    for (std::size_t vertex_index = 0;
         vertex_index < cells[cell_index].vertices.size(); ++vertex_index) {
      const auto point_index = cells[cell_index].vertices[vertex_index];
      const auto& point = configuration.points()[point_index];
      const double x = static_cast<double>(point.x / coordinate_scale);
      const double y = static_cast<double>(point.y / coordinate_scale);
      subdivision << cell_index << ',' << vertex_index << ',' << x << ',' << y
                  << '\n';
    }
  }
}

void write_stable_projection_plot_data(
    const std::filesystem::path& prefix,
    const PointConfiguration& configuration,
    const SolverResult& result) {
  if (!result.stable_projection_available ||
      result.stable_projection.size() !=
          static_cast<Eigen::Index>(configuration.size())) {
    throw std::invalid_argument(
        "Stable projection plot data requested without a stable projection.");
  }

  const RegularTriangulationOracle oracle(configuration);
  const Eigen::VectorXd& surface_values = result.stable_projection;
  const GkzVector surface_triangulation =
      oracle.minimize(surface_values, /*keep_faces=*/true);
  if (!surface_triangulation.triangulation) {
    throw std::runtime_error("The oracle did not return triangulation faces.");
  }
  const Eigen::VectorXd sigma_vee = lower_envelope_values(
      configuration, surface_values, *surface_triangulation.triangulation);

  const std::filesystem::path surface_path = prefix.string() + "_surface.csv";
  const std::filesystem::path triangles_path =
      prefix.string() + "_triangles.csv";
  const std::filesystem::path subdivision_path =
      prefix.string() + "_subdivision.csv";
  std::ofstream surface(surface_path);
  std::ofstream triangles(triangles_path);
  std::ofstream subdivision(subdivision_path);
  if (!surface || !triangles || !subdivision) {
    throw std::runtime_error(
        "Cannot open one of the stable projection plot files for prefix: " +
        prefix.string());
  }

  const long double coordinate_scale =
      configuration.is_polygon_level()
          ? static_cast<long double>(configuration.level())
          : 1.0L;
  const bool has_psi = configuration.is_polygon_level();
  const long double psi_factor =
      has_psi
          ? as_long_double(configuration.base_twice_area()) *
                static_cast<long double>(configuration.level()) *
                static_cast<long double>(configuration.level()) *
                static_cast<long double>(configuration.level())
          : 0.0L;

  surface << "x,y,sigma,sigma_vee,psi\n" << std::setprecision(17);
  for (std::size_t i = 0; i < configuration.size(); ++i) {
    const auto& point = configuration.points()[i];
    const double x = static_cast<double>(point.x / coordinate_scale);
    const double y = static_cast<double>(point.y / coordinate_scale);
    const double height = surface_values[static_cast<Eigen::Index>(i)];
    const double envelope = sigma_vee[static_cast<Eigen::Index>(i)];
    surface << x << ',' << y << ',' << height << ',' << envelope << ',';
    if (has_psi) {
      surface << static_cast<double>(psi_factor * envelope -
                                     2.0L * configuration.level());
    }
    surface << '\n';
  }

  triangles << "i,j,l\n";
  for (const auto& face : *surface_triangulation.triangulation) {
    triangles << face[0] << ',' << face[1] << ',' << face[2] << '\n';
  }

  const auto cells = subdivision_cells(configuration, surface_values,
                                       *surface_triangulation.triangulation);
  subdivision << "cell,vertex,x,y\n" << std::setprecision(17);
  for (std::size_t cell_index = 0; cell_index < cells.size(); ++cell_index) {
    for (std::size_t vertex_index = 0;
         vertex_index < cells[cell_index].vertices.size(); ++vertex_index) {
      const auto point_index = cells[cell_index].vertices[vertex_index];
      const auto& point = configuration.points()[point_index];
      const double x = static_cast<double>(point.x / coordinate_scale);
      const double y = static_cast<double>(point.y / coordinate_scale);
      subdivision << cell_index << ',' << vertex_index << ',' << x << ',' << y
                  << '\n';
    }
  }
}

}  // namespace gkz
