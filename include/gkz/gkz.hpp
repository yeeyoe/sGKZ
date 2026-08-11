#pragma once

#include <CGAL/Gmpq.h>
#include <Eigen/Core>

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace gkz {

using WideInt = __int128_t;

struct IntPoint {
  std::int64_t x = 0;
  std::int64_t y = 0;

  friend bool operator==(const IntPoint&, const IntPoint&) = default;
};

std::string to_string(WideInt value);

class PointConfiguration {
 public:
  static PointConfiguration from_points(std::vector<IntPoint> points);
  static PointConfiguration from_points_file(const std::filesystem::path& path);
  static PointConfiguration from_lattice_polygon(
      const std::vector<IntPoint>& polygon_vertices, std::int64_t k);
  static PointConfiguration from_polygon_file(
      const std::filesystem::path& path, std::int64_t k);

  const std::vector<IntPoint>& points() const { return points_; }
  const std::vector<IntPoint>& hull() const { return hull_; }
  WideInt twice_area() const { return twice_area_; }
  // level() is k for A_k = kP cap Z^2, and 0 for a general input set.
  std::int64_t level() const { return level_; }
  bool is_polygon_level() const { return level_ > 0; }
  // The doubled area of the unscaled base polygon P. It is 0 for --points.
  WideInt base_twice_area() const { return base_twice_area_; }
  std::size_t size() const { return points_.size(); }
  double centroid_x() const;
  double centroid_y() const;

 private:
  PointConfiguration(std::vector<IntPoint> points,
                     std::vector<IntPoint> hull, WideInt twice_area,
                     std::int64_t level = 0, WideInt base_twice_area = 0);

  std::vector<IntPoint> points_;
  std::vector<IntPoint> hull_;
  WideInt twice_area_ = 0;
  std::int64_t level_ = 0;
  WideInt base_twice_area_ = 0;
};

struct GkzVector {
  Eigen::VectorXd values;
  std::vector<WideInt> area_numerators;
  std::size_t visible_vertices = 0;
  std::size_t hidden_vertices = 0;
  std::size_t triangles = 0;
  // Point indices of the finite faces of the lower regular triangulation.
  // Shared because GkzVector instances are copied into the active set.
  std::shared_ptr<const std::vector<std::array<std::size_t, 3>>>
      triangulation;

  // GKZ vectors for one point configuration have the same denominator.
  bool has_same_area_numerators(const GkzVector& other) const {
    return area_numerators == other.area_numerators;
  }
};

struct AffineFunction {
  // ell_A(x,y) = coefficients[0] + coefficients[1] x + coefficients[2] y.
  // For A_k, x and y are the normalized coordinates in P, not in kP.
  std::array<double, 3> coefficients{};
  std::array<CGAL::Gmpq, 3> exact_coefficients{};
  Eigen::VectorXd values;
  std::vector<CGAL::Gmpq> exact_values;
};

AffineFunction compute_ell(const PointConfiguration& configuration);

class RegularTriangulationOracle {
 public:
  explicit RegularTriangulationOracle(const PointConfiguration& configuration);
  ~RegularTriangulationOracle();

  // keep_faces controls whether the returned vector retains the full
  // triangulation faces (needed for plot data). The Frank--Wolfe loop
  // passes false to avoid storing ~2n triangles per active vector.
  GkzVector minimize(const Eigen::VectorXd& heights,
                     bool keep_faces = true) const;
  GkzVector minimize_exact(const std::vector<CGAL::Gmpq>& heights,
                           bool keep_faces = true) const;
  GkzVector minimize_exact_integer(
      const std::vector<CGAL::Gmpz>& heights, bool keep_faces = true) const;

 private:
  const PointConfiguration& configuration_;
  struct Cache;
  std::unique_ptr<Cache> cache_;
};

struct SolverOptions {
  double tolerance = 1e-11;
  double absolute_tolerance = 1e-14;
  double correction_tolerance = 1e-14;
  double prune_tolerance = 1e-15;
  int max_iterations = 500;
  int max_correction_steps = 10000;
  int exact_max_active = 128;
  bool exact_certification = true;
  bool verbose = false;
  // Experimental projected-face probes. Disabled by default so the standard
  // fully-corrective Frank--Wolfe trajectory is unchanged.
  bool projection = false;
  int projection_window = 32;
  double projection_stall_ratio = 0.90;
  double projection_relative_gap = 1e-2;
  int projection_probe_period = 16;
  double projection_rank_tolerance = 1e-11;
};

struct ExactCertificate {
  bool certified = false;
  std::string message;
  std::vector<CGAL::Gmpq> sigma;
  std::vector<CGAL::Gmpq> coefficients;
  CGAL::Gmpq norm_squared = 0;
  CGAL::Gmpq support_value = 0;
  // When certification fails because the exact oracle finds a positive
  // Frank--Wolfe gap, witness holds the offending secondary-polytope
  // vertex so the Frank--Wolfe loop can continue with it.
  bool has_witness = false;
  GkzVector witness;
};

// Minimize the norm exactly over conv(active), then use the exact regular
// triangulation oracle to certify global optimality over Sigma(A).
ExactCertificate certify_active_set_exact(
    const PointConfiguration& configuration,
    const std::vector<GkzVector>& active);

struct SolverResult {
  bool converged = false;
  int iterations = 0;
  Eigen::VectorXd sigma;
  double norm_squared = 0.0;
  double gap = 0.0;
  double l2_error_bound = 0.0;
  std::vector<GkzVector> active_vectors;
  Eigen::VectorXd coefficients;
  ExactCertificate exact;
  bool projection_enabled = false;
  int projection_start_iteration = -1;
  int projection_probes = 0;
  int projection_new_vertices = 0;
  std::size_t projection_affine_rank = 0;
};

class ShortestGkzSolver {
 public:
  explicit ShortestGkzSolver(SolverOptions options = {})
      : options_(options) {}

  SolverResult solve(const PointConfiguration& configuration) const;

 private:
  SolverOptions options_;
};

void write_result_csv(const std::filesystem::path& path,
                      const PointConfiguration& configuration,
                      const SolverResult& result);

void write_plot_data(const std::filesystem::path& prefix,
                     const PointConfiguration& configuration,
                     const SolverResult& result,
                     const AffineFunction& ell);

}  // namespace gkz
