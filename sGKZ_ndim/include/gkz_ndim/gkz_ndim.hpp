#pragma once
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>
#include <Eigen/Core>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace sgkz {
using Rational = CGAL::Gmpq;
using Integer = CGAL::Gmpz;
using IntVector = std::vector<std::int64_t>;
using Wide = __int128_t;

std::string rational_string(const Rational& value);
std::string integer_string(const Integer& value);

class PointConfiguration {
 public:
  static PointConfiguration from_points(std::vector<IntVector> points);
  static PointConfiguration from_points_file(const std::filesystem::path& path);
  static PointConfiguration from_lattice_polytope(std::vector<IntVector> vertices,
                                                   std::int64_t k,
                                                   std::size_t max_points = 1000000);
  static PointConfiguration from_polytope_file(const std::filesystem::path& path,
                                                std::int64_t k,
                                                std::size_t max_points = 1000000);
  int dimension() const { return dimension_; }
  std::size_t size() const { return points_.size(); }
  const std::vector<IntVector>& points() const { return points_; }
  const std::vector<IntVector>& hull_vertices() const { return hull_vertices_; }
  const std::vector<std::vector<std::size_t>>& hull_simplices() const { return hull_simplices_; }
  const Integer& determinant_volume() const { return determinant_volume_; }
  const Integer& base_determinant_volume() const { return base_determinant_volume_; }
  std::int64_t level() const { return level_; }
  bool is_polytope_level() const { return level_ > 0; }
  Rational volume() const;
  Rational base_volume() const;
  std::vector<Rational> centroid() const;
 private:
  PointConfiguration(std::vector<IntVector> points,
                     std::vector<IntVector> hull_vertices,
                     std::vector<std::vector<std::size_t>> hull_simplices,
                     Integer determinant_volume,
                     std::int64_t level,
                     Integer base_determinant_volume);
  int dimension_ = 0;
  std::vector<IntVector> points_;
  std::vector<IntVector> hull_vertices_;
  std::vector<std::vector<std::size_t>> hull_simplices_;
  Integer determinant_volume_ = 0;
  std::int64_t level_ = 0;
  Integer base_determinant_volume_ = 0;
};

struct GkzVector {
  Eigen::VectorXd values;
  std::vector<Integer> numerators;
  std::size_t cells = 0;
  std::size_t visible_vertices = 0;
  std::size_t hidden_vertices = 0;
  std::shared_ptr<const std::vector<std::vector<std::size_t>>> triangulation;
  bool same_numerators(const GkzVector& other) const { return numerators == other.numerators; }
};

struct AffineFunction {
  std::vector<double> coefficients;
  std::vector<Rational> exact_coefficients;
  Eigen::VectorXd values;
  std::vector<Rational> exact_values;
};
AffineFunction compute_ell(const PointConfiguration& configuration);

class RegularTriangulationOracle {
 public:
  explicit RegularTriangulationOracle(const PointConfiguration& configuration);
  GkzVector minimize(const Eigen::VectorXd& heights, bool keep_cells = false) const;
  GkzVector minimize_exact(const std::vector<Rational>& heights, bool keep_cells = false) const;
  GkzVector minimize_exact_integer(const std::vector<Integer>& heights, bool keep_cells = false) const;
 private:
  const PointConfiguration& configuration_;
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
};

struct ExactCertificate {
  bool certified = false;
  std::string message;
  std::vector<Rational> sigma;
  std::vector<Rational> coefficients;
  Rational norm_squared = 0;
  bool has_witness = false;
  GkzVector witness;
};

struct SolverResult {
  bool converged = false;
  int iterations = 0;
  Eigen::VectorXd sigma;
  double norm_squared = 0;
  double gap = 0;
  double l2_error_bound = 0;
  std::vector<GkzVector> active_vectors;
  Eigen::VectorXd coefficients;
  ExactCertificate exact;
};

class ShortestGkzSolver {
 public:
  explicit ShortestGkzSolver(SolverOptions options = {}) : options_(options) {}
  SolverResult solve(const PointConfiguration& configuration) const;
 private:
  SolverOptions options_;
};

void write_result_csv(const std::filesystem::path& path,
                      const PointConfiguration& configuration,
                      const SolverResult& result,
                      bool include_ell);

}  // namespace sgkz
