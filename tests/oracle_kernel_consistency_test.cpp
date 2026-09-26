// Consistency test between the numerical (EPICK) and exact (EPECK)
// regular-triangulation oracles. The two kernels may legitimately return
// different refinements of the same lower hull, so the test does not compare
// triangulations. It compares the support value H(h) = <h, v> as an exact
// rational number, which every minimizing regular triangulation must attain
// regardless of the kernel's tie-breaking.
#include "gkz/gkz.hpp"

#include <Eigen/Core>

#include <cmath>
#include <cstdint>
#include <exception>
#include <iostream>
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

// Deterministic SplitMix64 generator; reproducible across platforms.
class SplitMix64 {
 public:
  explicit SplitMix64(std::uint64_t seed) : state_(seed) {}

  std::uint64_t next_u64() {
    std::uint64_t z = (state_ += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  }

  // Uniform double in (-magnitude, magnitude), excluding denormals.
  double next_double(double magnitude) {
    const double unit = static_cast<double>(next_u64() >> 11) * 0x1p-53;
    return (2.0 * unit - 1.0) * magnitude;
  }

 private:
  std::uint64_t state_;
};

CGAL::Gmpq support_value(const Eigen::VectorXd& heights,
                         const gkz::GkzVector& vertex,
                         const gkz::PointConfiguration& configuration) {
  CGAL::Gmpq value(0);
  const CGAL::Gmpq denominator(
      gkz::to_string(3 * configuration.twice_area()).c_str());
  for (std::size_t i = 0; i < configuration.size(); ++i) {
    // Gmpq(double) converts the binary64 value exactly.
    value += CGAL::Gmpq(heights[static_cast<Eigen::Index>(i)]) *
             CGAL::Gmpq(
                 gkz::to_string(vertex.area_numerators[i]).c_str()) /
             denominator;
  }
  return value;
}

void check_consistency(const gkz::PointConfiguration& configuration,
                       const Eigen::VectorXd& heights,
                       const std::string& label) {
  const gkz::RegularTriangulationOracle oracle(configuration);
  const gkz::GkzVector numeric_vertex = oracle.minimize(heights);

  std::vector<CGAL::Gmpq> exact_heights;
  exact_heights.reserve(configuration.size());
  for (Eigen::Index i = 0; i < heights.size(); ++i) {
    exact_heights.emplace_back(heights[i]);
  }
  const gkz::GkzVector exact_vertex = oracle.minimize_exact(exact_heights);

  const CGAL::Gmpq numeric_support =
      support_value(heights, numeric_vertex, configuration);
  const CGAL::Gmpq exact_support =
      support_value(heights, exact_vertex, configuration);
  std::ostringstream message;
  message << label << ": numerical and exact oracles disagree on H(h). "
          << "numeric=" << numeric_support << " exact=" << exact_support;
  require(numeric_support == exact_support, message.str());
}

std::vector<gkz::PointConfiguration> test_configurations() {
  std::vector<gkz::PointConfiguration> configurations;
  configurations.push_back(gkz::PointConfiguration::from_points(
      {{0, 0}, {1, 0}, {0, 1}, {1, 1}}));
  configurations.push_back(gkz::PointConfiguration::from_points(
      {{0, 0}, {0, 2}, {2, 0}, {2, 2}, {1, 1}, {6, 5}}));
  configurations.push_back(gkz::PointConfiguration::from_points(
      {{0, 0}, {2, 0}, {3, 1}, {2, 3}, {0, 3}, {-1, 1}}));
  configurations.push_back(gkz::PointConfiguration::from_lattice_polygon(
      {{0, 0}, {1, 0}, {1, 1}, {0, 1}}, 4));
  configurations.push_back(gkz::PointConfiguration::from_lattice_polygon(
      {{0, 0}, {4, 1}, {3, 2}, {0, 3}, {-1, 1}}, 3));
  return configurations;
}

void test_zero_and_constant_heights() {
  for (const auto& configuration : test_configurations()) {
    const auto n = static_cast<Eigen::Index>(configuration.size());
    check_consistency(configuration, Eigen::VectorXd::Zero(n),
                      "zero heights");
    check_consistency(configuration, Eigen::VectorXd::Constant(n, 2.5),
                      "constant heights");
  }
}

void test_random_heights() {
  const auto configurations = test_configurations();
  SplitMix64 generator(0x5eedu);
  const double magnitudes[] = {1e-9, 1e-3, 1.0, 1e3, 1e9};
  for (const auto& configuration : configurations) {
    const auto n = static_cast<Eigen::Index>(configuration.size());
    for (const double magnitude : magnitudes) {
      for (int trial = 0; trial < 8; ++trial) {
        Eigen::VectorXd heights(n);
        for (Eigen::Index i = 0; i < n; ++i) {
          heights[i] = generator.next_double(magnitude);
        }
        check_consistency(configuration, heights,
                          "random heights, magnitude=" +
                              std::to_string(magnitude));
      }
    }
  }
}

// Heights taken from the actual Frank--Wolfe iterates exercise the nearly
// degenerate direction: the candidate sigma approaches the normal fan of
// sigma_A, where the numerical kernel is most likely to flip a wall.
void test_frank_wolfe_iterate_heights() {
  for (const auto& configuration : test_configurations()) {
    gkz::SolverOptions options;
    options.max_iterations = 12;
    options.exact_certification = false;
    const gkz::SolverResult result =
        gkz::ShortestGkzSolver(options).solve(configuration);
    check_consistency(configuration, result.sigma,
                      "Frank--Wolfe iterate heights");
  }
}

}  // namespace

int main() {
  try {
    test_zero_and_constant_heights();
    test_random_heights();
    test_frank_wolfe_iterate_heights();
    std::cout << "All oracle kernel consistency tests passed.\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "TEST FAILURE: " << error.what() << '\n';
    return 1;
  }
}
