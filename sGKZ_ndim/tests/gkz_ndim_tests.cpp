#include "gkz_ndim/gkz_ndim.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

int main() {
  using sgkz::IntVector;

  const auto simplex = sgkz::PointConfiguration::from_points(
      {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
  assert(simplex.dimension() == 3 && simplex.size() == 4);
  assert(sgkz::rational_string(simplex.volume()) == "1/6");

  sgkz::RegularTriangulationOracle oracle(simplex);
  const auto initial = oracle.minimize(Eigen::VectorXd::Zero(4));
  assert(std::abs(initial.values.sum() - 1.0) < 1e-12);
  for (int i = 0; i < 4; ++i) assert(std::abs(initial.values[i] - .25) < 1e-12);

  const auto simplex_result = sgkz::ShortestGkzSolver().solve(simplex);
  assert(simplex_result.exact.certified);

  const auto six_points = sgkz::PointConfiguration::from_points(
      {{0, 0}, {0, 2}, {2, 0}, {2, 2}, {1, 1}, {6, 5}});
  assert(six_points.dimension() == 2 && six_points.size() == 6);
  assert(sgkz::rational_string(six_points.volume()) == "11");
  const auto six_result = sgkz::ShortestGkzSolver().solve(six_points);
  assert(six_result.exact.certified);
  assert(sgkz::rational_string(six_result.exact.norm_squared) == "52739/299475");
  const std::vector<std::string> expected_six_sigma = {
      "15/121", "5/33", "493/3025", "1511/9075", "1318/9075", "2267/9075"};
  for (std::size_t i = 0; i < expected_six_sigma.size(); ++i)
    assert(sgkz::rational_string(six_result.exact.sigma[i]) == expected_six_sigma[i]);

  const std::vector<IntVector> cube_vertices = {
      {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
      {1, 1, 0}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
  const auto cube = sgkz::PointConfiguration::from_lattice_polytope(cube_vertices, 1);
  assert(cube.size() == 8 && sgkz::rational_string(cube.volume()) == "1");

  const auto cube2 = sgkz::PointConfiguration::from_lattice_polytope(cube_vertices, 2);
  assert(cube2.size() == 27 && sgkz::rational_string(cube2.volume()) == "8");
  const auto cube2_result = sgkz::ShortestGkzSolver().solve(cube2);
  assert(cube2_result.converged && cube2_result.exact.certified);

  bool limited = false;
  try {
    (void)sgkz::PointConfiguration::from_lattice_polytope(cube_vertices, 2, 10);
  } catch (const std::runtime_error&) {
    limited = true;
  }
  assert(limited);

  std::cout << "ok\n";
}
