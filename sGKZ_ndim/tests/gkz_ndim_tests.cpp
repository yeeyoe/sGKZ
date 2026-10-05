#include "gkz_ndim/gkz_ndim.hpp"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void check(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

std::vector<std::string> split_csv(const std::string& line) {
  std::vector<std::string> fields;
  std::stringstream stream(line);
  std::string field;
  while (std::getline(stream, field, ',')) fields.push_back(field);
  return fields;
}

void check_six_points_fixture() {
  std::ifstream input(SGKZ_SIX_POINTS_FIXTURE);
  check(static_cast<bool>(input), "cannot open six_points.points fixture");

  std::vector<sgkz::IntVector> recorded_points;
  std::vector<std::string> recorded_sigma;
  std::string line;
  while (std::getline(input, line)) {
    const auto first = line.find_first_not_of(" \t");
    if (first == std::string::npos || line[first] != '#') continue;
    const auto fields = split_csv(line.substr(first + 1));
    if (fields.size() != 4) continue;
    try {
      recorded_points.push_back({std::stoll(fields[0]), std::stoll(fields[1])});
      recorded_sigma.push_back(fields[3]);
    } catch (const std::invalid_argument&) {
      // Other comments, including prose and headers, are not result rows.
    }
  }

  const auto configuration =
      sgkz::PointConfiguration::from_points_file(SGKZ_SIX_POINTS_FIXTURE);
  check(recorded_points.size() == configuration.size(),
        "six_points.points must record one result row per point");
  for (std::size_t i = 0; i < configuration.size(); ++i)
    check(recorded_points[i] == configuration.points()[i],
          "six_points.points result rows must follow point order");

  const auto result = sgkz::ShortestGkzSolver().solve(configuration);
  check(result.exact.certified, "six_points.points result is not exactly certified");
  for (std::size_t i = 0; i < recorded_sigma.size(); ++i)
    check(sgkz::rational_string(result.exact.sigma[i]) == recorded_sigma[i],
          "six_points.points recorded exact GKZ coordinate changed at row " +
              std::to_string(i));
}
}  // namespace

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

  check_six_points_fixture();

  std::cout << "ok\n";
}
