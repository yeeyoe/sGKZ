#include "chowk/chowk.hpp"

#include <cassert>
#include <fstream>
#include <iostream>
#include <cmath>

int main() {
  using chowk::IntVector;
  const std::vector<IntVector> triangle{{0, 0}, {1, 0}, {0, 1}};
  const auto polytope = kstab_highdim::build_polytope(triangle);
  assert(polytope.dimension == 2);
  assert(polytope.volume == kstab_highdim::Rational(1, 2));
  kstab_highdim::AffineFunction affine;
  affine.constant = 3;
  affine.coefficients = {2, -1};
  assert(kstab_highdim::evaluate_df_exact(polytope, affine) == 0);

  const std::vector<IntVector> six{{0, 0}, {0, 2}, {2, 0}, {2, 2}, {1, 1}, {6, 5}};
  const auto six_configuration = chowk_sgkz::PointConfiguration::from_points(six);
  auto six_result = chowk_sgkz::ShortestGkzSolver().solve(six_configuration);
  const double expected[] = {15.0 / 121.0, 5.0 / 33.0, 493.0 / 3025.0,
                             1511.0 / 9075.0, 1318.0 / 9075.0, 2267.0 / 9075.0};
  assert(six_result.sigma.size() == 6);
  for (int i = 0; i < 6; ++i) assert(std::abs(six_result.sigma[i] - expected[i]) < 1e-10);

  const auto root = std::filesystem::temp_directory_path() / "chowk-tests";
  std::filesystem::remove_all(root);
  const auto input = root / "triangle.points";
  std::filesystem::create_directories(root);
  {
    std::ofstream file(input);
    file << "# triangle\n0,0\n1 0\n0 1\n";
  }
  chowk::Options options;
  options.max_k = 2;
  options.max_points = 100;
  options.max_iterations = 10;
  const auto result = chowk::run(input, options, root / "output");
  assert(!result.levels.empty());
  assert(std::filesystem::exists(root / "output" / "summary.json"));
  assert(std::filesystem::exists(root / "output" / "k_0001" / "points.csv"));
  std::ifstream csv(root / "output" / "k_0001" / "points.csv");
  std::string header;
  std::getline(csv, header);
  assert(header.find("q1,q2,x1,x2,rho,envelope_value,contact_residual") != std::string::npos);
  std::filesystem::remove_all(root);
  std::cout << "chowk tests passed\n";
}
