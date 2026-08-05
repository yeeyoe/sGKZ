#include "k_stability_highdim.hpp"

#include <iostream>
#include <stdexcept>

namespace {
using namespace kstab_highdim;

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

std::vector<IntVector> cube(int d) {
  std::vector<IntVector> points;
  const std::size_t count = std::size_t{1} << d;
  for (std::size_t mask = 0; mask < count; ++mask) {
    IntVector point(d);
    for (int i = 0; i < d; ++i) point[i] = (mask >> i) & 1U;
    points.push_back(std::move(point));
  }
  return points;
}

std::vector<IntVector> standard_simplex(int d) {
  std::vector<IntVector> points(d + 1, IntVector(d, 0));
  for (int i = 0; i < d; ++i) points[i + 1][i] = 1;
  return points;
}

std::vector<IntVector> trapezoid_cube3() {
  const std::vector<IntVector> base = {{0, 0}, {2, 0}, {1, 1}, {0, 1}};
  std::vector<IntVector> points;
  for (const auto& vertex : base) {
    for (int mask = 0; mask < 8; ++mask) {
      IntVector point = {vertex[0], vertex[1], (mask >> 0) & 1,
                         (mask >> 1) & 1, (mask >> 2) & 1};
      points.push_back(std::move(point));
    }
  }
  return points;
}

void test_cube() {
  for (int d : {2, 3, 5}) {
    const auto polytope = build_polytope(cube(d));
    require(polytope.ell.size() == static_cast<std::size_t>(d + 1), "ell dimension");
    require(polytope.ell[0] == 2 * d, "cube ell constant");
    for (int i = 1; i <= d; ++i) require(polytope.ell[i] == 0, "cube ell slope");
    AffineFunction branch;
    branch.coefficients.assign(d, 0);
    branch.coefficients[0] = 1;
    branch.constant = Rational(-1, 2);
    const Rational value = evaluate_df_exact(polytope, ConvexPLFunction{{branch}});
    require(value == Rational(1, 4), "cube hinge value");
  }
}

void test_simplex() {
  for (int d : {2, 5}) {
    const auto polytope = build_polytope(standard_simplex(d));
    require(polytope.ell[0] == d * (d + 1), "simplex ell constant");
    for (int i = 1; i <= d; ++i) require(polytope.ell[i] == 0, "simplex ell slope");
  }
}

void test_pl_invariance() {
  const auto polytope = build_polytope(cube(2));
  AffineFunction first;
  first.coefficients = {Rational(1), Rational(0)};
  first.constant = Rational(-1, 2);
  AffineFunction duplicate = first;
  const auto one = evaluate_df_exact(polytope, ConvexPLFunction{{first}});
  const auto two = evaluate_df_exact(polytope, ConvexPLFunction{{first, duplicate}});
  require(one == two, "duplicate PL branch must not change value");
  AffineFunction second;
  second.coefficients = {Rational(0), Rational(1)};
  second.constant = Rational(-1, 2);
  const auto two_creases = evaluate_df_exact(polytope, ConvexPLFunction{{first, second}});
  require(two_creases == Rational(5, 12), "two-branch PL value");
}

void test_redundant_vertices_and_affine_identity() {
  const auto polytope = build_polytope({{0, 0}, {2, 0}, {2, 2}, {0, 2}, {1, 1}});
  require(polytope.points.size() == 4, "redundant interior point must be excluded from hull");
  require(polytope.volume == 4 && polytope.boundary_measure == 8,
          "redundant point must not change geometry");

  Rational integral_ell = polytope.ell[0] * polytope.volume;
  for (int i = 0; i < polytope.dimension; ++i) {
    integral_ell += polytope.ell[i + 1] * polytope.first_moment[i];
  }
  require(integral_ell == polytope.boundary_measure,
          "ell_P must satisfy the constant affine identity");
  for (int k = 0; k < polytope.dimension; ++k) {
    Rational integral_coordinate = polytope.ell[0] * polytope.first_moment[k];
    for (int i = 0; i < polytope.dimension; ++i) {
      integral_coordinate += polytope.ell[i + 1] * polytope.second_moment[k][i];
    }
    require(integral_coordinate == polytope.boundary_first_moment[k],
            "ell_P must satisfy a coordinate affine identity");
  }

  AffineFunction constant;
  constant.coefficients.assign(polytope.dimension, 0);
  constant.constant = 1;
  require(evaluate_df_exact(polytope, constant) == 0,
          "M_l(1) must vanish");
  for (int k = 0; k < polytope.dimension; ++k) {
    AffineFunction coordinate;
    coordinate.coefficients.assign(polytope.dimension, 0);
    coordinate.coefficients[k] = 1;
    require(evaluate_df_exact(polytope, coordinate) == 0,
            "M_l(x_i) must vanish");
  }
}

void test_product_polytope() {
  const auto polytope = build_polytope(trapezoid_cube3());
  require(polytope.ell[0] == Rational(132, 13), "product ell constant");
  require(polytope.ell[1] == 0 && polytope.ell[2] == Rational(-24, 13),
          "product ell trapezoid slope");
}

void test_rational_approximation() {
  require(approximate_rational(0.5 + 1e-9, 10) == Rational(1, 2), "rational approximation");
}

void test_negative_certificate() {
  auto polytope = build_polytope(cube(2));
  polytope.ell[0] = 8;
  AffineFunction branch;
  branch.coefficients = {Rational(1), Rational(0)};
  branch.constant = Rational(-1, 2);
  const auto result = certify_function(polytope, ConvexPLFunction{{branch}});
  require(result.certified && result.value == Rational(-1, 4), "negative witness certification");
}

void test_search_negative_witness() {
  auto polytope = build_polytope(cube(2));
  polytope.ell[0] = 8;
  SearchOptions options;
  options.pieces = 1;
  options.population = 16;
  options.generations = 8;
  options.quadrature_samples = 1024;
  options.seed = 0;
  const auto witness = search_pl_witness(polytope, options);
  require(witness.normalized < -1e-6, "search should find the mechanical negative witness");
  require(certify_witness(polytope, witness, 1000000).certified,
          "search witness should certify");
}

}  // namespace

int main() {
  try {
    test_cube();
    test_simplex();
    test_pl_invariance();
    test_redundant_vertices_and_affine_identity();
    test_product_polytope();
    test_rational_approximation();
    test_negative_certificate();
    test_search_negative_witness();
    std::cout << "highdim tests passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
