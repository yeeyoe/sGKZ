#pragma once

#include <CGAL/Gmpq.h>

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace kstab_highdim {

using Rational = CGAL::Gmpq;
using Vector = std::vector<Rational>;
using IntVector = std::vector<std::int64_t>;

struct AffineFunction {
  Vector coefficients;
  Rational constant = 0;
};

struct ConvexPLFunction {
  std::vector<AffineFunction> branches;
};

struct Simplex {
  std::vector<Vector> vertices;
};

struct Facet {
  std::vector<IntVector> vertices;
  Rational lattice_measure = 0;
};

struct LatticePolytope {
  int dimension = 0;
  std::vector<IntVector> input_points;
  std::vector<IntVector> points;
  std::vector<Simplex> simplices;
  std::vector<Facet> boundary_facets;
  Rational volume = 0;
  Rational boundary_measure = 0;
  Vector first_moment;
  std::vector<std::vector<Rational>> second_moment;
  Vector boundary_first_moment;
  Vector ell;
};

LatticePolytope build_polytope(const std::vector<IntVector>& points);
LatticePolytope parse_polytope_file(const std::filesystem::path& path);
std::vector<AffineFunction> parse_pl_file(const std::filesystem::path& path,
                                          int dimension);

Vector compute_ell_p(const LatticePolytope& polytope);
Rational evaluate_df_exact(const LatticePolytope& polytope,
                           const ConvexPLFunction& function);
Rational evaluate_df_exact(const LatticePolytope& polytope,
                           const AffineFunction& function);

struct SearchOptions {
  int pieces = 1;
  std::size_t population = 0;
  std::size_t generations = 80;
  std::size_t quadrature_samples = 20000;
  std::uint64_t seed = 0;
  unsigned threads = 1;
  bool verbose = false;
};

struct NumericalWitness {
  ConvexPLFunction function;
  std::vector<ConvexPLFunction> alternatives;
  double value = 0.0;
  double normalized = 0.0;
  long evaluations = 0;
};

NumericalWitness search_pl_witness(const LatticePolytope& polytope,
                                   const SearchOptions& options);

struct CertifyResult {
  bool certified = false;
  ConvexPLFunction function;
  Rational value = 0;
};

CertifyResult certify_witness(const LatticePolytope& polytope,
                              const NumericalWitness& witness,
                              std::int64_t max_denominator);

CertifyResult certify_function(const LatticePolytope& polytope,
                               const ConvexPLFunction& function);

Rational approximate_rational(double value, std::int64_t denominator_cap);
double rational_to_double(const Rational& value);
std::string rational_string(const Rational& value);
std::string format_affine(const AffineFunction& function);

}  // namespace kstab_highdim
