#pragma once

#include "chowk/highdim.hpp"
#include "chowk/sgkz.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace chowk {
using Rational = CGAL::Gmpq;
using IntVector = std::vector<std::int64_t>;

struct Options {
  int max_k = 1;
  std::size_t max_points = 1000000;
  double threshold = -1e-6;
  double tolerance = 1e-11;
  double absolute_tolerance = 1e-14;
  int max_iterations = 500;
  bool verbose = false;
};

struct EnvelopeBranch { std::vector<Rational> coefficients; };

struct Envelope {
  std::vector<EnvelopeBranch> branches;
  double max_violation = 0;
  double max_contact_residual = 0;
  bool valid = false;
};

struct LevelResult {
  int k = 0;
  std::vector<IntVector> integer_points;
  std::vector<std::vector<Rational>> points;
  std::vector<double> rho;
  Envelope envelope;
  Rational interior_integral = 0;
  Rational boundary_integral = 0;
  Rational df = 0;
  double normalized_df = 0;
  bool unstable = false;
  std::string error;
  chowk_sgkz::SolverResult solver;
};

struct RunResult {
  std::vector<LevelResult> levels;
  int candidate_k = 0;
  bool unstable = false;
  bool failed = false;
};

std::vector<IntVector> read_vertices(const std::filesystem::path& path);
RunResult run(const std::filesystem::path& path, const Options& options,
             const std::filesystem::path& output_dir);
void write_level(const std::filesystem::path& directory,
                 const kstab_highdim::LatticePolytope& polytope,
                 const LevelResult& result, double threshold);
}  // namespace chowk
