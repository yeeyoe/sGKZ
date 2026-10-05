#include "chowk/chowk.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>

namespace chowk {
namespace {
using Rat = Rational;
using Vec = std::vector<Rat>;

std::string rat_string(const Rat& x) { return kstab_highdim::rational_string(x); }
Rat rat(std::int64_t x) { return Rat(static_cast<long>(x)); }

std::string decimal_string(double value) {
  std::ostringstream out;
  out.setf(std::ios::fixed);
  out << std::setprecision(17) << value;
  return out.str();
}

Rat decimal_rational(const std::string& text) {
  if (text.empty()) throw std::invalid_argument("empty decimal");
  bool negative = false;
  std::size_t begin = 0;
  if (text[begin] == '+' || text[begin] == '-') {
    negative = text[begin] == '-';
    ++begin;
  }
  const auto dot = text.find('.', begin);
  const std::size_t fractional = dot == std::string::npos ? 0 : text.size() - dot - 1;
  std::string digits = dot == std::string::npos ? text.substr(begin)
                                                 : text.substr(begin, dot - begin) + text.substr(dot + 1);
  if (digits.empty()) digits = "0";
  for (char c : digits) if (c < '0' || c > '9') throw std::invalid_argument("invalid decimal: " + text);
  std::string denominator = "1" + std::string(fractional, '0');
  CGAL::Gmpz numerator(digits.c_str());
  if (negative) numerator = -numerator;
  return Rat(numerator, CGAL::Gmpz(denominator.c_str()));
}

std::vector<std::string> tokens(std::string line) {
  if (const auto p = line.find('#'); p != std::string::npos) line.resize(p);
  std::replace(line.begin(), line.end(), ',', ' ');
  std::istringstream in(line);
  std::vector<std::string> out;
  std::string token;
  while (in >> token) out.push_back(token);
  return out;
}

bool solve_linear(std::vector<std::vector<Rat>> a, std::vector<Rat> b,
                  std::vector<Rat>& x) {
  const std::size_t n = a.size();
  if (b.size() != n) return false;
  for (const auto& row : a) if (row.size() != n) return false;
  for (std::size_t c = 0; c < n; ++c) {
    std::size_t p = c;
    while (p < n && a[p][c] == 0) ++p;
    if (p == n) return false;
    if (p != c) { std::swap(a[p], a[c]); std::swap(b[p], b[c]); }
    const Rat q = a[c][c];
    for (std::size_t j = c; j < n; ++j) a[c][j] /= q;
    b[c] /= q;
    for (std::size_t r = 0; r < n; ++r) if (r != c && a[r][c] != 0) {
      const Rat f = a[r][c];
      for (std::size_t j = c; j < n; ++j) a[r][j] -= f * a[c][j];
      b[r] -= f * b[c];
    }
  }
  x = std::move(b);
  return true;
}

Rat affine_value(const std::vector<Rat>& a, const Vec& x) {
  Rat result = a[0];
  for (std::size_t i = 1; i < a.size(); ++i) result += a[i] * x[i - 1];
  return result;
}

Envelope make_envelope(const chowk_sgkz::PointConfiguration& configuration,
                       const std::vector<double>& rho, int k) {
  Envelope empty;
  const int d = configuration.dimension();
  if (rho.size() != configuration.size()) return empty;
  std::vector<Rat> heights;
  heights.reserve(rho.size());
  for (double value : rho) heights.push_back(decimal_rational(decimal_string(value)));
  const auto lower = chowk_sgkz::RegularTriangulationOracle(configuration)
      .minimize_exact(heights, true);
  if (!lower.triangulation) return empty;

  std::set<std::vector<Rat>> seen;
  for (const auto& cell : *lower.triangulation) {
    if (cell.size() != static_cast<std::size_t>(d + 1)) continue;
    std::vector<std::vector<Rat>> matrix(d + 1, std::vector<Rat>(d + 1));
    std::vector<Rat> rhs(d + 1);
    for (int row = 0; row <= d; ++row) {
      const auto index = cell[static_cast<std::size_t>(row)];
      for (int col = 0; col < d; ++col)
        matrix[row][col] = rat(configuration.points()[index][col]);
      matrix[row][d] = 1;
      rhs[row] = heights[index];
    }
    std::vector<Rat> plane;
    if (!solve_linear(matrix, rhs, plane)) continue;
    std::vector<Rat> branch(static_cast<std::size_t>(d + 1));
    branch[0] = plane[d];
    for (int i = 0; i < d; ++i) branch[i + 1] = plane[i] * k;
    if (seen.insert(branch).second) empty.branches.push_back({std::move(branch)});
  }
  if (empty.branches.empty()) return empty;
  for (std::size_t i = 0; i < configuration.size(); ++i) {
    Vec x;
    for (auto coordinate : configuration.points()[i]) x.push_back(rat(coordinate) / k);
    Rat value = affine_value(empty.branches.front().coefficients, x);
    for (const auto& branch : empty.branches)
      value = std::max(value, affine_value(branch.coefficients, x));
    const Rat difference = heights[i] - value;
    empty.max_violation = std::max(empty.max_violation,
                                   std::max(0.0, -CGAL::to_double(difference)));
    empty.max_contact_residual = std::max(empty.max_contact_residual,
                                          std::abs(CGAL::to_double(difference)));
  }
  empty.valid = empty.max_violation <= 1e-10;
  return empty;
}

kstab_highdim::ConvexPLFunction as_pl(const Envelope& envelope) {
  kstab_highdim::ConvexPLFunction function;
  for (const auto& branch : envelope.branches) {
    kstab_highdim::AffineFunction affine;
    affine.constant = branch.coefficients[0];
    affine.coefficients.assign(branch.coefficients.begin() + 1,
                               branch.coefficients.end());
    function.branches.push_back(std::move(affine));
  }
  return function;
}

double polytope_diameter(const kstab_highdim::LatticePolytope& polytope) {
  double diameter = 0.0;
  for (const auto& p : polytope.input_points) for (const auto& q : polytope.input_points) {
    double squared = 0.0;
    for (int i = 0; i < polytope.dimension; ++i) {
      const double delta = static_cast<double>(p[i] - q[i]);
      squared += delta * delta;
    }
    diameter = std::max(diameter, std::sqrt(squared));
  }
  return diameter;
}

double normalization_scale(const kstab_highdim::LatticePolytope& polytope) {
  const double diameter = polytope_diameter(polytope);
  double sup_ell = 0.0;
  for (const auto& p : polytope.input_points) {
    double value = CGAL::to_double(polytope.ell[0]);
    for (int i = 0; i < polytope.dimension; ++i)
      value += CGAL::to_double(polytope.ell[i + 1]) * static_cast<double>(p[i]);
    sup_ell = std::max(sup_ell, std::abs(value));
  }
  return std::max(1.0, diameter * (CGAL::to_double(polytope.boundary_measure) +
                                  CGAL::to_double(polytope.volume) * sup_ell));
}

void write_json_string(std::ostream& out, const std::string& value) {
  out << '"';
  for (char c : value) { if (c == '"' || c == '\\') out << '\\'; out << c; }
  out << '"';
}
}  // namespace

std::vector<IntVector> read_vertices(const std::filesystem::path& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open polytope file: " + path.string());
  std::vector<IntVector> result;
  std::string line;
  int dimension = 0;
  while (std::getline(input, line)) {
    const auto fields = tokens(line);
    if (fields.empty()) continue;
    if (dimension == 0) dimension = static_cast<int>(fields.size());
    if (static_cast<int>(fields.size()) != dimension)
      throw std::runtime_error("inconsistent coordinate count");
    IntVector point;
    for (const auto& field : fields) {
      std::size_t consumed = 0;
      const auto value = std::stoll(field, &consumed);
      if (consumed != field.size()) throw std::runtime_error("non-integer coordinate");
      point.push_back(static_cast<std::int64_t>(value));
    }
    result.push_back(std::move(point));
  }
  if (result.empty()) throw std::runtime_error("polytope has no vertices");
  return result;
}

RunResult run(const std::filesystem::path& path, const Options& options,
              const std::filesystem::path& output) {
  if (options.max_k < 1) throw std::invalid_argument("--max-k must be positive");
  std::filesystem::create_directories(output);
  const auto vertices = read_vertices(path);
  const auto polytope = kstab_highdim::build_polytope(vertices);
  RunResult result;
  for (int k = 1; k <= options.max_k; ++k) {
    LevelResult level;
    level.k = k;
    const auto configuration = chowk_sgkz::PointConfiguration::from_lattice_polytope(
        vertices, k, options.max_points);
    level.integer_points = configuration.points();
    for (const auto& point : level.integer_points) {
      Vec x;
      for (auto coordinate : point) x.push_back(rat(coordinate) / k);
      level.points.push_back(std::move(x));
    }
    chowk_sgkz::SolverOptions solver_options;
    solver_options.tolerance = options.tolerance;
    solver_options.absolute_tolerance = options.absolute_tolerance;
    solver_options.max_iterations = options.max_iterations;
    solver_options.exact_certification = false;
    solver_options.verbose = options.verbose;
    level.solver = chowk_sgkz::ShortestGkzSolver(solver_options).solve(configuration);
    level.rho.resize(static_cast<std::size_t>(level.solver.sigma.size()));
    for (Eigen::Index i = 0; i < level.solver.sigma.size(); ++i)
      level.rho[static_cast<std::size_t>(i)] = level.solver.sigma[i];
    try {
      level.envelope = make_envelope(configuration, level.rho, k);
    } catch (const std::exception& error) {
      level.error = std::string("lower envelope construction failed: ") + error.what();
    }
    if (!level.envelope.valid) level.error = "lower envelope validation failed";
    if (level.error.empty()) {
      const auto integrals = kstab_highdim::evaluate_pl_integrals_exact(polytope, as_pl(level.envelope));
      level.boundary_integral = integrals.boundary;
      level.interior_integral = integrals.interior;
      level.df = level.boundary_integral - level.interior_integral;
      level.normalized_df = CGAL::to_double(level.df) / normalization_scale(polytope);
      level.unstable = level.normalized_df <= options.threshold;
    }
    std::ostringstream name;
    name << "k_" << std::setw(4) << std::setfill('0') << k;
    write_level(output / name.str(), polytope, level, options.threshold);
    result.levels.push_back(level);
    if (!level.error.empty()) result.failed = true;
    if (level.unstable) {
      result.unstable = true;
      result.candidate_k = k;
      write_level(output / "candidate", polytope, level, options.threshold);
      break;
    }
  }
  std::ofstream summary(output / "summary.json");
  summary << "{\n  \"input\": "; write_json_string(summary, path.string());
  summary << ",\n  \"dimension\": " << polytope.dimension
          << ",\n  \"max_k\": " << options.max_k
          << ",\n  \"volume\": "; write_json_string(summary, rat_string(polytope.volume));
  summary << ",\n  \"boundary_measure\": "; write_json_string(summary, rat_string(polytope.boundary_measure));
  summary << ",\n  \"ell_P\": [";
  for (std::size_t i = 0; i < polytope.ell.size(); ++i) { if (i) summary << ','; write_json_string(summary, rat_string(polytope.ell[i])); }
  summary << "]"
          << ",\n  \"diameter\": " << polytope_diameter(polytope)
          << ",\n  \"candidate_k\": " << result.candidate_k
          << ",\n  \"levels_written\": " << result.levels.size()
          << ",\n  \"output_dir\": "; write_json_string(summary, output.string());
  summary << ",\n  \"candidate_dir\": ";
  write_json_string(summary, result.unstable ? (output / "candidate").string() : "");
  summary << ",\n  \"relative_K_status\": ";
  write_json_string(summary, result.unstable ? "unstable" : (result.failed ? "error" : "no_counterexample_found"));
  summary << "\n}\n";
  return result;
}

void write_level(const std::filesystem::path& directory,
                 const kstab_highdim::LatticePolytope& polytope,
                 const LevelResult& result, double threshold) {
  std::filesystem::create_directories(directory);
  std::ofstream csv(directory / "points.csv");
  for (int i = 0; i < polytope.dimension; ++i) { if (i) csv << ','; csv << 'q' << i + 1; }
  for (int i = 0; i < polytope.dimension; ++i) csv << ",x" << i + 1;
  csv << ",rho,envelope_value,contact_residual\n";
  for (std::size_t i = 0; i < result.integer_points.size(); ++i) {
    for (auto coordinate : result.integer_points[i]) csv << coordinate << ',';
    for (std::size_t j = 0; j < result.points[i].size(); ++j) {
      if (j) csv << ','; csv << rat_string(result.points[i][j]);
    }
    Vec envelope_point = result.points[i];
    Rat value = 0;
    if (!result.envelope.branches.empty()) {
      value = affine_value(result.envelope.branches.front().coefficients, envelope_point);
      for (const auto& branch : result.envelope.branches)
        value = std::max(value, affine_value(branch.coefficients, envelope_point));
    }
    const double rho = i < result.rho.size() ? result.rho[i] : 0.0;
    const auto rho_text = decimal_string(rho);
    csv << ',' << rho_text << ',' << rat_string(value)
        << ',' << rat_string(decimal_rational(rho_text) - value) << '\n';
  }
  std::ofstream pl(directory / "envelope.pl");
  for (const auto& branch : result.envelope.branches) {
    for (std::size_t i = 0; i < branch.coefficients.size(); ++i) {
      if (i) pl << ',';
      pl << rat_string(branch.coefficients[i]);
    }
    pl << '\n';
  }
  std::ofstream json(directory / "summary.json");
  json << "{\n  \"k\": " << result.k
       << ",\n  \"dimension\": " << polytope.dimension
       << ",\n  \"volume\": "; write_json_string(json, rat_string(polytope.volume));
  json << ",\n  \"boundary_measure\": "; write_json_string(json, rat_string(polytope.boundary_measure));
  json << ",\n  \"diameter\": " << polytope_diameter(polytope)
       << ",\n  \"ell_P\": [";
  for (std::size_t i = 0; i < polytope.ell.size(); ++i) { if (i) json << ','; write_json_string(json, rat_string(polytope.ell[i])); }
  json << "]"
       << ",\n  \"point_count\": " << result.integer_points.size()
       << ",\n  \"solver_converged\": " << (result.solver.converged ? "true" : "false")
       << ",\n  \"iterations\": " << result.solver.iterations
       << ",\n  \"active_size\": " << result.solver.active_vectors.size()
       << ",\n  \"norm_squared\": " << std::setprecision(17) << result.solver.norm_squared
       << ",\n  \"gap\": " << result.solver.gap
       << ",\n  \"l2_error_bound\": " << result.solver.l2_error_bound
       << ",\n  \"branch_count\": " << result.envelope.branches.size()
       << ",\n  \"boundary_integral\": "; write_json_string(json, rat_string(result.boundary_integral));
  json << ",\n  \"interior_integral\": "; write_json_string(json, rat_string(result.interior_integral));
  json << ",\n  \"M_rel\": "; write_json_string(json, rat_string(result.df));
  json << ",\n  \"normalized_df\": " << result.normalized_df
       << ",\n  \"threshold\": " << threshold
       << ",\n  \"relative_K_status\": ";
  write_json_string(json, !result.error.empty() ? "error" : (result.unstable ? "unstable" : "no_counterexample_found"));
  json << ",\n  \"max_envelope_violation\": " << result.envelope.max_violation
       << ",\n  \"max_contact_residual\": " << result.envelope.max_contact_residual
       << ",\n  \"unstable\": " << (result.unstable ? "true" : "false");
  json << ",\n  \"points_file\": "; write_json_string(json, (directory / "points.csv").string());
  json << ",\n  \"envelope_file\": "; write_json_string(json, (directory / "envelope.pl").string());
  if (!result.error.empty()) { json << ",\n  \"error\": "; write_json_string(json, result.error); }
  json << "\n}\n";
}
}  // namespace chowk
