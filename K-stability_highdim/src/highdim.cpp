#include "k_stability_highdim.hpp"
#include "delaunay_adapter.hpp"

#include <CGAL/Gmpz.h>

#include <boost/math/optimization/differential_evolution.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace kstab_highdim {
namespace {

Rational rat(std::int64_t value) { return Rational(static_cast<long>(value)); }

Rational abs_rat(const Rational& value) { return value < 0 ? -value : value; }

Rational factorial(int n) {
  Rational result = 1;
  for (int i = 2; i <= n; ++i) result *= i;
  return result;
}

bool equal_vector(const Vector& a, const Vector& b) {
  return a == b;
}

Rational determinant(std::vector<std::vector<Rational>> matrix) {
  const int n = static_cast<int>(matrix.size());
  if (n == 0) return 1;
  for (const auto& row : matrix) {
    if (static_cast<int>(row.size()) != n) {
      throw std::invalid_argument("determinant requires a square matrix");
    }
  }
  Rational result = 1;
  int sign = 1;
  for (int column = 0; column < n; ++column) {
    int pivot = column;
    while (pivot < n && matrix[pivot][column] == 0) ++pivot;
    if (pivot == n) return 0;
    if (pivot != column) {
      std::swap(matrix[pivot], matrix[column]);
      sign = -sign;
    }
    const Rational pivot_value = matrix[column][column];
    result *= pivot_value;
    for (int row = column + 1; row < n; ++row) {
      const Rational factor = matrix[row][column] / pivot_value;
      for (int j = column + 1; j < n; ++j) {
        matrix[row][j] -= factor * matrix[column][j];
      }
    }
  }
  return sign < 0 ? -result : result;
}

Vector solve_linear(std::vector<std::vector<Rational>> matrix,
                    Vector rhs) {
  const int n = static_cast<int>(matrix.size());
  if (static_cast<int>(rhs.size()) != n) {
    throw std::invalid_argument("linear system size mismatch");
  }
  for (const auto& row : matrix) {
    if (static_cast<int>(row.size()) != n) {
      throw std::invalid_argument("linear system must be square");
    }
  }
  for (int column = 0; column < n; ++column) {
    int pivot = column;
    while (pivot < n && matrix[pivot][column] == 0) ++pivot;
    if (pivot == n) throw std::runtime_error("singular exact linear system");
    if (pivot != column) {
      std::swap(matrix[pivot], matrix[column]);
      std::swap(rhs[pivot], rhs[column]);
    }
    const Rational pivot_value = matrix[column][column];
    for (int j = column; j < n; ++j) matrix[column][j] /= pivot_value;
    rhs[column] /= pivot_value;
    for (int row = 0; row < n; ++row) {
      if (row == column) continue;
      const Rational factor = matrix[row][column];
      if (factor == 0) continue;
      for (int j = column; j < n; ++j) {
        matrix[row][j] -= factor * matrix[column][j];
      }
      rhs[row] -= factor * rhs[column];
    }
  }
  return rhs;
}

int matrix_rank(std::vector<std::vector<Rational>> matrix) {
  if (matrix.empty()) return 0;
  const int rows = static_cast<int>(matrix.size());
  const int columns = static_cast<int>(matrix.front().size());
  int rank = 0;
  for (int column = 0; column < columns && rank < rows; ++column) {
    int pivot = rank;
    while (pivot < rows && matrix[pivot][column] == 0) ++pivot;
    if (pivot == rows) continue;
    std::swap(matrix[rank], matrix[pivot]);
    const Rational value = matrix[rank][column];
    for (int j = column; j < columns; ++j) matrix[rank][j] /= value;
    for (int row = rank + 1; row < rows; ++row) {
      const Rational factor = matrix[row][column];
      for (int j = column; j < columns; ++j) {
        matrix[row][j] -= factor * matrix[rank][j];
      }
    }
    ++rank;
  }
  return rank;
}

struct SimplexMoments {
  Rational volume = 0;
  Vector first;
  std::vector<std::vector<Rational>> second;
};

SimplexMoments simplex_moments(const std::vector<Vector>& vertices) {
  const int d = static_cast<int>(vertices.front().size());
  if (static_cast<int>(vertices.size()) != d + 1) {
    throw std::invalid_argument("simplex has wrong number of vertices");
  }
  std::vector<std::vector<Rational>> edges(d, Vector(d));
  for (int row = 0; row < d; ++row) {
    for (int column = 0; column < d; ++column) {
      edges[row][column] = vertices[row + 1][column] - vertices[0][column];
    }
  }
  const Rational det = abs_rat(determinant(edges));
  SimplexMoments result;
  result.volume = det / factorial(d);
  result.first.assign(d, 0);
  result.second.assign(d, Vector(d, 0));
  Vector sums(d, 0);
  for (const auto& vertex : vertices) {
    for (int i = 0; i < d; ++i) sums[i] += vertex[i];
  }
  for (int i = 0; i < d; ++i) {
    result.first[i] = result.volume * sums[i] / (d + 1);
    for (int j = 0; j < d; ++j) {
      Rational diagonal = 0;
      for (const auto& vertex : vertices) diagonal += vertex[i] * vertex[j];
      result.second[i][j] = result.volume * (sums[i] * sums[j] + diagonal) /
                            ((d + 1) * (d + 2));
    }
  }
  return result;
}

Rational product_integral(const AffineFunction& left,
                          const AffineFunction& right,
                          const SimplexMoments& moments) {
  Rational value = left.constant * right.constant * moments.volume;
  for (std::size_t i = 0; i < left.coefficients.size(); ++i) {
    value += (left.coefficients[i] * right.constant +
              right.coefficients[i] * left.constant) * moments.first[i];
    for (std::size_t j = 0; j < left.coefficients.size(); ++j) {
      value += left.coefficients[i] * right.coefficients[j] *
               moments.second[i][j];
    }
  }
  return value;
}

std::int64_t integer_coordinate(const Rational& value) {
  if (value.denominator() != 1) {
    throw std::runtime_error("expected an integer boundary coordinate");
  }
  const CGAL::Gmpz numerator = value.numerator();
  if (numerator < CGAL::Gmpz(static_cast<long>(std::numeric_limits<std::int64_t>::min())) ||
      numerator > CGAL::Gmpz(static_cast<long>(std::numeric_limits<std::int64_t>::max()))) {
    throw std::runtime_error("integer coordinate exceeds int64 range");
  }
  return static_cast<std::int64_t>(CGAL::to_double(numerator));
}

std::int64_t gcd_abs(std::int64_t a, std::int64_t b) {
  if (a < 0) a = -a;
  if (b < 0) b = -b;
  return std::gcd(a, b);
}

std::int64_t gcd_vector(const std::vector<std::int64_t>& values) {
  std::int64_t result = 0;
  for (const auto value : values) result = gcd_abs(result, value);
  return result;
}

Rational facet_measure(const std::vector<IntVector>& vertices) {
  const int d = static_cast<int>(vertices.front().size());
  const int k = d - 1;
  std::vector<std::vector<std::int64_t>> edge(k, std::vector<std::int64_t>(d));
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < d; ++j) edge[i][j] = vertices[i + 1][j] - vertices[0][j];
  }
  std::vector<std::int64_t> cofactors(d, 0);
  for (int omitted = 0; omitted < d; ++omitted) {
    std::vector<std::vector<Rational>> minor(k, Vector(k));
    for (int i = 0; i < k; ++i) {
      int column = 0;
      for (int j = 0; j < d; ++j) {
        if (j == omitted) continue;
        minor[i][column++] = rat(edge[i][j]);
      }
    }
    Rational value = determinant(minor);
    if (omitted % 2 != 0) value = -value;
    cofactors[omitted] = integer_coordinate(value);
  }
  const std::int64_t primitive_scale = gcd_vector(cofactors);
  if (primitive_scale == 0) throw std::runtime_error("degenerate boundary facet");
  return rat(primitive_scale) / factorial(k);
}

void add_simplex_moments(const SimplexMoments& source, Rational& volume,
                         Vector& first,
                         std::vector<std::vector<Rational>>& second) {
  volume += source.volume;
  for (std::size_t i = 0; i < first.size(); ++i) {
    first[i] += source.first[i];
    for (std::size_t j = 0; j < first.size(); ++j) {
      second[i][j] += source.second[i][j];
    }
  }
}

std::vector<IntVector> extract_integer_vertices(const std::vector<Vector>& vertices) {
  std::vector<IntVector> result;
  result.reserve(vertices.size());
  for (const auto& vertex : vertices) {
    IntVector integer;
    integer.reserve(vertex.size());
    for (const auto& coordinate : vertex) integer.push_back(integer_coordinate(coordinate));
    result.push_back(std::move(integer));
  }
  return result;
}

}  // namespace

LatticePolytope build_polytope(const std::vector<IntVector>& input) {
  if (input.empty()) throw std::invalid_argument("polytope has no points");
  const int d = static_cast<int>(input.front().size());
  if (d < 2) throw std::invalid_argument("polytope dimension must be at least 2");
  std::set<IntVector> unique;
  for (const auto& point : input) {
    if (static_cast<int>(point.size()) != d) {
      throw std::invalid_argument("polytope rows have different dimensions");
    }
    if (!unique.insert(point).second) throw std::invalid_argument("duplicate polytope point");
  }
  std::vector<Vector> rational_points;
  rational_points.reserve(input.size());
  for (const auto& point : input) {
    Vector rational;
    rational.reserve(point.size());
    for (const auto coordinate : point) rational.push_back(rat(coordinate));
    rational_points.push_back(std::move(rational));
  }
  std::vector<std::vector<Rational>> differences;
  for (std::size_t i = 1; i < rational_points.size(); ++i) {
    Vector difference(d);
    for (int j = 0; j < d; ++j) difference[j] = rational_points[i][j] - rational_points[0][j];
    differences.push_back(std::move(difference));
  }
  if (matrix_rank(differences) != d) {
    throw std::invalid_argument("polytope points are not full dimensional");
  }

  LatticePolytope polytope;
  polytope.dimension = d;
  polytope.input_points = input;
  polytope.first_moment.assign(d, 0);
  polytope.second_moment.assign(d, Vector(d, 0));
  polytope.boundary_first_moment.assign(d, 0);
  const auto exact = detail::triangulate_exact(rational_points, d);
  polytope.simplices = exact.cells;
  polytope.points = extract_integer_vertices(exact.hull_points);
  if (polytope.points.empty()) polytope.points = input;
  if (polytope.simplices.empty()) throw std::runtime_error("cannot triangulate polytope");

  for (const auto& simplex : polytope.simplices) {
    add_simplex_moments(simplex_moments(simplex.vertices), polytope.volume,
                        polytope.first_moment, polytope.second_moment);
  }

  for (const auto& facet_points : exact.boundary_facets) {
      Facet facet;
      facet.vertices = extract_integer_vertices(facet_points);
      facet.lattice_measure = facet_measure(facet.vertices);
      polytope.boundary_measure += facet.lattice_measure;
      for (const auto& point : facet.vertices) {
        for (int j = 0; j < d; ++j) polytope.boundary_first_moment[j] += facet.lattice_measure * rat(point[j]) / d;
      }
      polytope.boundary_facets.push_back(std::move(facet));
  }
  if (polytope.boundary_facets.empty()) throw std::runtime_error("cannot extract boundary facets");
  polytope.ell = compute_ell_p(polytope);
  return polytope;
}

Vector compute_ell_p(const LatticePolytope& polytope) {
  const int d = polytope.dimension;
  if (d < 2 || static_cast<int>(polytope.first_moment.size()) != d ||
      static_cast<int>(polytope.boundary_first_moment.size()) != d ||
      static_cast<int>(polytope.second_moment.size()) != d ||
      static_cast<int>(polytope.second_moment.front().size()) != d) {
    throw std::invalid_argument("invalid polytope moments for ell_P");
  }
  for (const auto& row : polytope.second_moment) {
    if (static_cast<int>(row.size()) != d) {
      throw std::invalid_argument("invalid second moment dimensions for ell_P");
    }
  }
  std::vector<std::vector<Rational>> matrix(d + 1, Vector(d + 1, 0));
  matrix[0][0] = polytope.volume;
  for (int i = 0; i < d; ++i) {
    matrix[0][i + 1] = polytope.first_moment[i];
    matrix[i + 1][0] = polytope.first_moment[i];
    for (int j = 0; j < d; ++j) matrix[i + 1][j + 1] = polytope.second_moment[i][j];
  }
  Vector rhs(d + 1, 0);
  rhs[0] = polytope.boundary_measure;
  for (int i = 0; i < d; ++i) rhs[i + 1] = polytope.boundary_first_moment[i];
  const Vector ell = solve_linear(matrix, rhs);

  // Recheck the defining identity on the affine basis exactly. This catches
  // malformed externally constructed moment data instead of returning a
  // plausible-looking but invalid affine function.
  Rational interior_constant = ell[0] * polytope.volume;
  for (int i = 0; i < d; ++i) interior_constant += ell[i + 1] * polytope.first_moment[i];
  if (polytope.boundary_measure - interior_constant != 0) {
    throw std::runtime_error("ell_P failed the constant affine identity");
  }
  for (int k = 0; k < d; ++k) {
    Rational interior_coordinate = ell[0] * polytope.first_moment[k];
    for (int i = 0; i < d; ++i) interior_coordinate += ell[i + 1] * polytope.second_moment[k][i];
    if (polytope.boundary_first_moment[k] - interior_coordinate != 0) {
      throw std::runtime_error("ell_P failed an affine-coordinate identity");
    }
  }
  return ell;
}

std::string rational_string(const Rational& value) {
  std::ostringstream stream;
  if (value.denominator() == 1) {
    stream << value.numerator();
  } else {
    stream << value;
  }
  return stream.str();
}

double rational_to_double(const Rational& value) { return CGAL::to_double(value); }

Rational parse_rational_token(const std::string& token) {
  const auto slash = token.find('/');
  if (slash == std::string::npos) {
    std::size_t consumed = 0;
    const long long value = std::stoll(token, &consumed);
    if (consumed != token.size()) throw std::invalid_argument("expected integer or p/q rational");
    return Rational(static_cast<long>(value));
  }
  if (slash == 0 || slash + 1 == token.size() || token.find('/', slash + 1) != std::string::npos) {
    throw std::invalid_argument("malformed rational: " + token);
  }
  const long long numerator = std::stoll(token.substr(0, slash));
  const long long denominator = std::stoll(token.substr(slash + 1));
  if (denominator == 0) throw std::invalid_argument("zero rational denominator");
  return Rational(CGAL::Gmpz(static_cast<long>(numerator)),
                  CGAL::Gmpz(static_cast<long>(denominator)));
}

std::vector<std::string> tokens_without_comment(std::string line) {
  if (const auto comment = line.find('#'); comment != std::string::npos) line.erase(comment);
  std::replace(line.begin(), line.end(), ',', ' ');
  std::istringstream stream(line);
  std::vector<std::string> tokens;
  std::string token;
  while (stream >> token) tokens.push_back(token);
  return tokens;
}

std::vector<IntVector> parse_integer_rows(const std::filesystem::path& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open polytope file: " + path.string());
  std::vector<IntVector> rows;
  std::string line;
  std::size_t line_number = 0;
  int dimension = 0;
  while (std::getline(input, line)) {
    ++line_number;
    const auto tokens = tokens_without_comment(line);
    if (tokens.empty()) continue;
    if (dimension == 0) dimension = static_cast<int>(tokens.size());
    if (static_cast<int>(tokens.size()) != dimension) {
      throw std::runtime_error("inconsistent coordinate count at line " + std::to_string(line_number));
    }
    IntVector row;
    for (const auto& token : tokens) {
      std::size_t consumed = 0;
      const long long value = std::stoll(token, &consumed);
      if (consumed != token.size()) throw std::runtime_error("non-integer coordinate at line " + std::to_string(line_number));
      row.push_back(static_cast<std::int64_t>(value));
    }
    rows.push_back(std::move(row));
  }
  return rows;
}

LatticePolytope parse_polytope_file(const std::filesystem::path& path) {
  return build_polytope(parse_integer_rows(path));
}

std::vector<AffineFunction> parse_pl_file(const std::filesystem::path& path, int dimension) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open PL file: " + path.string());
  std::vector<AffineFunction> branches;
  std::string line;
  std::size_t line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    const auto tokens = tokens_without_comment(line);
    if (tokens.empty()) continue;
    if (static_cast<int>(tokens.size()) != dimension + 1) {
      throw std::runtime_error("PL line " + std::to_string(line_number) + " must have d+1 coefficients");
    }
    AffineFunction branch;
    branch.coefficients.reserve(dimension);
    for (int i = 0; i < dimension; ++i) branch.coefficients.push_back(parse_rational_token(tokens[i]));
    branch.constant = parse_rational_token(tokens.back());
    branches.push_back(std::move(branch));
  }
  if (branches.empty() || branches.size() > 8) throw std::invalid_argument("PL branch count must be in [1,8]");
  return branches;
}

AffineFunction ell_function(const LatticePolytope& polytope) {
  AffineFunction result;
  result.constant = polytope.ell[0];
  result.coefficients.assign(polytope.ell.begin() + 1, polytope.ell.end());
  return result;
}

Rational evaluate_affine(const AffineFunction& function, const Vector& point) {
  Rational value = function.constant;
  if (function.coefficients.size() != point.size()) throw std::invalid_argument("affine dimension mismatch");
  for (std::size_t i = 0; i < point.size(); ++i) value += function.coefficients[i] * point[i];
  return value;
}

struct Inequality {
  Vector coefficients;
  Rational constant = 0;
};

bool satisfies(const Inequality& inequality, const Vector& point) {
  Rational value = inequality.constant;
  for (std::size_t i = 0; i < point.size(); ++i) value += inequality.coefficients[i] * point[i];
  return value >= 0;
}

std::vector<Vector> region_vertices(const std::vector<Vector>& simplex,
                                    const std::vector<Inequality>& inequalities,
                                    int k) {
  if (static_cast<int>(simplex.size()) != k + 1) throw std::invalid_argument("region source is not a simplex");
  const int ambient_dimension = static_cast<int>(simplex.front().size());
  if (static_cast<int>(simplex.size()) != k + 1) {
    throw std::invalid_argument("region source has wrong intrinsic dimension");
  }
  std::vector<std::vector<Rational>> basis(k, Vector(ambient_dimension));
  for (int row = 0; row < k; ++row) {
    for (int column = 0; column < ambient_dimension; ++column) basis[row][column] = simplex[row + 1][column] - simplex[0][column];
  }
  std::vector<Inequality> local;
  local.reserve(inequalities.size());
  for (std::size_t inequality_index = 0; inequality_index < inequalities.size(); ++inequality_index) {
    const auto& inequality = inequalities[inequality_index];
    Inequality transformed;
    if (inequality_index < static_cast<std::size_t>(k + 1)) {
      transformed = inequality;
    } else {
      transformed.coefficients.assign(k, 0);
      transformed.constant = inequality.constant;
      for (int i = 0; i < ambient_dimension; ++i) transformed.constant += inequality.coefficients[i] * simplex[0][i];
      for (int j = 0; j < k; ++j) {
        for (int i = 0; i < ambient_dimension; ++i) transformed.coefficients[j] += inequality.coefficients[i] * basis[j][i];
      }
    }
    local.push_back(std::move(transformed));
  }
  std::vector<Vector> result;
  if (local.size() < static_cast<std::size_t>(k)) return result;
  std::vector<int> choice(k);
  std::function<void(int, int)> enumerate = [&](int depth, int next) {
    if (depth == k) {
      std::vector<std::vector<Rational>> matrix(k, Vector(k));
      Vector rhs(k, 0);
      for (int row = 0; row < k; ++row) {
        matrix[row] = local[choice[row]].coefficients;
        rhs[row] = -local[choice[row]].constant;
      }
      Rational det = determinant(matrix);
      if (det == 0) return;
      Vector point = solve_linear(std::move(matrix), std::move(rhs));
      for (const auto& inequality : local) {
        if (!satisfies(inequality, point)) return;
      }
      Vector global(k == 0 ? 0 : simplex[0].size());
      for (std::size_t i = 0; i < global.size(); ++i) {
        global[i] = simplex[0][i];
        for (int j = 0; j < k; ++j) global[i] += basis[j][i] * point[j];
      }
      for (const auto& existing : result) if (equal_vector(existing, global)) return;
      result.push_back(std::move(global));
      return;
    }
    for (int i = next; i <= static_cast<int>(local.size()) - (k - depth); ++i) {
      choice[depth] = i;
      enumerate(depth + 1, i + 1);
    }
  };
  if (k == 0) return {simplex.front()};
  enumerate(0, 0);
  return result;
}

std::vector<Inequality> simplex_inequalities(int k) {
  std::vector<Inequality> result;
  for (int i = 0; i < k; ++i) {
    Inequality inequality;
    inequality.coefficients.assign(k, 0);
    inequality.coefficients[i] = 1;
    result.push_back(std::move(inequality));
  }
  Inequality last;
  last.coefficients.assign(k, -1);
  last.constant = 1;
  result.push_back(std::move(last));
  return result;
}

std::vector<Inequality> dominance_inequalities(const std::vector<Vector>& simplex,
                                               const std::vector<AffineFunction>& branches,
                                               int branch_index) {
  const int d = static_cast<int>(simplex.front().size());
  const int intrinsic_dimension = static_cast<int>(simplex.size()) - 1;
  std::vector<Inequality> result = simplex_inequalities(intrinsic_dimension);
  AffineFunction zero;
  zero.coefficients.assign(d, 0);
  std::vector<AffineFunction> all;
  all.push_back(zero);
  all.insert(all.end(), branches.begin(), branches.end());
  const AffineFunction& branch = all[branch_index + 1];
  for (int q = 0; q < static_cast<int>(all.size()); ++q) {
    if (q == branch_index + 1) continue;
    Inequality inequality;
    inequality.coefficients.resize(d);
    for (int i = 0; i < d; ++i) inequality.coefficients[i] = branch.coefficients[i] - all[q].coefficients[i];
    inequality.constant = branch.constant - all[q].constant;
    result.push_back(std::move(inequality));
  }
  return result;
}

Rational integrate_region(const std::vector<Vector>& source,
                          const std::vector<Inequality>& inequalities,
                          const AffineFunction& branch,
                          const AffineFunction* ell,
                          const Rational& boundary_scale = Rational(0)) {
  const int k = static_cast<int>(source.size()) - 1;
  const auto vertices = region_vertices(source, inequalities, k);
  if (vertices.size() < static_cast<std::size_t>(k + 1)) return 0;
  const auto triangulation = detail::triangulate_exact(vertices, k).cells;
  Rational result = 0;
  for (const auto& simplex : triangulation) {
    if (ell != nullptr) {
      const auto moments = simplex_moments(simplex.vertices);
      result += product_integral(branch, *ell, moments);
      continue;
    }
    const int ambient = static_cast<int>(source.front().size());
    std::vector<std::vector<Rational>> basis(k, Vector(ambient));
    for (int i = 0; i < k; ++i) {
      for (int j = 0; j < ambient; ++j) basis[i][j] = source[i + 1][j] - source[0][j];
    }
    std::vector<int> columns;
    std::function<bool(int)> find_columns = [&](int next) {
      if (static_cast<int>(columns.size()) == k) {
        std::vector<std::vector<Rational>> minor(k, Vector(k));
        for (int i = 0; i < k; ++i) for (int j = 0; j < k; ++j) minor[i][j] = basis[i][columns[j]];
        return determinant(minor) != 0;
      }
      for (int column = next; column < ambient; ++column) {
        columns.push_back(column);
        if (find_columns(column + 1)) return true;
        columns.pop_back();
      }
      return false;
    };
    find_columns(0);
    if (static_cast<int>(columns.size()) != k) throw std::runtime_error("cannot choose boundary coordinates");
    std::vector<std::vector<Rational>> matrix(k, Vector(k));
    for (int i = 0; i < k; ++i) for (int j = 0; j < k; ++j) matrix[i][j] = basis[i][columns[j]];
    std::vector<Vector> local;
    local.reserve(simplex.vertices.size());
    for (const auto& vertex : simplex.vertices) {
      Vector rhs(k, 0);
      for (int j = 0; j < k; ++j) rhs[j] = vertex[columns[j]] - source[0][columns[j]];
      local.push_back(solve_linear(matrix, std::move(rhs)));
    }
    std::vector<std::vector<Rational>> edges(k, Vector(k));
    for (int i = 0; i < k; ++i) for (int j = 0; j < k; ++j) edges[i][j] = local[i + 1][j] - local[0][j];
    const Rational normalized_volume = abs_rat(determinant(edges));
    Rational sum = 0;
    for (const auto& vertex : simplex.vertices) sum += evaluate_affine(branch, vertex);
    result += boundary_scale * normalized_volume * sum / (k + 1);
  }
  return result;
}

std::vector<Vector> rational_facet_points(const Facet& facet) {
  std::vector<Vector> result;
  for (const auto& point : facet.vertices) {
    Vector rational;
    for (const auto coordinate : point) rational.push_back(rat(coordinate));
    result.push_back(std::move(rational));
  }
  return result;
}

Rational exact_pl_value(const LatticePolytope& polytope,
                        const ConvexPLFunction& function) {
  if (function.branches.empty() || function.branches.size() > 8) {
    throw std::invalid_argument("PL branch count must be in [1,8]");
  }
  const int d = polytope.dimension;
  std::vector<AffineFunction> branches;
  for (const auto& branch : function.branches) {
    if (static_cast<int>(branch.coefficients.size()) != d) throw std::invalid_argument("PL dimension mismatch");
    bool duplicate = false;
    for (const auto& existing : branches) {
      if (existing.coefficients == branch.coefficients && existing.constant == branch.constant) {
        duplicate = true;
        break;
      }
    }
    if (!duplicate) branches.push_back(branch);
  }
  const AffineFunction ell = ell_function(polytope);
  Rational interior = 0;
  Rational boundary = 0;
  for (const auto& simplex : polytope.simplices) {
    for (std::size_t j = 0; j < branches.size(); ++j) {
      const auto inequalities = dominance_inequalities(simplex.vertices, branches, static_cast<int>(j));
      interior += integrate_region(simplex.vertices, inequalities, branches[j], &ell);
    }
  }
  for (const auto& facet : polytope.boundary_facets) {
    const auto source = rational_facet_points(facet);
    for (std::size_t j = 0; j < branches.size(); ++j) {
      const auto inequalities = dominance_inequalities(source, branches, static_cast<int>(j));
      boundary += integrate_region(source, inequalities, branches[j], nullptr, facet.lattice_measure);
    }
  }
  return boundary - interior;
}

struct DoubleSample {
  std::vector<double> point;
  double weight = 0.0;
  bool boundary = false;
};

std::vector<int> first_primes(int count) {
  std::vector<int> primes;
  for (int candidate = 2; static_cast<int>(primes.size()) < count; ++candidate) {
    bool prime = true;
    for (int p : primes) {
      if (p * p > candidate) break;
      if (candidate % p == 0) { prime = false; break; }
    }
    if (prime) primes.push_back(candidate);
  }
  return primes;
}

double halton(std::size_t index, int base) {
  double result = 0.0;
  double factor = 1.0 / static_cast<double>(base);
  while (index != 0) {
    result += factor * static_cast<double>(index % static_cast<std::size_t>(base));
    index /= static_cast<std::size_t>(base);
    factor /= static_cast<double>(base);
  }
  return result;
}

std::vector<double> simplex_sample(const std::vector<Vector>& simplex,
                                   std::size_t index, int prime_offset) {
  const int k = static_cast<int>(simplex.size()) - 1;
  const auto primes = first_primes(k + prime_offset);
  std::vector<double> cuts;
  cuts.reserve(k + 2);
  cuts.push_back(0.0);
  for (int i = 0; i < k; ++i) cuts.push_back(halton(index + 1, primes[i + prime_offset]));
  cuts.push_back(1.0);
  std::sort(cuts.begin(), cuts.end());
  std::vector<double> bary(k + 1);
  for (int i = 0; i <= k; ++i) bary[i] = cuts[i + 1] - cuts[i];
  std::vector<double> result(simplex.front().size(), 0.0);
  for (int i = 0; i <= k; ++i) for (std::size_t j = 0; j < result.size(); ++j) result[j] += bary[i] * rational_to_double(simplex[i][j]);
  return result;
}

std::vector<DoubleSample> make_samples(const LatticePolytope& polytope,
                                        std::size_t count) {
  if (count < 16) count = 16;
  std::vector<std::vector<Vector>> cells;
  cells.reserve(polytope.simplices.size());
  std::vector<double> cell_weights;
  for (const auto& simplex : polytope.simplices) {
    cells.push_back(simplex.vertices);
    cell_weights.push_back(rational_to_double(simplex_moments(simplex.vertices).volume));
  }
  const double total_volume = rational_to_double(polytope.volume);
  std::vector<DoubleSample> samples;
  samples.reserve(2 * count);
  for (std::size_t s = 0; s < count; ++s) {
    const double pick = halton(s + 1, 2) * total_volume;
    double cumulative = 0.0;
    std::size_t cell = cells.size() - 1;
    for (std::size_t i = 0; i < cell_weights.size(); ++i) {
      cumulative += cell_weights[i];
      if (pick <= cumulative) { cell = i; break; }
    }
    samples.push_back({simplex_sample(cells[cell], s + 1, 1), total_volume / static_cast<double>(count), false});
  }
  std::vector<double> facet_weights;
  double total_boundary = rational_to_double(polytope.boundary_measure);
  double cumulative = 0.0;
  for (std::size_t s = 0; s < count; ++s) {
    const double pick = halton(s + 1, 3) * total_boundary;
    cumulative = 0.0;
    std::size_t facet_index = polytope.boundary_facets.size() - 1;
    for (std::size_t i = 0; i < polytope.boundary_facets.size(); ++i) {
      cumulative += rational_to_double(polytope.boundary_facets[i].lattice_measure);
      if (pick <= cumulative) { facet_index = i; break; }
    }
    std::vector<Vector> facet;
    for (const auto& point : polytope.boundary_facets[facet_index].vertices) {
      Vector rational;
      for (const auto coordinate : point) rational.push_back(rat(coordinate));
      facet.push_back(std::move(rational));
    }
    samples.push_back({simplex_sample(facet, s + 1, 2), total_boundary / static_cast<double>(count), true});
  }
  return samples;
}

struct NumericFunction {
  std::vector<std::vector<double>> coefficients;
  std::vector<double> constants;
};

double numeric_value(const NumericFunction& function, const std::vector<double>& point) {
  double value = 0.0;
  for (std::size_t j = 0; j < function.coefficients.size(); ++j) {
    double branch = function.constants[j];
    for (std::size_t i = 0; i < point.size(); ++i) branch += function.coefficients[j][i] * point[i];
    value = std::max(value, branch);
  }
  return value;
}

double numeric_df(const LatticePolytope& polytope,
                  const std::vector<DoubleSample>& samples,
                  const NumericFunction& function) {
  double interior = 0.0;
  double boundary = 0.0;
  std::vector<double> ell_coefficients;
  for (std::size_t i = 1; i < polytope.ell.size(); ++i) ell_coefficients.push_back(rational_to_double(polytope.ell[i]));
  const double ell_constant = rational_to_double(polytope.ell[0]);
  for (const auto& sample : samples) {
    const double value = numeric_value(function, sample.point);
    if (sample.boundary) boundary += sample.weight * value;
    else interior += sample.weight * value * (ell_constant + std::inner_product(ell_coefficients.begin(), ell_coefficients.end(), sample.point.begin(), 0.0));
  }
  return boundary - interior;
}

NumericFunction decode_parameters(const LatticePolytope& polytope,
                                  const std::vector<double>& parameters,
                                  int pieces) {
  const int d = polytope.dimension;
  const int block = d;
  NumericFunction result;
  result.coefficients.assign(pieces, std::vector<double>(d));
  result.constants.assign(pieces, 0.0);
  std::vector<double> logits(pieces, 0.0);
  for (int j = 1; j < pieces; ++j) logits[j] = parameters[pieces * block + j - 1];
  double max_logit = *std::max_element(logits.begin(), logits.end());
  std::vector<double> weights(pieces);
  double weight_sum = 0.0;
  for (int j = 0; j < pieces; ++j) { weights[j] = std::exp(logits[j] - max_logit); weight_sum += weights[j]; }
  for (double& weight : weights) weight /= weight_sum;
  for (int j = 0; j < pieces; ++j) {
    const int offset = j * block;
    std::vector<double> direction(d, 0.0);
    double product = 1.0;
    for (int i = 0; i < d - 1; ++i) {
      const double angle = parameters[offset + i];
      direction[i] = product * std::cos(angle);
      product *= std::sin(angle);
    }
    direction[d - 1] = product;
    double lo = 0.0, hi = 0.0;
    for (const auto& point : polytope.input_points) {
      double projection = 0.0;
      for (int i = 0; i < d; ++i) projection += direction[i] * static_cast<double>(point[i]);
      lo = std::min(lo, projection); hi = std::max(hi, projection);
    }
    const double t = lo + parameters[offset + d - 1] * (hi - lo);
    for (int i = 0; i < d; ++i) result.coefficients[j][i] = weights[j] * direction[i];
    result.constants[j] = -weights[j] * t;
  }
  return result;
}

ConvexPLFunction rationalize_numeric(const NumericFunction& numeric,
                                     std::int64_t cap) {
  ConvexPLFunction result;
  for (std::size_t j = 0; j < numeric.coefficients.size(); ++j) {
    AffineFunction branch;
    for (const double coefficient : numeric.coefficients[j]) branch.coefficients.push_back(approximate_rational(coefficient, cap));
    branch.constant = approximate_rational(numeric.constants[j], cap);
    result.branches.push_back(std::move(branch));
  }
  return result;
}

NumericFunction numeric_from_function(const ConvexPLFunction& function) {
  NumericFunction result;
  for (const auto& branch : function.branches) {
    std::vector<double> coefficients;
    for (const auto& coefficient : branch.coefficients) coefficients.push_back(rational_to_double(coefficient));
    result.coefficients.push_back(std::move(coefficients));
    result.constants.push_back(rational_to_double(branch.constant));
  }
  return result;
}

double polytope_diameter(const LatticePolytope& polytope) {
  double diameter = 0.0;
  for (const auto& p : polytope.input_points) {
    for (const auto& q : polytope.input_points) {
      double squared = 0.0;
      for (int i = 0; i < polytope.dimension; ++i) {
        const double difference = static_cast<double>(p[i] - q[i]);
        squared += difference * difference;
      }
      diameter = std::max(diameter, std::sqrt(squared));
    }
  }
  return diameter;
}

double normalization_scale(const LatticePolytope& polytope) {
  double sup_ell = 0.0;
  for (const auto& point : polytope.input_points) {
    double value = rational_to_double(polytope.ell[0]);
    for (int i = 0; i < polytope.dimension; ++i) value += rational_to_double(polytope.ell[i + 1]) * static_cast<double>(point[i]);
    sup_ell = std::max(sup_ell, std::abs(value));
  }
  return std::max(1.0, polytope_diameter(polytope) *
                           (rational_to_double(polytope.boundary_measure) +
                            rational_to_double(polytope.volume) * sup_ell));
}

Rational evaluate_df_exact(const LatticePolytope& polytope,
                           const ConvexPLFunction& function) {
  return exact_pl_value(polytope, function);
}

Rational evaluate_df_exact(const LatticePolytope& polytope,
                           const AffineFunction& function) {
  if (static_cast<int>(function.coefficients.size()) != polytope.dimension) {
    throw std::invalid_argument("affine dimension mismatch");
  }
  Rational boundary = function.constant * polytope.boundary_measure;
  for (int i = 0; i < polytope.dimension; ++i) {
    boundary += function.coefficients[i] * polytope.boundary_first_moment[i];
  }
  const AffineFunction ell = ell_function(polytope);
  const auto moments = SimplexMoments{polytope.volume, polytope.first_moment,
                                     polytope.second_moment};
  return boundary - product_integral(function, ell, moments);
}

Rational approximate_rational(double value, std::int64_t cap) {
  if (cap < 1 || !std::isfinite(value) || std::abs(value) > 1e15) {
    throw std::invalid_argument("invalid rational approximation request");
  }
  const bool negative = value < 0.0;
  const double target = std::abs(value);
  std::int64_t h_m2 = 0, h_m1 = 1, k_m2 = 1, k_m1 = 0;
  double remainder = target;
  for (int iteration = 0; iteration < 128; ++iteration) {
    const double integer_part = std::floor(remainder);
    if (integer_part > 4.0e18) break;
    const std::int64_t a = static_cast<std::int64_t>(integer_part);
    const __int128 h_wide = static_cast<__int128>(a) * h_m1 + h_m2;
    const __int128 k_wide = static_cast<__int128>(a) * k_m1 + k_m2;
    if (k_wide > cap || h_wide > std::numeric_limits<std::int64_t>::max()) {
      const std::int64_t t = k_m1 > 0 ? (cap - k_m2) / k_m1 : 0;
      Rational candidate(CGAL::Gmpz(static_cast<long>(h_m1)), CGAL::Gmpz(static_cast<long>(k_m1)));
      if (t >= 1) {
        const std::int64_t hs = t * h_m1 + h_m2;
        const std::int64_t ks = t * k_m1 + k_m2;
        const double semi_error = std::abs(target - static_cast<double>(hs) / ks);
        const double previous_error = std::abs(target - static_cast<double>(h_m1) / k_m1);
        if (semi_error <= previous_error) candidate = Rational(CGAL::Gmpz(static_cast<long>(hs)), CGAL::Gmpz(static_cast<long>(ks)));
      }
      return negative ? -candidate : candidate;
    }
    const std::int64_t h = static_cast<std::int64_t>(h_wide);
    const std::int64_t k = static_cast<std::int64_t>(k_wide);
    h_m2 = h_m1; h_m1 = h; k_m2 = k_m1; k_m1 = k;
    if (remainder == integer_part) break;
    remainder = 1.0 / (remainder - integer_part);
  }
  Rational result(CGAL::Gmpz(static_cast<long>(h_m1)), CGAL::Gmpz(static_cast<long>(k_m1)));
  return negative ? -result : result;
}

std::string format_affine(const AffineFunction& function) {
  std::ostringstream stream;
  for (std::size_t i = 0; i < function.coefficients.size(); ++i) {
    if (i != 0) stream << ' ';
    stream << rational_string(function.coefficients[i]);
  }
  if (!function.coefficients.empty()) stream << ' ';
  stream << rational_string(function.constant);
  return stream.str();
}

NumericalWitness search_pl_witness(const LatticePolytope& polytope,
                                   const SearchOptions& options) {
  if (options.pieces < 1 || options.pieces > 8) throw std::invalid_argument("pieces must be in [1,8]");
  const int d = polytope.dimension;
  const std::size_t parameter_dimension = static_cast<std::size_t>(options.pieces * d + options.pieces - 1);
  const std::size_t population = options.population == 0
                                     ? std::max<std::size_t>(32, 8 * parameter_dimension)
                                     : options.population;
  auto samples = make_samples(polytope, options.quadrature_samples);
  std::vector<double> lower(parameter_dimension), upper(parameter_dimension);
  constexpr double pi = 3.141592653589793238462643383279502884;
  for (int branch = 0; branch < options.pieces; ++branch) {
    const std::size_t offset = static_cast<std::size_t>(branch * d);
    for (int angle = 0; angle < d - 1; ++angle) {
      lower[offset + angle] = 0.0;
      upper[offset + angle] = (angle == d - 2 ? 2.0 * pi : pi);
    }
    lower[offset + d - 1] = 0.0;
    upper[offset + d - 1] = 1.0;
  }
  for (int branch = 1; branch < options.pieces; ++branch) {
    lower[static_cast<std::size_t>(options.pieces * d + branch - 1)] = -6.0;
    upper[static_cast<std::size_t>(options.pieces * d + branch - 1)] = 6.0;
  }
  std::atomic<long> evaluations{0};
  auto objective = [&](const std::vector<double>& parameters) {
    ++evaluations;
    const auto function = decode_parameters(polytope, parameters, options.pieces);
    return numeric_df(polytope, samples, function);
  };
  boost::math::optimization::differential_evolution_parameters<std::vector<double>> parameters;
  parameters.lower_bounds = lower;
  parameters.upper_bounds = upper;
  parameters.NP = population;
  parameters.max_generations = options.generations;
  parameters.mutation_factor = 0.7;
  parameters.crossover_probability = 0.8;
  parameters.threads = std::max(1u, options.threads);
  std::mt19937_64 generator(options.seed);
  std::vector<std::pair<std::vector<double>, double>> queries;
  const auto best_parameters = boost::math::optimization::differential_evolution(
      objective, parameters, generator, std::numeric_limits<double>::quiet_NaN(),
      nullptr, nullptr, &queries);
  const auto best_numeric = decode_parameters(polytope, best_parameters, options.pieces);
  NumericalWitness result;
  result.function = rationalize_numeric(best_numeric, 1000000);
  result.value = numeric_df(polytope, samples, best_numeric);
  result.normalized = result.value / normalization_scale(polytope);
  result.evaluations = evaluations.load();
  std::sort(queries.begin(), queries.end(),
            [](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
  auto signature = [](const ConvexPLFunction& function) {
    std::vector<std::string> branches;
    branches.reserve(function.branches.size());
    for (const auto& branch : function.branches) branches.push_back(format_affine(branch));
    std::sort(branches.begin(), branches.end());
    std::ostringstream stream;
    for (const auto& branch : branches) stream << branch << ';';
    return stream.str();
  };
  std::set<std::string> seen;
  seen.insert(signature(result.function));
  for (const auto& query : queries) {
    if (result.alternatives.size() >= 15) break;
    const auto candidate = decode_parameters(polytope, query.first, options.pieces);
    auto rational_candidate = rationalize_numeric(candidate, 1000000);
    if (seen.insert(signature(rational_candidate)).second) {
      result.alternatives.push_back(std::move(rational_candidate));
    }
  }
  if (options.verbose) {
    std::cerr << "search pieces=" << options.pieces << " evaluations=" << result.evaluations
              << " M=" << std::setprecision(17) << result.value << '\n';
  }
  return result;
}

CertifyResult certify_function(const LatticePolytope& polytope,
                               const ConvexPLFunction& function) {
  CertifyResult result;
  result.function = function;
  result.value = exact_pl_value(polytope, function);
  result.certified = result.value < 0;
  return result;
}

CertifyResult certify_witness(const LatticePolytope& polytope,
                              const NumericalWitness& witness,
                              std::int64_t max_denominator) {
  if (max_denominator < 1) throw std::invalid_argument("certification denominator must be positive");
  const NumericFunction numeric = numeric_from_function(witness.function);
  CertifyResult result;
  std::vector<ConvexPLFunction> candidates;
  candidates.push_back(witness.function);
  candidates.insert(candidates.end(), witness.alternatives.begin(), witness.alternatives.end());
  std::int64_t cap = 10;
  while (true) {
    for (const auto& seed : candidates) {
      const auto candidate = certify_function(polytope, rationalize_numeric(numeric_from_function(seed), cap));
      if (candidate.certified) return candidate;
    }
    if (cap >= max_denominator) break;
    if (cap > max_denominator / 10) cap = max_denominator;
    else cap *= 10;
  }
  result.function = rationalize_numeric(numeric, max_denominator);
  result.value = exact_pl_value(polytope, result.function);
  result.certified = result.value < 0;
  if (!result.certified) {
    for (const auto& seed : candidates) {
      result.function = rationalize_numeric(numeric_from_function(seed), max_denominator);
      result.value = exact_pl_value(polytope, result.function);
      if (result.value < 0) {
        result.certified = true;
        break;
      }
    }
  }
  return result;
}

}  // namespace kstab_highdim
