#include "decompose.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <queue>

namespace kstab {
namespace {

Rational rat64(std::int64_t value) { return Rational(static_cast<long>(value)); }

void validate(const MeasuredPolygon& polygon) {
  if (polygon.vertices_ccw.size() < 3 ||
      polygon.vertices_ccw.size() != polygon.edge_measures.size()) {
    throw std::invalid_argument(
        "MeasuredPolygon requires at least three vertices and one measure per edge.");
  }
  if (polygon_moments(polygon.vertices_ccw).area <= 0) {
    throw std::invalid_argument("MeasuredPolygon must be counterclockwise with positive area.");
  }
  for (std::size_t i = 0; i < polygon.vertices_ccw.size(); ++i) {
    const QPoint& p = polygon.vertices_ccw[i];
    const QPoint& q = polygon.vertices_ccw[(i + 1) % polygon.vertices_ccw.size()];
    const QPoint& r = polygon.vertices_ccw[(i + 2) % polygon.vertices_ccw.size()];
    const Rational turn = (q.x - p.x) * (r.y - p.y) -
                          (q.y - p.y) * (r.x - p.x);
    if (turn <= 0) {
      throw std::invalid_argument(
          "MeasuredPolygon must be strictly convex with no collinear vertices.");
    }
  }
  for (const Rational& measure : polygon.edge_measures) {
    if (measure < 0) throw std::invalid_argument("Boundary measures must be nonnegative.");
  }
}

Rational det3(const std::array<std::array<Rational, 3>, 3>& m) {
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
         m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
         m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

std::array<Rational, 3> solve3(
    std::array<std::array<Rational, 3>, 3> matrix,
    const std::array<Rational, 3>& rhs) {
  const Rational determinant = det3(matrix);
  if (determinant == 0) throw std::runtime_error("The 3x3 moment system is singular.");
  std::array<Rational, 3> result;
  for (int column = 0; column < 3; ++column) {
    auto replaced = matrix;
    for (int row = 0; row < 3; ++row) replaced[row][column] = rhs[row];
    result[column] = det3(replaced) / determinant;
  }
  return result;
}

QPoint interpolate(const QPoint& p, const QPoint& q, const Rational& t) {
  return {p.x + t * (q.x - p.x), p.y + t * (q.y - p.y)};
}

Rational value(const std::array<Rational, 3>& ell, const QPoint& point) {
  return ell[0] + ell[1] * point.x + ell[2] * point.y;
}

std::vector<QPoint> clip_exact(const std::vector<QPoint>& polygon,
                               const Rational& a, const Rational& b,
                               const Rational& c) {
  std::vector<QPoint> output;
  for (std::size_t i = 0; i < polygon.size(); ++i) {
    const QPoint& p = polygon[i];
    const QPoint& q = polygon[(i + 1) % polygon.size()];
    const Rational sp = a * p.x + b * p.y + c;
    const Rational sq = a * q.x + b * q.y + c;
    if (sp >= 0) output.push_back(p);
    if ((sp >= 0) != (sq >= 0)) output.push_back(interpolate(p, q, sp / (sp - sq)));
  }
  return output;
}

struct DMeasuredPolygon {
  std::vector<DPoint> vertices;
  std::vector<double> measures;
};

DMeasuredPolygon to_double(const MeasuredPolygon& polygon) {
  DMeasuredPolygon result;
  result.vertices.reserve(polygon.vertices_ccw.size());
  result.measures.reserve(polygon.edge_measures.size());
  for (const QPoint& point : polygon.vertices_ccw) {
    result.vertices.push_back({rational_to_double(point.x), rational_to_double(point.y)});
  }
  for (const Rational& measure : polygon.edge_measures) {
    result.measures.push_back(rational_to_double(measure));
  }
  return result;
}

std::vector<DPoint> clip_double(const std::vector<DPoint>& polygon, double a,
                                double b, double c) {
  std::vector<DPoint> result;
  for (std::size_t i = 0; i < polygon.size(); ++i) {
    const DPoint& p = polygon[i];
    const DPoint& q = polygon[(i + 1) % polygon.size()];
    const double sp = a * p.x + b * p.y + c;
    const double sq = a * q.x + b * q.y + c;
    if (sp >= 0.0) result.push_back(p);
    if ((sp >= 0.0) != (sq >= 0.0)) {
      const double t = std::clamp(sp / (sp - sq), 0.0, 1.0);
      result.push_back({p.x + t * (q.x - p.x), p.y + t * (q.y - p.y)});
    }
  }
  return result;
}

double df_double(const DMeasuredPolygon& polygon,
                 const std::array<double, 3>& ell, double a, double b,
                 double c, long* evaluations) {
  if (evaluations != nullptr) ++*evaluations;
  const DMoments moments = polygon_moments_double(clip_double(polygon.vertices, a, b, c));
  const double interior =
      a * ell[1] * moments.ixx + b * ell[2] * moments.iyy +
      (a * ell[2] + b * ell[1]) * moments.ixy +
      (a * ell[0] + c * ell[1]) * moments.ix +
      (b * ell[0] + c * ell[2]) * moments.iy + c * ell[0] * moments.area;
  double boundary = 0.0;
  for (std::size_t i = 0; i < polygon.vertices.size(); ++i) {
    if (polygon.measures[i] == 0.0) continue;
    const DPoint& p = polygon.vertices[i];
    const DPoint& q = polygon.vertices[(i + 1) % polygon.vertices.size()];
    const double sp = a * p.x + b * p.y + c;
    const double sq = a * q.x + b * q.y + c;
    if (sp <= 0.0 && sq <= 0.0) continue;
    double lo = 0.0, hi = 1.0;
    if (sp < 0.0) lo = std::clamp(sp / (sp - sq), 0.0, 1.0);
    else if (sq < 0.0) hi = std::clamp(sp / (sp - sq), 0.0, 1.0);
    const double g0 = sp + lo * (sq - sp);
    const double g1 = sp + hi * (sq - sp);
    boundary += polygon.measures[i] * (hi - lo) * (g0 + g1) / 2.0;
  }
  return boundary - interior;
}

double affine_at(const std::array<Rational, 3>& ell, double x, double y) {
  return rational_to_double(ell[0]) + rational_to_double(ell[1]) * x +
         rational_to_double(ell[2]) * y;
}

std::array<double, 2> residual(const MeasuredPolygon& parent, std::size_t i,
                               double s, std::size_t j, double t) {
  try {
    const ChordSplit split = split_by_chord(parent, i, Rational(s), j, Rational(t));
    const auto first = compute_ell_p(split.first_piece);
    const auto second = compute_ell_p(split.second_piece);
    return {affine_at(first, rational_to_double(split.first.x), rational_to_double(split.first.y)) -
                affine_at(second, rational_to_double(split.first.x), rational_to_double(split.first.y)),
            affine_at(first, rational_to_double(split.second.x), rational_to_double(split.second.y)) -
                affine_at(second, rational_to_double(split.second.x), rational_to_double(split.second.y))};
  } catch (const std::exception&) {
    return {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
  }
}

double norm2(const std::array<double, 2>& r) { return std::hypot(r[0], r[1]); }

bool exact_continuity(const ChordSplit& split,
                      const std::array<Rational, 3>& first,
                      const std::array<Rational, 3>& second) {
  return value(first, split.first) == value(second, split.first) &&
         value(first, split.second) == value(second, split.second);
}

bool solve_root(const MeasuredPolygon& parent, std::size_t i, std::size_t j,
                double& s, double& t, const DecomposeOptions& options) {
  constexpr double kMin = 1e-7;
  for (int iteration = 0; iteration < options.root_iterations; ++iteration) {
    const auto f = residual(parent, i, s, j, t);
    const double base = norm2(f);
    if (!std::isfinite(base)) return false;
    if (base <= options.root_tolerance) return true;
    const double h = 1e-5;
    const auto fs = residual(parent, i, std::clamp(s + h, kMin, 1.0 - kMin), j, t);
    const auto ft = residual(parent, i, s, j, std::clamp(t + h, kMin, 1.0 - kMin));
    const double j00 = (fs[0] - f[0]) / h;
    const double j10 = (fs[1] - f[1]) / h;
    const double j01 = (ft[0] - f[0]) / h;
    const double j11 = (ft[1] - f[1]) / h;
    const double determinant = j00 * j11 - j01 * j10;
    if (std::fabs(determinant) < 1e-14 || !std::isfinite(determinant)) return false;
    const double ds = (-f[0] * j11 + j01 * f[1]) / determinant;
    const double dt = (j10 * f[0] - j00 * f[1]) / determinant;
    bool improved = false;
    for (double scale = 1.0; scale >= 1.0 / 256.0; scale *= 0.5) {
      const double next_s = std::clamp(s + scale * ds, kMin, 1.0 - kMin);
      const double next_t = std::clamp(t + scale * dt, kMin, 1.0 - kMin);
      if (norm2(residual(parent, i, next_s, j, next_t)) < base) {
        s = next_s;
        t = next_t;
        improved = true;
        break;
      }
    }
    if (!improved) return false;
  }
  return norm2(residual(parent, i, s, j, t)) <= options.root_tolerance;
}

bool solve_root_one_variable(const MeasuredPolygon& parent, std::size_t i,
                             double fixed_s, std::size_t j, double& t,
                             bool first_fixed,
                             const DecomposeOptions& options) {
  constexpr double kMin = 1e-7;
  for (int iteration = 0; iteration < options.root_iterations; ++iteration) {
    const auto f = first_fixed ? residual(parent, i, fixed_s, j, t)
                               : residual(parent, i, t, j, fixed_s);
    const double base = norm2(f);
    if (!std::isfinite(base)) return false;
    if (base <= options.root_tolerance) return true;
    const double h = 1e-5;
    const double next = std::clamp(t + h, kMin, 1.0 - kMin);
    const auto g = first_fixed ? residual(parent, i, fixed_s, j, next)
                               : residual(parent, i, next, j, fixed_s);
    const double derivative = ((g[0] - f[0]) * f[0] + (g[1] - f[1]) * f[1]) / h;
    const double curvature = ((g[0] - f[0]) * (g[0] - f[0]) +
                              (g[1] - f[1]) * (g[1] - f[1])) / (h * h);
    if (!(curvature > 1e-18)) return false;
    const double step = -derivative / curvature;
    bool improved = false;
    for (double scale = 1.0; scale >= 1.0 / 256.0; scale *= 0.5) {
      const double candidate = std::clamp(t + scale * step, kMin, 1.0 - kMin);
      const auto r = first_fixed ? residual(parent, i, fixed_s, j, candidate)
                                 : residual(parent, i, candidate, j, fixed_s);
      if (norm2(r) < base) { t = candidate; improved = true; break; }
    }
    if (!improved) return false;
  }
  const auto f = first_fixed ? residual(parent, i, fixed_s, j, t)
                             : residual(parent, i, t, j, fixed_s);
  return norm2(f) <= options.root_tolerance;
}

bool concave(const ChordSplit& split, const std::array<Rational, 3>& first,
             const std::array<Rational, 3>& second) {
  const auto centroid = [](const MeasuredPolygon& polygon) {
    double x = 0.0, y = 0.0;
    for (const QPoint& point : polygon.vertices_ccw) {
      x += rational_to_double(point.x);
      y += rational_to_double(point.y);
    }
    const double n = static_cast<double>(polygon.vertices_ccw.size());
    return DPoint{x / n, y / n};
  };
  const DPoint first_probe = centroid(split.first_piece);
  const DPoint second_probe = centroid(split.second_piece);
  constexpr double tolerance = 1e-8;
  return affine_at(first, first_probe.x, first_probe.y) <=
             affine_at(second, first_probe.x, first_probe.y) + tolerance &&
         affine_at(second, second_probe.x, second_probe.y) <=
             affine_at(first, second_probe.x, second_probe.y) + tolerance;
}

std::optional<DecompositionCandidate> make_candidate(
    const MeasuredPolygon& parent, std::size_t i, double s, std::size_t j,
    double t, const DecomposeOptions& options, bool require_exact) {
  try {
    const ChordSplit split = split_by_chord(parent, i, Rational(s), j, Rational(t));
    const auto first = compute_ell_p(split.first_piece);
    const auto second = compute_ell_p(split.second_piece);
    const auto r = residual(parent, i, s, j, t);
    if (!std::isfinite(norm2(r)) || norm2(r) > options.root_tolerance ||
        !concave(split, first, second)) return std::nullopt;
    DecompositionCandidate candidate;
    candidate.first_edge = i;
    candidate.second_edge = j;
    candidate.first_parameter = s;
    candidate.second_parameter = t;
    candidate.continuity_residual = norm2(r);
    candidate.concave = true;
    candidate.first_ell = first;
    candidate.second_ell = second;
    candidate.first_piece = split.first_piece;
    candidate.second_piece = split.second_piece;
    candidate.first_endpoint = split.first;
    candidate.second_endpoint = split.second;
    const Rational sr = approximate_rational(s, options.rational_max_denominator);
    const Rational tr = approximate_rational(t, options.rational_max_denominator);
    const ChordSplit exact_split = split_by_chord(parent, i, sr, j, tr);
    const auto exact_first = compute_ell_p(exact_split.first_piece);
    const auto exact_second = compute_ell_p(exact_split.second_piece);
    if (exact_continuity(exact_split, exact_first, exact_second) &&
        concave(exact_split, exact_first, exact_second)) {
      candidate.certified_rational = true;
      candidate.first_parameter_rational = sr;
      candidate.second_parameter_rational = tr;
      candidate.first_piece = exact_split.first_piece;
      candidate.second_piece = exact_split.second_piece;
      candidate.first_endpoint = exact_split.first;
      candidate.second_endpoint = exact_split.second;
      candidate.first_ell = exact_first;
      candidate.second_ell = exact_second;
    }
    if (require_exact && !candidate.certified_rational) return std::nullopt;
    candidate.first_search = search_witness(candidate.first_piece, candidate.first_ell,
                                             options.witness_search);
    candidate.second_search = search_witness(candidate.second_piece, candidate.second_ell,
                                              options.witness_search);
    if (candidate.first_search.unstable || candidate.second_search.unstable) {
      return std::nullopt;
    }
    return candidate;
  } catch (const std::exception&) {
    return std::nullopt;
  }
}

}  // namespace

std::string rational_string(const Rational& value) {
  std::ostringstream stream;
  if (value.denominator() == 1) stream << value.numerator();
  else stream << value;
  return stream.str();
}

MeasuredPolygon make_measured_polygon(const std::vector<IntPoint>& vertices,
                                      const std::vector<bool>& null_edges) {
  if (!null_edges.empty() && null_edges.size() != vertices.size()) {
    throw std::invalid_argument("null measure mask size mismatch.");
  }
  MeasuredPolygon result;
  for (const IntPoint& point : vertices) result.vertices_ccw.push_back({rat64(point.x), rat64(point.y)});
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    const IntPoint& p = vertices[i];
    const IntPoint& q = vertices[(i + 1) % vertices.size()];
    const auto length = std::gcd<std::int64_t>(p.x - q.x, p.y - q.y);
    result.edge_measures.push_back(!null_edges.empty() && null_edges[i] ? 0 : rat64(length));
  }
  validate(result);
  return result;
}

BoundaryMoments boundary_moments(const MeasuredPolygon& polygon) {
  validate(polygon);
  BoundaryMoments result;
  for (std::size_t i = 0; i < polygon.vertices_ccw.size(); ++i) {
    const QPoint& p = polygon.vertices_ccw[i];
    const QPoint& q = polygon.vertices_ccw[(i + 1) % polygon.vertices_ccw.size()];
    const Rational& measure = polygon.edge_measures[i];
    result.length += measure;
    result.ix += measure * (p.x + q.x) / 2;
    result.iy += measure * (p.y + q.y) / 2;
  }
  return result;
}

std::array<Rational, 3> compute_ell_p(const MeasuredPolygon& polygon) {
  validate(polygon);
  const Moments moments = polygon_moments(polygon.vertices_ccw);
  const BoundaryMoments boundary = boundary_moments(polygon);
  return solve3({{{moments.area, moments.ix, moments.iy},
                  {moments.ix, moments.ixx, moments.ixy},
                  {moments.iy, moments.ixy, moments.iyy}}},
                {boundary.length, boundary.ix, boundary.iy});
}

Rational df_simple_exact(const MeasuredPolygon& polygon,
                         const std::array<Rational, 3>& ell,
                         const Rational& a, const Rational& b,
                         const Rational& c) {
  validate(polygon);
  const Moments moments = polygon_moments(clip_exact(polygon.vertices_ccw, a, b, c));
  const Rational interior =
      a * ell[1] * moments.ixx + b * ell[2] * moments.iyy +
      (a * ell[2] + b * ell[1]) * moments.ixy +
      (a * ell[0] + c * ell[1]) * moments.ix +
      (b * ell[0] + c * ell[2]) * moments.iy + c * ell[0] * moments.area;
  Rational boundary = 0;
  for (std::size_t i = 0; i < polygon.vertices_ccw.size(); ++i) {
    const QPoint& p = polygon.vertices_ccw[i];
    const QPoint& q = polygon.vertices_ccw[(i + 1) % polygon.vertices_ccw.size()];
    const Rational sp = a * p.x + b * p.y + c;
    const Rational sq = a * q.x + b * q.y + c;
    if (sp <= 0 && sq <= 0) continue;
    Rational lo = 0, hi = 1;
    if (sp < 0) lo = sp / (sp - sq);
    else if (sq < 0) hi = sp / (sp - sq);
    const Rational g0 = sp + lo * (sq - sp);
    const Rational g1 = sp + hi * (sq - sp);
    boundary += polygon.edge_measures[i] * (hi - lo) * (g0 + g1) / 2;
  }
  return boundary - interior;
}

MeasuredSearchResult search_witness(const MeasuredPolygon& polygon,
                                    const std::array<Rational, 3>& ell,
                                    const SearchOptions& options) {
  validate(polygon);
  const DMeasuredPolygon doubles = to_double(polygon);
  const auto ell_d = ell_to_double(ell);
  double min_x = doubles.vertices.front().x, min_y = doubles.vertices.front().y;
  double boundary_length = 0.0, diameter = 0.0, sup_ell = 0.0;
  for (std::size_t i = 0; i < doubles.vertices.size(); ++i) {
    boundary_length += doubles.measures[i];
    min_x = std::min(min_x, doubles.vertices[i].x);
    min_y = std::min(min_y, doubles.vertices[i].y);
    sup_ell = std::max(sup_ell, std::fabs(affine_at(ell, doubles.vertices[i].x, doubles.vertices[i].y)));
    for (const DPoint& other : doubles.vertices) diameter = std::max(diameter, std::hypot(doubles.vertices[i].x - other.x, doubles.vertices[i].y - other.y));
  }
  MeasuredSearchResult result;
  result.normalization = boundary_length * std::max(sup_ell, 1e-300) * diameter;
  double best = std::numeric_limits<double>::infinity(), best_theta = 0.0, best_t = 0.0;
  constexpr double kTwoPi = 6.2831853071795864769;
  for (int k = 0; k < options.theta_steps; ++k) {
    const double theta = kTwoPi * k / options.theta_steps;
    const double a = std::cos(theta), b = std::sin(theta);
    std::vector<double> values;
    for (const DPoint& point : doubles.vertices) values.push_back(a * point.x + b * point.y);
    const auto [lo_it, hi_it] = std::minmax_element(values.begin(), values.end());
    const double lo = *lo_it, hi = *hi_it;
    for (int sample = 0; sample <= options.t_steps; ++sample) {
      const double t = lo + (hi - lo) * sample / options.t_steps;
      const double candidate = df_double(doubles, ell_d, a, b, -t, &result.evaluations);
      if (candidate < best) { best = candidate; best_theta = theta; best_t = t; }
    }
  }
  result.witness = {std::cos(best_theta), std::sin(best_theta), best_t, best,
                    result.normalization > 0 ? best / result.normalization : best};
  result.unstable = result.witness.normalized < -1e-6;
  return result;
}

ChordSplit split_by_chord(const MeasuredPolygon& parent, std::size_t edge_i,
                          const Rational& s, std::size_t edge_j,
                          const Rational& t) {
  validate(parent);
  const std::size_t n = parent.vertices_ccw.size();
  if (edge_i >= n || edge_j >= n || edge_i == edge_j || s < 0 || s > 1 || t < 0 || t > 1) {
    throw std::invalid_argument("Invalid chord parameters.");
  }
  const QPoint first = interpolate(parent.vertices_ccw[edge_i], parent.vertices_ccw[(edge_i + 1) % n], s);
  const QPoint second = interpolate(parent.vertices_ccw[edge_j], parent.vertices_ccw[(edge_j + 1) % n], t);
  if (first.x == second.x && first.y == second.y) throw std::invalid_argument("Degenerate chord.");
  const auto make_piece = [&](std::size_t start_edge, const Rational& start_t,
                              const QPoint& start, std::size_t end_edge,
                              const Rational& end_t, const QPoint& end) {
    MeasuredPolygon piece;
    piece.vertices_ccw.push_back(start);
    std::size_t edge = start_edge;
    Rational offset = start_t;
    while (true) {
      const QPoint& edge_end = parent.vertices_ccw[(edge + 1) % n];
      if (edge == end_edge) {
        piece.vertices_ccw.push_back(end);
        piece.edge_measures.push_back(parent.edge_measures[edge] * (end_t - offset));
        break;
      }
      piece.vertices_ccw.push_back(edge_end);
      piece.edge_measures.push_back(parent.edge_measures[edge] * (1 - offset));
      edge = (edge + 1) % n;
      offset = 0;
    }
    piece.edge_measures.push_back(0);  // New chord end -> start.
    validate(piece);
    return piece;
  };
  ChordSplit result;
  result.first = first;
  result.second = second;
  result.first_piece = make_piece(edge_i, s, first, edge_j, t, second);
  result.second_piece = make_piece(edge_j, t, second, edge_i, s, first);
  return result;
}

std::vector<DecompositionCandidate> find_chord_decompositions(
    const MeasuredPolygon& parent, const DecomposeOptions& options) {
  validate(parent);
  std::vector<DecompositionCandidate> result;
  std::set<std::tuple<std::size_t, std::size_t, long long, long long>> seen;
  const std::size_t n = parent.vertices_ccw.size();
  const int starts = std::max(2, options.root_starts);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 2; j < n; ++j) {
      if (i == 0 && j + 1 == n) continue;
      for (int si = 1; si < starts; ++si) for (int tj = 1; tj < starts; ++tj) {
        double s = static_cast<double>(si) / starts;
        double t = static_cast<double>(tj) / starts;
        if (!solve_root(parent, i, j, s, t, options)) continue;
        const auto key = std::make_tuple(i, j, std::llround(s * 1e8), std::llround(t * 1e8));
        if (!seen.insert(key).second) continue;
        if (auto candidate = make_candidate(parent, i, s, j, t, options, false)) {
          result.push_back(std::move(*candidate));
        }
      }

      // Endpoint cases are lower-dimensional: fix one parameter at an edge
      // endpoint and solve the remaining least-squares condition.
      // A vertex is represented uniquely as parameter 0 on its outgoing
      // boundary edge; using parameter 1 would duplicate that vertex.
      for (const double fixed : {0.0}) {
        for (int sample = 1; sample < 2 * starts; ++sample) {
          const double moving = static_cast<double>(sample) / (2 * starts);
          double t = moving;
          if (solve_root_one_variable(parent, i, fixed, j, t, true, options)) {
            if (auto candidate = make_candidate(parent, i, fixed, j, t, options, true)) {
              result.push_back(std::move(*candidate));
            }
          }
          double s = moving;
          if (solve_root_one_variable(parent, i, fixed, j, s, false, options)) {
            if (auto candidate = make_candidate(parent, i, s, j, fixed, options, true)) {
              result.push_back(std::move(*candidate));
            }
          }
        }
      }
      if (auto candidate = make_candidate(parent, i, 0.0, j, 0.0, options, true)) {
        result.push_back(std::move(*candidate));
      }
    }
  }
  return result;
}

TreeTopologySummary validate_tree_topology(
    const DecompositionTreeTopology& topology, std::size_t parent_edge_count) {
  if (topology.nodes.empty() || parent_edge_count == 0) {
    throw std::invalid_argument("Tree topology needs nodes and a parent boundary.");
  }
  if (topology.edges.size() + 1 != topology.nodes.size()) {
    throw std::invalid_argument("Decomposition topology must have |E| = |V|-1.");
  }
  std::vector<std::vector<std::size_t>> incident(topology.nodes.size());
  for (std::size_t edge = 0; edge < topology.edges.size(); ++edge) {
    const auto [first, second] = topology.edges[edge];
    if (first >= topology.nodes.size() || second >= topology.nodes.size() ||
        first == second) {
      throw std::invalid_argument("Tree edge has invalid endpoints.");
    }
    incident[first].push_back(edge);
    incident[second].push_back(edge);
  }
  std::vector<bool> visited(topology.nodes.size(), false);
  std::queue<std::size_t> queue;
  queue.push(0);
  visited[0] = true;
  while (!queue.empty()) {
    const std::size_t node = queue.front();
    queue.pop();
    for (const std::size_t edge : incident[node]) {
      const auto [first, second] = topology.edges[edge];
      const std::size_t neighbor = first == node ? second : first;
      if (!visited[neighbor]) {
        visited[neighbor] = true;
        queue.push(neighbor);
      }
    }
  }
  if (std::find(visited.begin(), visited.end(), false) != visited.end()) {
    throw std::invalid_argument("Decomposition topology must be connected.");
  }

  TreeTopologySummary summary;
  for (std::size_t node = 0; node < topology.nodes.size(); ++node) {
    const DecompositionTreeNode& descriptor = topology.nodes[node];
    if (descriptor.kind == DecompositionNodeKind::interior) {
      ++summary.interior_nodes;
      if (incident[node].size() != 3 || descriptor.cyclic_edges.size() != 3) {
        throw std::invalid_argument("Interior tree nodes must have degree 3 and a cyclic order.");
      }
      auto expected = incident[node];
      auto actual = descriptor.cyclic_edges;
      std::sort(expected.begin(), expected.end());
      std::sort(actual.begin(), actual.end());
      if (expected != actual) {
        throw std::invalid_argument("Interior cyclic order must list exactly its incident edges.");
      }
      if (descriptor.boundary_edge || descriptor.boundary_vertex) {
        throw std::invalid_argument("Interior tree nodes cannot have boundary placement.");
      }
    } else {
      ++summary.boundary_leaves;
      if (incident[node].size() != 1 || !descriptor.cyclic_edges.empty()) {
        throw std::invalid_argument("Boundary tree nodes must be degree-1 leaves.");
      }
      if (descriptor.boundary_edge.has_value() == descriptor.boundary_vertex.has_value()) {
        throw std::invalid_argument("Boundary leaf must specify exactly one parent edge or vertex.");
      }
      if ((descriptor.boundary_edge && *descriptor.boundary_edge >= parent_edge_count) ||
          (descriptor.boundary_vertex && *descriptor.boundary_vertex >= parent_edge_count)) {
        throw std::invalid_argument("Boundary leaf placement is outside the parent polygon.");
      }
      if (descriptor.boundary_edge) ++summary.free_variables;
    }
  }
  if (summary.boundary_leaves != summary.interior_nodes + 2) {
    throw std::invalid_argument("A degree-3 decomposition tree must have I+2 boundary leaves.");
  }
  summary.free_variables += 2 * summary.interior_nodes;
  summary.continuity_equations = 2 * summary.interior_nodes + summary.boundary_leaves;
  return summary;
}

}  // namespace kstab
