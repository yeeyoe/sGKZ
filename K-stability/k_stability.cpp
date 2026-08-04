#include "k_stability.hpp"

#include <CGAL/Gmpz.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace kstab {
namespace {

using WideInt = __int128_t;

Rational rat64(std::int64_t value) {
  return Rational(static_cast<long>(value));
}

Rational make_rational(std::int64_t numerator, std::int64_t denominator) {
  if (denominator <= 0) {
    throw std::invalid_argument("make_rational: denominator must be positive.");
  }
  return Rational(CGAL::Gmpz(static_cast<long>(numerator)),
                  CGAL::Gmpz(static_cast<long>(denominator)));
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

WideInt cross(const IntPoint& origin, const IntPoint& p, const IntPoint& q) {
  const WideInt px = static_cast<WideInt>(p.x) - origin.x;
  const WideInt py = static_cast<WideInt>(p.y) - origin.y;
  const WideInt qx = static_cast<WideInt>(q.x) - origin.x;
  const WideInt qy = static_cast<WideInt>(q.y) - origin.y;
  return px * qy - py * qx;
}

std::int64_t edge_lattice_length(const IntPoint& p, const IntPoint& q) {
  return std::gcd<std::int64_t>(p.x - q.x, p.y - q.y);
}

// 3x3 行列式，按第一行展开。
Rational det3(const std::array<std::array<Rational, 3>, 3>& m) {
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
         m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
         m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

// Cramer 法则解 3x3 系统；奇异时抛异常。
std::array<Rational, 3> solve3(
    std::array<std::array<Rational, 3>, 3> matrix,
    const std::array<Rational, 3>& rhs) {
  const Rational determinant = det3(matrix);
  if (determinant == 0) {
    throw std::runtime_error("The 3x3 moment system is singular.");
  }
  std::array<Rational, 3> solution;
  for (int column = 0; column < 3; ++column) {
    auto replaced = matrix;
    for (int row = 0; row < 3; ++row) {
      replaced[row][column] = rhs[row];
    }
    solution[column] = det3(replaced) / determinant;
  }
  return solution;
}

// 半平面 {a x + b y + c >= 0} 裁剪（精确）。
std::vector<QPoint> clip_halfplane_exact(const std::vector<QPoint>& polygon,
                                         const Rational& a, const Rational& b,
                                         const Rational& c) {
  std::vector<QPoint> output;
  const std::size_t n = polygon.size();
  for (std::size_t i = 0; i < n; ++i) {
    const QPoint& p = polygon[i];
    const QPoint& q = polygon[(i + 1) % n];
    const Rational sp = a * p.x + b * p.y + c;
    const Rational sq = a * q.x + b * q.y + c;
    const bool p_in = sp >= 0;
    const bool q_in = sq >= 0;
    if (p_in) {
      output.push_back(p);
    }
    if (p_in != q_in) {
      const Rational tau = sp / (sp - sq);
      output.push_back({p.x + tau * (q.x - p.x), p.y + tau * (q.y - p.y)});
    }
  }
  return output;
}

// 半平面 {a x + b y + c >= 0} 裁剪（double，带 snap 容差）。
std::vector<DPoint> clip_halfplane_double(const std::vector<DPoint>& polygon,
                                          double a, double b, double c) {
  const std::size_t n = polygon.size();
  std::vector<double> values(n);
  double scale = 1.0;
  for (std::size_t i = 0; i < n; ++i) {
    values[i] = a * polygon[i].x + b * polygon[i].y + c;
    scale = std::max(scale, std::fabs(values[i]));
  }
  const double tolerance = 1e-12 * scale;
  for (double& value : values) {
    if (std::fabs(value) <= tolerance) {
      value = 0.0;
    }
  }
  std::vector<DPoint> output;
  for (std::size_t i = 0; i < n; ++i) {
    const DPoint& p = polygon[i];
    const DPoint& q = polygon[(i + 1) % n];
    const double sp = values[i];
    const double sq = values[(i + 1) % n];
    const bool p_in = sp >= 0.0;
    const bool q_in = sq >= 0.0;
    if (p_in) {
      output.push_back(p);
    }
    if (p_in != q_in) {
      double tau = sp / (sp - sq);
      tau = std::clamp(tau, 0.0, 1.0);
      output.push_back({p.x + tau * (q.x - p.x), p.y + tau * (q.y - p.y)});
    }
  }
  return output;
}

// 平移后的双精度求值管线。顶点与 ell 都已按平移调整过。
struct PreparedPolygon {
  std::vector<DPoint> vertices;         // 平移后的顶点
  std::vector<std::int64_t> gcds;       // 每条边的格点长度（平移不变）
  std::array<double, 3> ell;            // 平移后的 ell_P
  double boundary_length = 0.0;         // |∂P|_{dσ}
};

double df_prepared(const PreparedPolygon& prepared, double a, double b,
                   double c, long* evaluation_counter) {
  if (evaluation_counter != nullptr) {
    ++(*evaluation_counter);
  }
  const std::vector<DPoint> clipped =
      clip_halfplane_double(prepared.vertices, a, b, c);
  const DMoments moments = polygon_moments_double(clipped);
  const auto& ell = prepared.ell;
  const double interior =
      a * ell[1] * moments.ixx + b * ell[2] * moments.iyy +
      (a * ell[2] + b * ell[1]) * moments.ixy +
      (a * ell[0] + c * ell[1]) * moments.ix +
      (b * ell[0] + c * ell[2]) * moments.iy + c * ell[0] * moments.area;

  double boundary = 0.0;
  const std::size_t n = prepared.vertices.size();
  double scale = 1.0;
  std::vector<double> values(n);
  for (std::size_t i = 0; i < n; ++i) {
    values[i] =
        a * prepared.vertices[i].x + b * prepared.vertices[i].y + c;
    scale = std::max(scale, std::fabs(values[i]));
  }
  const double tolerance = 1e-12 * scale;
  for (double& value : values) {
    if (std::fabs(value) <= tolerance) {
      value = 0.0;
    }
  }
  for (std::size_t i = 0; i < n; ++i) {
    const double sp = values[i];
    const double sq = values[(i + 1) % n];
    if (sp <= 0.0 && sq <= 0.0) {
      continue;
    }
    double tau0 = 0.0;
    double tau1 = 1.0;
    if (sp < 0.0) {
      tau0 = std::clamp(sp / (sp - sq), 0.0, 1.0);
    } else if (sq < 0.0) {
      tau1 = std::clamp(sp / (sp - sq), 0.0, 1.0);
    }
    const double s0 = sp + tau0 * (sq - sp);
    const double s1 = sp + tau1 * (sq - sp);
    boundary += static_cast<double>(prepared.gcds[i]) * (tau1 - tau0) *
                (s0 + s1) / 2.0;
  }
  return boundary - interior;
}

// golden-section 最小化；返回 (argmin, min)。
template <typename F>
std::pair<double, double> golden_minimize(F&& function, double lo, double hi,
                                          int iterations,
                                          long* evaluation_counter) {
  const double ratio = (std::sqrt(5.0) - 1.0) / 2.0;
  double c = hi - ratio * (hi - lo);
  double d = lo + ratio * (hi - lo);
  double fc = function(c);
  double fd = function(d);
  if (evaluation_counter != nullptr) {
    *evaluation_counter += 2;
  }
  for (int i = 0; i < iterations; ++i) {
    if (fc < fd) {
      hi = d;
      d = c;
      fd = fc;
      c = hi - ratio * (hi - lo);
      fc = function(c);
    } else {
      lo = c;
      c = d;
      fc = fd;
      d = lo + ratio * (hi - lo);
      fd = function(d);
    }
    if (evaluation_counter != nullptr) {
      ++(*evaluation_counter);
    }
  }
  if (fc < fd) {
    return {c, fc};
  }
  return {d, fd};
}

PreparedPolygon prepare_polygon(const std::vector<IntPoint>& vertices,
                                const std::array<double, 3>& ell,
                                double shift_x, double shift_y) {
  PreparedPolygon prepared;
  prepared.vertices.reserve(vertices.size());
  for (const IntPoint& vertex : vertices) {
    prepared.vertices.push_back({static_cast<double>(vertex.x) - shift_x,
                                 static_cast<double>(vertex.y) - shift_y});
  }
  prepared.gcds.reserve(vertices.size());
  double boundary = 0.0;
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    const std::int64_t g =
        edge_lattice_length(vertices[i], vertices[(i + 1) % vertices.size()]);
    prepared.gcds.push_back(g);
    boundary += static_cast<double>(g);
  }
  prepared.ell = {ell[0] + ell[1] * shift_x + ell[2] * shift_y, ell[1],
                  ell[2]};
  prepared.boundary_length = boundary;
  return prepared;
}

}  // namespace

std::vector<IntPoint> normalize_polygon(std::vector<IntPoint> vertices) {
  if (vertices.size() >= 2 && vertices.front().x == vertices.back().x &&
      vertices.front().y == vertices.back().y) {
    vertices.pop_back();  // 接受首尾重复的闭合写法
  }
  if (vertices.size() < 3) {
    throw std::runtime_error("A polygon needs at least 3 distinct vertices.");
  }
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    const IntPoint& p = vertices[i];
    const IntPoint& q = vertices[(i + 1) % vertices.size()];
    if (p.x == q.x && p.y == q.y) {
      throw std::runtime_error("Consecutive duplicate polygon vertices.");
    }
  }

  // 合并共线连续顶点；拒绝回溯与非凸。
  bool merged = true;
  while (merged) {
    merged = false;
    const std::size_t n = vertices.size();
    for (std::size_t i = 0; i < n; ++i) {
      const IntPoint& p = vertices[i];
      const IntPoint& q = vertices[(i + 1) % n];
      const IntPoint& r = vertices[(i + 2) % n];
      if (cross(p, q, r) != 0) {
        continue;
      }
      const WideInt dot =
          (static_cast<WideInt>(q.x) - p.x) *
              (static_cast<WideInt>(r.x) - q.x) +
          (static_cast<WideInt>(q.y) - p.y) *
              (static_cast<WideInt>(r.y) - q.y);
      if (dot < 0) {
        throw std::runtime_error("Polygon boundary backtracks at a vertex.");
      }
      vertices.erase(vertices.begin() +
                     static_cast<std::ptrdiff_t>((i + 1) % n));
      merged = true;
      break;
    }
  }
  if (vertices.size() < 3) {
    throw std::runtime_error("Polygon degenerates to fewer than 3 vertices.");
  }

  int sign = 0;
  WideInt twice_area = 0;
  const std::size_t n = vertices.size();
  for (std::size_t i = 0; i < n; ++i) {
    const IntPoint& p = vertices[i];
    const IntPoint& q = vertices[(i + 1) % n];
    const IntPoint& r = vertices[(i + 2) % n];
    const WideInt turn = cross(p, q, r);
    if (turn == 0) {
      throw std::runtime_error("Collinear vertex survived merging.");
    }
    const int turn_sign = turn > 0 ? 1 : -1;
    if (sign == 0) {
      sign = turn_sign;
    } else if (turn_sign != sign) {
      throw std::runtime_error("Polygon is not convex.");
    }
    twice_area += static_cast<WideInt>(p.x) * q.y -
                  static_cast<WideInt>(p.y) * q.x;
  }
  if (twice_area == 0) {
    throw std::runtime_error("Polygon has zero area.");
  }
  if (sign < 0) {
    std::reverse(vertices.begin(), vertices.end());
  }
  return vertices;
}

std::vector<IntPoint> parse_polygon_file(const std::filesystem::path& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("Cannot open polygon file: " + path.string());
  }
  std::vector<IntPoint> vertices;
  std::string line;
  std::int64_t line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    if (const auto comment = line.find('#'); comment != std::string::npos) {
      line.erase(comment);
    }
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream parser(line);
    long long x = 0;
    long long y = 0;
    if (!(parser >> x)) {
      continue;
    }
    if (!(parser >> y)) {
      throw std::runtime_error("Missing y coordinate at " + path.string() +
                               ":" + std::to_string(line_number));
    }
    std::string extra;
    if (parser >> extra) {
      throw std::runtime_error("Unexpected third field at " + path.string() +
                               ":" + std::to_string(line_number));
    }
    vertices.push_back({static_cast<std::int64_t>(x),
                        static_cast<std::int64_t>(y)});
  }
  return normalize_polygon(std::move(vertices));
}

Moments polygon_moments(const std::vector<QPoint>& vertices_ccw) {
  Moments moments;
  const std::size_t n = vertices_ccw.size();
  for (std::size_t i = 0; i < n; ++i) {
    const QPoint& u = vertices_ccw[i];
    const QPoint& v = vertices_ccw[(i + 1) % n];
    const Rational delta = u.x * v.y - u.y * v.x;
    moments.area += delta / 2;
    moments.ix += delta * (u.x + v.x) / 6;
    moments.iy += delta * (u.y + v.y) / 6;
    moments.ixx += delta * (u.x * u.x + u.x * v.x + v.x * v.x) / 12;
    moments.iyy += delta * (u.y * u.y + u.y * v.y + v.y * v.y) / 12;
    moments.ixy += delta *
                   (2 * u.x * u.y + u.x * v.y + v.x * u.y + 2 * v.x * v.y) / 24;
  }
  return moments;
}

DMoments polygon_moments_double(const std::vector<DPoint>& vertices_ccw) {
  DMoments moments;
  const std::size_t n = vertices_ccw.size();
  for (std::size_t i = 0; i < n; ++i) {
    const DPoint& u = vertices_ccw[i];
    const DPoint& v = vertices_ccw[(i + 1) % n];
    const double delta = u.x * v.y - u.y * v.x;
    moments.area += delta / 2.0;
    moments.ix += delta * (u.x + v.x) / 6.0;
    moments.iy += delta * (u.y + v.y) / 6.0;
    moments.ixx += delta * (u.x * u.x + u.x * v.x + v.x * v.x) / 12.0;
    moments.iyy += delta * (u.y * u.y + u.y * v.y + v.y * v.y) / 12.0;
    moments.ixy += delta *
                   (2.0 * u.x * u.y + u.x * v.y + v.x * u.y +
                    2.0 * v.x * v.y) /
                   24.0;
  }
  return moments;
}

BoundaryMoments boundary_moments(
    const std::vector<IntPoint>& vertices_ccw) {
  BoundaryMoments moments;
  const std::size_t n = vertices_ccw.size();
  for (std::size_t i = 0; i < n; ++i) {
    const IntPoint& p = vertices_ccw[i];
    const IntPoint& q = vertices_ccw[(i + 1) % n];
    const Rational g = rat64(edge_lattice_length(p, q));
    moments.length += g;
    moments.ix += g * (rat64(p.x) + rat64(q.x)) / 2;
    moments.iy += g * (rat64(p.y) + rat64(q.y)) / 2;
  }
  return moments;
}

std::array<Rational, 3> compute_ell_p(
    const std::vector<IntPoint>& vertices_ccw) {
  std::vector<QPoint> rational_vertices;
  rational_vertices.reserve(vertices_ccw.size());
  for (const IntPoint& vertex : vertices_ccw) {
    rational_vertices.push_back({rat64(vertex.x), rat64(vertex.y)});
  }
  const Moments moments = polygon_moments(rational_vertices);
  const BoundaryMoments boundary = boundary_moments(vertices_ccw);
  return solve3({{{moments.area, moments.ix, moments.iy},
                  {moments.ix, moments.ixx, moments.ixy},
                  {moments.iy, moments.ixy, moments.iyy}}},
                {boundary.length, boundary.ix, boundary.iy});
}

Rational df_simple_exact(const std::vector<IntPoint>& vertices_ccw,
                         const std::array<Rational, 3>& ell,
                         const Rational& a, const Rational& b,
                         const Rational& c) {
  std::vector<QPoint> polygon;
  polygon.reserve(vertices_ccw.size());
  for (const IntPoint& vertex : vertices_ccw) {
    polygon.push_back({rat64(vertex.x), rat64(vertex.y)});
  }
  const std::vector<QPoint> clipped = clip_halfplane_exact(polygon, a, b, c);
  const Moments moments = polygon_moments(clipped);
  const Rational interior =
      a * ell[1] * moments.ixx + b * ell[2] * moments.iyy +
      (a * ell[2] + b * ell[1]) * moments.ixy +
      (a * ell[0] + c * ell[1]) * moments.ix +
      (b * ell[0] + c * ell[2]) * moments.iy + c * ell[0] * moments.area;

  Rational boundary = 0;
  const std::size_t n = vertices_ccw.size();
  for (std::size_t i = 0; i < n; ++i) {
    const IntPoint& p = vertices_ccw[i];
    const IntPoint& q = vertices_ccw[(i + 1) % n];
    const Rational sp = a * rat64(p.x) + b * rat64(p.y) + c;
    const Rational sq = a * rat64(q.x) + b * rat64(q.y) + c;
    if (sp <= 0 && sq <= 0) {
      continue;
    }
    Rational tau0 = 0;
    Rational tau1 = 1;
    if (sp < 0) {
      tau0 = sp / (sp - sq);
    } else if (sq < 0) {
      tau1 = sp / (sp - sq);
    }
    const Rational s0 = sp + tau0 * (sq - sp);
    const Rational s1 = sp + tau1 * (sq - sp);
    boundary += rat64(edge_lattice_length(p, q)) * (tau1 - tau0) *
                (s0 + s1) / 2;
  }
  return boundary - interior;
}

double df_simple_double(const std::vector<IntPoint>& vertices_ccw,
                        const std::array<double, 3>& ell, double a, double b,
                        double c) {
  double shift_x = 0.0;
  double shift_y = 0.0;
  if (!vertices_ccw.empty()) {
    std::int64_t min_x = vertices_ccw.front().x;
    std::int64_t min_y = vertices_ccw.front().y;
    for (const IntPoint& vertex : vertices_ccw) {
      min_x = std::min(min_x, vertex.x);
      min_y = std::min(min_y, vertex.y);
    }
    shift_x = static_cast<double>(min_x);
    shift_y = static_cast<double>(min_y);
  }
  const PreparedPolygon prepared =
      prepare_polygon(vertices_ccw, ell, shift_x, shift_y);
  return df_prepared(prepared, a, b, c + a * shift_x + b * shift_y, nullptr);
}

Rational approximate_rational(double value, std::int64_t cap) {
  if (cap < 1) {
    throw std::invalid_argument("approximate_rational: cap must be >= 1.");
  }
  if (!std::isfinite(value) || std::fabs(value) > 1e15) {
    throw std::invalid_argument("approximate_rational: value out of range.");
  }
  const bool negative = value < 0.0;
  const double x = std::fabs(value);
  std::int64_t hm2 = 0, hm1 = 1, km2 = 1, km1 = 0;
  double y = x;
  for (int iteration = 0; iteration < 64; ++iteration) {
    const double ay = std::floor(y);
    if (ay > 4.0e18) {
      break;
    }
    const std::int64_t a = static_cast<std::int64_t>(ay);
    const std::int64_t h = a * hm1 + hm2;
    const std::int64_t k = a * km1 + km2;
    if (k > cap) {
      // 半收敛子：t = (cap - km2) / km1
      const std::int64_t t = km1 > 0 ? (cap - km2) / km1 : 0;
      if (t >= 1) {
        const std::int64_t hs = t * hm1 + hm2;
        const std::int64_t ks = t * km1 + km2;
        const double error_semi =
            std::fabs(x - static_cast<double>(hs) / static_cast<double>(ks));
        const double error_prev =
            std::fabs(x - static_cast<double>(hm1) / static_cast<double>(km1));
        const Rational result = error_semi <= error_prev
                                    ? make_rational(hs, ks)
                                    : make_rational(hm1, km1);
        return negative ? -result : result;
      }
      const Rational result = make_rational(hm1, km1);
      return negative ? -result : result;
    }
    hm2 = hm1;
    hm1 = h;
    km2 = km1;
    km1 = k;
    if (y == ay) {
      const Rational result = make_rational(hm1, km1);
      return negative ? -result : result;
    }
    y = 1.0 / (y - ay);
  }
  const Rational result = make_rational(hm1, km1);
  return negative ? -result : result;
}

SearchResult search_witness(const std::vector<IntPoint>& vertices_ccw,
                            const std::array<Rational, 3>& ell,
                            const SearchOptions& options) {
  double min_x = static_cast<double>(vertices_ccw.front().x);
  double min_y = static_cast<double>(vertices_ccw.front().y);
  for (const IntPoint& vertex : vertices_ccw) {
    min_x = std::min(min_x, static_cast<double>(vertex.x));
    min_y = std::min(min_y, static_cast<double>(vertex.y));
  }
  const auto ell_d = ell_to_double(ell);
  const PreparedPolygon prepared =
      prepare_polygon(vertices_ccw, ell_d, min_x, min_y);

  // 归一化除数：|∂P| · sup|ell_P| · diam。
  double sup_ell = 0.0;
  double diameter = 0.0;
  for (const IntPoint& vertex : vertices_ccw) {
    sup_ell = std::max(
        sup_ell, std::fabs(ell_d[0] + ell_d[1] * vertex.x +
                           ell_d[2] * vertex.y));
    for (const IntPoint& other : vertices_ccw) {
      diameter = std::max(
          diameter, std::hypot(static_cast<double>(vertex.x - other.x),
                               static_cast<double>(vertex.y - other.y)));
    }
  }
  const double normalization =
      prepared.boundary_length * std::max(sup_ell, 1e-300) * diameter;

  // 方向集合：均匀网格 + 边法向 + 顶点对法向。
  constexpr double kTwoPi = 6.283185307179586476925286766559;
  std::vector<double> thetas;
  thetas.reserve(static_cast<std::size_t>(options.theta_steps) + 8);
  for (int j = 0; j < options.theta_steps; ++j) {
    thetas.push_back(kTwoPi * static_cast<double>(j) /
                     static_cast<double>(options.theta_steps));
  }
  const auto add_normal = [&](double dx, double dy) {
    const double phi = std::atan2(dy, dx);
    for (const double offset : {phi + kTwoPi / 4, phi - kTwoPi / 4}) {
      double angle = std::fmod(offset, kTwoPi);
      if (angle < 0.0) {
        angle += kTwoPi;
      }
      thetas.push_back(angle);
    }
  };
  const std::size_t n = vertices_ccw.size();
  for (std::size_t i = 0; i < n; ++i) {
    const IntPoint& p = vertices_ccw[i];
    const IntPoint& q = vertices_ccw[(i + 1) % n];
    add_normal(static_cast<double>(q.x - p.x),
               static_cast<double>(q.y - p.y));
  }
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 2; j < n; ++j) {
      if (i == 0 && j == n - 1) {
        continue;  // 相邻边已覆盖
      }
      add_normal(
          static_cast<double>(vertices_ccw[j].x - vertices_ccw[i].x),
          static_cast<double>(vertices_ccw[j].y - vertices_ccw[i].y));
    }
  }

  long evaluations = 0;
  struct Candidate {
    double value;
    double theta;
    double t;
  };
  std::vector<Candidate> best_candidates;
  const auto consider = [&](double value, double theta, double t) {
    best_candidates.push_back({value, theta, t});
    std::push_heap(best_candidates.begin(), best_candidates.end(),
                   [](const Candidate& lhs, const Candidate& rhs) {
                     return lhs.value > rhs.value;
                   });
    if (best_candidates.size() > 5) {
      std::pop_heap(best_candidates.begin(), best_candidates.end(),
                    [](const Candidate& lhs, const Candidate& rhs) {
                      return lhs.value > rhs.value;
                    });
      best_candidates.pop_back();
    }
  };

  for (const double theta : thetas) {
    const double ux = std::cos(theta);
    const double uy = std::sin(theta);
    double w_min = 0.0;
    double w_max = 0.0;
    std::vector<double> projections(n);
    for (std::size_t i = 0; i < n; ++i) {
      projections[i] =
          ux * prepared.vertices[i].x + uy * prepared.vertices[i].y;
    }
    const auto [min_it, max_it] =
        std::minmax_element(projections.begin(), projections.end());
    w_min = *min_it;
    w_max = *max_it;
    if (w_max <= w_min) {
      continue;
    }
    const double width = w_max - w_min;

    const auto evaluate = [&](double t) {
      return df_prepared(prepared, ux, uy, -t, &evaluations);
    };

    double best_t = w_min + width / 2.0;
    double best_value = evaluate(best_t);
    // 均匀网格
    for (int j = 0; j < options.t_steps; ++j) {
      const double t = w_min + width * (static_cast<double>(j) + 0.5) /
                                   static_cast<double>(options.t_steps);
      const double value = evaluate(t);
      if (value < best_value) {
        best_value = value;
        best_t = t;
      }
    }
    // 顶点投影 breakpoint（kink）
    std::sort(projections.begin(), projections.end());
    for (const double w : projections) {
      if (w <= w_min || w >= w_max) {
        continue;
      }
      const double value = evaluate(w);
      if (value < best_value) {
        best_value = value;
        best_t = w;
      }
    }
    // 当前方向的局部细化（括号必须夹在 [w_min, w_max] 内：
    // 区间外 g 在多边形上是仿射函数，对非真实 ell 没有下界）
    if (options.refine) {
      const double bracket = width / static_cast<double>(options.t_steps);
      const double refine_lo = std::max(best_t - bracket, w_min);
      const double refine_hi = std::min(best_t + bracket, w_max);
      if (refine_hi - refine_lo > 1e-15) {
        const auto [refined_t, refined_value] = golden_minimize(
            evaluate, refine_lo, refine_hi, 40, &evaluations);
        if (refined_value < best_value) {
          best_value = refined_value;
          best_t = refined_t;
        }
      }
    }
    if (options.verbose) {
      std::cerr << "theta=" << theta << " t=" << best_t
                << " M=" << best_value << '\n';
    }
    consider(best_value, theta, best_t);
  }

  // 交替 (theta, t) 细化最优候选
  Candidate best{0.0, 0.0, 0.0};
  bool have_candidate = false;
  for (const Candidate& candidate : best_candidates) {
    double theta = candidate.theta;
    double t = candidate.t;
    double value = candidate.value;
    if (options.refine) {
      double width_theta =
          4.0 * kTwoPi / static_cast<double>(options.theta_steps);
      for (int round = 0; round < 3; ++round) {
        const auto t_function = [&](double tt) {
          return df_prepared(prepared, std::cos(theta), std::sin(theta), -tt,
                             &evaluations);
        };
        // 重新计算该方向的投影范围作为 t 括号
        const double ux = std::cos(theta);
        const double uy = std::sin(theta);
        double w_lo = ux * prepared.vertices[0].x + uy * prepared.vertices[0].y;
        double w_hi = w_lo;
        for (const DPoint& v : prepared.vertices) {
          const double w = ux * v.x + uy * v.y;
          w_lo = std::min(w_lo, w);
          w_hi = std::max(w_hi, w);
        }
        const double w_width =
            (w_hi - w_lo) /
            static_cast<double>(std::max(options.t_steps, 1)) * 4.0;
        const double refine_lo = std::max(t - w_width, w_lo);
        const double refine_hi = std::min(t + w_width, w_hi);
        if (refine_hi - refine_lo > 1e-15) {
          const auto [new_t, value_t] = golden_minimize(
              t_function, refine_lo, refine_hi, 40, &evaluations);
          if (value_t < value) {
            t = new_t;
            value = value_t;
          }
        }
        const auto theta_function = [&](double th) {
          return df_prepared(prepared, std::cos(th), std::sin(th), -t,
                             &evaluations);
        };
        const auto [new_theta, value_theta] =
            golden_minimize(theta_function, theta - width_theta,
                            theta + width_theta, 40, &evaluations);
        if (value_theta < value) {
          theta = new_theta;
          value = value_theta;
        }
        width_theta /= 8.0;
      }
    }
    if (!have_candidate || value < best.value) {
      best = {value, theta, t};
      have_candidate = true;
    }
  }

  SearchResult result;
  result.evaluations = evaluations;
  result.normalization = normalization;
  if (have_candidate) {
    const double ux = std::cos(best.theta);
    const double uy = std::sin(best.theta);
    result.witness.ux = ux;
    result.witness.uy = uy;
    // 平移回原坐标系：t_orig = t' + u·shift
    result.witness.t = best.t + ux * min_x + uy * min_y;
    result.witness.value = best.value;
    result.witness.normalized =
        normalization > 0.0 ? best.value / normalization : best.value;
    result.unstable = result.witness.normalized < -1e-6;
  }
  return result;
}

CertifyResult certify_witness(const std::vector<IntPoint>& vertices_ccw,
                              const std::array<Rational, 3>& ell,
                              const Witness& witness,
                              std::int64_t max_denominator) {
  CertifyResult result;
  const double target_a = witness.ux;
  const double target_b = witness.uy;
  const double target_c = -witness.t;
  for (std::int64_t cap = 10; cap <= max_denominator; cap *= 10) {
    const Rational a = approximate_rational(target_a, cap);
    const Rational b = approximate_rational(target_b, cap);
    const Rational c = approximate_rational(target_c, cap);
    const Rational value =
        df_simple_exact(vertices_ccw, ell, a, b, c);
    if (value < 0) {
      result.certified = true;
      result.coefficients = {a, b, c};
      result.value = value;
      return result;
    }
  }
  return result;
}

bool write_svg(const std::filesystem::path& path,
               const std::vector<IntPoint>& vertices_ccw,
               const std::array<Rational, 3>& ell, const Witness* witness) {
  std::ofstream out(path);
  if (!out) {
    return false;
  }
  double min_x = static_cast<double>(vertices_ccw.front().x);
  double max_x = min_x;
  double min_y = static_cast<double>(vertices_ccw.front().y);
  double max_y = min_y;
  for (const IntPoint& vertex : vertices_ccw) {
    min_x = std::min(min_x, static_cast<double>(vertex.x));
    max_x = std::max(max_x, static_cast<double>(vertex.x));
    min_y = std::min(min_y, static_cast<double>(vertex.y));
    max_y = std::max(max_y, static_cast<double>(vertex.y));
  }
  constexpr double kSize = 640.0;
  constexpr double kPad = 48.0;
  const double range_x = std::max(max_x - min_x, 1e-9);
  const double range_y = std::max(max_y - min_y, 1e-9);
  const double scale =
      std::min((kSize - 2 * kPad) / range_x, (kSize - 2 * kPad) / range_y);
  const double tx = kPad - scale * min_x;
  const double ty = kPad + scale * max_y;

  out << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 "
      << static_cast<int>(kSize) << ' ' << static_cast<int>(kSize)
      << "\" width=\"" << static_cast<int>(kSize) << "\" height=\""
      << static_cast<int>(kSize) << "\">\n";
  out << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
  out << "<g transform=\"matrix(" << scale << " 0 0 " << -scale << ' ' << tx
      << ' ' << ty << ")\" fill=\"none\" vector-effect=\"non-scaling-stroke\">\n";

  out << "<polygon points=\"";
  for (const IntPoint& vertex : vertices_ccw) {
    out << vertex.x << ',' << vertex.y << ' ';
  }
  out << "\" fill=\"#eff6ff\" stroke=\"#1f2937\" stroke-width=\"1.5\"/>\n";

  if (witness != nullptr) {
    const double ux = witness->ux;
    const double uy = witness->uy;
    const double t = witness->t;
    // P ∩ {s ≥ 0} 阴影
    std::vector<DPoint> polygon_d;
    polygon_d.reserve(vertices_ccw.size());
    for (const IntPoint& vertex : vertices_ccw) {
      polygon_d.push_back({static_cast<double>(vertex.x),
                           static_cast<double>(vertex.y)});
    }
    const std::vector<DPoint> wedge =
        clip_halfplane_double(polygon_d, ux, uy, -t);
    if (wedge.size() >= 3) {
      out << "<polygon points=\"";
      for (const DPoint& vertex : wedge) {
        out << vertex.x << ',' << vertex.y << ' ';
      }
      out << "\" fill=\"#fca5a5\" fill-opacity=\"0.4\" stroke=\"none\"/>\n";
    }
    // 折痕弦：与每条边求交（含折痕恰过顶点的退化情形），
    // 取两个端点并外延 10%
    std::vector<DPoint> crossings;
    const std::size_t n = polygon_d.size();
    double coord_scale = 1.0;
    for (const DPoint& vertex : polygon_d) {
      coord_scale = std::max(
          coord_scale, std::max(std::fabs(vertex.x), std::fabs(vertex.y)));
    }
    const double on_line_tolerance = 1e-9 * coord_scale;
    for (std::size_t i = 0; i < n; ++i) {
      const DPoint& p = polygon_d[i];
      const DPoint& q = polygon_d[(i + 1) % n];
      const double sp = ux * p.x + uy * p.y - t;
      const double sq = ux * q.x + uy * q.y - t;
      if ((sp < 0.0 && sq > 0.0) || (sp > 0.0 && sq < 0.0)) {
        const double tau = sp / (sp - sq);
        crossings.push_back(
            {p.x + tau * (q.x - p.x), p.y + tau * (q.y - p.y)});
      } else if (std::fabs(sp) <= on_line_tolerance) {
        crossings.push_back(p);
      }
    }
    // 按坐标排序后去重（折痕过顶点时相邻边会各贡献一次）
    std::sort(crossings.begin(), crossings.end(),
              [](const DPoint& lhs, const DPoint& rhs) {
                if (lhs.x != rhs.x) {
                  return lhs.x < rhs.x;
                }
                return lhs.y < rhs.y;
              });
    crossings.erase(
        std::unique(crossings.begin(), crossings.end(),
                    [&](const DPoint& lhs, const DPoint& rhs) {
                      return std::hypot(lhs.x - rhs.x, lhs.y - rhs.y) <=
                             on_line_tolerance;
                    }),
        crossings.end());
    if (crossings.size() >= 2) {
      DPoint a = crossings.front();
      DPoint b = crossings.front();
      double best = -1.0;
      for (const DPoint& c1 : crossings) {
        for (const DPoint& c2 : crossings) {
          const double distance = std::hypot(c1.x - c2.x, c1.y - c2.y);
          if (distance > best) {
            best = distance;
            a = c1;
            b = c2;
          }
        }
      }
      const double extend = 0.1;
      const double mx = (a.x + b.x) / 2.0;
      const double my = (a.y + b.y) / 2.0;
      const double hx = (b.x - a.x) * (0.5 + extend);
      const double hy = (b.y - a.y) * (0.5 + extend);
      out << "<line x1=\"" << mx - hx << "\" y1=\"" << my - hy << "\" x2=\""
          << mx + hx << "\" y2=\"" << my + hy
          << "\" stroke=\"#dc2626\" stroke-width=\"2\" "
             "stroke-dasharray=\"6 3\"/>\n";
    }
  }
  out << "</g>\n";

  // 文本放在翻转组外，使用屏幕坐标。
  out << "<text x=\"16\" y=\"24\" font-family=\"monospace\" font-size=\"14\" "
         "fill=\"#111827\">ell_P(x,y) = "
      << rational_string(ell[0]) << " + (" << rational_string(ell[1])
      << ") x + (" << rational_string(ell[2]) << ") y</text>\n";
  if (witness != nullptr) {
    out << "<text x=\"16\" y=\"46\" font-family=\"monospace\" "
           "font-size=\"14\" fill=\"#b91c1c\">M_l = "
        << witness->value << "  (relatively K-unstable)</text>\n";
  }
  out << "</svg>\n";
  return static_cast<bool>(out);
}

std::array<double, 3> ell_to_double(const std::array<Rational, 3>& ell) {
  return {rational_to_double(ell[0]), rational_to_double(ell[1]),
          rational_to_double(ell[2])};
}

double rational_to_double(const Rational& value) {
  return CGAL::to_double(value);
}

}  // namespace kstab
