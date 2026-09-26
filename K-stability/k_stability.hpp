#pragma once

// 连续 Donaldson 函数 ell_P 的精确计算，以及基于 Donaldson 事实的
// 相对 K-不稳定性检测。数学定义见 paper/K-stability.tex。
//
// 本模块完全自包含：只依赖 CGAL::Gmpq，不链接 gkz_core。

#include <CGAL/Gmpq.h>

#include <array>
#include <cstdint>
#include <filesystem>
#include <vector>

namespace kstab {

using Rational = CGAL::Gmpq;

struct IntPoint {
  std::int64_t x = 0;
  std::int64_t y = 0;
};

// 规范化后的多边形及其边界测度。null_measure_edges[i] 对应
// vertices_ccw[i] -> vertices_ccw[(i + 1) % n]；true 表示该边 dσ = 0。
// 空掩码表示所有边都使用标准格点测度。
struct PolygonInput {
  std::vector<IntPoint> vertices_ccw;
  std::vector<bool> null_measure_edges;
};

// 有理坐标点（裁剪后的多边形顶点）。
struct QPoint {
  Rational x;
  Rational y;
};

// double 坐标点（数值搜索路径）。
struct DPoint {
  double x = 0.0;
  double y = 0.0;
};

// 多边形 P 上的矩：对 1, x, y, x^2, xy, y^2 的面积分。
struct Moments {
  Rational area = 0;
  Rational ix = 0;
  Rational iy = 0;
  Rational ixx = 0;
  Rational ixy = 0;
  Rational iyy = 0;
};

struct DMoments {
  double area = 0.0;
  double ix = 0.0;
  double iy = 0.0;
  double ixx = 0.0;
  double ixy = 0.0;
  double iyy = 0.0;
};

// 边界格点测度 dσ 下的矩：length = |∂P|_{dσ}，ix = ∫_{∂P} x dσ 等。
struct BoundaryMoments {
  Rational length = 0;
  Rational ix = 0;
  Rational iy = 0;
};

// 解析多边形文件。顶点部分每行 "x y"；可选的单独一行
// "null measure edges" 之后，每行 "x1 y1 x2 y2" 指定一条 dσ = 0
// 的无向边。'#' 之后为注释，逗号视为空白。零测度边必须匹配规范化后的
// 真实边；重复或不存在的边报错。
PolygonInput parse_polygon_input_file(const std::filesystem::path& path);

// 只返回规范化顶点的兼容接口；忽略文件中的零测度边信息。
std::vector<IntPoint> parse_polygon_file(const std::filesystem::path& path);

// 对内存中的顶点做同样的校验与 CCW 归一化。
std::vector<IntPoint> normalize_polygon(std::vector<IntPoint> vertices);

// CCW 顶点的带符号扇形矩（从原点）。两版结果数学上相同。
Moments polygon_moments(const std::vector<QPoint>& vertices_ccw);
DMoments polygon_moments_double(const std::vector<DPoint>& vertices_ccw);

// 边界格点测度矩。边 p->q 的格点长度 g = gcd(|dx|,|dy|)，
// 仿射 h 在该边上 ∫ h dσ = g (h(p)+h(q))/2。
BoundaryMoments boundary_moments(
    const std::vector<IntPoint>& vertices_ccw,
    const std::vector<bool>& null_measure_edges = {});

// ell_P(x,y) = coefficients[0] + coefficients[1] x + coefficients[2] y，
// 精确满足：对所有仿射 h，∫_{∂P} h dσ = ∫_P h ell_P dx。
std::array<Rational, 3> compute_ell_p(
    const std::vector<IntPoint>& vertices_ccw,
    const std::vector<bool>& null_measure_edges = {});

// M_ℓ(max{a x + b y + c, 0}) 的精确值。
Rational df_simple_exact(const std::vector<IntPoint>& vertices_ccw,
                         const std::array<Rational, 3>& ell,
                         const Rational& a, const Rational& b,
                         const Rational& c,
                         const std::vector<bool>& null_measure_edges = {});

// M_ℓ(max{a x + b y + c, 0}) 的 double 值（搜索用）。
double df_simple_double(const std::vector<IntPoint>& vertices_ccw,
                        const std::array<double, 3>& ell, double a, double b,
                        double c,
                        const std::vector<bool>& null_measure_edges = {});

struct SearchOptions {
  int theta_steps = 720;   // [0, 2π) 上均匀方向数
  int t_steps = 512;       // 每个方向上 t 的均匀采样数
  bool refine = true;      // 是否对最优候选做局部细化
  bool verbose = false;    // 逐方向输出最小值到 stderr
};

// 简单凸函数 witness：g = max{<x,u> - t, 0}，u 为单位法向。
struct Witness {
  double ux = 0.0;
  double uy = 0.0;
  double t = 0.0;
  double value = 0.0;       // M_ℓ(g)（double）
  double normalized = 0.0;  // value / (|∂P| · sup|ℓ_P| · diam)
};

struct SearchResult {
  bool unstable = false;      // normalized < -1e-6
  Witness witness;            // 最优（最负）候选
  double normalization = 1.0; // 归一化除数
  long evaluations = 0;       // df_simple_double 调用次数
};

SearchResult search_witness(const std::vector<IntPoint>& vertices_ccw,
                            const std::array<Rational, 3>& ell,
                            const SearchOptions& options,
                            const std::vector<bool>& null_measure_edges = {});

struct CertifyResult {
  bool certified = false;
  std::array<Rational, 3> coefficients{};  // 认证 witness 的 a, b, c
  Rational value = 0;                      // 精确 M_ℓ < 0
};

// 把数值 witness 的有理化逼近逐档（分母上限 10, 100, ..., max_denominator）
// 送入精确路径，首个 M_ℓ < 0 者即为认证 witness。
CertifyResult certify_witness(const std::vector<IntPoint>& vertices_ccw,
                              const std::array<Rational, 3>& ell,
                              const Witness& witness,
                              std::int64_t max_denominator,
                              const std::vector<bool>& null_measure_edges = {});

// 连分数有理逼近：返回分母不超过 cap 的最佳逼近。
Rational approximate_rational(double value, std::int64_t cap);

// 写出多边形 + witness 折痕线（及 P∩{s≥0} 阴影）的 SVG。
// witness 为空指针时只画多边形。成功返回 true。
bool write_svg(const std::filesystem::path& path,
               const std::vector<IntPoint>& vertices_ccw,
               const std::array<Rational, 3>& ell,
               const Witness* witness);

std::array<double, 3> ell_to_double(const std::array<Rational, 3>& ell);
double rational_to_double(const Rational& value);

}  // namespace kstab
