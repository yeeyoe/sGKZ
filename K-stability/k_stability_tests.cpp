// K-stability 模块的无框架回归测试。
// 手算基准（均已独立验证）：
//   单位正方形      ell_P ≡ 4，M_l(max{x-1/2,0}) = 1/4
//   三角形          (0,0),(1,0),(0,1)  ell_P ≡ 6
//   梯形            (0,0),(2,0),(1,1),(0,1)  ell_P = (54-24y)/13

#include "k_stability.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using kstab::IntPoint;
using kstab::Rational;

void require(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void require_close(double actual, double expected, double tolerance,
                   const std::string& message) {
  if (std::abs(actual - expected) > tolerance) {
    throw std::runtime_error(message + ": actual=" + std::to_string(actual) +
                             ", expected=" + std::to_string(expected));
  }
}

Rational make_rational(std::int64_t numerator, std::int64_t denominator = 1) {
  return Rational(CGAL::Gmpz(static_cast<long>(numerator)),
                  CGAL::Gmpz(static_cast<long>(denominator)));
}

std::vector<IntPoint> unit_square() {
  return {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
}

std::vector<IntPoint> right_triangle() {
  return {{0, 0}, {1, 0}, {0, 1}};
}

std::vector<IntPoint> trapezoid() {
  return {{0, 0}, {2, 0}, {1, 1}, {0, 1}};
}

void test_ell_p_exact() {
  // 单位正方形：ell_P ≡ 4
  {
    const auto ell = kstab::compute_ell_p(unit_square());
    require(ell[0] == 4 && ell[1] == 0 && ell[2] == 0,
            "unit square ell_P should be constant 4");
  }
  // 三角形：ell_P ≡ 6
  {
    const auto ell = kstab::compute_ell_p(right_triangle());
    require(ell[0] == 6 && ell[1] == 0 && ell[2] == 0,
            "right triangle ell_P should be constant 6");
  }
  // 梯形：ell_P = (54 - 24 y)/13
  {
    const auto ell = kstab::compute_ell_p(trapezoid());
    require(ell[0] == make_rational(54, 13) && ell[1] == 0 &&
                ell[2] == make_rational(-24, 13),
            "trapezoid ell_P should be (54 - 24 y)/13");
  }
}

void test_polygon_moments() {
  // 标准单形：∫1 = 1/2, ∫x = ∫y = 1/6, ∫x² = ∫y² = 1/12, ∫xy = 1/24
  const std::vector<kstab::QPoint> simplex = {
      {Rational(0), Rational(0)},
      {Rational(1), Rational(0)},
      {Rational(0), Rational(1)}};
  const kstab::Moments moments = kstab::polygon_moments(simplex);
  require(moments.area == make_rational(1, 2), "simplex area");
  require(moments.ix == make_rational(1, 6), "simplex ∫x");
  require(moments.iy == make_rational(1, 6), "simplex ∫y");
  require(moments.ixx == make_rational(1, 12), "simplex ∫x²");
  require(moments.iyy == make_rational(1, 12), "simplex ∫y²");
  require(moments.ixy == make_rational(1, 24), "simplex ∫xy");

  // 定向反转变号
  const std::vector<kstab::QPoint> reversed = {simplex[0], simplex[2],
                                               simplex[1]};
  const kstab::Moments negative = kstab::polygon_moments(reversed);
  require(negative.area == make_rational(-1, 2),
          "reversed simplex area should flip sign");

  // 平移三角形（原点在外部）：(2,1),(4,1),(2,3)，面积 2
  const std::vector<kstab::QPoint> shifted = {
      {Rational(2), Rational(1)},
      {Rational(4), Rational(1)},
      {Rational(2), Rational(3)}};
  const kstab::Moments sm = kstab::polygon_moments(shifted);
  require(sm.area == 2, "shifted triangle area");
  // 质心 (8/3, 5/3)：∫x = 2·8/3 = 16/3
  require(sm.ix == make_rational(16, 3), "shifted triangle ∫x");
  require(sm.iy == make_rational(10, 3), "shifted triangle ∫y");
}

void test_boundary_moments() {
  // 边 (0,0)->(2,2) 的格点长度为 2
  const std::vector<IntPoint> skinny = {{0, 0}, {2, 2}, {0, 2}};
  const auto moments = kstab::boundary_moments(skinny);
  // 三条边：gcd(2,2)=2, gcd(2,0)=2, gcd(0,2)=2 → 总长 6
  require(moments.length == 6, "skinny triangle boundary length");

  require(kstab::boundary_moments(unit_square()).length == 4,
          "unit square boundary length");
  require(kstab::boundary_moments(right_triangle()).length == 3,
          "right triangle boundary length");
  require(kstab::boundary_moments(trapezoid()).length == 5,
          "trapezoid boundary length");

  // 正方形 ∫_{∂P} x dσ = 2（右边贡献 1，上下边各 1/2）
  require(kstab::boundary_moments(unit_square()).ix == 2,
          "unit square boundary ∫x");
}

void test_defining_property() {
  // ell_P 定义性质：任意仿射 h 精确满足 M_l(h) = 0。
  const auto vertices = trapezoid();
  const auto ell = kstab::compute_ell_p(vertices);
  const Rational affine_h[][3] = {
      {Rational(1), Rational(0), Rational(0)},
      {Rational(0), Rational(1), Rational(0)},
      {Rational(0), Rational(0), Rational(1)},
      {make_rational(3, 7), make_rational(-2, 5), make_rational(11, 13)},
  };
  for (const auto& h : affine_h) {
    const Rational value =
        kstab::df_simple_exact(vertices, ell, h[0], h[1], h[2]);
    require(value == 0,
            "M_l must vanish on affine functions (defining property)");
  }
}

void test_unit_square_hand_computed() {
  // 手算：单位正方形 g = max{x - 1/2, 0} → M_l = 1/4
  const auto ell = kstab::compute_ell_p(unit_square());
  const Rational value = kstab::df_simple_exact(
      unit_square(), ell, Rational(1), Rational(0), make_rational(-1, 2));
  require(value == make_rational(1, 4),
          "unit square M_l(max{x-1/2,0}) should be exactly 1/4");
}

void test_out_of_range_zero() {
  // t 在投影范围之外时 M = 0（double 路径）
  const auto vertices = trapezoid();
  const auto ell = kstab::ell_to_double(kstab::compute_ell_p(vertices));
  // 方向 (0,1)：投影范围 [0,1]
  require_close(kstab::df_simple_double(vertices, ell, 0.0, 1.0, 1.0), 0.0,
                1e-12, "M should be 0 above w_max");
  require_close(kstab::df_simple_double(vertices, ell, 0.0, 1.0, -1.0), 0.0,
                1e-12, "M should be 0 below w_min (affine on P)");
}

void test_double_exact_consistency() {
  const auto vertices = trapezoid();
  const auto ell_exact = kstab::compute_ell_p(vertices);
  const auto ell_d = kstab::ell_to_double(ell_exact);
  const struct {
    double a, b, c;
  } lines[] = {
      {0.0, 1.0, -0.4}, {1.0, 0.0, -0.7}, {0.6, 0.8, -0.9}, {-1.0, 0.25, 1.3},
  };
  for (const auto& line : lines) {
    const Rational exact_value = kstab::df_simple_exact(
        vertices, ell_exact, kstab::approximate_rational(line.a, 1000000),
        kstab::approximate_rational(line.b, 1000000),
        kstab::approximate_rational(line.c, 1000000));
    const double double_value =
        kstab::df_simple_double(vertices, ell_d, line.a, line.b, line.c);
    require_close(double_value, kstab::rational_to_double(exact_value), 1e-9,
                  "double and exact M_l should agree");
  }
}

void test_degenerate_creases() {
  const auto vertices = unit_square();
  const auto ell = kstab::ell_to_double(kstab::compute_ell_p(vertices));
  // 折痕过两顶点（对角线 x+y=1）
  const double diagonal =
      kstab::df_simple_double(vertices, ell, 1.0, 1.0, -1.0);
  require(std::isfinite(diagonal), "crease through vertices must be finite");
  // 折痕沿边（x=0）：g = max{x,0} 是仿射 → M = 0
  require_close(kstab::df_simple_double(vertices, ell, 1.0, 0.0, 0.0), 0.0,
                1e-9, "crease along edge: g affine on P, M = 0");
  // 裁剪为空
  require_close(kstab::df_simple_double(vertices, ell, 1.0, 0.0, -5.0), 0.0,
                1e-12, "empty clip should give M = 0");
  // 精确路径：折痕沿边
  const auto ell_exact = kstab::compute_ell_p(vertices);
  require(kstab::df_simple_exact(vertices, ell_exact, Rational(1),
                                 Rational(0), Rational(0)) == 0,
          "exact crease along edge should give M = 0");
}

void test_parser() {
  // 顺时针输入归一化后 ell_P 相同
  const std::vector<IntPoint> clockwise = {{0, 0}, {0, 1}, {1, 1}, {1, 0}};
  const auto ccw = kstab::normalize_polygon(clockwise);
  const auto ell = kstab::compute_ell_p(ccw);
  require(ell[0] == 4 && ell[1] == 0 && ell[2] == 0,
          "clockwise input should normalize to CCW");

  // 首尾重复闭合写法
  const std::vector<IntPoint> closed = {
      {0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0}};
  require(kstab::normalize_polygon(closed).size() == 4,
          "closing duplicate vertex should be dropped");

  // 共线连续顶点被合并，边界测度不变
  const std::vector<IntPoint> with_collinear = {
      {0, 0}, {1, 0}, {2, 0}, {1, 1}, {0, 1}};
  const auto merged = kstab::normalize_polygon(with_collinear);
  require(merged.size() == 4, "collinear vertex should be merged");
  require(kstab::boundary_moments(merged).length == 5,
          "merged polygon boundary length preserved");

  // 非凸输入报错
  bool thrown = false;
  try {
    kstab::normalize_polygon({{0, 0}, {2, 0}, {1, 1}, {2, 2}, {0, 2}});
  } catch (const std::runtime_error&) {
    thrown = true;
  }
  require(thrown, "non-convex polygon should be rejected");

  // 面积为零报错
  thrown = false;
  try {
    kstab::normalize_polygon({{0, 0}, {1, 1}, {2, 2}});
  } catch (const std::runtime_error&) {
    thrown = true;
  }
  require(thrown, "zero-area polygon should be rejected");
}

void test_approximate_rational() {
  require(kstab::approximate_rational(0.5 + 1e-9, 10) == make_rational(1, 2),
          "0.5+1e-9 should round to 1/2 at cap 10");
  require(kstab::approximate_rational(-2.25, 100) == make_rational(-9, 4),
          "-2.25 should be -9/4");
  const Rational pi_approx =
      kstab::approximate_rational(3.141592653589793, 100);
  // 分母 ≤ 100 的最佳逼近是 311/99（|π−311/99| < |π−22/7|）
  require(pi_approx == make_rational(311, 99),
          "pi at cap 100 should be 311/99");
  const Rational pi_better =
      kstab::approximate_rational(3.141592653589793, 1000);
  require(pi_better == make_rational(355, 113),
          "pi at cap 1000 should be 355/113");
}

void test_semistable_square_sweep() {
  // 单位正方形全扫：不得找到反例
  const auto vertices = unit_square();
  const auto ell = kstab::compute_ell_p(vertices);
  kstab::SearchOptions options;
  options.theta_steps = 180;
  options.t_steps = 128;
  const kstab::SearchResult result =
      kstab::search_witness(vertices, ell, options);
  require(!result.unstable, "unit square must be no_counterexample_found");
  require(result.witness.normalized > -1e-8,
          "unit square sweep min should be ~0");
}

void test_semistable_integration() {
  // 更多真实多边形（相对 K-半稳定，与 Wang--Zhou 定理一致）：
  // 搜索不得找到反例。
  const std::vector<IntPoint> my_pentagon = {
      {0, 0}, {4, 1}, {3, 2}, {0, 3}, {-1, 1}};
  const std::vector<IntPoint> wang_zhou_a1 = {
      {0, 21}, {1, 21}, {2, 20}, {3, 18}, {4, 15},
      {5, 11}, {7, 1}, {7, -1}, {0, -1}};
  for (const auto& vertices : {trapezoid(), my_pentagon, wang_zhou_a1}) {
    const auto ell = kstab::compute_ell_p(vertices);
    kstab::SearchOptions options;
    options.theta_steps = 120;
    options.t_steps = 128;
    const kstab::SearchResult result =
        kstab::search_witness(vertices, ell, options);
    require(!result.unstable,
            "semistable polygon must be no_counterexample_found");
  }
}

void test_unstable_machinery() {
  // 机械测试：故意传入错误的 ell（ℓ̃ ≡ 8 而非真正的 ℓ_P ≡ 4），
  // 验证当 M_ℓ 确实存在负值时，搜索能找到、认证能确认。
  // 单位正方形、ℓ̃ ≡ 8：g = max{x-1/2, 0} 时
  //   M = ∫_{∂P} g dσ - 8 ∫_P g = 3/4 - 8·(1/8) = -1/4。
  const auto vertices = unit_square();
  const std::array<Rational, 3> fake_ell = {Rational(8), Rational(0),
                                            Rational(0)};
  kstab::SearchOptions options;
  options.theta_steps = 180;
  options.t_steps = 128;
  const kstab::SearchResult result =
      kstab::search_witness(vertices, fake_ell, options);
  require(result.unstable,
          "machinery test: search must detect the negative dip");
  require(result.witness.value <= -0.2,
          "machinery test: witness value should be <= -0.2");
  const kstab::CertifyResult certification =
      kstab::certify_witness(vertices, fake_ell, result.witness, 1000000);
  require(certification.certified,
          "machinery test: certification must succeed");
  require(certification.value < 0,
          "machinery test: certified M_l must be negative");

  // 带 witness 的 SVG：应包含折痕线与状态文本（真实输入全是半稳定，
  // CLI 路径画不到 witness，借此机械用例覆盖折痕/阴影绘制）。
  // 手工构造折痕 x = 1/2，保证弦与多边形内部相交。
  kstab::Witness mid_witness;
  mid_witness.ux = 1.0;
  mid_witness.uy = 0.0;
  mid_witness.t = 0.5;
  mid_witness.value = -0.25;
  const std::filesystem::path svg_path =
      std::filesystem::temp_directory_path() / "kstab_machinery_witness.svg";
  const auto ell = kstab::compute_ell_p(vertices);
  require(kstab::write_svg(svg_path, vertices, ell, &mid_witness),
          "write_svg with witness should succeed");
  std::ifstream svg(svg_path);
  std::stringstream buffer;
  buffer << svg.rdbuf();
  const std::string content = buffer.str();
  require(content.find("<svg") != std::string::npos &&
              content.find("</svg>") != std::string::npos,
          "SVG should be well-formed");
  require(content.find("<line") != std::string::npos,
          "SVG with witness should contain the crease line");
  require(content.find("relatively K-unstable") != std::string::npos,
          "SVG with witness should contain the status text");
  std::filesystem::remove(svg_path);
}

}  // namespace

int main() {
  const std::pair<std::string, void (*)()> tests[] = {
      {"ell_p_exact", test_ell_p_exact},
      {"polygon_moments", test_polygon_moments},
      {"boundary_moments", test_boundary_moments},
      {"defining_property", test_defining_property},
      {"unit_square_hand_computed", test_unit_square_hand_computed},
      {"out_of_range_zero", test_out_of_range_zero},
      {"double_exact_consistency", test_double_exact_consistency},
      {"degenerate_creases", test_degenerate_creases},
      {"parser", test_parser},
      {"approximate_rational", test_approximate_rational},
      {"semistable_square_sweep", test_semistable_square_sweep},
      {"semistable_integration", test_semistable_integration},
      {"unstable_machinery", test_unstable_machinery},
  };
  int failed = 0;
  for (const auto& [name, test] : tests) {
    try {
      test();
      std::cout << "PASS " << name << '\n';
    } catch (const std::exception& error) {
      ++failed;
      std::cerr << "FAIL " << name << ": " << error.what() << '\n';
    }
  }
  if (failed != 0) {
    std::cerr << failed << " test(s) failed.\n";
    return 1;
  }
  std::cout << "All " << std::size(tests) << " tests passed.\n";
  return 0;
}
