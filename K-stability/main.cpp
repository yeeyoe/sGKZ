// k_stability：精确计算格点多边形的 ell_P，并用 Donaldson 事实
// 搜索简单凸函数 witness g = max{<x,u> - t, 0} 使 M_l(g) < 0。
// 数学定义见 paper/K-stability.tex。

#include "k_stability.hpp"

#include <CGAL/Gmpz.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

struct Arguments {
  std::optional<std::filesystem::path> polygon;
  std::optional<std::filesystem::path> svg;
  std::optional<std::string> check_line;
  kstab::SearchOptions search;
  bool certify = false;
  std::int64_t certify_max_denom = 1048576;
};

void print_usage(std::ostream& stream) {
  stream
      << "Usage:\n"
      << "  k_stability --polygon FILE [options]\n\n"
      << "Options:\n"
      << "  --theta-steps N        Crease directions in [0, 2*pi) (default 720).\n"
      << "  --t-steps N            Offset samples per direction (default 512).\n"
      << "  --no-refine            Skip local (theta,t) refinement.\n"
      << "  --certify              Certify the witness with rational arithmetic.\n"
      << "  --certify-max-denom N  Denominator cap for certification.\n"
      << "  --svg FILE             Write polygon + witness crease SVG.\n"
      << "  --check-line \"a b c\"   Evaluate M_l(max{a x + b y + c, 0}).\n"
      << "  Polygon files may end with a 'null measure edges' section; each\n"
      << "  following line 'x1 y1 x2 y2' sets that polygon edge's dσ to zero.\n"
      << "  --verbose              Print per-direction minima.\n"
      << "  --help                 Show this message.\n";
}

std::string require_value(int argc, char** argv, int& index,
                          const std::string& option) {
  if (index + 1 >= argc) {
    throw std::invalid_argument("Missing value after " + option + ".");
  }
  return argv[++index];
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments arguments;
  for (int i = 1; i < argc; ++i) {
    const std::string option = argv[i];
    if (option == "--help") {
      print_usage(std::cout);
      std::exit(0);
    } else if (option == "--polygon") {
      arguments.polygon = require_value(argc, argv, i, option);
    } else if (option == "--svg") {
      arguments.svg = require_value(argc, argv, i, option);
    } else if (option == "--check-line") {
      arguments.check_line = require_value(argc, argv, i, option);
    } else if (option == "--theta-steps") {
      arguments.search.theta_steps =
          std::stoi(require_value(argc, argv, i, option));
    } else if (option == "--t-steps") {
      arguments.search.t_steps =
          std::stoi(require_value(argc, argv, i, option));
    } else if (option == "--no-refine") {
      arguments.search.refine = false;
    } else if (option == "--certify") {
      arguments.certify = true;
    } else if (option == "--certify-max-denom") {
      arguments.certify_max_denom =
          std::stoll(require_value(argc, argv, i, option));
    } else if (option == "--verbose") {
      arguments.search.verbose = true;
    } else {
      throw std::invalid_argument("Unknown option: " + option);
    }
  }
  if (!arguments.polygon.has_value()) {
    throw std::invalid_argument("Missing required --polygon FILE.");
  }
  if (arguments.search.theta_steps <= 0 || arguments.search.t_steps <= 0 ||
      arguments.certify_max_denom < 1) {
    throw std::invalid_argument(
        "theta-steps/t-steps must be positive and certify-max-denom >= 1.");
  }
  return arguments;
}

std::string rational_string(const kstab::Rational& value) {
  std::ostringstream stream;
  if (value.denominator() == 1) {
    stream << value.numerator();
  } else {
    stream << value;
  }
  return stream.str();
}

// 解析 "p/q"、整数或浮点 token。前两种走精确路径。
kstab::Rational parse_coefficient(const std::string& token, bool& exact) {
  if (const auto slash = token.find('/'); slash != std::string::npos) {
    const long long numerator = std::stoll(token.substr(0, slash));
    const long long denominator = std::stoll(token.substr(slash + 1));
    if (denominator == 0) {
      throw std::invalid_argument("Zero denominator in --check-line.");
    }
    return kstab::Rational(CGAL::Gmpz(static_cast<long>(numerator)),
                           CGAL::Gmpz(static_cast<long>(denominator)));
  }
  if (token.find_first_not_of("+-0123456789") == std::string::npos) {
    return kstab::Rational(CGAL::Gmpz(static_cast<long>(std::stoll(token))));
  }
  exact = false;
  return kstab::Rational(std::stod(token));
}

}  // namespace

int main(int argc, char** argv) {
  Arguments arguments;
  try {
    arguments = parse_arguments(argc, argv);
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    print_usage(std::cerr);
    return 1;
  }

  try {
    const kstab::PolygonInput polygon =
        kstab::parse_polygon_input_file(*arguments.polygon);
    const std::vector<kstab::IntPoint>& vertices = polygon.vertices_ccw;
    const auto ell = kstab::compute_ell_p(vertices, polygon.null_measure_edges);
    const auto boundary =
        kstab::boundary_moments(vertices, polygon.null_measure_edges);

    __int128_t twice_area = 0;
    for (std::size_t i = 0; i < vertices.size(); ++i) {
      const auto& p = vertices[i];
      const auto& q = vertices[(i + 1) % vertices.size()];
      twice_area += static_cast<__int128_t>(p.x) * q.y -
                    static_cast<__int128_t>(p.y) * q.x;
    }
    std::ostringstream twice_area_stream;
    {
      __int128_t value = twice_area;
      std::string digits;
      if (value < 0) {
        twice_area_stream << '-';
        value = -value;
      }
      do {
        digits.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
      } while (value != 0);
      for (auto it = digits.rbegin(); it != digits.rend(); ++it) {
        twice_area_stream << *it;
      }
    }

    std::cout << "vertices=" << vertices.size() << '\n';
    std::cout << "twice_area=" << twice_area_stream.str() << '\n';
    std::cout << "boundary_length_dsigma=" << rational_string(boundary.length)
              << '\n';
    std::cout << "null_measure_edge_count="
              << std::count(polygon.null_measure_edges.begin(),
                            polygon.null_measure_edges.end(), true)
              << '\n';
    std::cout << "ell_P(x,y) = " << rational_string(ell[0]) << " + ("
              << rational_string(ell[1]) << ") x + ("
              << rational_string(ell[2]) << ") y\n";
    std::cout << std::boolalpha
              << "ell_P_constant=" << (ell[1] == 0 && ell[2] == 0) << '\n';

    if (arguments.check_line.has_value()) {
      std::istringstream parser(*arguments.check_line);
      std::string tokens[3];
      if (!(parser >> tokens[0] >> tokens[1] >> tokens[2])) {
        throw std::invalid_argument(
            "--check-line expects three coefficients \"a b c\".");
      }
      bool exact = true;
      const kstab::Rational a = parse_coefficient(tokens[0], exact);
      const kstab::Rational b = parse_coefficient(tokens[1], exact);
      const kstab::Rational c = parse_coefficient(tokens[2], exact);
      if (exact) {
        const kstab::Rational value =
            kstab::df_simple_exact(vertices, ell, a, b, c,
                                   polygon.null_measure_edges);
        std::cout << "check_line_exact=true\n";
        std::cout << "M_l=" << rational_string(value) << '\n';
      } else {
        const double value = kstab::df_simple_double(
            vertices, kstab::ell_to_double(ell), kstab::rational_to_double(a),
            kstab::rational_to_double(b), kstab::rational_to_double(c),
            polygon.null_measure_edges);
        std::cout << "check_line_exact=false\n";
        std::cout << "M_l=" << std::setprecision(17) << value << '\n';
      }
      return 0;
    }

    const kstab::SearchResult result =
        kstab::search_witness(vertices, ell, arguments.search,
                              polygon.null_measure_edges);
    std::cout << std::setprecision(17);
    std::cout << "search_evaluations=" << result.evaluations << '\n';
    std::cout << "sweep_min_M_l=" << result.witness.value << '\n';
    std::cout << "sweep_min_M_l_normalized=" << result.witness.normalized
              << '\n';
    std::cout << "witness_theta=" << std::atan2(result.witness.uy,
                                                 result.witness.ux)
              << '\n';
    std::cout << "witness_t=" << result.witness.t << '\n';
    std::cout << "witness_g=max{(" << result.witness.ux << ") x + ("
              << result.witness.uy << ") y + (" << -result.witness.t
              << "), 0}\n";
    std::cout << "relative_K_status="
              << (result.unstable ? "unstable" : "no_counterexample_found")
              << '\n';

    bool certified = false;
    if (arguments.certify && result.unstable) {
      const kstab::CertifyResult certification = kstab::certify_witness(
          vertices, ell, result.witness, arguments.certify_max_denom,
          polygon.null_measure_edges);
      certified = certification.certified;
      std::cout << "certified=" << certified << '\n';
      if (certified) {
        std::cout << "certified_g=max{("
                  << rational_string(certification.coefficients[0])
                  << ") x + ("
                  << rational_string(certification.coefficients[1])
                  << ") y + ("
                  << rational_string(certification.coefficients[2])
                  << "), 0}\n";
        std::cout << "certified_M_l=" << rational_string(certification.value)
                  << '\n';
      }
    } else if (arguments.certify) {
      std::cout << "certified=false\n";
    }

    if (arguments.svg.has_value()) {
      const kstab::Witness* witness =
          result.unstable ? &result.witness : nullptr;
      if (!kstab::write_svg(*arguments.svg, vertices, ell, witness)) {
        throw std::runtime_error("Cannot write SVG: " +
                                 arguments.svg->string());
      }
      std::cout << "svg=" << arguments.svg->string() << '\n';
    }

    if (!result.unstable) {
      return 2;
    }
    if (arguments.certify && !certified) {
      return 2;
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
