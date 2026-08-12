#include "gkz/gkz.hpp"

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
  std::optional<std::filesystem::path> points;
  std::optional<std::filesystem::path> polygon;
  std::optional<std::filesystem::path> output;
  std::optional<std::filesystem::path> plot_prefix;
  std::int64_t k = 0;
  gkz::SolverOptions options;
};

void print_usage(std::ostream& stream) {
  stream
      << "Usage:\n"
      << "  shortest_gkz --points FILE [options]\n"
      << "  shortest_gkz --polygon FILE --k INTEGER [options]\n\n"
      << "Options:\n"
      << "  --output FILE            Write sigma and ell_A values to CSV.\n"
      << "  --plot-prefix PREFIX     Write plot CSVs for plot_results.py.\n"
      << "  --tolerance VALUE        Relative Frank--Wolfe gap tolerance.\n"
      << "  --absolute-tolerance V   Absolute gap tolerance.\n"
      << "  --prune-tolerance V      Drop QP coefficients at or below V (default: 1e-15).\n"
      << "  --max-iterations N       Maximum active-set expansion steps.\n"
      << "  --exact-max-active N     Exact QP variable limit; 0 means no limit.\n"
      << "  --no-exact               Skip exact rational certification.\n"
      << "  --projection             Enable experimental stable-face projection.\n"
      << "  --projection-window N    Gap history window (default: 32).\n"
      << "  --projection-stall-ratio R  Minimum window gap ratio (default: 0.90).\n"
      << "  --projection-relative-gap R  Activation gap / initial gap (default: 1e-2).\n"
      << "  --projection-rank-stall-window N  Rank-stable observations (default: 64).\n"
      << "  --projection-rank-tolerance R  Relative affine-rank tolerance (default: 1e-11).\n"
      << "  --verbose                Print one line per iteration.\n"
      << "  --help                   Show this message.\n";
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
    } else if (option == "--points") {
      arguments.points = require_value(argc, argv, i, option);
    } else if (option == "--polygon") {
      arguments.polygon = require_value(argc, argv, i, option);
    } else if (option == "--k") {
      arguments.k = std::stoll(require_value(argc, argv, i, option));
    } else if (option == "--output") {
      arguments.output = require_value(argc, argv, i, option);
    } else if (option == "--plot-prefix") {
      arguments.plot_prefix = require_value(argc, argv, i, option);
    } else if (option == "--tolerance") {
      arguments.options.tolerance =
          std::stod(require_value(argc, argv, i, option));
    } else if (option == "--absolute-tolerance") {
      arguments.options.absolute_tolerance =
          std::stod(require_value(argc, argv, i, option));
    } else if (option == "--prune-tolerance") {
      arguments.options.prune_tolerance =
          std::stod(require_value(argc, argv, i, option));
    } else if (option == "--max-iterations") {
      arguments.options.max_iterations =
          std::stoi(require_value(argc, argv, i, option));
    } else if (option == "--exact-max-active") {
      arguments.options.exact_max_active =
          std::stoi(require_value(argc, argv, i, option));
    } else if (option == "--no-exact") {
      arguments.options.exact_certification = false;
    } else if (option == "--projection") {
      arguments.options.projection = true;
    } else if (option == "--projection-window") {
      arguments.options.projection_window =
          std::stoi(require_value(argc, argv, i, option));
    } else if (option == "--projection-stall-ratio") {
      arguments.options.projection_stall_ratio =
          std::stod(require_value(argc, argv, i, option));
    } else if (option == "--projection-relative-gap") {
      arguments.options.projection_relative_gap =
          std::stod(require_value(argc, argv, i, option));
    } else if (option == "--projection-rank-stall-window") {
      arguments.options.projection_rank_stall_window =
          std::stoi(require_value(argc, argv, i, option));
    } else if (option == "--projection-rank-tolerance") {
      arguments.options.projection_rank_tolerance =
          std::stod(require_value(argc, argv, i, option));
    } else if (option == "--verbose") {
      arguments.options.verbose = true;
    } else {
      throw std::invalid_argument("Unknown option: " + option);
    }
  }

  if (arguments.points.has_value() == arguments.polygon.has_value()) {
    throw std::invalid_argument(
        "Specify exactly one of --points and --polygon.");
  }
  if (arguments.polygon && arguments.k <= 0) {
    throw std::invalid_argument("--polygon requires a positive --k.");
  }
  if (arguments.options.tolerance < 0.0 ||
      arguments.options.absolute_tolerance < 0.0 ||
      arguments.options.prune_tolerance < 0.0 ||
      arguments.options.max_iterations < 0 ||
      arguments.options.exact_max_active < 0 ||
      arguments.options.projection_window <= 0 ||
      arguments.options.projection_rank_stall_window <= 0 ||
      arguments.options.projection_stall_ratio <= 0.0 ||
      arguments.options.projection_stall_ratio > 1.0 ||
      arguments.options.projection_relative_gap <= 0.0 ||
      arguments.options.projection_relative_gap > 1.0 ||
      arguments.options.projection_rank_tolerance <= 0.0) {
    throw std::invalid_argument(
        "Invalid tolerance, iteration count, or projection parameter.");
  }
  return arguments;
}

std::string rational_string(const CGAL::Gmpq& value) {
  std::ostringstream stream;
  if (value.denominator() == 1) {
    stream << value.numerator();
  } else {
    stream << value;
  }
  return stream.str();
}

std::string coordinate_string(std::int64_t value,
                              const gkz::PointConfiguration& configuration) {
  if (!configuration.is_polygon_level()) {
    return std::to_string(value);
  }
  return rational_string(CGAL::Gmpq(value) /
                         CGAL::Gmpq(configuration.level()));
}

std::string affine_expression(const gkz::AffineFunction& ell,
                              bool exact) {
  std::ostringstream stream;
  stream << (exact ? "ell_A(x,y) = " : "ell_A(x,y) ~= ");
  if (exact) {
    stream << rational_string(ell.exact_coefficients[0]) << " + ("
           << rational_string(ell.exact_coefficients[1]) << ") x + ("
           << rational_string(ell.exact_coefficients[2]) << ") y";
  } else {
    stream << std::setprecision(17) << ell.coefficients[0] << " + ("
           << ell.coefficients[1] << ") x + (" << ell.coefficients[2]
           << ") y";
  }
  return stream.str();
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    const gkz::PointConfiguration configuration = arguments.points
        ? gkz::PointConfiguration::from_points_file(*arguments.points)
        : gkz::PointConfiguration::from_polygon_file(*arguments.polygon,
                                                       arguments.k);

    std::cout << "points=" << configuration.size() << '\n' << std::flush;
    const gkz::ShortestGkzSolver solver(arguments.options);
    const gkz::SolverResult result = solver.solve(configuration);
    const gkz::AffineFunction ell = gkz::compute_ell(configuration);

    bool ell_equals_sigma = false;
    if (result.exact.certified) {
      ell_equals_sigma = true;
      for (std::size_t i = 0; i < configuration.size(); ++i) {
        if (ell.exact_values[i] != result.exact.sigma[i]) {
          ell_equals_sigma = false;
          break;
        }
      }
    } else {
      const double comparison_tolerance = 1e-9;
      ell_equals_sigma =
          (ell.values - result.sigma).lpNorm<Eigen::Infinity>() <=
          comparison_tolerance;
    }
    const bool ell_is_constant = ell.exact_coefficients[1] == 0 &&
                                  ell.exact_coefficients[2] == 0;

    std::cout << std::setprecision(17);
    std::cout << "twice_area="
              << gkz::to_string(configuration.twice_area()) << '\n';
    if (configuration.is_polygon_level()) {
      std::cout << "level=" << configuration.level() << '\n';
      std::cout << "base_twice_area="
                << gkz::to_string(configuration.base_twice_area()) << '\n';
    }
    std::cout << "converged=" << std::boolalpha << result.converged << '\n';
    std::cout << "iterations=" << result.iterations << '\n';
    std::cout << "active_size=" << result.active_vectors.size() << '\n';
    std::cout << "norm_squared=" << result.norm_squared << '\n';
    std::cout << "gap=" << result.gap << '\n';
    std::cout << "l2_error_bound=" << result.l2_error_bound << '\n';
    std::cout << "exact_certified=" << result.exact.certified << '\n';
    std::cout << "projection_enabled=" << result.projection_enabled << '\n';
    std::cout << "projection_start_iteration="
              << result.projection_start_iteration << '\n';
    std::cout << "stable_projection_available="
              << result.stable_projection_available << '\n';
    std::cout << "stable_projection_stop_iteration="
              << result.stable_projection_stop_iteration << '\n';
    std::cout << "stable_projection_rank="
              << result.stable_projection_rank << '\n';
    std::cout << "stable_projection_observations="
              << result.stable_projection_observations << '\n';
    std::cout << "stable_projection_stop_reason="
              << (result.stable_projection_stop_reason.empty()
                      ? "none"
                      : result.stable_projection_stop_reason)
              << '\n';
    if (result.stable_projection_available) {
      std::cout << "stable_projection_norm_squared="
                << result.stable_projection_norm_squared << '\n';
    }
    std::cout << "final_qp_performed=" << result.final_qp_performed << '\n';
    if (result.final_qp_performed) {
      std::cout << "final_qp_norm_squared=" << result.final_qp_norm_squared
                << '\n';
      std::cout << "final_qp_gap=" << result.final_qp_gap << '\n';
    }
    if (arguments.options.exact_certification) {
      std::cout << "exact_message=" << result.exact.message << '\n';
      if (result.exact.certified) {
        std::cout << "exact_norm_squared="
                  << rational_string(result.exact.norm_squared) << '\n';
      }
    }

    std::cout << affine_expression(ell, true) << '\n';
    std::cout << "ell_A_constant=" << ell_is_constant << '\n';
    std::cout << "ell_A_sigma_comparison="
              << (result.exact.certified ? "exact" : "numerical") << '\n';
    std::cout << "ell_A_equals_shortest_GKZ=" << ell_equals_sigma << '\n';
    if (result.exact.certified) {
      std::cout << "relative_Chow_semistable=" << ell_equals_sigma << '\n';
    } else {
      std::cout << "relative_Chow_semistable=unverified\n";
    }
    std::cout << "ell_A_on_A=";
    if (configuration.size() <= 256 ||
        (!arguments.output && !arguments.plot_prefix)) {
      for (std::size_t i = 0; i < configuration.size(); ++i) {
        if (i != 0) {
          std::cout << ';';
        }
        std::cout << '('
                  << coordinate_string(configuration.points()[i].x,
                                       configuration)
                  << ','
                  << coordinate_string(configuration.points()[i].y,
                                       configuration)
                  << ")="
                  << rational_string(ell.exact_values[i]);
      }
      std::cout << '\n';
    } else {
      std::cout << "(written to "
                << (arguments.output ? arguments.output->string()
                                     : arguments.plot_prefix->string() +
                                           "_ell.csv")
                << "; " << configuration.size() << " values)\n";
    }

    if (arguments.output) {
      gkz::write_result_csv(*arguments.output, configuration, result);
      std::cout << "output=" << arguments.output->string() << '\n';
    }
    if (arguments.plot_prefix) {
      gkz::write_plot_data(*arguments.plot_prefix, configuration, result, ell);
      std::cout << "plot_surface="
                << arguments.plot_prefix->string() + "_surface.csv" << '\n';
      std::cout << "plot_triangles="
                << arguments.plot_prefix->string() + "_triangles.csv" << '\n';
      std::cout << "plot_subdivision="
                << arguments.plot_prefix->string() + "_subdivision.csv"
                << '\n';
      std::cout << "plot_ell=" << arguments.plot_prefix->string() + "_ell.csv"
                << '\n';
      if (result.stable_projection_available) {
        const std::filesystem::path stable_prefix =
            arguments.plot_prefix->string() + "_stable";
        gkz::write_stable_projection_plot_data(stable_prefix, configuration,
                                               result);
        std::cout << "stable_plot_surface="
                  << stable_prefix.string() + "_surface.csv" << '\n';
        std::cout << "stable_plot_triangles="
                  << stable_prefix.string() + "_triangles.csv" << '\n';
        std::cout << "stable_plot_subdivision="
                  << stable_prefix.string() + "_subdivision.csv" << '\n';
      }
    }
    return result.converged ? 0 : 2;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    print_usage(std::cerr);
    return 1;
  }
}
