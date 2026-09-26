#include "unstable_delzant.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

namespace {

std::string value(int argc, char** argv, int& i, const std::string& option) {
  if (++i >= argc) throw std::invalid_argument("missing value after " + option);
  return argv[i];
}

void usage(std::ostream& out) {
  out << "Usage: unstable_Delzant --d D --time-limit SEC [options]\n\n"
      << "Search relatively K-unstable smooth lattice polygons.\n"
      << "Options:\n"
      << "  --N N                   Fan coordinate bound (default 12)\n"
      << "  --M M                   Initial edge-length bound (default 12)\n"
      << "  --max-diameter D       Required max L-infinity diameter after scale normalization\n"
      << "  --seed S                Random seed (default 1)\n"
      << "  --database FILE         SQLite database\n"
      << "  --output-dir DIR        Result directory\n"
      << "  --theta-steps N         Numerical direction steps (default 360)\n"
      << "  --t-steps N             Numerical offset steps (default 256)\n"
      << "  --certify-max-denom N   Rational witness cap (default 1048576)\n"
      << "  --help                  Show this help\n";
}

}  // namespace

int main(int argc, char** argv) {
  delzant::SearchOptions options;
  bool have_d = false;
  bool have_time_limit = false;
  bool have_max_diameter = false;
  try {
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "--help") { usage(std::cout); return 0; }
      if (arg == "--d") { options.d = std::stoi(value(argc, argv, i, arg)); have_d = true; }
      else if (arg == "--time-limit") {
        options.time_limit_seconds = std::stod(value(argc, argv, i, arg));
        have_time_limit = true;
      } else if (arg == "--N") options.N = std::stoi(value(argc, argv, i, arg));
      else if (arg == "--M") options.M = std::stoi(value(argc, argv, i, arg));
      else if (arg == "--max-diameter") {
        options.max_diameter = std::stoll(value(argc, argv, i, arg));
        have_max_diameter = true;
      }
      else if (arg == "--seed") options.seed = std::stoull(value(argc, argv, i, arg));
      else if (arg == "--database") options.database = value(argc, argv, i, arg);
      else if (arg == "--output-dir") options.output_directory = value(argc, argv, i, arg);
      else if (arg == "--theta-steps") options.theta_steps = std::stoi(value(argc, argv, i, arg));
      else if (arg == "--t-steps") options.t_steps = std::stoi(value(argc, argv, i, arg));
      else if (arg == "--certify-max-denom")
        options.certify_max_denominator = std::stoll(value(argc, argv, i, arg));
      else throw std::invalid_argument("unknown option: " + arg);
    }
    if (!have_d || !have_time_limit || !have_max_diameter)
      throw std::invalid_argument("--d, --time-limit, and --max-diameter are required");
    if (options.d < 3 || options.N < 1 || options.M < 1 ||
        options.max_diameter < 1 || options.time_limit_seconds <= 0 || options.theta_steps < 8 ||
        options.t_steps < 8 || options.certify_max_denominator < 10)
      throw std::invalid_argument("invalid search parameter");

    const auto summary = delzant::run_search(options);
    std::cout << "database=" << options.database.string() << '\n'
              << "max_diameter=" << options.max_diameter << '\n'
              << "generated=" << summary.generated << '\n'
              << "oversized=" << summary.oversized << '\n'
              << "duplicates=" << summary.duplicates << '\n'
              << "tested=" << summary.tested << '\n'
              << "probed=" << summary.probed << '\n'
              << "confirmed=" << summary.confirmed << '\n'
              << "finalized=" << summary.finalized << '\n'
              << "df_evaluations=" << summary.df_evaluations << '\n'
              << "certified_candidates=" << summary.certified_candidates << '\n'
              << "probe_seconds=" << summary.probe_seconds << '\n'
              << "confirm_seconds=" << summary.confirm_seconds << '\n'
              << "final_seconds=" << summary.final_seconds << '\n'
              << "generation_seconds=" << summary.generation_seconds << '\n'
              << "total_seconds=" << summary.total_seconds << '\n'
              << "tested_per_second="
              << (summary.total_seconds > 0.0
                      ? static_cast<double>(summary.tested) / summary.total_seconds
                      : 0.0)
              << '\n'
              << "found=" << std::boolalpha << summary.found << '\n';
    if (summary.found) {
      std::cout << "candidate_key=" << summary.candidate.key << '\n'
                << "vertex_count=" << summary.candidate.vertices.size() << '\n'
                << "twice_area=" << summary.candidate.twice_area << '\n'
                << "ell_P=" << summary.candidate.ell[0] << " + ("
                << summary.candidate.ell[1] << ")x + ("
                << summary.candidate.ell[2] << ")y\n"
                << "witness=" << summary.certificate.coefficients[0] << ' '
                << summary.certificate.coefficients[1] << ' '
                << summary.certificate.coefficients[2] << '\n'
                << "certified_M_l=" << summary.certificate.value << '\n'
                << "polygon_file="
                << (options.output_directory /
                    ("unstable_d" + std::to_string(options.d) + ".polygon")).string()
                << '\n';
      for (const auto& p : summary.candidate.vertices)
        std::cout << "vertex=" << p.x << ' ' << p.y << '\n';
      return 0;
    }
    std::cout << "status=no_certified_candidate_before_timeout\n";
    return 2;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    usage(std::cerr);
    return 1;
  }
}
