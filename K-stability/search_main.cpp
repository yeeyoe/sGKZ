#include "search.hpp"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void usage(std::ostream& out) {
  out << "Usage: k_stability_search --d D --N N --M M --time-limit SEC [options]\n\n"
      << "Options:\n"
      << "  --backfill-q           Fill exact Q_P(g)^2 for existing certified witnesses\n"
      << "  --database FILE       SQLite state file (default K-stability/k_stability_search.sqlite)\n"
      << "  --output-dir DIR      Result report directory (default .)\n"
      << "  --shell-seconds SEC   Per-shell time slice (default 60)\n"
      << "  --beam-width N        Candidates per generation batch (default 48)\n"
      << "  --seed N              Deterministic random seed\n"
      << "  --stop-on-first       Stop after the first exact certification\n"
      << "  --smooth-only         Search only polygons smooth at every vertex\n"
      << "  --certify-max-denom N Rational denominator cap\n"
      << "  --verbose             Print progress\n";
}

std::string value(int argc, char** argv, int& i, const std::string& option) {
  if (++i >= argc) throw std::invalid_argument("missing value after " + option);
  return argv[i];
}

}  // namespace

int main(int argc, char** argv) {
  kstab::AreaSearchOptions options;
  bool backfill_q = false;
  try {
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "--help") { usage(std::cout); return 0; }
      if (arg == "--backfill-q") { backfill_q = true; continue; }
      if (arg == "--d") options.d = std::stoi(value(argc, argv, i, arg));
      else if (arg == "--N") options.initial_N = std::stoi(value(argc, argv, i, arg));
      else if (arg == "--M") options.initial_M = std::stoi(value(argc, argv, i, arg));
      else if (arg == "--time-limit") options.time_limit_seconds = std::stod(value(argc, argv, i, arg));
      else if (arg == "--shell-seconds") options.shell_seconds = std::stod(value(argc, argv, i, arg));
      else if (arg == "--beam-width") options.beam_width = std::stoi(value(argc, argv, i, arg));
      else if (arg == "--seed") options.seed = std::stoull(value(argc, argv, i, arg));
      else if (arg == "--database") options.database = value(argc, argv, i, arg);
      else if (arg == "--output-dir") options.output_directory = value(argc, argv, i, arg);
      else if (arg == "--certify-max-denom") options.certify_max_denominator = std::stoll(value(argc, argv, i, arg));
      else if (arg == "--stop-on-first") options.stop_on_first = true;
      else if (arg == "--smooth-only") options.smooth_only = true;
      else if (arg == "--verbose") options.verbose = true;
      else throw std::invalid_argument("unknown option: " + arg);
    }
    if (backfill_q) {
      kstab::SearchDatabase database(options.database);
      const auto summary = database.backfill_q_squared();
      std::cout << "backfill_scanned=" << summary.scanned << '\n'
                << "backfill_success=" << summary.backfilled << '\n'
                << "backfill_skipped=" << summary.skipped << '\n'
                << "backfill_errors=" << summary.errors << '\n';
      return summary.errors == 0 ? 0 : 1;
    }
    if (options.d < 3 || options.initial_N < 1 || options.initial_M < 1 ||
        options.beam_width < 1 || options.time_limit_seconds <= 0 ||
        options.shell_seconds <= 0) {
      throw std::invalid_argument("d, N, M, beam-width and time limits must be positive");
    }
    const kstab::SearchSummary summary = kstab::run_search(options);
    std::cout << "database=" << options.database.string() << '\n'
              << "report_file="
              << (options.output_directory / "k_stability_search_result.txt").string() << '\n'
              << "total tested: " << summary.total_tested << '\n'
              << "new tested: " << summary.new_tested << '\n'
              << "total verified unstable: " << summary.total_verified_unstable << '\n'
              << "new verified unstable: " << summary.new_verified_unstable << '\n'
              << (summary.smaller_volume_found ? "smaller volume found" :
                                                   "same least volume") << '\n'
              << "The top 5 with least volume:\n";
    for (const auto& verified : summary.top_verified) {
      std::cout << "key=" << verified.key << " twice_area=" << verified.twice_area
                << " q_squared_exact=" << verified.q_squared_exact
                << " q_squared_value=" << verified.q_squared_value << '\n';
    }
    std::cout << "generated=" << summary.generated << '\n'
              << "rejected=" << summary.rejected << '\n'
              << "probes=" << summary.probes << '\n'
              << "confirms=" << summary.confirms << '\n'
              << "finals=" << summary.finals << '\n'
              << "verified=" << summary.verified << '\n'
              << "unverified=" << summary.unverified << '\n'
              << "skipped=" << summary.skipped << '\n'
              << "have_verified=" << std::boolalpha << summary.have_verified << '\n';
    if (summary.have_verified) {
      std::cout << "best_twice_area=" << summary.best_twice_area << '\n'
                << "first_verified_key=" << summary.first_verified_key << '\n'
                << "best_verified_key=" << summary.best_verified_key << '\n';
    }
    return summary.have_verified ? 0 : 2;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    usage(std::cerr);
    return 1;
  }
}
