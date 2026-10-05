#include "chowk/chowk.hpp"

#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>

namespace {
void usage() {
  std::cerr << "Usage: chowk --polytope FILE --max-k K [options]\n"
            << "  --output-dir DIR              output directory (default: output)\n"
            << "  --threshold VALUE             normalized DF threshold\n"
            << "  --tolerance VALUE             solver relative tolerance\n"
            << "  --absolute-tolerance VALUE    solver absolute tolerance\n"
            << "  --max-iterations N            solver iteration limit\n"
            << "  --max-points N                lattice point limit\n"
            << "  --verbose                     print solver progress\n";
}

bool need_value(int index, int argc, const char* option) {
  if (index + 1 >= argc) {
    std::cerr << option << " requires a value\n";
    return false;
  }
  return true;
}
}  // namespace

int main(int argc, char** argv) {
  try {
    std::filesystem::path input;
    std::filesystem::path output = "output";
    chowk::Options options;
    bool have_input = false;
    bool have_k = false;
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "--polytope") {
        if (!need_value(i, argc, "--polytope")) return 1;
        input = argv[++i]; have_input = true;
      } else if (arg == "--max-k") {
        if (!need_value(i, argc, "--max-k")) return 1;
        options.max_k = std::stoi(argv[++i]); have_k = true;
      } else if (arg == "--output-dir") {
        if (!need_value(i, argc, "--output-dir")) return 1;
        output = argv[++i];
      } else if (arg == "--threshold") {
        if (!need_value(i, argc, "--threshold")) return 1;
        options.threshold = std::stod(argv[++i]);
      } else if (arg == "--tolerance") {
        if (!need_value(i, argc, "--tolerance")) return 1;
        options.tolerance = std::stod(argv[++i]);
      } else if (arg == "--absolute-tolerance") {
        if (!need_value(i, argc, "--absolute-tolerance")) return 1;
        options.absolute_tolerance = std::stod(argv[++i]);
      } else if (arg == "--max-iterations") {
        if (!need_value(i, argc, "--max-iterations")) return 1;
        options.max_iterations = std::stoi(argv[++i]);
      } else if (arg == "--max-points") {
        if (!need_value(i, argc, "--max-points")) return 1;
        options.max_points = static_cast<std::size_t>(std::stoull(argv[++i]));
      } else if (arg == "--verbose") {
        options.verbose = true;
      } else if (arg == "--help" || arg == "-h") {
        usage(); return 0;
      } else {
        std::cerr << "unknown option: " << arg << "\n";
        usage(); return 1;
      }
    }
    if (!have_input || !have_k) {
      usage(); return 1;
    }
    const auto result = chowk::run(input, options, output);
    const auto& last = result.levels.back();
    std::cout << "relative_K_status="
              << (result.unstable ? "unstable" : (result.failed ? "error" : "no_counterexample_found")) << '\n';
    std::cout << "levels=" << result.levels.size() << '\n';
    if (result.unstable) std::cout << "candidate_k=" << result.candidate_k << '\n';
    std::cout << std::setprecision(17) << "normalized_df=" << last.normalized_df << '\n';
    if (result.failed) return 1;
    return result.unstable ? 0 : 2;
  } catch (const std::exception& error) {
    std::cerr << "chowk: " << error.what() << '\n';
    return 1;
  }
}
