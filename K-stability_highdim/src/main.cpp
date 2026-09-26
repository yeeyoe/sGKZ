#include "k_stability_highdim.hpp"

#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace {

struct Arguments {
  std::optional<std::filesystem::path> polytope;
  std::optional<std::filesystem::path> check_pl;
  bool ell_only = false;
  std::optional<int> pieces;
  kstab_highdim::SearchOptions search;
  std::int64_t certify_max_denominator = 1048576;
};

void usage(std::ostream& stream) {
  stream << "Usage:\n"
         << "  k_stability_highdim --polytope FILE --ell-only\n"
         << "  k_stability_highdim --polytope FILE --pieces M [options]\n"
         << "  k_stability_highdim --polytope FILE --check-pl FILE\n\n"
         << "Options:\n"
         << "  --pieces M              Number of nonzero affine branches (1..8).\n"
         << "  --population N          Differential-evolution population.\n"
         << "  --generations N         Differential-evolution generations (default 80).\n"
         << "  --quadrature-samples N  Halton samples per integral (default 20000).\n"
         << "  --seed N                Deterministic search seed (default 0).\n"
         << "  --threads N             Search worker threads (default 1).\n"
         << "  --certify-max-denom N   Rational certification cap (default 1048576).\n"
         << "  --verbose               Print search summary to stderr.\n"
         << "  --help                  Show this message.\n";
}

std::string require_value(int argc, char** argv, int& index, const std::string& option) {
  if (index + 1 >= argc) throw std::invalid_argument("missing value after " + option);
  return argv[++index];
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments result;
  for (int i = 1; i < argc; ++i) {
    const std::string option = argv[i];
    if (option == "--help") {
      usage(std::cout);
      std::exit(0);
    } else if (option == "--polytope") {
      result.polytope = require_value(argc, argv, i, option);
    } else if (option == "--ell-only") {
      result.ell_only = true;
    } else if (option == "--pieces") {
      result.pieces = std::stoi(require_value(argc, argv, i, option));
    } else if (option == "--check-pl") {
      result.check_pl = require_value(argc, argv, i, option);
    } else if (option == "--population") {
      result.search.population = std::stoull(require_value(argc, argv, i, option));
    } else if (option == "--generations") {
      result.search.generations = std::stoull(require_value(argc, argv, i, option));
    } else if (option == "--quadrature-samples") {
      result.search.quadrature_samples = std::stoull(require_value(argc, argv, i, option));
    } else if (option == "--seed") {
      result.search.seed = std::stoull(require_value(argc, argv, i, option));
    } else if (option == "--threads") {
      result.search.threads = static_cast<unsigned>(std::stoul(require_value(argc, argv, i, option)));
    } else if (option == "--certify-max-denom") {
      result.certify_max_denominator = std::stoll(require_value(argc, argv, i, option));
    } else if (option == "--verbose") {
      result.search.verbose = true;
    } else {
      throw std::invalid_argument("unknown option: " + option);
    }
  }
  if (!result.polytope.has_value()) throw std::invalid_argument("missing --polytope FILE");
  const int modes = static_cast<int>(result.ell_only) + static_cast<int>(result.check_pl.has_value()) + static_cast<int>(result.pieces.has_value());
  if (modes != 1) throw std::invalid_argument("choose exactly one of --ell-only, --pieces, or --check-pl");
  if (result.pieces.has_value() && (*result.pieces < 1 || *result.pieces > 8)) throw std::invalid_argument("pieces must be in [1,8]");
  if (result.search.threads == 0 || result.search.generations == 0 || result.search.quadrature_samples == 0 || result.certify_max_denominator < 1) {
    throw std::invalid_argument("search counts and certification cap must be positive");
  }
  return result;
}

void print_ell(const kstab_highdim::LatticePolytope& polytope) {
  std::cout << "dimension=" << polytope.dimension << '\n'
            << "input_points=" << polytope.input_points.size() << '\n'
            << "hull_points=" << polytope.points.size() << '\n'
            << "redundant_points=" << (polytope.input_points.size() - polytope.points.size()) << '\n'
            << "simplices=" << polytope.simplices.size() << '\n'
            << "boundary_facets=" << polytope.boundary_facets.size() << '\n'
            << "volume=" << kstab_highdim::rational_string(polytope.volume) << '\n'
            << "boundary_measure=" << kstab_highdim::rational_string(polytope.boundary_measure) << '\n';
  std::cout << "ell_P=" << kstab_highdim::rational_string(polytope.ell[0]);
  for (int i = 0; i < polytope.dimension; ++i) {
    std::cout << " + (" << kstab_highdim::rational_string(polytope.ell[i + 1]) << ")x" << (i + 1);
  }
  std::cout << '\n';
}

}  // namespace

int main(int argc, char** argv) {
  Arguments arguments;
  try {
    arguments = parse_arguments(argc, argv);
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    usage(std::cerr);
    return 1;
  }
  try {
    const auto polytope = kstab_highdim::parse_polytope_file(*arguments.polytope);
    print_ell(polytope);
    if (arguments.ell_only) return 0;

    if (arguments.check_pl.has_value()) {
      const auto branches = kstab_highdim::parse_pl_file(*arguments.check_pl, polytope.dimension);
      kstab_highdim::ConvexPLFunction function{branches};
      const auto result = kstab_highdim::certify_function(polytope, function);
      std::cout << "check_pl_branches=" << branches.size() << '\n'
                << "check_pl_M_l=" << kstab_highdim::rational_string(result.value) << '\n'
                << "certified=" << std::boolalpha << result.certified << '\n'
                << "relative_K_status=" << (result.certified ? "unstable" : "no_counterexample_found") << '\n';
      return result.certified ? 0 : 2;
    }

    arguments.search.pieces = *arguments.pieces;
    const auto witness = kstab_highdim::search_pl_witness(polytope, arguments.search);
    const std::size_t parameter_dimension =
        static_cast<std::size_t>(arguments.search.pieces) *
            static_cast<std::size_t>(polytope.dimension) +
        static_cast<std::size_t>(arguments.search.pieces) - 1;
    const std::size_t effective_population =
        arguments.search.population == 0
            ? std::max<std::size_t>(32, 8 * parameter_dimension)
            : arguments.search.population;
    std::cout << std::setprecision(17)
              << "search_pieces=" << arguments.search.pieces << '\n'
              << "search_parameter_dimension=" << parameter_dimension << '\n'
              << "search_population=" << effective_population << '\n'
              << "search_generations=" << arguments.search.generations << '\n'
              << "search_quadrature_samples=" << arguments.search.quadrature_samples << '\n'
              << "search_seed=" << arguments.search.seed << '\n'
              << "search_threads=" << arguments.search.threads << '\n'
              << "certify_max_denom=" << arguments.certify_max_denominator << '\n'
              << "search_evaluations=" << witness.evaluations << '\n'
              << "sweep_min_M_l=" << witness.value << '\n'
              << "sweep_min_M_l_normalized=" << witness.normalized << '\n';
    for (std::size_t i = 0; i < witness.function.branches.size(); ++i) {
      std::cout << "witness_L" << (i + 1) << '=' << kstab_highdim::format_affine(witness.function.branches[i]) << '\n';
    }
    if (witness.normalized >= -1e-6) {
      std::cout << "relative_K_status=no_counterexample_found\n";
      return 2;
    }
    const auto certification = kstab_highdim::certify_witness(polytope, witness, arguments.certify_max_denominator);
    std::cout << "certified=" << std::boolalpha << certification.certified << '\n'
              << "certified_M_l=" << kstab_highdim::rational_string(certification.value) << '\n';
    if (certification.certified) {
      for (std::size_t i = 0; i < certification.function.branches.size(); ++i) {
        std::cout << "certified_L" << (i + 1) << '=' << kstab_highdim::format_affine(certification.function.branches[i]) << '\n';
      }
      std::cout << "relative_K_status=unstable\n";
      return 0;
    }
    std::cout << "relative_K_status=unverified_candidate\n";
    return 3;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
