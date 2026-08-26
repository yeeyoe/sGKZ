#include "decompose.hpp"

#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>

namespace {

struct Arguments {
  std::optional<std::filesystem::path> polygon;
  kstab::DecomposeOptions options;
  std::int64_t certify_max_denom = 1048576;
};

void usage(std::ostream& output) {
  output << "Usage: k_stability_decompose --polygon FILE [options]\n\n"
         << "Search one-chord decompositions of a certified relatively K-unstable polygon.\n"
         << "The chord has dσ=0; inherited boundary segments retain their input measure.\n\n"
         << "Options:\n"
         << "  --theta-steps N            Witness directions (default 720).\n"
         << "  --t-steps N                Witness offsets per direction (default 512).\n"
         << "  --no-refine                Skip parent witness local refinement.\n"
         << "  --root-starts N            Root starts per edge pair (default 7).\n"
         << "  --root-iterations N        Newton iterations per start (default 60).\n"
         << "  --root-tolerance X         Continuity residual tolerance (default 1e-10).\n"
         << "  --rational-max-denom N     Rational reconstruction cap (default 1048576).\n"
         << "  --certify-max-denom N      Parent witness certification cap (default 1048576).\n"
         << "  --help                     Show this message.\n";
}

std::string value(int argc, char** argv, int& index, const std::string& option) {
  if (++index >= argc) throw std::invalid_argument("Missing value after " + option + ".");
  return argv[index];
}

Arguments parse(int argc, char** argv) {
  Arguments result;
  for (int i = 1; i < argc; ++i) {
    const std::string option = argv[i];
    if (option == "--help") { usage(std::cout); std::exit(0); }
    else if (option == "--polygon") result.polygon = value(argc, argv, i, option);
    else if (option == "--theta-steps") result.options.witness_search.theta_steps = std::stoi(value(argc, argv, i, option));
    else if (option == "--t-steps") result.options.witness_search.t_steps = std::stoi(value(argc, argv, i, option));
    else if (option == "--no-refine") result.options.witness_search.refine = false;
    else if (option == "--root-starts") result.options.root_starts = std::stoi(value(argc, argv, i, option));
    else if (option == "--root-iterations") result.options.root_iterations = std::stoi(value(argc, argv, i, option));
    else if (option == "--root-tolerance") result.options.root_tolerance = std::stod(value(argc, argv, i, option));
    else if (option == "--rational-max-denom") result.options.rational_max_denominator = std::stoll(value(argc, argv, i, option));
    else if (option == "--certify-max-denom") result.certify_max_denom = std::stoll(value(argc, argv, i, option));
    else throw std::invalid_argument("Unknown option: " + option);
  }
  if (!result.polygon) throw std::invalid_argument("Missing required --polygon FILE.");
  if (result.options.witness_search.theta_steps <= 0 || result.options.witness_search.t_steps <= 0 ||
      result.options.root_starts < 2 || result.options.root_iterations <= 0 ||
      result.options.root_tolerance <= 0 || result.options.rational_max_denominator < 1 ||
      result.certify_max_denom < 1) {
    throw std::invalid_argument("All counts, tolerances, and denominator limits must be positive.");
  }
  return result;
}

void print_point(const kstab::QPoint& point) {
  std::cout << '(' << kstab::rational_string(point.x) << ','
            << kstab::rational_string(point.y) << ')';
}

void print_ell(const std::array<kstab::Rational, 3>& ell) {
  std::cout << kstab::rational_string(ell[0]) << " + ("
            << kstab::rational_string(ell[1]) << ") x + ("
            << kstab::rational_string(ell[2]) << ") y";
}

void print_piece(const char* name, const kstab::MeasuredPolygon& polygon,
                 const std::array<kstab::Rational, 3>& ell,
                 const kstab::MeasuredSearchResult& search,
                 bool exact) {
  std::cout << name << "_vertices=";
  for (const auto& point : polygon.vertices_ccw) {
    if (exact) print_point(point);
    else std::cout << '(' << kstab::rational_to_double(point.x) << ','
                   << kstab::rational_to_double(point.y) << ')';
    std::cout << ' ';
  }
  std::cout << '\n' << name << "_edge_measures=";
  for (const auto& measure : polygon.edge_measures) {
    if (exact) std::cout << kstab::rational_string(measure);
    else std::cout << kstab::rational_to_double(measure);
    std::cout << ' ';
  }
  std::cout << '\n' << name << "_ell=";
  if (exact) print_ell(ell);
  else std::cout << kstab::rational_to_double(ell[0]) << " + ("
                 << kstab::rational_to_double(ell[1]) << ") x + ("
                 << kstab::rational_to_double(ell[2]) << ") y";
  std::cout << '\n' << name << "_relative_K_status="
            << (search.unstable ? "unstable" : "no_counterexample_found") << '\n';
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse(argc, argv);
    const kstab::PolygonInput input = kstab::parse_polygon_input_file(*arguments.polygon);
    const auto parent_ell = kstab::compute_ell_p(input.vertices_ccw, input.null_measure_edges);
    const kstab::SearchResult parent_search =
        kstab::search_witness(input.vertices_ccw, parent_ell, arguments.options.witness_search,
                              input.null_measure_edges);
    std::cout << std::setprecision(17);
    std::cout << "parent_relative_K_status="
              << (parent_search.unstable ? "unstable" : "no_counterexample_found") << '\n';
    if (!parent_search.unstable) {
      std::cerr << "error: parent polygon has no numerical negative witness.\n";
      return 2;
    }
    const kstab::CertifyResult certificate = kstab::certify_witness(
        input.vertices_ccw, parent_ell, parent_search.witness,
        arguments.certify_max_denom, input.null_measure_edges);
    std::cout << "parent_certified=" << std::boolalpha << certificate.certified << '\n';
    if (!certificate.certified) {
      std::cerr << "error: parent negative witness could not be certified.\n";
      return 2;
    }
    std::cout << "parent_certified_M_l=" << kstab::rational_string(certificate.value) << '\n';
    const kstab::MeasuredPolygon parent =
        kstab::make_measured_polygon(input.vertices_ccw, input.null_measure_edges);
    const auto candidates = kstab::find_chord_decompositions(parent, arguments.options);
    std::cout << "decomposition_candidates=" << candidates.size() << '\n';
    for (std::size_t index = 0; index < candidates.size(); ++index) {
      const auto& candidate = candidates[index];
      std::cout << "candidate=" << index << '\n'
                << "candidate_edges=" << candidate.first_edge << ',' << candidate.second_edge << '\n'
                << "candidate_parameters=" << candidate.first_parameter << ',' << candidate.second_parameter << '\n'
                << "candidate_continuity_residual=" << candidate.continuity_residual << '\n'
                << "candidate_concave=" << candidate.concave << '\n'
                << "candidate_certified_rational=" << candidate.certified_rational << '\n';
      std::cout << "candidate_endpoints=";
      if (candidate.certified_rational) { print_point(candidate.first_endpoint); std::cout << ','; print_point(candidate.second_endpoint); }
      else std::cout << '(' << kstab::rational_to_double(candidate.first_endpoint.x) << ',' << kstab::rational_to_double(candidate.first_endpoint.y) << "),("
                     << kstab::rational_to_double(candidate.second_endpoint.x) << ',' << kstab::rational_to_double(candidate.second_endpoint.y) << ')';
      std::cout << '\n';
      if (candidate.certified_rational) {
        std::cout << "candidate_parameters_rational="
                  << kstab::rational_string(*candidate.first_parameter_rational) << ','
                  << kstab::rational_string(*candidate.second_parameter_rational) << '\n';
      }
      print_piece("first_piece", candidate.first_piece, candidate.first_ell, candidate.first_search, candidate.certified_rational);
      print_piece("second_piece", candidate.second_piece, candidate.second_ell, candidate.second_search, candidate.certified_rational);
    }
    return candidates.empty() ? 2 : 0;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
