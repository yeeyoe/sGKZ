#include "search.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

void test_geometry() {
  kstab::PolygonCandidate candidate;
  const bool ok = kstab::build_candidate(
      4, {{1, 0}, {0, 1}, {-1, 1}}, {2, 2, 2}, candidate);
  require(ok, "valid convex candidate should build");
  require(candidate.steps[0] == 1 && candidate.steps[1] == 1 &&
              candidate.steps[2] == 1,
          "common step gcd should be normalized");
  require(candidate.vertices.front().x == 0 && candidate.vertices.front().y == 0,
          "candidate must start at origin");
  for (const auto& vertex : candidate.vertices)
    require(vertex.y >= 0, "candidate vertices must be integral and upper-facing");
  require(candidate.vertices.back().y > 0, "last vertex must be above origin");
  require(candidate.directions.size() == 4 && candidate.steps.size() == 4,
          "closing direction and step must be recorded");
  require(candidate.twice_area == 3, "exact twice area should be recorded");
  require(candidate.facet_normals.size() == 4, "all facet normals must be recorded");
}

void test_rejection() {
  kstab::PolygonCandidate candidate;
  require(!kstab::build_candidate(4, {{1, 0}, {0, 1}, {1, 0}}, {1, 1, 1}, candidate),
          "repeated direction must be rejected");
}

void test_profile_and_database() {
  const auto path = std::filesystem::temp_directory_path() / "kstab_search_test.sqlite";
  std::filesystem::remove(path);
  kstab::PolygonCandidate candidate;
  require(kstab::build_candidate(4, {{1, 0}, {0, 1}, {-1, 1}}, {1, 1, 1}, candidate),
          "database fixture candidate");
  kstab::SearchDatabase database(path);
  database.save_candidate(candidate, "unverified");
  kstab::DetectorOutcome outcome;
  outcome.profile = kstab::DetectorProfile{kstab::DetectorTier::probe, 4}.fingerprint();
  database.save_attempt(candidate, outcome);
  require(database.has_attempt(candidate.key, outcome.profile), "same profile must be skipped");
  require(!database.has_attempt(candidate.key,
                               kstab::DetectorProfile{kstab::DetectorTier::confirm, 8}.fingerprint()),
          "changed profile must be eligible");
  kstab::DetectorProfile changed_cap{kstab::DetectorTier::probe, 4};
  changed_cap.certification_cap = 17;
  require(!database.has_attempt(candidate.key, changed_cap.fingerprint()),
          "changed certification cap must be eligible");
  require(database.load_candidates().size() == 1, "candidate should survive SQLite round trip");
  const auto report = std::filesystem::temp_directory_path() / "kstab_search_report.txt";
  kstab::SearchSummary none;
  database.write_report(report, none);
  std::ifstream no_result(report);
  std::stringstream no_result_text;
  no_result_text << no_result.rdbuf();
  require(no_result_text.str() == "没找到\n", "empty search report should say not found");

  kstab::DetectorOutcome verified;
  verified.profile = kstab::DetectorProfile{kstab::DetectorTier::confirm, 8}.fingerprint();
  verified.verified_unstable = true;
  verified.certification.certified = true;
  verified.certification.coefficients = {kstab::Rational(1), kstab::Rational(0),
                                         kstab::Rational(0)};
  verified.certification.value = kstab::Rational(-1);
  database.save_attempt(candidate, verified);
  kstab::SearchSummary found;
  found.have_verified = true;
  found.best_twice_area = candidate.twice_area;
  found.best_verified_key = candidate.key;
  found.first_verified_key = candidate.key;
  database.write_report(report, found);
  std::ifstream found_result(report);
  std::stringstream found_result_text;
  found_result_text << found_result.rdbuf();
  require(found_result_text.str().find("vertices=0:0;1:0;1:1;0:2") !=
              std::string::npos,
          "verified report should include candidate vertices");
  require(found_result_text.str().find("attempt_exact_value=-1") !=
              std::string::npos,
          "verified report should include exact witness");
  std::filesystem::remove(report);
  std::filesystem::remove(path);
}

}  // namespace

int main() {
  try {
    test_geometry();
    test_rejection();
    test_profile_and_database();
    std::cout << "All search tests passed.\n";
  } catch (const std::exception& error) {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
