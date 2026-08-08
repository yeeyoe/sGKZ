#include "search.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <sqlite3.h>

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
  require(!kstab::build_candidate(
              6, {{1, 0}, {1, 3}, {-1, -2}, {1, 1}, {1, 4}},
              {2, 4, 3, 2, 2}, candidate),
          "interior vertex must be rejected by the closing wedge condition");
}

void test_vertex_singularity_flags() {
  kstab::PolygonCandidate smooth;
  require(kstab::build_candidate(4, {{1, 0}, {0, 1}, {-1, 0}}, {1, 1, 1}, smooth),
          "smooth square candidate should build");
  const auto smooth_flags = kstab::compute_vertex_singularity_flags(smooth);
  require(smooth_flags == std::vector<int>({0, 0, 0, 0}),
          "unit square should be smooth at every vertex");

  kstab::PolygonCandidate singular;
  require(kstab::build_candidate(3, {{1, 0}, {-1, 1}}, {2, 1}, singular),
          "singular triangle should build");
  const auto singular_flags = kstab::compute_vertex_singularity_flags(singular);
  require(singular_flags == std::vector<int>({0, 0, 1}) &&
              std::accumulate(singular_flags.begin(), singular_flags.end(), 0) == 1,
          "singular triangle should have one singular vertex");

  const auto path = std::filesystem::temp_directory_path() / "kstab_singularity_test.sqlite";
  std::filesystem::remove(path);
  kstab::SearchDatabase database(path);
  database.save_candidate(singular, "pending");
  database.save_validation(singular, "test-profile", "verified_unstable", "final");
  const auto loaded = database.load_candidates();
  require(loaded.size() == 1 && loaded.front().vertex_singularity_flags == singular_flags &&
              loaded.front().singular_vertex_count == 1,
          "verified singularity metadata should persist");
  std::filesystem::remove(path);
  std::filesystem::remove(path.string() + "-wal");
  std::filesystem::remove(path.string() + "-shm");
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

  kstab::AreaSearchOptions profile_options;
  const std::string full_profile =
      kstab::current_validation_profile(profile_options);
  database.save_validation(candidate, full_profile, "unverified", "complete");
  const auto unverified = database.get_validation(candidate.key, full_profile);
  require(unverified && unverified->status == "unverified" &&
              unverified->last_stage == "complete",
          "candidate-level unverified validation should round trip");
  require(!database.has_verified_candidate(candidate.key),
          "unverified validation must not be considered verified");
  auto changed_options = profile_options;
  changed_options.certify_max_denominator += 1;
  const std::string changed_profile =
      kstab::current_validation_profile(changed_options);
  require(changed_profile != full_profile &&
              !database.get_validation(candidate.key, changed_profile),
          "changed complete profile must be eligible for revalidation");

  kstab::DetectorOutcome pending_attempt;
  pending_attempt.profile = full_profile + "|stage=confirm";
  pending_attempt.numerical_negative = true;
  pending_attempt.numerical.witness.value = -0.25;
  database.save_validation(candidate, full_profile, "pending", "confirm");
  database.save_attempt(candidate, pending_attempt);
  const auto attempt = database.get_attempt(candidate.key, pending_attempt.profile);
  require(attempt && attempt->status == "unverified" &&
              attempt->numerical_negative && attempt->value == -0.25,
          "pending stage attempt should be resumable without re-running it");

  database.save_validation(candidate, full_profile, "verified_unstable", "confirm");
  require(database.has_verified_candidate(candidate.key),
          "verified validation must be skipped permanently");
  require(database.get_validation(candidate.key, full_profile)->status ==
              "verified_unstable",
          "verified validation status should be retained");
  const auto loaded_verified = database.load_candidates();
  require(loaded_verified.size() == 1 && loaded_verified.front().singular_vertex_count == 0 &&
              loaded_verified.front().vertex_singularity_flags == std::vector<int>({0, 0, 0, 0}),
          "verified candidate singularity metadata should round trip");

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
  require(found_result_text.str().find("vertex_singularity_flags=0,0,0,0") !=
              std::string::npos &&
              found_result_text.str().find("singular_vertex_count=0") !=
                  std::string::npos,
          "verified report should include vertex singularity metadata");
  std::filesystem::remove(report);
  std::filesystem::remove(path);
}

void test_large_rational_load() {
  const auto path = std::filesystem::temp_directory_path() / "kstab_large_rational.sqlite";
  std::filesystem::remove(path);
  std::filesystem::remove(path.string() + "-wal");
  std::filesystem::remove(path.string() + "-shm");

  kstab::PolygonCandidate candidate;
  require(kstab::build_candidate(4, {{1, 0}, {0, 1}, {-1, 1}}, {1, 1, 1}, candidate),
          "large-rational fixture candidate");
  {
    kstab::SearchDatabase database(path);
    database.save_candidate(candidate, "pending");
  }

  sqlite3* db = nullptr;
  require(sqlite3_open(path.string().c_str(), &db) == SQLITE_OK,
          "open large-rational fixture database");
  const char* sql =
      "UPDATE candidates SET ell0='4335276737992173/15688438743253116926'";
  require(sqlite3_exec(db, sql, nullptr, nullptr, nullptr) == SQLITE_OK,
          "write large-rational fixture");
  sqlite3_close(db);

  kstab::SearchDatabase database(path);
  const auto loaded = database.load_candidates();
  require(loaded.size() == 1, "large-rational candidate should load");
  std::ostringstream coefficient;
  coefficient << loaded.front().ell[0];
  require(coefficient.str() == "4335276737992173/15688438743253116926",
          "large-rational ell coefficient must round trip exactly");

  std::filesystem::remove(path);
  std::filesystem::remove(path.string() + "-wal");
  std::filesystem::remove(path.string() + "-shm");
}

void test_dimension_scoping() {
  const auto path = std::filesystem::temp_directory_path() / "kstab_dimension_scope.sqlite";
  std::filesystem::remove(path);
  std::filesystem::remove(path.string() + "-wal");
  std::filesystem::remove(path.string() + "-shm");

  kstab::PolygonCandidate d4;
  kstab::PolygonCandidate d3;
  require(kstab::build_candidate(4, {{1, 0}, {0, 1}, {-1, 1}}, {1, 1, 1}, d4),
          "d4 dimension fixture");
  require(kstab::build_candidate(3, {{1, 0}, {-1, 1}}, {2, 1}, d3),
          "d3 dimension fixture");
  {
    kstab::SearchDatabase database(path);
    database.save_candidate(d4, "unverified");
    database.save_candidate(d3, "pending");
    kstab::DetectorOutcome d4_probe;
    d4_probe.profile = "dimension-test|probe";
    kstab::DetectorOutcome d4_confirm;
    d4_confirm.profile = "dimension-test|confirm";
    kstab::DetectorOutcome d3_probe;
    d3_probe.profile = "dimension-test|probe";
    database.save_attempt(d4, d4_probe);
    database.save_attempt(d4, d4_confirm);
    database.save_attempt(d3, d3_probe);
    database.save_state("d4|shell", "0");
    database.save_state("d4|rng", "");

    require(database.count_candidates() == 2, "all-dimension count should include both fixtures");
    require(database.count_candidates(3) == 1 && database.count_candidates(4) == 1,
            "dimension-scoped candidate counts should be isolated");
    require(database.count_status("unverified", 3) == 0 &&
                database.count_status("unverified", 4) == 1,
            "dimension-scoped status counts should be isolated");
    require(database.count_tested(3) == 1 && database.count_tested(4) == 1,
            "tested count should deduplicate stages and isolate dimensions");
    require(database.load_candidates(3).size() == 1 && database.load_candidates(3).front().d == 3,
            "d3 load must exclude d4 candidates");
    require(database.load_candidates(4).size() == 1 && database.load_candidates(4).front().d == 4,
            "d4 load must exclude d3 candidates");
    require(database.load_state("d4|shell") == "0" &&
                database.load_state("d4|rng").empty(),
            "dimension-scoped state should round trip");

    database.save_state("generator_revision", "incremental-convex-v1");
  }

  kstab::AreaSearchOptions options;
  options.d = 4;
  options.initial_N = 1;
  options.initial_M = 1;
  options.beam_width = 1;
  options.time_limit_seconds = 0.2;
  options.shell_seconds = 0.2;
  options.database = path;
  options.output_directory = std::filesystem::temp_directory_path() /
                             "kstab_dimension_scope_output";
  std::filesystem::remove_all(options.output_directory);
  kstab::run_search(options);

  kstab::SearchDatabase checked(path);
  const std::string probe_profile =
      kstab::current_validation_profile(options) + "|stage=probe";
  require(!checked.has_attempt(d3.key, probe_profile),
          "d4 search must not probe a d3 candidate");
  std::filesystem::remove_all(options.output_directory);

  std::filesystem::remove(path);
  std::filesystem::remove(path.string() + "-wal");
  std::filesystem::remove(path.string() + "-shm");
}

void test_search_statistics() {
  const auto path = std::filesystem::temp_directory_path() / "kstab_search_statistics.sqlite";
  const auto output = std::filesystem::temp_directory_path() / "kstab_search_statistics_output";
  std::filesystem::remove(path);
  std::filesystem::remove(path.string() + "-wal");
  std::filesystem::remove(path.string() + "-shm");
  std::filesystem::remove_all(output);

  kstab::PolygonCandidate candidate;
  require(kstab::build_candidate(4, {{1, 0}, {0, 1}, {-1, 0}}, {1, 1, 1}, candidate),
          "statistics unit-square fixture");
  candidate.ell = {kstab::Rational(8), kstab::Rational(0), kstab::Rational(0)};
  {
    kstab::SearchDatabase database(path);
    database.save_candidate(candidate, "pending");
    database.save_state("generator_revision", "incremental-convex-v1");
  }

  kstab::AreaSearchOptions options;
  options.d = 4;
  options.initial_N = 1;
  options.initial_M = 1;
  options.beam_width = 1;
  options.time_limit_seconds = 0.5;
  options.shell_seconds = 0.5;
  options.stop_on_first = true;
  options.database = path;
  options.output_directory = output;

  const auto first = kstab::run_search(options);
  require(first.total_tested == 1 && first.new_tested == 1,
          "first statistics run should test one unique candidate");
  require(first.total_verified_unstable == 1 && first.new_verified_unstable == 1,
          "first statistics run should add one verified candidate");
  require(first.smaller_volume_found && first.top_verified.size() == 1 &&
              first.top_verified.front().key == candidate.key &&
              first.top_verified.front().twice_area == candidate.twice_area,
          "first statistics run should report the new least-volume candidate");

  const auto second = kstab::run_search(options);
  require(second.total_tested == 1 && second.new_tested == 0 &&
              second.total_verified_unstable == 1 && second.new_verified_unstable == 0,
          "restart statistics run should skip the existing verified candidate");
  require(!second.smaller_volume_found && second.top_verified.size() == 1,
          "restart statistics run should keep the same least volume");

  kstab::PolygonCandidate larger_a;
  kstab::PolygonCandidate larger_b;
  require(kstab::build_candidate(4, {{1, 0}, {0, 1}, {-1, 2}}, {1, 1, 1}, larger_a),
          "top-five area fixture A");
  require(kstab::build_candidate(4, {{1, 0}, {0, 1}, {-1, 3}}, {1, 1, 1}, larger_b),
          "top-five area fixture B");
  {
    kstab::SearchDatabase database(path);
    const auto save_verified = [&](const kstab::PolygonCandidate& fixture) {
      database.save_candidate(fixture, "pending");
      kstab::DetectorOutcome outcome;
      outcome.profile = "statistics|final";
      outcome.verified_unstable = true;
      outcome.certification.certified = true;
      outcome.certification.coefficients = {kstab::Rational(1), kstab::Rational(0),
                                             kstab::Rational(0)};
      outcome.certification.value = kstab::Rational(-1);
      database.save_attempt(fixture, outcome);
    };
    save_verified(larger_a);
    save_verified(larger_b);
    const auto top = database.top_verified(4, 5);
    require(top.size() == 3 && top[0].key == candidate.key &&
                top[0].twice_area < top[1].twice_area &&
                top[1].twice_area < top[2].twice_area,
            "verified candidates should be ranked by twice_area");
  }

  std::filesystem::remove_all(output);
  std::filesystem::remove(path);
  std::filesystem::remove(path.string() + "-wal");
  std::filesystem::remove(path.string() + "-shm");
}

void test_smooth_only_search() {
  const auto path = std::filesystem::temp_directory_path() / "kstab_smooth_only.sqlite";
  const auto output = std::filesystem::temp_directory_path() / "kstab_smooth_only_output";
  std::filesystem::remove(path);
  std::filesystem::remove(path.string() + "-wal");
  std::filesystem::remove(path.string() + "-shm");
  std::filesystem::remove_all(output);

  kstab::AreaSearchOptions options;
  options.d = 4;
  options.initial_N = 3;
  options.initial_M = 3;
  options.beam_width = 8;
  options.time_limit_seconds = 0.2;
  options.shell_seconds = 0.2;
  options.smooth_only = true;
  options.database = path;
  options.output_directory = output;

  kstab::PolygonCandidate singular;
  require(kstab::build_candidate(4, {{1, 0}, {0, 1}, {-1, 1}}, {1, 1, 2}, singular),
          "singular candidate should build for smooth-only filtering");
  {
    kstab::SearchDatabase seed_database(path);
    seed_database.save_candidate(singular, "pending");
    seed_database.save_state("generator_revision", "incremental-convex-v1");
  }
  kstab::run_search(options);

  kstab::SearchDatabase database(path);
  const auto candidates = database.load_candidates();
  require(!candidates.empty(), "smooth-only smoke search should generate candidates");
  bool checked_generated = false;
  for (const auto& candidate : candidates) {
    if (candidate.key == singular.key) continue;
    checked_generated = true;
    const auto flags = kstab::compute_vertex_singularity_flags(candidate);
    require(std::all_of(flags.begin(), flags.end(), [](int flag) { return flag == 0; }),
            "smooth-only search must persist only smooth generated candidates");
  }
  require(checked_generated, "smooth-only smoke search should add a generated candidate");
  const std::string probe_profile =
      kstab::current_validation_profile(options) + "|stage=probe";
  require(!database.has_attempt(singular.key, probe_profile),
          "smooth-only search must skip a loaded singular candidate");
  std::filesystem::remove_all(output);
  std::filesystem::remove(path);
  std::filesystem::remove(path.string() + "-wal");
  std::filesystem::remove(path.string() + "-shm");
}

}  // namespace

int main() {
  try {
    test_geometry();
    test_rejection();
    test_vertex_singularity_flags();
    test_profile_and_database();
    test_large_rational_load();
    test_dimension_scoping();
    test_search_statistics();
    test_smooth_only_search();
    std::cout << "All search tests passed.\n";
  } catch (const std::exception& error) {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
