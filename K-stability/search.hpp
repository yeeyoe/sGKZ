#pragma once

#include "k_stability.hpp"

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace kstab {

struct Direction {
  std::int64_t x = 0;
  std::int64_t y = 0;
};

struct PolygonCandidate {
  int d = 0;
  std::vector<Direction> directions;
  std::vector<std::int64_t> steps;
  std::vector<IntPoint> vertices;
  std::vector<Direction> facet_normals;
  std::array<Rational, 3> ell{};
  std::int64_t twice_area = 0;
  std::string key;
  double probe_score = 0.0;
  std::vector<int> vertex_singularity_flags;
  int singular_vertex_count = -1;
};

enum class DetectorTier { probe, confirm, final };

struct DetectorProfile {
  DetectorTier tier = DetectorTier::probe;
  int integer_normal_bound = 4;
  std::int64_t certification_cap = 1048576;
  std::string revision = "area-search-detector-v2";
  std::string fingerprint() const;
};

struct DetectorOutcome {
  bool verified_unstable = false;
  bool numerical_negative = false;
  SearchResult numerical;
  CertifyResult certification;
  std::string profile;
};

struct ValidationRecord {
  std::string profile;
  std::string status;
  std::string last_stage;
};

struct AttemptRecord {
  std::string status;
  double value = 0.0;
  bool numerical_negative = false;
  bool has_exact_witness = false;
  std::string q_squared_exact;
  double q_squared_value = 0.0;
  bool has_q_squared = false;
};

// Exact Q_P(g)^2 for an authenticated simple convex witness.
Rational q_squared_exact(const PolygonCandidate& candidate,
                         const std::array<Rational, 3>& coefficients,
                         const Rational& exact_value);

struct BackfillSummary {
  std::uint64_t scanned = 0;
  std::uint64_t backfilled = 0;
  std::uint64_t skipped = 0;
  std::uint64_t errors = 0;
};

struct VerifiedCandidateSummary {
  std::string key;
  std::int64_t twice_area = 0;
  std::string q_squared_exact;
  double q_squared_value = 0.0;
};

struct AreaSearchOptions {
  int d = 6;
  int initial_N = 4;
  int initial_M = 4;
  int beam_width = 48;
  std::uint64_t seed = 1;
  double time_limit_seconds = 60.0;
  double shell_seconds = 60.0;
  bool stop_on_first = false;
  bool smooth_only = false;
  bool verbose = false;
  std::filesystem::path database = "K-stability/k_stability_search.sqlite";
  std::filesystem::path output_directory = ".";
  std::int64_t certify_max_denominator = 1048576;
};

struct SearchSummary {
  std::uint64_t generated = 0;
  std::uint64_t rejected = 0;
  std::uint64_t probes = 0;
  std::uint64_t confirms = 0;
  std::uint64_t finals = 0;
  std::uint64_t verified = 0;
  std::uint64_t unverified = 0;
  std::uint64_t skipped = 0;
  std::uint64_t total_tested = 0;
  std::uint64_t new_tested = 0;
  std::uint64_t total_verified_unstable = 0;
  std::uint64_t new_verified_unstable = 0;
  bool smaller_volume_found = false;
  std::vector<VerifiedCandidateSummary> top_verified;
  std::string first_verified_key;
  std::string best_verified_key;
  std::int64_t best_twice_area = 0;
  bool have_verified = false;
};

// Build one strictly convex lattice candidate from the first d-1 facets.
// The closing direction and step are derived and are deliberately unbounded.
bool build_candidate(int d, const std::vector<Direction>& first_directions,
                     const std::vector<std::int64_t>& first_steps,
                     PolygonCandidate& result, std::string* reason = nullptr);

// Convert a normalized CCW polygon to the canonical candidate representation
// used by the area search. Returns false when the polygon is not in the
// searcher's origin-based canonical form.
bool candidate_from_vertices(const std::vector<IntPoint>& vertices,
                             PolygonCandidate& result,
                             std::string* reason = nullptr);

std::vector<int> compute_vertex_singularity_flags(
    const PolygonCandidate& candidate);

std::vector<Direction> primitive_directions(int bound);
std::string candidate_key(const PolygonCandidate& candidate);
std::string current_validation_profile(const AreaSearchOptions& options);

DetectorOutcome detect_candidate(const PolygonCandidate& candidate,
                                 const DetectorProfile& profile,
                                 std::int64_t max_denominator);

class SearchDatabase {
 public:
  explicit SearchDatabase(const std::filesystem::path& path);
  ~SearchDatabase();
  SearchDatabase(const SearchDatabase&) = delete;
  SearchDatabase& operator=(const SearchDatabase&) = delete;

  void save_candidate(const PolygonCandidate& candidate, const std::string& status);
  void update_probe_score(const std::string& key, double score);
  void ensure_generator_revision(const std::string& revision);
  void save_vertex_singularity(const PolygonCandidate& candidate);
  bool has_attempt(const std::string& key, const std::string& profile) const;
  void save_attempt(const PolygonCandidate& candidate, const DetectorOutcome& outcome);
  BackfillSummary backfill_q_squared();
  std::optional<AttemptRecord> get_attempt(
      const std::string& candidate_key, const std::string& profile) const;
  std::optional<ValidationRecord> get_validation(
      const std::string& candidate_key, const std::string& profile) const;
  bool has_verified_candidate(const std::string& candidate_key) const;
  void save_validation(const PolygonCandidate& candidate,
                       const std::string& profile,
                       const std::string& status,
                       const std::string& last_stage);
  // If dimension is nonnegative, load only candidates with that d. The
  // default keeps the administrative all-dimensions view for database tools.
  std::vector<PolygonCandidate> load_candidates(int dimension = -1) const;
  bool candidate_is_verified(const std::string& key) const;
  void save_state(const std::string& name, const std::string& value);
  std::string load_state(const std::string& name) const;
  std::uint64_t count_candidates(int dimension = -1) const;
  std::uint64_t count_status(const std::string& status, int dimension = -1) const;
  std::uint64_t count_tested(int dimension) const;
  std::uint64_t count_verified(int dimension) const;
  std::vector<VerifiedCandidateSummary> top_verified(int dimension,
                                                      std::size_t limit) const;
  std::optional<std::int64_t> min_verified_twice_area(int dimension) const;
  void write_report(const std::filesystem::path& path,
                    const SearchSummary& summary, int dimension = -1) const;

 private:
  struct Impl;
  Impl* impl_ = nullptr;
};

SearchSummary run_search(const AreaSearchOptions& options);

}  // namespace kstab
