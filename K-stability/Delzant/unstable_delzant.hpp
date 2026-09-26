#pragma once

#include "../k_stability.hpp"

#include <cstdint>
#include <filesystem>
#include <random>
#include <string>
#include <vector>

namespace delzant {

struct Direction {
  std::int64_t x = 0;
  std::int64_t y = 0;
};

struct Candidate {
  std::vector<Direction> fan;
  std::vector<std::int64_t> edge_lengths;
  std::vector<kstab::IntPoint> vertices;
  std::array<kstab::Rational, 3> ell{};
  std::int64_t twice_area = 0;
  std::string key;
};

// Makes one smooth complete fan with d rays and coordinates bounded by N.
bool generate_fan(int d, int N, std::mt19937_64& rng,
                  std::vector<Direction>& fan);

// Samples positive integer lengths, closes the polygon exactly, and normalizes
// its translation and common edge-length factor.
bool build_candidate(const std::vector<Direction>& fan, int M,
                     std::mt19937_64& rng, Candidate& candidate,
                     std::int64_t max_diameter, bool* exceeded_max_diameter = nullptr);

bool is_smooth_complete_fan(const std::vector<Direction>& fan);
bool is_smooth_lattice_polygon(const Candidate& candidate);
std::string canonical_candidate_key(const std::vector<Direction>& fan,
                                    const std::vector<std::int64_t>& edge_lengths);
bool should_audit_candidate(const std::string& canonical_key);

struct SearchOptions {
  int d = 0;
  int N = 12;
  int M = 12;
  std::int64_t max_diameter = 0;
  std::uint64_t seed = 1;
  double time_limit_seconds = 600.0;
  std::int64_t certify_max_denominator = 1048576;
  std::filesystem::path database = "K-stability/Delzant/unstable_delzant.sqlite";
  std::filesystem::path output_directory = "K-stability/Delzant/results";
  int theta_steps = 360;
  int t_steps = 256;
};

struct SearchSummary {
  std::uint64_t generated = 0;
  std::uint64_t oversized = 0;
  std::uint64_t duplicates = 0;
  std::uint64_t tested = 0;
  std::uint64_t probed = 0;
  std::uint64_t confirmed = 0;
  std::uint64_t finalized = 0;
  std::uint64_t df_evaluations = 0;
  std::uint64_t certified_candidates = 0;
  double probe_seconds = 0.0;
  double confirm_seconds = 0.0;
  double final_seconds = 0.0;
  double generation_seconds = 0.0;
  double total_seconds = 0.0;
  bool found = false;
  Candidate candidate;
  kstab::CertifyResult certificate;
};

SearchSummary run_search(const SearchOptions& options);

}  // namespace delzant
