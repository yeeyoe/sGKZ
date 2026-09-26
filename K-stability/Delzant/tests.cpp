#include "unstable_delzant.hpp"

#include <CGAL/Gmpz.h>
#include <sqlite3.h>

#include <algorithm>
#include <cmath>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

void test_fan_and_polygon_generation() {
  std::mt19937_64 rng(7);
  std::vector<delzant::Direction> fan;
  bool built = false;
  delzant::Candidate candidate;
  for (int attempt = 0; attempt < 10000 && !built; ++attempt) {
    require(delzant::generate_fan(8, 24, rng, fan),
            "bounded smooth fan generation failed");
    require(fan.size() == 8 && delzant::is_smooth_complete_fan(fan),
            "generated fan is not smooth and complete");
    built = delzant::build_candidate(fan, 8, rng, candidate, 1000000);
  }
  require(built, "could not close a generated smooth fan with positive lengths");
  require(delzant::is_smooth_lattice_polygon(candidate),
          "generated polygon failed the Delzant checks");
  require(candidate.twice_area > 0 && candidate.vertices.size() == 8,
          "generated polygon has invalid area or vertex count");
  for (std::size_t i = 0; i < fan.size(); ++i) {
    const auto& p = candidate.vertices[i];
    const auto& q = candidate.vertices[(i + 1) % fan.size()];
    require(static_cast<__int128_t>(q.x) - p.x ==
                static_cast<__int128_t>(candidate.edge_lengths[i]) * candidate.fan[i].x,
            "stored edge vector disagrees with its fan ray");
  }

  const std::vector<delzant::Direction> invalid{{1, 0}, {2, 1}, {0, 1}, {-1, 0}, {0, -1}};
  require(!delzant::is_smooth_complete_fan(invalid),
          "fan with a non-unimodular adjacent pair was accepted");
}

void test_scale_dedup_and_size_limit() {
  const std::vector<delzant::Direction> fan{{1, 0}, {0, 1}, {-1, 2}, {0, -1}};
  const std::vector<std::int64_t> lengths{1, 1, 1, 3};
  std::vector<std::int64_t> scaled{3, 3, 3, 9};
  require(delzant::canonical_candidate_key(fan, lengths) ==
              delzant::canonical_candidate_key(fan, scaled),
          "overall integer scaling produced a distinct candidate key");
  auto rotated_fan = fan;
  auto rotated_lengths = lengths;
  std::rotate(rotated_fan.begin(), rotated_fan.begin() + 2, rotated_fan.end());
  std::rotate(rotated_lengths.begin(), rotated_lengths.begin() + 2, rotated_lengths.end());
  require(delzant::canonical_candidate_key(fan, lengths) ==
              delzant::canonical_candidate_key(rotated_fan, rotated_lengths),
          "cyclic polygon representation produced a distinct candidate key");
  std::size_t audits = 0;
  for (int i = 0; i < 10000; ++i) {
    const std::string key = "candidate-" + std::to_string(i);
    const bool selected = delzant::should_audit_candidate(key);
    require(selected == delzant::should_audit_candidate(key),
            "audit selection is not reproducible");
    audits += selected;
  }
  require(audits > 350 && audits < 650,
          "stable-hash audit sample is not close to five percent");

  std::mt19937_64 rng(11);
  delzant::Candidate candidate;
  bool oversized = false;
  require(delzant::build_candidate(fan, 1, rng, candidate, 3, &oversized),
          "polygon exactly at max-diameter was rejected");
  require(!oversized && candidate.ell[0] != 0,
          "boundary-size candidate did not finish construction");
  rng.seed(11);
  require(!delzant::build_candidate(fan, 1, rng, candidate, 2, &oversized) &&
              oversized,
          "polygon above max-diameter was not rejected early");
}

void test_exact_relative_df_certificate() {
  // The existing d5_a74 witness is an exact negative DF regression fixture.
  std::vector<kstab::IntPoint> vertices{{0, 0}, {1, 0}, {-17, 21},
                                        {-15, 16}, {-10, 10}};
  vertices = kstab::normalize_polygon(std::move(vertices));
  const auto ell = kstab::compute_ell_p(vertices);
  const kstab::Rational a(CGAL::Gmpz(-5), CGAL::Gmpz(6));
  const kstab::Rational b(CGAL::Gmpz(-5), CGAL::Gmpz(9));
  const kstab::Rational c(CGAL::Gmpz(-7), CGAL::Gmpz(3));
  const auto exact = kstab::df_simple_exact(vertices, ell, a, b, c);
  require(exact < 0, "known exact witness is not relatively destabilizing");
  std::vector<kstab::IntPoint> scaled_vertices;
  for (const auto& point : vertices)
    scaled_vertices.push_back({3 * point.x, 3 * point.y});
  const auto scaled_ell = kstab::compute_ell_p(scaled_vertices);
  require(kstab::df_simple_exact(scaled_vertices, scaled_ell, a, b, 3 * c) < 0,
          "positive integral scaling changed the exact instability sign");

  kstab::Witness witness;
  const double ad = kstab::rational_to_double(a);
  const double bd = kstab::rational_to_double(b);
  const double cd = kstab::rational_to_double(c);
  const double norm = std::hypot(ad, bd);
  witness.ux = ad / norm;
  witness.uy = bd / norm;
  witness.t = -cd / norm;
  const auto certificate = kstab::certify_witness(vertices, ell, witness, 1000000);
  require(certificate.certified && certificate.value < 0,
          "numerical witness did not recover an exact negative DF certificate");
  kstab::SearchOptions quick;
  quick.theta_steps = 32;
  quick.t_steps = 16;
  quick.refine = false;
  const auto numerical = kstab::search_witness(vertices, ell, quick);
  const auto ell_double = kstab::ell_to_double(ell);
  const double reevaluated = kstab::df_simple_double(
      vertices, ell_double, numerical.witness.ux, numerical.witness.uy,
      -numerical.witness.t);
  require(std::fabs(reevaluated - numerical.witness.value) < 1e-9,
          "scratch-buffer witness value differs from independent DF evaluation");
  require(numerical.evaluations > 0,
          "numerical witness search did not report DF evaluation count");
}

void test_timeout_and_resume() {
  const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
  const auto base = std::filesystem::temp_directory_path() /
                    ("unstable_delzant_test_" + std::to_string(stamp));
  delzant::SearchOptions options;
  options.d = 3;
  options.N = 1;
  options.M = 1;
  options.max_diameter = 1;
  options.seed = 19;
  options.time_limit_seconds = 0.02;
  options.theta_steps = 32;
  options.t_steps = 16;
  options.database = base.string() + ".sqlite";
  options.output_directory = base.string() + "_output";
  const auto first = delzant::run_search(options);
  require(!first.found && first.tested == 1,
          "the stable three-ray fan should be tested without a certificate");
  const auto resumed = delzant::run_search(options);
  require(!resumed.found && resumed.tested == 0,
          "resume should skip the already tested canonical candidate");
  sqlite3* db = nullptr;
  require(sqlite3_open(options.database.c_str(), &db) == SQLITE_OK,
          "could not inspect stage profile in search database");
  sqlite3_stmt* stmt = nullptr;
  require(sqlite3_prepare_v2(db,
      "SELECT detector_profile,stage,probe_evaluations FROM candidates LIMIT 1",
      -1, &stmt, nullptr) == SQLITE_OK,
      "could not read persisted detector stage");
  require(sqlite3_step(stmt) == SQLITE_ROW,
          "search did not persist the tested candidate stage");
  const auto* profile = sqlite3_column_text(stmt, 0);
  const auto* stage = sqlite3_column_text(stmt, 1);
  require(profile && std::string(reinterpret_cast<const char*>(profile)).find("staged-v1") != std::string::npos,
          "search database did not store detector profile");
  require(stage && (std::string(reinterpret_cast<const char*>(stage)) == "probe_complete" ||
                    std::string(reinterpret_cast<const char*>(stage)) == "confirm_complete" ||
                    std::string(reinterpret_cast<const char*>(stage)) == "final_complete"),
          "screening stage was not persisted");
  require(sqlite3_column_int64(stmt, 2) > 0,
          "probe DF evaluation count was not persisted");
  sqlite3_finalize(stmt);
  sqlite3_close(db);
  std::filesystem::remove(options.database);
  std::filesystem::remove(options.database.string() + "-wal");
  std::filesystem::remove(options.database.string() + "-shm");
}

void test_restore_certified_result_and_write_files() {
  const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
  const auto base = std::filesystem::temp_directory_path() /
                    ("unstable_delzant_result_test_" + std::to_string(stamp));
  const auto database_path = base.string() + ".sqlite";
  const auto output_directory = base.string() + "_output";

  std::vector<kstab::IntPoint> vertices{{0, 0}, {1, 0}, {-17, 21},
                                        {-15, 16}, {-10, 10}};
  vertices = kstab::normalize_polygon(std::move(vertices));
  std::vector<delzant::Direction> fan;
  std::vector<std::int64_t> lengths;
  std::string vertices_text;
  std::string fan_text;
  std::string steps_text;
  std::int64_t twice_area = 0;
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    const auto& p = vertices[i];
    const auto& q = vertices[(i + 1) % vertices.size()];
    const auto dx = q.x - p.x;
    const auto dy = q.y - p.y;
    const auto length = std::gcd(std::llabs(dx), std::llabs(dy));
    fan.push_back({dx / length, dy / length});
    lengths.push_back(length);
    twice_area += p.x * q.y - p.y * q.x;
    if (i) {
      vertices_text += ';';
      fan_text += ';';
      steps_text += ',';
    }
    vertices_text += std::to_string(p.x) + ':' + std::to_string(p.y);
    fan_text += std::to_string(fan.back().x) + ':' + std::to_string(fan.back().y);
    steps_text += std::to_string(length);
  }
  const auto ell = kstab::compute_ell_p(vertices);
  const kstab::Rational a(CGAL::Gmpz(-5), CGAL::Gmpz(6));
  const kstab::Rational b(CGAL::Gmpz(-5), CGAL::Gmpz(9));
  const kstab::Rational c(CGAL::Gmpz(-7), CGAL::Gmpz(3));
  kstab::Witness numerical_witness;
  const double ad = kstab::rational_to_double(a);
  const double bd = kstab::rational_to_double(b);
  const double norm = std::hypot(ad, bd);
  numerical_witness.ux = ad / norm;
  numerical_witness.uy = bd / norm;
  numerical_witness.t = -kstab::rational_to_double(c) / norm;
  const auto cert = kstab::certify_witness(vertices, ell, numerical_witness, 1000000);
  require(cert.certified && cert.value < 0,
          "fixture for stored result is not exactly certified unstable");

  sqlite3* db = nullptr;
  require(sqlite3_open(database_path.c_str(), &db) == SQLITE_OK,
          "could not create result-test database");
  const char* schema =
      "CREATE TABLE candidates(key TEXT PRIMARY KEY,d INTEGER,fan TEXT,steps TEXT,"
      "vertices TEXT,twice_area INTEGER,ell0 TEXT,ell1 TEXT,ell2 TEXT,status TEXT,"
      "witness_a TEXT,witness_b TEXT,witness_c TEXT,df_value TEXT);"
      "CREATE TABLE state(name TEXT PRIMARY KEY,value TEXT);";
  require(sqlite3_exec(db, schema, nullptr, nullptr, nullptr) == SQLITE_OK,
          "could not create result-test schema");
  const char* insert =
      "INSERT INTO candidates VALUES('fixture',5,?,?,?,?,?,?,?,'verified_unstable',?,?,?,?)";
  sqlite3_stmt* stmt = nullptr;
  require(sqlite3_prepare_v2(db, insert, -1, &stmt, nullptr) == SQLITE_OK,
          "could not prepare certified fixture row");
  const std::string ell0 = [&] { std::ostringstream s; s << ell[0]; return s.str(); }();
  const std::string ell1 = [&] { std::ostringstream s; s << ell[1]; return s.str(); }();
  const std::string ell2 = [&] { std::ostringstream s; s << ell[2]; return s.str(); }();
  const std::string wa = [&] { std::ostringstream s; s << cert.coefficients[0]; return s.str(); }();
  const std::string wb = [&] { std::ostringstream s; s << cert.coefficients[1]; return s.str(); }();
  const std::string wc = [&] { std::ostringstream s; s << cert.coefficients[2]; return s.str(); }();
  const std::string df = [&] { std::ostringstream s; s << cert.value; return s.str(); }();
  sqlite3_bind_text(stmt, 1, fan_text.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt, 2, steps_text.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt, 3, vertices_text.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_int64(stmt, 4, twice_area);
  sqlite3_bind_text(stmt, 5, ell0.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt, 6, ell1.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt, 7, ell2.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt, 8, wa.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt, 9, wb.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt, 10, wc.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt, 11, df.c_str(), -1, SQLITE_TRANSIENT);
  require(sqlite3_step(stmt) == SQLITE_DONE, "could not insert certified fixture row");
  sqlite3_finalize(stmt);
  sqlite3_close(db);

  delzant::SearchOptions options;
  options.d = 5;
  options.max_diameter = 100;
  options.database = database_path;
  options.output_directory = output_directory;
  options.time_limit_seconds = 0.02;
  const auto restored = delzant::run_search(options);
  require(restored.found && restored.certificate.value < 0,
          "search did not restore its existing certified candidate");
  for (const auto& extension : {".polygon", ".svg", ".txt"})
    require(std::filesystem::exists(output_directory + "/unstable_d5" + extension),
            "restoring a certified result did not write all result files");
  std::filesystem::remove(database_path);
  std::filesystem::remove(database_path + "-wal");
  std::filesystem::remove(database_path + "-shm");
  std::filesystem::remove_all(output_directory);
}

}  // namespace

int main() {
  try {
    test_fan_and_polygon_generation();
    test_scale_dedup_and_size_limit();
    test_exact_relative_df_certificate();
    test_timeout_and_resume();
    test_restore_certified_result_and_write_files();
    std::cout << "unstable_Delzant tests passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "test failure: " << error.what() << '\n';
    return 1;
  }
}
