#include "search.hpp"

#include <sqlite3.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <limits>
#include <optional>
#include <queue>
#include <random>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <tuple>

namespace kstab {
namespace {

using WideInt = __int128_t;

std::string rat_string(const Rational& value) {
  std::ostringstream out;
  out << value;
  return out.str();
}

std::string sql_quote(const std::string& text) {
  std::string result = "'";
  for (char c : text) {
    if (c == '\'') result += "''";
    else result += c;
  }
  result += '\'';
  return result;
}

std::string join_dirs(const std::vector<Direction>& dirs) {
  std::ostringstream out;
  for (std::size_t i = 0; i < dirs.size(); ++i) {
    if (i) out << ';';
    out << dirs[i].x << ':' << dirs[i].y;
  }
  return out.str();
}

std::vector<Direction> split_dirs(const std::string& text) {
  std::vector<Direction> result;
  std::istringstream input(text);
  std::string item;
  while (std::getline(input, item, ';')) {
    const std::size_t colon = item.find(':');
    if (colon == std::string::npos) continue;
    result.push_back({std::stoll(item.substr(0, colon)),
                      std::stoll(item.substr(colon + 1))});
  }
  return result;
}

std::string join_steps(const std::vector<std::int64_t>& steps) {
  std::ostringstream out;
  for (std::size_t i = 0; i < steps.size(); ++i) {
    if (i) out << ',';
    out << steps[i];
  }
  return out.str();
}

std::string join_points(const std::vector<IntPoint>& points) {
  std::ostringstream out;
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (i) out << ';';
    out << points[i].x << ':' << points[i].y;
  }
  return out.str();
}

std::vector<std::int64_t> split_steps(const std::string& text) {
  std::vector<std::int64_t> result;
  std::istringstream input(text);
  std::string item;
  while (std::getline(input, item, ',')) result.push_back(std::stoll(item));
  return result;
}

bool primitive(Direction p) {
  return (p.x != 0 || p.y != 0) &&
         std::gcd(std::llabs(p.x), std::llabs(p.y)) == 1;
}

WideInt cross(Direction a, Direction b) {
  return static_cast<WideInt>(a.x) * b.y - static_cast<WideInt>(a.y) * b.x;
}

std::string profile_name(DetectorTier tier) {
  switch (tier) {
    case DetectorTier::probe: return "probe:32x16:norefine:b4";
    case DetectorTier::confirm: return "confirm:128x64:refine:b8";
    case DetectorTier::final: return "final:720x512:refine:b16";
  }
  return "unknown";
}

SearchOptions detector_options(DetectorTier tier) {
  SearchOptions result;
  if (tier == DetectorTier::probe) {
    result.theta_steps = 32;
    result.t_steps = 16;
    result.refine = false;
  } else if (tier == DetectorTier::confirm) {
    result.theta_steps = 128;
    result.t_steps = 64;
    result.refine = true;
  } else {
    result.theta_steps = 720;
    result.t_steps = 512;
    result.refine = true;
  }
  return result;
}

void set_witness_from_line(Witness& witness, double a, double b, double c,
                           double value) {
  const double length = std::hypot(a, b);
  witness.ux = length > 0 ? a / length : 0;
  witness.uy = length > 0 ? b / length : 0;
  witness.t = length > 0 ? -c / length : 0;
  witness.value = value;
}

bool exact_integer_fallback(const PolygonCandidate& candidate, int bound,
                            DetectorOutcome& outcome) {
  std::vector<Direction> normals = candidate.facet_normals;
  for (int a = -bound; a <= bound; ++a) {
    for (int b = -bound; b <= bound; ++b) {
      if ((a == 0 && b == 0) || std::gcd(std::abs(a), std::abs(b)) != 1) continue;
      normals.push_back({a, b});
    }
  }
  std::unordered_set<std::string> seen;
  for (const Direction n : normals) {
    std::vector<Rational> offsets;
    offsets.reserve(candidate.vertices.size() * 2);
    for (std::size_t i = 0; i < candidate.vertices.size(); ++i) {
      const auto& p = candidate.vertices[i];
      const auto& q = candidate.vertices[(i + 1) % candidate.vertices.size()];
      offsets.emplace_back(-(Rational(CGAL::Gmpz(static_cast<long>(n.x))) * p.x +
                             Rational(CGAL::Gmpz(static_cast<long>(n.y))) * p.y));
      offsets.emplace_back(-(Rational(CGAL::Gmpz(static_cast<long>(n.x))) *
                             (Rational(CGAL::Gmpz(static_cast<long>(p.x))) +
                              Rational(CGAL::Gmpz(static_cast<long>(q.x)))) / 2 +
                             Rational(CGAL::Gmpz(static_cast<long>(n.y))) *
                                 (Rational(CGAL::Gmpz(static_cast<long>(p.y))) +
                                  Rational(CGAL::Gmpz(static_cast<long>(q.y)))) / 2));
    }
    for (const Rational& c : offsets) {
      const std::string marker = std::to_string(n.x) + ":" +
                                  std::to_string(n.y) + ":" + rat_string(c);
      if (!seen.insert(marker).second) continue;
      const Rational value = df_simple_exact(
          candidate.vertices, candidate.ell,
          Rational(CGAL::Gmpz(static_cast<long>(n.x))),
          Rational(CGAL::Gmpz(static_cast<long>(n.y))), c);
      if (value < 0) {
        outcome.verified_unstable = true;
        outcome.certification.certified = true;
        outcome.certification.coefficients = {
            Rational(CGAL::Gmpz(static_cast<long>(n.x))),
            Rational(CGAL::Gmpz(static_cast<long>(n.y))), c};
        outcome.certification.value = value;
        outcome.numerical_negative = true;
        set_witness_from_line(outcome.numerical.witness,
                              static_cast<double>(n.x), static_cast<double>(n.y),
                              rational_to_double(c), rational_to_double(value));
        return true;
      }
    }
  }
  return false;
}

std::string serialize_rng(const std::mt19937_64& rng) {
  std::ostringstream out;
  out << rng;
  return out.str();
}

void restore_rng(std::mt19937_64& rng, const std::string& state) {
  if (!state.empty()) {
    std::istringstream in(state);
    in >> rng;
  }
}

std::int64_t shell_value(int base, std::uint64_t shell, bool n_value) {
  const std::uint64_t exponent = n_value ? shell / 2 : (shell + 1) / 2;
  const std::uint64_t multiplier = std::uint64_t{1} << std::min<std::uint64_t>(exponent, 30);
  if (base > std::numeric_limits<int>::max() / static_cast<int>(multiplier))
    return std::numeric_limits<int>::max();
  return static_cast<std::int64_t>(base * multiplier);
}

}  // namespace

std::string DetectorProfile::fingerprint() const {
  return revision + ":" + profile_name(tier) + ":bound" +
         std::to_string(integer_normal_bound) + ":cap" +
         std::to_string(certification_cap);
}

std::vector<Direction> primitive_directions(int bound) {
  std::vector<Direction> result;
  for (int x = -bound; x <= bound; ++x) {
    for (int y = -bound; y <= bound; ++y) {
      if (primitive({x, y})) result.push_back({x, y});
    }
  }
  std::sort(result.begin(), result.end(), [](Direction a, Direction b) {
    const double aa = std::atan2(static_cast<double>(a.y), static_cast<double>(a.x));
    const double ab = std::atan2(static_cast<double>(b.y), static_cast<double>(b.x));
    const double xa = aa < 0 ? aa + 2 * M_PI : aa;
    const double xb = ab < 0 ? ab + 2 * M_PI : ab;
    return xa == xb ? (a.x == b.x ? a.y < b.y : a.x < b.x) : xa < xb;
  });
  return result;
}

bool build_candidate(int d, const std::vector<Direction>& first_directions,
                     const std::vector<std::int64_t>& first_steps,
                     PolygonCandidate& result, std::string* reason) {
  auto fail = [&](const std::string& message) {
    if (reason) *reason = message;
    return false;
  };
  if (d < 3 || first_directions.size() != static_cast<std::size_t>(d - 1) ||
      first_steps.size() != static_cast<std::size_t>(d - 1))
    return fail("wrong candidate arity");
  if (first_directions.front().x != 1 || first_directions.front().y != 0)
    return fail("p1 must be (1,0)");
  for (std::size_t i = 0; i < first_directions.size(); ++i) {
    if (!primitive(first_directions[i]) || first_steps[i] <= 0)
      return fail("nonprimitive direction or nonpositive step");
    if (i && cross(first_directions[i - 1], first_directions[i]) <= 0)
      return fail("directions are not strictly increasing");
  }
  std::int64_t common = 0;
  for (const auto k : first_steps) common = std::gcd(common, k);
  std::vector<std::int64_t> steps = first_steps;
  if (common > 1) for (auto& k : steps) k /= common;

  std::vector<IntPoint> vertices{{0, 0}};
  for (std::size_t i = 0; i < first_directions.size(); ++i) {
    const WideInt x = static_cast<WideInt>(vertices.back().x) +
                      static_cast<WideInt>(steps[i]) * first_directions[i].x;
    const WideInt y = static_cast<WideInt>(vertices.back().y) +
                      static_cast<WideInt>(steps[i]) * first_directions[i].y;
    if (x < std::numeric_limits<std::int64_t>::min() ||
        x > std::numeric_limits<std::int64_t>::max() ||
        y < std::numeric_limits<std::int64_t>::min() ||
        y > std::numeric_limits<std::int64_t>::max())
      return fail("coordinate overflow");
    vertices.push_back({static_cast<std::int64_t>(x), static_cast<std::int64_t>(y)});
    if (i > 0 && vertices.back().y <= 0) return fail("vertex not in upper half-plane");
  }
  if (vertices.back().x == std::numeric_limits<std::int64_t>::min() ||
      vertices.back().y == std::numeric_limits<std::int64_t>::min())
    return fail("closing edge overflow");
  const Direction closing{-vertices.back().x, -vertices.back().y};
  const std::int64_t closing_step = std::gcd(std::llabs(closing.x), std::llabs(closing.y));
  if (closing_step <= 0) return fail("zero closing edge");
  const Direction closing_primitive{closing.x / closing_step, closing.y / closing_step};
  if (cross(first_directions.back(), closing_primitive) <= 0 ||
      cross(closing_primitive, first_directions.front()) <= 0)
    return fail("closing turn is not strictly convex");
  if (vertices.back().y <= 0) return fail("last vertex not in upper half-plane");
  WideInt area = 0;
  for (std::size_t i = 0; i < vertices.size(); ++i)
    area += static_cast<WideInt>(vertices[i].x) * vertices[(i + 1) % vertices.size()].y -
            static_cast<WideInt>(vertices[i].y) * vertices[(i + 1) % vertices.size()].x;
  if (area <= 0 || area > std::numeric_limits<std::int64_t>::max()) return fail("invalid area");
  try {
    const auto normalized = normalize_polygon(vertices);
    if (normalized.size() != static_cast<std::size_t>(d))
      return fail("self-intersection or non-strict polygon");
    for (std::size_t i = 0; i < vertices.size(); ++i)
      if (normalized[i].x != vertices[i].x || normalized[i].y != vertices[i].y)
        return fail("polygon normalization changed candidate");
  } catch (const std::exception& error) {
    return fail(error.what());
  }
  result = {};
  result.d = d;
  result.directions = first_directions;
  result.directions.push_back(closing_primitive);
  result.steps = steps;
  result.steps.push_back(closing_step);
  result.vertices = vertices;
  for (const Direction p : result.directions) result.facet_normals.push_back({-p.y, p.x});
  result.ell = compute_ell_p(vertices);
  result.twice_area = static_cast<std::int64_t>(area);
  result.key = candidate_key(result);
  return true;
}

std::string candidate_key(const PolygonCandidate& candidate) {
  return "d" + std::to_string(candidate.d) + "|p=" + join_dirs(candidate.directions) +
         "|k=" + join_steps(candidate.steps);
}

std::string current_validation_profile(const AreaSearchOptions& options) {
  // This is the cache key for the complete probe -> confirm -> final pipeline.
  // Bump both revisions when geometry validation or detector mathematics changes.
  return "validation-v2|geometry=strict-convex-hull-v2|detector=area-search-detector-v1|"
         "probe=32x16:norefine|confirm=128x64:refine|final=720x512:refine|"
         "integer-normals=4,8,16|certify-cap=" +
         std::to_string(options.certify_max_denominator);
}

DetectorOutcome detect_candidate(const PolygonCandidate& candidate,
                                 const DetectorProfile& profile,
                                 std::int64_t max_denominator) {
  DetectorOutcome outcome;
  outcome.profile = profile.fingerprint();
  const SearchResult numeric = search_witness(candidate.vertices, candidate.ell,
                                               detector_options(profile.tier));
  outcome.numerical = numeric;
  outcome.numerical_negative = numeric.unstable;
  if (numeric.unstable) {
    outcome.certification = certify_witness(candidate.vertices, candidate.ell,
                                            numeric.witness, max_denominator);
    outcome.verified_unstable = outcome.certification.certified;
  }
  if (!outcome.verified_unstable) {
    const int bound = profile.tier == DetectorTier::probe ? 4 :
                      profile.tier == DetectorTier::confirm ? 8 : 16;
    exact_integer_fallback(candidate, std::min(bound, profile.integer_normal_bound), outcome);
  }
  return outcome;
}

struct SearchDatabase::Impl {
  sqlite3* db = nullptr;
};

SearchDatabase::SearchDatabase(const std::filesystem::path& path) : impl_(new Impl) {
  if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
  if (sqlite3_open(path.string().c_str(), &impl_->db) != SQLITE_OK)
    throw std::runtime_error("cannot open SQLite database: " + path.string());
  const char* schema =
      "PRAGMA journal_mode=WAL;"
      "CREATE TABLE IF NOT EXISTS candidates(key TEXT PRIMARY KEY,d INTEGER,directions TEXT,steps TEXT,twice_area TEXT,probe_score REAL,status TEXT,ell0 TEXT,ell1 TEXT,ell2 TEXT,vertices TEXT,normals TEXT);"
      "CREATE TABLE IF NOT EXISTS attempts(candidate_key TEXT,profile TEXT,status TEXT,value REAL,witness_ux REAL,witness_uy REAL,witness_t REAL,exact_a TEXT,exact_b TEXT,exact_c TEXT,exact_value TEXT,numerical_negative INTEGER NOT NULL DEFAULT 0,PRIMARY KEY(candidate_key,profile));"
      "CREATE TABLE IF NOT EXISTS candidate_validations(candidate_key TEXT,validation_profile TEXT,status TEXT,last_stage TEXT,updated_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,PRIMARY KEY(candidate_key,validation_profile));"
      "CREATE TABLE IF NOT EXISTS state(name TEXT PRIMARY KEY,value TEXT);";
  char* error = nullptr;
  if (sqlite3_exec(impl_->db, schema, nullptr, nullptr, &error) != SQLITE_OK) {
    const std::string message = error ? error : "schema error";
    sqlite3_free(error);
    throw std::runtime_error(message);
  }
  // Keep databases made by an earlier development revision usable.
  sqlite3_exec(impl_->db, "ALTER TABLE attempts ADD COLUMN numerical_negative INTEGER NOT NULL DEFAULT 0;", nullptr, nullptr, nullptr);
  sqlite3_exec(impl_->db, "ALTER TABLE candidates ADD COLUMN vertices TEXT;", nullptr, nullptr, nullptr);
  sqlite3_exec(impl_->db, "ALTER TABLE candidates ADD COLUMN normals TEXT;", nullptr, nullptr, nullptr);
}

SearchDatabase::~SearchDatabase() {
  if (impl_) {
    if (impl_->db) sqlite3_close(impl_->db);
    delete impl_;
  }
}

void SearchDatabase::save_candidate(const PolygonCandidate& c, const std::string& status) {
  const std::string sql = "INSERT OR IGNORE INTO candidates(key,d,directions,steps,twice_area,probe_score,status,ell0,ell1,ell2,vertices,normals) VALUES(" +
      sql_quote(c.key) + "," + std::to_string(c.d) + "," + sql_quote(join_dirs(c.directions)) + "," +
      sql_quote(join_steps(c.steps)) + "," + std::to_string(c.twice_area) + "," +
      std::to_string(c.probe_score) + "," + sql_quote(status) + "," + sql_quote(rat_string(c.ell[0])) + "," +
      sql_quote(rat_string(c.ell[1])) + "," + sql_quote(rat_string(c.ell[2]) ) + "," +
      sql_quote(join_points(c.vertices)) + "," + sql_quote(join_dirs(c.facet_normals)) + ");";
  char* error = nullptr;
  if (sqlite3_exec(impl_->db, sql.c_str(), nullptr, nullptr, &error) != SQLITE_OK) {
    const std::string message = error ? error : "candidate insert failed";
    sqlite3_free(error);
    throw std::runtime_error(message);
  }
}

void SearchDatabase::update_probe_score(const std::string& key, double score) {
  const std::string sql = "UPDATE candidates SET probe_score=" + std::to_string(score) +
                          " WHERE key=" + sql_quote(key) + ";";
  sqlite3_exec(impl_->db, sql.c_str(), nullptr, nullptr, nullptr);
}

bool SearchDatabase::has_attempt(const std::string& key, const std::string& profile) const {
  const std::string sql = "SELECT 1 FROM attempts WHERE candidate_key=" + sql_quote(key) +
                          " AND profile=" + sql_quote(profile) + " LIMIT 1;";
  sqlite3_stmt* statement = nullptr;
  sqlite3_prepare_v2(impl_->db, sql.c_str(), -1, &statement, nullptr);
  const bool found = sqlite3_step(statement) == SQLITE_ROW;
  sqlite3_finalize(statement);
  return found;
}

void SearchDatabase::save_attempt(const PolygonCandidate& c, const DetectorOutcome& o) {
  const auto& w = o.numerical.witness;
  std::string sql = "INSERT OR REPLACE INTO attempts "
      "(candidate_key,profile,status,value,witness_ux,witness_uy,witness_t,"
      "exact_a,exact_b,exact_c,exact_value,numerical_negative) VALUES(" +
      sql_quote(c.key) + "," +
      sql_quote(o.profile) + "," + sql_quote(o.verified_unstable ? "verified_unstable" : "unverified") + "," +
      std::to_string(w.value) + ",";
  if (o.verified_unstable) {
    sql += std::to_string(w.ux) + "," + std::to_string(w.uy) + "," + std::to_string(w.t) + "," +
           sql_quote(rat_string(o.certification.coefficients[0])) + "," +
           sql_quote(rat_string(o.certification.coefficients[1])) + "," +
           sql_quote(rat_string(o.certification.coefficients[2])) + "," +
           sql_quote(rat_string(o.certification.value));
  } else {
    sql += "NULL,NULL,NULL,NULL,NULL,NULL,NULL";
  }
  sql += "," + std::to_string(o.numerical_negative ? 1 : 0) + ");";
  char* error = nullptr;
  if (sqlite3_exec(impl_->db, sql.c_str(), nullptr, nullptr, &error) != SQLITE_OK) {
    const std::string message = error ? error : "attempt insert failed";
    sqlite3_free(error);
    throw std::runtime_error(message);
  }
  if (o.verified_unstable) {
    const std::string update = "UPDATE candidates SET status='verified_unstable' WHERE key=" +
                               sql_quote(c.key) + ";";
    sqlite3_exec(impl_->db, update.c_str(), nullptr, nullptr, nullptr);
  }
}

std::optional<AttemptRecord> SearchDatabase::get_attempt(
    const std::string& key, const std::string& profile) const {
  const std::string sql =
      "SELECT status,value,numerical_negative,exact_a,exact_b,exact_c,exact_value "
      "FROM attempts WHERE candidate_key=" + sql_quote(key) +
      " AND profile=" + sql_quote(profile) + " LIMIT 1;";
  sqlite3_stmt* statement = nullptr;
  if (sqlite3_prepare_v2(impl_->db, sql.c_str(), -1, &statement, nullptr) != SQLITE_OK)
    return std::nullopt;
  std::optional<AttemptRecord> result;
  if (sqlite3_step(statement) == SQLITE_ROW) {
    const auto text_column = [&](int column) {
      const auto* text = sqlite3_column_text(statement, column);
      return text ? std::string(reinterpret_cast<const char*>(text)) : std::string();
    };
    const bool exact = sqlite3_column_type(statement, 3) != SQLITE_NULL &&
                       sqlite3_column_type(statement, 4) != SQLITE_NULL &&
                       sqlite3_column_type(statement, 5) != SQLITE_NULL &&
                       sqlite3_column_type(statement, 6) != SQLITE_NULL &&
                       !text_column(6).empty();
    result = AttemptRecord{text_column(0), sqlite3_column_double(statement, 1),
                           sqlite3_column_int(statement, 2) != 0, exact};
  }
  sqlite3_finalize(statement);
  return result;
}

std::optional<ValidationRecord> SearchDatabase::get_validation(
    const std::string& key, const std::string& profile) const {
  const std::string sql =
      "SELECT validation_profile,status,last_stage FROM candidate_validations "
      "WHERE candidate_key=" + sql_quote(key) +
      " AND validation_profile=" + sql_quote(profile) + " LIMIT 1;";
  sqlite3_stmt* statement = nullptr;
  sqlite3_prepare_v2(impl_->db, sql.c_str(), -1, &statement, nullptr);
  std::optional<ValidationRecord> result;
  if (sqlite3_step(statement) == SQLITE_ROW) {
    const auto text_column = [&](int column) {
      const auto* text = sqlite3_column_text(statement, column);
      return text ? std::string(reinterpret_cast<const char*>(text)) : std::string();
    };
    result = ValidationRecord{text_column(0), text_column(1), text_column(2)};
  }
  sqlite3_finalize(statement);
  return result;
}

bool SearchDatabase::has_verified_candidate(const std::string& key) const {
  const std::string sql =
      "SELECT 1 FROM candidate_validations WHERE candidate_key=" +
      sql_quote(key) + " AND status='verified_unstable' "
      "UNION SELECT 1 FROM candidates AS c JOIN attempts AS a "
      "ON a.candidate_key=c.key WHERE c.key=" + sql_quote(key) +
      " AND c.status='verified_unstable' AND a.status='verified_unstable' "
      "AND a.exact_a IS NOT NULL AND a.exact_b IS NOT NULL "
      "AND a.exact_c IS NOT NULL AND a.exact_value IS NOT NULL "
      "AND length(trim(a.exact_a)) > 0 AND length(trim(a.exact_b)) > 0 "
      "AND length(trim(a.exact_c)) > 0 AND length(trim(a.exact_value)) > 0 LIMIT 1;";
  sqlite3_stmt* statement = nullptr;
  sqlite3_prepare_v2(impl_->db, sql.c_str(), -1, &statement, nullptr);
  const bool result = sqlite3_step(statement) == SQLITE_ROW;
  sqlite3_finalize(statement);
  return result;
}

void SearchDatabase::save_validation(const PolygonCandidate& c,
                                     const std::string& profile,
                                     const std::string& status,
                                     const std::string& last_stage) {
  const std::string sql =
      "INSERT OR REPLACE INTO candidate_validations "
      "(candidate_key,validation_profile,status,last_stage,updated_at) VALUES(" +
      sql_quote(c.key) + "," + sql_quote(profile) + "," + sql_quote(status) +
      "," + sql_quote(last_stage) + ",CURRENT_TIMESTAMP);";
  char* error = nullptr;
  if (sqlite3_exec(impl_->db, sql.c_str(), nullptr, nullptr, &error) != SQLITE_OK) {
    const std::string message = error ? error : "validation insert failed";
    sqlite3_free(error);
    throw std::runtime_error(message);
  }
  if (status == "verified_unstable") {
    const std::string update = "UPDATE candidates SET status='verified_unstable' WHERE key=" +
                               sql_quote(c.key) + ";";
    sqlite3_exec(impl_->db, update.c_str(), nullptr, nullptr, nullptr);
  } else if (!has_verified_candidate(c.key)) {
    const std::string update = "UPDATE candidates SET status=" + sql_quote(status) +
                               " WHERE key=" + sql_quote(c.key) + ";";
    sqlite3_exec(impl_->db, update.c_str(), nullptr, nullptr, nullptr);
  }
}

std::vector<PolygonCandidate> SearchDatabase::load_candidates() const {
  std::vector<PolygonCandidate> result;
  sqlite3_stmt* statement = nullptr;
  sqlite3_prepare_v2(impl_->db, "SELECT d,directions,steps,probe_score FROM candidates;", -1,
                     &statement, nullptr);
  while (sqlite3_step(statement) == SQLITE_ROW) {
    const int d = sqlite3_column_int(statement, 0);
    const auto dirs = split_dirs(reinterpret_cast<const char*>(sqlite3_column_text(statement, 1)));
    const auto steps = split_steps(reinterpret_cast<const char*>(sqlite3_column_text(statement, 2)));
    PolygonCandidate candidate;
    if (build_candidate(d, std::vector<Direction>(dirs.begin(), dirs.begin() + std::min<int>(d - 1, dirs.size())),
                        std::vector<std::int64_t>(steps.begin(), steps.begin() + std::min<int>(d - 1, steps.size())), candidate)) {
      candidate.probe_score = sqlite3_column_double(statement, 3);
      result.push_back(std::move(candidate));
    }
  }
  sqlite3_finalize(statement);
  return result;
}

bool SearchDatabase::candidate_is_verified(const std::string& key) const {
  return has_verified_candidate(key);
}

void SearchDatabase::save_state(const std::string& name, const std::string& value) {
  const std::string sql = "INSERT OR REPLACE INTO state VALUES(" + sql_quote(name) + "," + sql_quote(value) + ");";
  sqlite3_exec(impl_->db, sql.c_str(), nullptr, nullptr, nullptr);
}

std::string SearchDatabase::load_state(const std::string& name) const {
  sqlite3_stmt* statement = nullptr;
  const std::string sql = "SELECT value FROM state WHERE name=" + sql_quote(name) + ";";
  sqlite3_prepare_v2(impl_->db, sql.c_str(), -1, &statement, nullptr);
  std::string result;
  if (sqlite3_step(statement) == SQLITE_ROW && sqlite3_column_text(statement, 0))
    result = reinterpret_cast<const char*>(sqlite3_column_text(statement, 0));
  sqlite3_finalize(statement);
  return result;
}

std::uint64_t SearchDatabase::count_candidates() const {
  sqlite3_stmt* statement = nullptr;
  sqlite3_prepare_v2(impl_->db, "SELECT count(*) FROM candidates;", -1, &statement, nullptr);
  sqlite3_step(statement);
  const auto result = static_cast<std::uint64_t>(sqlite3_column_int64(statement, 0));
  sqlite3_finalize(statement);
  return result;
}

std::uint64_t SearchDatabase::count_status(const std::string& status) const {
  sqlite3_stmt* statement = nullptr;
  const std::string sql = "SELECT count(*) FROM candidates WHERE status=" + sql_quote(status) + ";";
  sqlite3_prepare_v2(impl_->db, sql.c_str(), -1, &statement, nullptr);
  sqlite3_step(statement);
  const auto result = static_cast<std::uint64_t>(sqlite3_column_int64(statement, 0));
  sqlite3_finalize(statement);
  return result;
}

void SearchDatabase::write_report(const std::filesystem::path& path,
                                  const SearchSummary& summary) const {
  if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
  std::ofstream output(path);
  if (!output) throw std::runtime_error("cannot write search report: " + path.string());
  if (!summary.have_verified || summary.best_verified_key.empty()) {
    output << "没找到\n";
    return;
  }

  PolygonCandidate selected;
  bool found = false;
  for (const auto& candidate : load_candidates()) {
    if (candidate.key == summary.best_verified_key) {
      selected = candidate;
      found = true;
      break;
    }
  }
  if (!found) {
    output << "没找到\n";
    return;
  }
  const BoundaryMoments boundary = boundary_moments(selected.vertices);
  output << "status=verified_unstable\n"
         << "candidate_key=" << selected.key << '\n'
         << "d=" << selected.d << '\n'
         << "directions=" << join_dirs(selected.directions) << '\n'
         << "steps=" << join_steps(selected.steps) << '\n'
         << "vertices=" << join_points(selected.vertices) << '\n'
         << "facet_normals=" << join_dirs(selected.facet_normals) << '\n'
         << "twice_area=" << selected.twice_area << '\n'
         << "area=" << selected.twice_area << "/2\n"
         << "boundary_length_dsigma=" << rat_string(boundary.length) << '\n'
         << "boundary_ix=" << rat_string(boundary.ix) << '\n'
         << "boundary_iy=" << rat_string(boundary.iy) << '\n'
         << "ell_P=" << rat_string(selected.ell[0]) << " + ("
         << rat_string(selected.ell[1]) << ")x + ("
         << rat_string(selected.ell[2]) << ")y\n"
         << "first_verified_key=" << summary.first_verified_key << '\n'
         << "best_twice_area=" << summary.best_twice_area << '\n'
         << "database_candidates=" << count_candidates() << '\n'
         << "database_verified=" << count_status("verified_unstable") << '\n'
         << "database_unverified=" << count_status("unverified") << '\n';

  const std::string sql =
      "SELECT profile,status,value,witness_ux,witness_uy,witness_t,"
      "exact_a,exact_b,exact_c,exact_value FROM attempts WHERE candidate_key=" +
      sql_quote(selected.key) + " ORDER BY rowid;";
  sqlite3_stmt* statement = nullptr;
  sqlite3_prepare_v2(impl_->db, sql.c_str(), -1, &statement, nullptr);
  while (sqlite3_step(statement) == SQLITE_ROW) {
    const auto text_column = [&](int column) {
      const auto* text = sqlite3_column_text(statement, column);
      return text ? std::string(reinterpret_cast<const char*>(text)) : std::string();
    };
    output << "attempt_profile=" << text_column(0) << '\n'
           << "attempt_status=" << text_column(1) << '\n'
           << "attempt_numeric_value=" << sqlite3_column_double(statement, 2) << '\n';
    if (sqlite3_column_type(statement, 3) != SQLITE_NULL) {
      output << "attempt_witness_ux=" << sqlite3_column_double(statement, 3) << '\n'
             << "attempt_witness_uy=" << sqlite3_column_double(statement, 4) << '\n'
             << "attempt_witness_t=" << sqlite3_column_double(statement, 5) << '\n';
    }
    if (sqlite3_column_type(statement, 6) != SQLITE_NULL) {
      output << "attempt_exact_a=" << text_column(6) << '\n'
             << "attempt_exact_b=" << text_column(7) << '\n'
             << "attempt_exact_c=" << text_column(8) << '\n'
             << "attempt_exact_value=" << text_column(9) << '\n';
    }
  }
  sqlite3_finalize(statement);
}

SearchSummary run_search(const AreaSearchOptions& options) {
  if (options.d < 3 || options.initial_N < 1 || options.initial_M < 1 ||
      options.time_limit_seconds <= 0 || options.shell_seconds <= 0)
    throw std::invalid_argument("invalid search options");
  SearchDatabase database(options.database);
  const std::string validation_profile = current_validation_profile(options);
  std::mt19937_64 rng(options.seed);
  restore_rng(rng, database.load_state("rng"));
  std::unordered_set<std::string> queued;
  std::vector<PolygonCandidate> area_frontier;
  std::vector<PolygonCandidate> score_frontier;
  SearchSummary summary;
  const auto started = std::chrono::steady_clock::now();
  const auto deadline = started + std::chrono::duration<double>(options.time_limit_seconds);
  std::uint64_t shell = 0;
  if (!database.load_state("shell").empty()) shell = std::stoull(database.load_state("shell"));
  std::unordered_set<std::string> generated;
  bool stop_requested = false;

  const auto stage_name = [](DetectorTier tier) {
    switch (tier) {
      case DetectorTier::probe: return std::string("probe");
      case DetectorTier::confirm: return std::string("confirm");
      case DetectorTier::final: return std::string("final");
    }
    return std::string("unknown");
  };
  const auto stage_profile = [&](DetectorTier tier) {
    return validation_profile + "|stage=" + stage_name(tier);
  };
  const auto mark_verified = [&](const PolygonCandidate& candidate) {
    if (summary.first_verified_key.empty()) summary.first_verified_key = candidate.key;
    if (!summary.have_verified || candidate.twice_area < summary.best_twice_area) {
      summary.have_verified = true;
      summary.best_twice_area = candidate.twice_area;
      summary.best_verified_key = candidate.key;
    }
  };
  const auto enqueue = [&](PolygonCandidate candidate) {
    if (queued.insert(candidate.key).second) {
      area_frontier.push_back(candidate);
      score_frontier.push_back(std::move(candidate));
    }
  };
  const auto finish = [&]() {
    summary.unverified = database.count_status("unverified");
    summary.verified = database.count_status("verified_unstable");
    database.write_report(options.output_directory / "k_stability_search_result.txt",
                          summary);
    return summary;
  };

  const auto run_probe = [&](PolygonCandidate& candidate, bool resume_existing) {
    // A crash after the probe attempt was committed but before the validation
    // row was advanced must resume at confirm without evaluating the probe again.
    if (resume_existing) {
      const auto previous = database.get_attempt(candidate.key,
                                                 stage_profile(DetectorTier::probe));
      if (previous && previous->status == "verified_unstable" &&
          previous->has_exact_witness) {
        database.save_validation(candidate, validation_profile,
                                 "verified_unstable", "probe");
        mark_verified(candidate);
        if (options.stop_on_first) stop_requested = true;
        return false;
      }
      if (previous && previous->status == "unverified") {
        candidate.probe_score = previous->value;
        database.update_probe_score(candidate.key, candidate.probe_score);
        database.save_validation(candidate, validation_profile, "pending", "confirm");
        return true;
      }
    }
    database.save_validation(candidate, validation_profile, "pending", "probe");
    DetectorProfile profile{DetectorTier::probe, 4, options.certify_max_denominator};
    DetectorOutcome outcome = detect_candidate(candidate, profile,
                                                options.certify_max_denominator);
    outcome.profile = stage_profile(DetectorTier::probe);
    database.save_attempt(candidate, outcome);
    ++summary.probes;
    candidate.probe_score = outcome.numerical.witness.value;
    database.update_probe_score(candidate.key, candidate.probe_score);
    if (outcome.verified_unstable) {
      database.save_validation(candidate, validation_profile,
                               "verified_unstable", "probe");
      mark_verified(candidate);
      if (options.stop_on_first) stop_requested = true;
      return false;
    }
    // Probe is only a ranking stage. Even a nonnegative probe must reach the
    // area/score frontier so that confirm periodically samples shallow cases.
    database.save_validation(candidate, validation_profile, "pending", "confirm");
    return true;
  };

  const auto run_final = [&](PolygonCandidate& candidate) {
    const std::string profile_key = stage_profile(DetectorTier::final);
    if (const auto previous = database.get_attempt(candidate.key, profile_key)) {
      if (previous->status == "verified_unstable" && previous->has_exact_witness) {
        database.save_validation(candidate, validation_profile,
                                 "verified_unstable", "final");
        mark_verified(candidate);
        if (options.stop_on_first) stop_requested = true;
      } else if (previous->status == "unverified") {
        database.save_validation(candidate, validation_profile,
                                 "unverified", "complete");
        return;
      }
    }
    database.save_validation(candidate, validation_profile, "pending", "final");
    DetectorProfile profile{DetectorTier::final, 16, options.certify_max_denominator};
    DetectorOutcome outcome = detect_candidate(candidate, profile,
                                                options.certify_max_denominator);
    outcome.profile = profile_key;
    database.save_attempt(candidate, outcome);
    ++summary.finals;
    if (outcome.verified_unstable) {
      database.save_validation(candidate, validation_profile,
                               "verified_unstable", "final");
      mark_verified(candidate);
      if (options.stop_on_first) stop_requested = true;
    } else {
      database.save_validation(candidate, validation_profile,
                               "unverified", "complete");
    }
  };

  const auto prepare_loaded_candidate = [&](PolygonCandidate candidate,
                                             bool count_as_generated) {
    generated.insert(candidate.key);
    if (count_as_generated) ++summary.generated;
    if (database.has_verified_candidate(candidate.key)) {
      mark_verified(candidate);
      if (!database.get_validation(candidate.key, validation_profile)) {
        // A legacy verified row is reusable only after load_candidates has
        // rebuilt and strictly validated the polygon and an exact witness was
        // found by has_verified_candidate().
        database.save_validation(candidate, validation_profile,
                                 "verified_unstable", "legacy");
      }
      ++summary.skipped;
      return;
    }
    const auto validation = database.get_validation(candidate.key, validation_profile);
    if (validation && validation->status == "unverified") {
      ++summary.skipped;
      return;
    }
    if (!validation) {
      database.save_candidate(candidate, "pending");
      if (std::chrono::steady_clock::now() >= deadline) {
        stop_requested = true;
        return;
      }
      if (run_probe(candidate, false)) enqueue(std::move(candidate));
      return;
    }
    if (validation->status == "pending" && validation->last_stage == "probe") {
      if (std::chrono::steady_clock::now() >= deadline) {
        stop_requested = true;
        return;
      }
      if (run_probe(candidate, true)) enqueue(std::move(candidate));
      return;
    }
    if (validation->status == "pending" &&
        (validation->last_stage == "confirm" || validation->last_stage == "final")) {
      enqueue(std::move(candidate));
      return;
    }
    ++summary.skipped;
  };

  for (auto candidate : database.load_candidates()) {
    prepare_loaded_candidate(std::move(candidate), false);
  }
  if (options.stop_on_first && summary.have_verified) stop_requested = true;

  while (!stop_requested && std::chrono::steady_clock::now() < deadline) {
    const int n_bound = static_cast<int>(shell_value(options.initial_N, shell, true));
    const int m_bound = static_cast<int>(shell_value(options.initial_M, shell, false));
    if (options.verbose)
      std::cerr << "search shell=" << shell << " N=" << n_bound
                << " M=" << m_bound << '\n';
    const auto shell_deadline = std::min(deadline, std::chrono::steady_clock::now() +
        std::chrono::duration<double>(options.shell_seconds));
    const auto directions = primitive_directions(n_bound);
    std::uniform_int_distribution<std::size_t> direction_pick(0, directions.size() - 1);
    auto sample_step = [&]() {
      const int small = std::min(m_bound, 4);
      if (std::uniform_int_distribution<int>(0, 99)(rng) < 75)
        return static_cast<std::int64_t>(std::uniform_int_distribution<int>(1, small)(rng));
      return static_cast<std::int64_t>(std::uniform_int_distribution<int>(1, m_bound)(rng));
    };
    while (!stop_requested && std::chrono::steady_clock::now() < shell_deadline &&
           std::chrono::steady_clock::now() < deadline) {
      for (int batch = 0; batch < options.beam_width && !stop_requested &&
           std::chrono::steady_clock::now() < shell_deadline; ++batch) {
        std::vector<Direction> dirs{{1, 0}};
        std::vector<std::int64_t> steps{sample_step()};
        for (int i = 1; i < options.d - 1; ++i) {
          bool found = false;
          for (int attempt = 0; attempt < 24 && !found; ++attempt) {
            const Direction p = directions[direction_pick(rng)];
            if (cross(dirs.back(), p) > 0) {
              dirs.push_back(p);
              steps.push_back(sample_step());
              found = true;
            }
          }
          if (!found) break;
        }
        if (dirs.size() != static_cast<std::size_t>(options.d - 1) ||
            steps.size() != static_cast<std::size_t>(options.d - 1)) {
          ++summary.rejected;
          continue;
        }
        PolygonCandidate candidate;
        if (build_candidate(options.d, dirs, steps, candidate)) {
          if (generated.insert(candidate.key).second) {
            prepare_loaded_candidate(std::move(candidate), true);
          }
        }
        else ++summary.rejected;
      }
      auto pop_area = [&]() -> std::optional<PolygonCandidate> {
        if (area_frontier.empty()) return std::nullopt;
        auto it = std::min_element(area_frontier.begin(), area_frontier.end(), [](const auto& a, const auto& b) {
          const auto rank = [](const PolygonCandidate& candidate) {
            std::int64_t direction_norm = 0;
            for (int i = 0; i < candidate.d - 1; ++i)
              direction_norm = std::max(direction_norm, std::max(std::llabs(candidate.directions[i].x),
                                                                 std::llabs(candidate.directions[i].y)));
            std::int64_t max_step = 0;
            for (int i = 0; i < candidate.d - 1; ++i) max_step = std::max(max_step, candidate.steps[i]);
            std::int64_t perimeter = 0;
            for (const auto step : candidate.steps) perimeter += step;
            return std::tuple(candidate.twice_area, candidate.d, direction_norm,
                              max_step, perimeter, candidate.key);
          };
          return rank(a) < rank(b);
        });
        PolygonCandidate c = std::move(*it); area_frontier.erase(it); return c;
      };
      auto pop_score = [&]() -> std::optional<PolygonCandidate> {
        if (score_frontier.empty()) return std::nullopt;
        auto it = std::min_element(score_frontier.begin(), score_frontier.end(), [](const auto& a, const auto& b) {
          return std::tie(a.probe_score, a.key) < std::tie(b.probe_score, b.key);
        });
        PolygonCandidate c = std::move(*it); score_frontier.erase(it); return c;
      };
      for (int take = 0; take < 4 && !stop_requested &&
           std::chrono::steady_clock::now() < deadline; ++take) {
        const bool use_score = (take % 4) == 3;
        auto candidate = use_score ? pop_score() : pop_area();
        if (!candidate) break;
        const auto validation = database.get_validation(candidate->key, validation_profile);
        if (!validation || validation->status != "pending") {
          ++summary.skipped;
          continue;
        }
        if (validation->last_stage == "final") {
          run_final(*candidate);
          continue;
        }
        if (validation->last_stage == "probe") {
          if (run_probe(*candidate, true)) enqueue(std::move(*candidate));
          continue;
        }
        database.save_validation(*candidate, validation_profile, "pending", "confirm");
        const std::string confirm_key = stage_profile(DetectorTier::confirm);
        const auto previous_confirm = database.get_attempt(candidate->key, confirm_key);
        bool need_final = false;
        if (previous_confirm && previous_confirm->status == "verified_unstable" &&
            previous_confirm->has_exact_witness) {
          database.save_validation(*candidate, validation_profile,
                                   "verified_unstable", "confirm");
          mark_verified(*candidate);
          if (options.stop_on_first) stop_requested = true;
          continue;
        }
        if (previous_confirm && previous_confirm->status == "unverified") {
          need_final = previous_confirm->numerical_negative ||
                       previous_confirm->value < -1e-8 ||
                       candidate->probe_score < -1e-8;
        } else {
          DetectorProfile profile{DetectorTier::confirm, 8,
                                  options.certify_max_denominator};
          DetectorOutcome outcome = detect_candidate(*candidate, profile,
                                                     options.certify_max_denominator);
          outcome.profile = confirm_key;
          database.save_attempt(*candidate, outcome);
          ++summary.confirms;
          if (outcome.verified_unstable) {
            database.save_validation(*candidate, validation_profile,
                                     "verified_unstable", "confirm");
            mark_verified(*candidate);
            if (options.stop_on_first) stop_requested = true;
            continue;
          }
          need_final = outcome.numerical_negative || candidate->probe_score < -1e-8;
        }
        if (need_final) {
          database.save_validation(*candidate, validation_profile, "pending", "final");
          run_final(*candidate);
        } else {
          database.save_validation(*candidate, validation_profile,
                                   "unverified", "complete");
        }
      }
      database.save_state("rng", serialize_rng(rng));
      database.save_state("shell", std::to_string(shell));
    }
    ++shell;
    database.save_state("shell", std::to_string(shell));
    if (shell > 30) break;
  }
  return finish();
}

}  // namespace kstab
