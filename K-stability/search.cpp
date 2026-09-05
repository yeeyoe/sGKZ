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
#include <unordered_map>
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

std::vector<IntPoint> split_points(const std::string& text) {
  std::vector<IntPoint> result;
  std::istringstream input(text);
  std::string item;
  while (std::getline(input, item, ';')) {
    const std::size_t colon = item.find(':');
    if (colon == std::string::npos) throw std::invalid_argument("invalid point list");
    result.push_back({std::stoll(item.substr(0, colon)),
                      std::stoll(item.substr(colon + 1))});
  }
  return result;
}

std::vector<int> split_flags(const std::string& text) {
  std::vector<int> result;
  std::istringstream input(text);
  std::string item;
  while (std::getline(input, item, ',')) {
    const int flag = std::stoi(item);
    if (flag != 0 && flag != 1) throw std::invalid_argument("invalid singularity flag");
    result.push_back(flag);
  }
  return result;
}

Rational parse_rational(const std::string& text) {
  const std::size_t slash = text.find('/');
  if (slash == std::string::npos) return Rational(CGAL::Gmpz(text));

  // ell_P coefficients can have numerators and denominators larger than
  // int64_t.  Parse both parts directly as arbitrary-precision integers.
  const CGAL::Gmpz numerator(text.substr(0, slash));
  const CGAL::Gmpz denominator(text.substr(slash + 1));
  if (denominator == 0) throw std::invalid_argument("zero rational denominator");
  return Rational(numerator, denominator);
}

std::string join_flags(const std::vector<int>& flags) {
  std::ostringstream out;
  for (std::size_t i = 0; i < flags.size(); ++i) {
    if (i) out << ',';
    out << flags[i];
  }
  return out.str();
}

bool primitive(Direction p) {
  return (p.x != 0 || p.y != 0) &&
         std::gcd(std::llabs(p.x), std::llabs(p.y)) == 1;
}

WideInt cross(Direction a, Direction b) {
  return static_cast<WideInt>(a.x) * b.y - static_cast<WideInt>(a.y) * b.x;
}

WideInt abs_wide(WideInt value) {
  return value < 0 ? -value : value;
}

WideInt cross_to_origin(Direction direction, const IntPoint& vertex) {
  return -static_cast<WideInt>(direction.x) * vertex.y +
         static_cast<WideInt>(direction.y) * vertex.x;
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

std::optional<int> shell_value(int base, std::uint64_t shell, bool n_value) {
  // Keep each range for several equal time slices. Random sampling does not
  // exhaust a range, so this gives the smallest feasible polygons more
  // opportunities before enlarging the direction or step bound.
  //
  // shell : 0..3  4..7  8..11  12..15 ...
  // N     :   1      2       2        3 ...
  // M     :   1      1       2        2 ...
  constexpr std::uint64_t k_shells_per_range = 4;
  const std::uint64_t stage = shell / k_shells_per_range;
  const std::uint64_t multiplier =
      n_value ? ((stage + 1) / 2 + 1) : (stage / 2 + 1);
  if (multiplier > static_cast<std::uint64_t>(std::numeric_limits<int>::max()) /
                       static_cast<std::uint64_t>(base)) {
    return std::nullopt;
  }
  return static_cast<int>(static_cast<std::uint64_t>(base) * multiplier);
}

std::string dimension_state_name(int dimension, const std::string& name) {
  return "d" + std::to_string(dimension) + "|" + name;
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
    if (i) {
      if (cross(first_directions[i - 1], first_directions[i]) <= 0)
        return fail("directions are not strictly increasing");
    }
  }
  std::int64_t common = 0;
  for (const auto k : first_steps) common = std::gcd(common, k);
  std::vector<std::int64_t> steps = first_steps;
  if (common > 1) for (auto& k : steps) k /= common;

  std::vector<IntPoint> vertices{{0, 0}};
  for (std::size_t i = 0; i < first_directions.size(); ++i) {
    if (i > 0) {
      if (cross_to_origin(first_directions[i], vertices.back()) <= 0)
        return fail("direction leaves the current closing wedge");
    }
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

bool candidate_from_vertices(const std::vector<IntPoint>& vertices,
                             PolygonCandidate& result, std::string* reason) {
  auto fail = [&](const std::string& message) {
    if (reason) *reason = message;
    return false;
  };
  const int d = static_cast<int>(vertices.size());
  if (d < 3) return fail("polygon has fewer than three vertices");
  if (vertices.front().x != 0 || vertices.front().y != 0)
    return fail("first vertex is not the origin");
  std::vector<Direction> directions;
  std::vector<std::int64_t> steps;
  directions.reserve(d - 1);
  steps.reserve(d - 1);
  for (int i = 0; i < d - 1; ++i) {
    const WideInt dx = static_cast<WideInt>(vertices[i + 1].x) -
                       static_cast<WideInt>(vertices[i].x);
    const WideInt dy = static_cast<WideInt>(vertices[i + 1].y) -
                       static_cast<WideInt>(vertices[i].y);
    const WideInt abs_dx = dx < 0 ? -dx : dx;
    const WideInt abs_dy = dy < 0 ? -dy : dy;
    const WideInt length = std::gcd(abs_dx, abs_dy);
    if (length <= 0 || length > std::numeric_limits<std::int64_t>::max())
      return fail("polygon has an invalid edge length");
    const WideInt primitive_x = dx / length;
    const WideInt primitive_y = dy / length;
    if (primitive_x < std::numeric_limits<std::int64_t>::min() ||
        primitive_x > std::numeric_limits<std::int64_t>::max() ||
        primitive_y < std::numeric_limits<std::int64_t>::min() ||
        primitive_y > std::numeric_limits<std::int64_t>::max())
      return fail("polygon direction overflow");
    directions.push_back({static_cast<std::int64_t>(primitive_x),
                          static_cast<std::int64_t>(primitive_y)});
    steps.push_back(static_cast<std::int64_t>(length));
  }
  if (directions.front().x != 1 || directions.front().y != 0)
    return fail("first edge is not the canonical (1,0) direction");
  if (!build_candidate(d, directions, steps, result, reason)) return false;
  if (result.vertices.size() != vertices.size())
    return fail("polygon is not in canonical search geometry");
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    if (result.vertices[i].x != vertices[i].x || result.vertices[i].y != vertices[i].y)
      return fail("polygon is not in canonical search geometry");
  }
  return true;
}

std::string candidate_key(const PolygonCandidate& candidate) {
  return "d" + std::to_string(candidate.d) + "|p=" + join_dirs(candidate.directions) +
         "|k=" + join_steps(candidate.steps);
}

std::vector<int> compute_vertex_singularity_flags(
    const PolygonCandidate& candidate) {
  if (candidate.d < 3 || candidate.directions.size() != static_cast<std::size_t>(candidate.d))
    throw std::invalid_argument("candidate has invalid direction count");
  std::vector<int> flags(candidate.d, 0);
  for (int i = 0; i < candidate.d; ++i) {
    const Direction& incoming = candidate.directions[(i + candidate.d - 1) % candidate.d];
    const Direction& outgoing = candidate.directions[i];
    flags[i] = abs_wide(cross(incoming, outgoing)) == 1 ? 0 : 1;
  }
  return flags;
}

bool candidate_is_smooth(const PolygonCandidate& candidate) {
  const auto flags = compute_vertex_singularity_flags(candidate);
  return std::all_of(flags.begin(), flags.end(), [](int flag) { return flag == 0; });
}

std::string current_validation_profile(const AreaSearchOptions& options) {
  // This is the cache key for the complete probe -> confirm -> final pipeline.
  // Bump both revisions when geometry validation or detector mathematics changes.
  return "validation-v2|geometry=incremental-convex-v1|detector=area-search-detector-v2|"
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

Rational q_squared_exact(const PolygonCandidate& candidate,
                         const std::array<Rational, 3>& coefficients,
                         const Rational& exact_value) {
  if (candidate.twice_area <= 0)
    throw std::invalid_argument("Q squared requires positive polygon area");
  Rational supremum = 0;
  for (const IntPoint& vertex : candidate.vertices) {
    const Rational value = coefficients[0] * vertex.x + coefficients[1] * vertex.y +
                           coefficients[2];
    if (value > supremum) supremum = value;
  }
  if (supremum <= 0)
    throw std::invalid_argument("Q squared requires positive witness supremum");
  const Rational area(CGAL::Gmpz(std::to_string(candidate.twice_area)), CGAL::Gmpz("2"));
  return exact_value * exact_value / (area * supremum * supremum);
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
      "CREATE TABLE IF NOT EXISTS candidates(key TEXT PRIMARY KEY,d INTEGER,directions TEXT,steps TEXT,twice_area TEXT,probe_score REAL,status TEXT,ell0 TEXT,ell1 TEXT,ell2 TEXT,vertices TEXT,normals TEXT,vertex_singularity_flags TEXT,singular_vertex_count INTEGER,q_squared_exact TEXT,q_squared_value REAL);"
      "CREATE TABLE IF NOT EXISTS attempts(candidate_key TEXT,profile TEXT,status TEXT,value REAL,witness_ux REAL,witness_uy REAL,witness_t REAL,exact_a TEXT,exact_b TEXT,exact_c TEXT,exact_value TEXT,numerical_negative INTEGER NOT NULL DEFAULT 0,q_squared_exact TEXT,q_squared_value REAL,PRIMARY KEY(candidate_key,profile));"
      "CREATE TABLE IF NOT EXISTS candidate_validations(candidate_key TEXT,validation_profile TEXT,status TEXT,last_stage TEXT,updated_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,PRIMARY KEY(candidate_key,validation_profile));"
      "CREATE TABLE IF NOT EXISTS state(name TEXT PRIMARY KEY,value TEXT);"
      "CREATE INDEX IF NOT EXISTS idx_candidates_d_status ON candidates(d,status);";
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
  sqlite3_exec(impl_->db, "ALTER TABLE candidates ADD COLUMN vertex_singularity_flags TEXT;", nullptr, nullptr, nullptr);
  sqlite3_exec(impl_->db, "ALTER TABLE candidates ADD COLUMN singular_vertex_count INTEGER;", nullptr, nullptr, nullptr);
  sqlite3_exec(impl_->db, "ALTER TABLE candidates ADD COLUMN q_squared_exact TEXT;", nullptr, nullptr, nullptr);
  sqlite3_exec(impl_->db, "ALTER TABLE candidates ADD COLUMN q_squared_value REAL;", nullptr, nullptr, nullptr);
  sqlite3_exec(impl_->db, "ALTER TABLE attempts ADD COLUMN q_squared_exact TEXT;", nullptr, nullptr, nullptr);
  sqlite3_exec(impl_->db, "ALTER TABLE attempts ADD COLUMN q_squared_value REAL;", nullptr, nullptr, nullptr);
  sqlite3_exec(impl_->db, "CREATE INDEX IF NOT EXISTS idx_candidates_q_squared_value ON candidates(q_squared_value DESC);", nullptr, nullptr, nullptr);
}

SearchDatabase::~SearchDatabase() {
  if (impl_) {
    if (impl_->db) sqlite3_close(impl_->db);
    delete impl_;
  }
}

void SearchDatabase::save_candidate(const PolygonCandidate& c, const std::string& status) {
  std::string singularity_sql = "NULL,NULL";
  if (c.singular_vertex_count >= 0 &&
      c.vertex_singularity_flags.size() == static_cast<std::size_t>(c.d)) {
    singularity_sql = sql_quote(join_flags(c.vertex_singularity_flags)) + "," +
                      std::to_string(c.singular_vertex_count);
  }
  const std::string sql = "INSERT OR IGNORE INTO candidates(key,d,directions,steps,twice_area,probe_score,status,ell0,ell1,ell2,vertices,normals,vertex_singularity_flags,singular_vertex_count) VALUES(" +
      sql_quote(c.key) + "," + std::to_string(c.d) + "," + sql_quote(join_dirs(c.directions)) + "," +
      sql_quote(join_steps(c.steps)) + "," + std::to_string(c.twice_area) + "," +
      std::to_string(c.probe_score) + "," + sql_quote(status) + "," + sql_quote(rat_string(c.ell[0])) + "," +
      sql_quote(rat_string(c.ell[1])) + "," + sql_quote(rat_string(c.ell[2]) ) + "," +
      sql_quote(join_points(c.vertices)) + "," + sql_quote(join_dirs(c.facet_normals)) + "," +
      singularity_sql + ");";
  char* error = nullptr;
  if (sqlite3_exec(impl_->db, sql.c_str(), nullptr, nullptr, &error) != SQLITE_OK) {
    const std::string message = error ? error : "candidate insert failed";
    sqlite3_free(error);
    throw std::runtime_error(message);
  }
}

void SearchDatabase::ensure_generator_revision(const std::string& revision) {
  const std::string existing = load_state("generator_revision");
  if (!existing.empty()) {
    if (existing != revision)
      throw std::runtime_error("database generator revision mismatch: " + existing +
                               " (expected " + revision + ")");
    return;
  }
  if (count_candidates() != 0)
    throw std::runtime_error("database has candidates but no generator revision; refusing legacy database");
  save_state("generator_revision", revision);
}

void SearchDatabase::save_vertex_singularity(const PolygonCandidate& candidate) {
  const std::vector<int> flags = compute_vertex_singularity_flags(candidate);
  const int count = std::accumulate(flags.begin(), flags.end(), 0);
  const std::string sql = "UPDATE candidates SET vertex_singularity_flags=" +
                          sql_quote(join_flags(flags)) + ", singular_vertex_count=" +
                          std::to_string(count) + " WHERE key=" + sql_quote(candidate.key) + ";";
  char* error = nullptr;
  if (sqlite3_exec(impl_->db, sql.c_str(), nullptr, nullptr, &error) != SQLITE_OK) {
    const std::string message = error ? error : "singularity update failed";
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
  std::optional<Rational> q_squared;
  if (o.verified_unstable) {
    try {
      q_squared = q_squared_exact(c, o.certification.coefficients,
                                  o.certification.value);
    } catch (const std::exception&) {
      q_squared.reset();
    }
  }
  std::string sql = "INSERT OR REPLACE INTO attempts "
      "(candidate_key,profile,status,value,witness_ux,witness_uy,witness_t,"
      "exact_a,exact_b,exact_c,exact_value,numerical_negative,q_squared_exact,q_squared_value) VALUES(" +
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
  sql += "," + std::to_string(o.numerical_negative ? 1 : 0) + ",";
  if (q_squared) {
    sql += sql_quote(rat_string(*q_squared)) + "," +
           std::to_string(rational_to_double(*q_squared));
  } else {
    sql += "NULL,NULL";
  }
  sql += ");";
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
    if (q_squared) {
      bool replace_candidate_q = true;
      sqlite3_stmt* current = nullptr;
      const std::string current_query =
          "SELECT q_squared_exact FROM candidates WHERE key=" + sql_quote(c.key) + ";";
      if (sqlite3_prepare_v2(impl_->db, current_query.c_str(), -1, &current, nullptr) == SQLITE_OK &&
          sqlite3_step(current) == SQLITE_ROW && sqlite3_column_type(current, 0) != SQLITE_NULL) {
        const auto* text = sqlite3_column_text(current, 0);
        if (text && *text) replace_candidate_q = *q_squared > parse_rational(reinterpret_cast<const char*>(text));
      }
      sqlite3_finalize(current);
      if (replace_candidate_q) {
        const std::string q_update =
            "UPDATE candidates SET q_squared_exact=" + sql_quote(rat_string(*q_squared)) +
            ", q_squared_value=" + std::to_string(rational_to_double(*q_squared)) +
            " WHERE key=" + sql_quote(c.key) + ";";
        sqlite3_exec(impl_->db, q_update.c_str(), nullptr, nullptr, nullptr);
      }
      // Recompute from all authenticated attempts so replacing a profile with
      // a smaller witness cannot leave a stale candidate maximum.
      std::optional<Rational> maximum;
      sqlite3_stmt* all_attempts = nullptr;
      const std::string all_attempts_query =
          "SELECT q_squared_exact FROM attempts WHERE candidate_key=" +
          sql_quote(c.key) + " AND status='verified_unstable' AND q_squared_exact IS NOT NULL;";
      if (sqlite3_prepare_v2(impl_->db, all_attempts_query.c_str(), -1, &all_attempts, nullptr) == SQLITE_OK) {
        while (sqlite3_step(all_attempts) == SQLITE_ROW) {
          const auto* text = sqlite3_column_text(all_attempts, 0);
          if (!text || !*text) continue;
          const Rational value = parse_rational(reinterpret_cast<const char*>(text));
          if (!maximum || value > *maximum) maximum = value;
        }
      }
      sqlite3_finalize(all_attempts);
      if (maximum) {
        const std::string maximum_update =
            "UPDATE candidates SET q_squared_exact=" + sql_quote(rat_string(*maximum)) +
            ",q_squared_value=" + std::to_string(rational_to_double(*maximum)) +
            " WHERE key=" + sql_quote(c.key) + ";";
        sqlite3_exec(impl_->db, maximum_update.c_str(), nullptr, nullptr, nullptr);
      }
    }
  }
}

BackfillSummary SearchDatabase::backfill_q_squared() {
  BackfillSummary summary;
  const auto candidates = load_candidates();
  std::unordered_map<std::string, PolygonCandidate> by_key;
  by_key.reserve(candidates.size());
  for (const auto& candidate : candidates) by_key.emplace(candidate.key, candidate);

  char* error = nullptr;
  if (sqlite3_exec(impl_->db, "BEGIN IMMEDIATE;", nullptr, nullptr, &error) != SQLITE_OK) {
    const std::string message = error ? error : "cannot begin Q squared backfill";
    sqlite3_free(error);
    throw std::runtime_error(message);
  }
  try {
    const char* query =
        "SELECT candidate_key,profile,exact_a,exact_b,exact_c,exact_value "
        "FROM attempts WHERE status='verified_unstable';";
    sqlite3_stmt* statement = nullptr;
    if (sqlite3_prepare_v2(impl_->db, query, -1, &statement, nullptr) != SQLITE_OK)
      throw std::runtime_error("cannot prepare Q squared backfill query");
    while (sqlite3_step(statement) == SQLITE_ROW) {
      ++summary.scanned;
      const auto text_column = [&](int column) {
        const auto* text = sqlite3_column_text(statement, column);
        return text ? std::string(reinterpret_cast<const char*>(text)) : std::string();
      };
      const std::string key = text_column(0);
      const std::string profile = text_column(1);
      if (sqlite3_column_type(statement, 2) == SQLITE_NULL ||
          sqlite3_column_type(statement, 3) == SQLITE_NULL ||
          sqlite3_column_type(statement, 4) == SQLITE_NULL ||
          sqlite3_column_type(statement, 5) == SQLITE_NULL) {
        ++summary.skipped;
        continue;
      }
      const auto candidate_it = by_key.find(key);
      if (candidate_it == by_key.end()) {
        ++summary.skipped;
        continue;
      }
      try {
        const std::array<Rational, 3> coefficients = {
            parse_rational(text_column(2)), parse_rational(text_column(3)),
            parse_rational(text_column(4))};
        const Rational exact_value = parse_rational(text_column(5));
        const Rational q = q_squared_exact(candidate_it->second, coefficients, exact_value);
        const std::string update =
            "UPDATE attempts SET q_squared_exact=" + sql_quote(rat_string(q)) +
            ",q_squared_value=" + std::to_string(rational_to_double(q)) +
            " WHERE candidate_key=" + sql_quote(key) + " AND profile=" +
            sql_quote(profile) + ";";
        if (sqlite3_exec(impl_->db, update.c_str(), nullptr, nullptr, &error) != SQLITE_OK) {
          const std::string message = error ? error : "attempt Q squared update failed";
          sqlite3_free(error);
          throw std::runtime_error(message);
        }

        sqlite3_stmt* current = nullptr;
        const std::string current_query =
            "SELECT q_squared_exact FROM candidates WHERE key=" + sql_quote(key) + ";";
        if (sqlite3_prepare_v2(impl_->db, current_query.c_str(), -1, &current, nullptr) != SQLITE_OK)
          throw std::runtime_error("cannot read candidate Q squared");
        std::optional<Rational> current_q;
        if (sqlite3_step(current) == SQLITE_ROW && sqlite3_column_type(current, 0) != SQLITE_NULL) {
          const auto* text = sqlite3_column_text(current, 0);
          if (text && *text) current_q = parse_rational(reinterpret_cast<const char*>(text));
        }
        sqlite3_finalize(current);
        if (!current_q || q > *current_q) {
          const std::string candidate_update =
              "UPDATE candidates SET q_squared_exact=" + sql_quote(rat_string(q)) +
              ",q_squared_value=" + std::to_string(rational_to_double(q)) +
              " WHERE key=" + sql_quote(key) + ";";
          if (sqlite3_exec(impl_->db, candidate_update.c_str(), nullptr, nullptr, &error) != SQLITE_OK) {
            const std::string message = error ? error : "candidate Q squared update failed";
            sqlite3_free(error);
            throw std::runtime_error(message);
          }
        }
        ++summary.backfilled;
      } catch (const std::exception&) {
        ++summary.errors;
      }
    }
    sqlite3_finalize(statement);
    if (sqlite3_exec(impl_->db, "COMMIT;", nullptr, nullptr, &error) != SQLITE_OK) {
      const std::string message = error ? error : "cannot commit Q squared backfill";
      sqlite3_free(error);
      throw std::runtime_error(message);
    }
  } catch (...) {
    sqlite3_exec(impl_->db, "ROLLBACK;", nullptr, nullptr, nullptr);
    throw;
  }
  return summary;
}

std::optional<AttemptRecord> SearchDatabase::get_attempt(
    const std::string& key, const std::string& profile) const {
  const std::string sql =
      "SELECT status,value,numerical_negative,exact_a,exact_b,exact_c,exact_value,q_squared_exact,q_squared_value "
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
    AttemptRecord record{text_column(0), sqlite3_column_double(statement, 1),
                         sqlite3_column_int(statement, 2) != 0, exact, {}, 0.0, false};
    if (sqlite3_column_type(statement, 7) != SQLITE_NULL && !text_column(7).empty()) {
      record.q_squared_exact = text_column(7);
      record.has_q_squared = true;
      record.q_squared_value = sqlite3_column_double(statement, 8);
    }
    result = record;
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
    save_vertex_singularity(c);
  } else if (!has_verified_candidate(c.key)) {
    const std::string update = "UPDATE candidates SET status=" + sql_quote(status) +
                               " WHERE key=" + sql_quote(c.key) + ";";
    sqlite3_exec(impl_->db, update.c_str(), nullptr, nullptr, nullptr);
  }
}

std::vector<PolygonCandidate> SearchDatabase::load_candidates(int dimension) const {
  if (dimension < -1 || dimension == 0)
    throw std::invalid_argument("candidate dimension must be positive or -1");
  std::vector<PolygonCandidate> result;
  sqlite3_stmt* statement = nullptr;
  std::string query =
      "SELECT key,d,directions,steps,twice_area,probe_score,status,ell0,ell1,ell2,"
      "vertices,normals,vertex_singularity_flags,singular_vertex_count FROM candidates";
  if (dimension >= 1) query += " WHERE d=" + std::to_string(dimension);
  query += ";";
  if (sqlite3_prepare_v2(impl_->db, query.c_str(), -1, &statement, nullptr) != SQLITE_OK)
    throw std::runtime_error("cannot prepare candidate load query");
  while (sqlite3_step(statement) == SQLITE_ROW) {
    const auto text_column = [&](int column) {
      const auto* text = sqlite3_column_text(statement, column);
      if (!text) throw std::runtime_error("NULL candidate field");
      return std::string(reinterpret_cast<const char*>(text));
    };
    const auto optional_text_column = [&](int column) -> std::optional<std::string> {
      if (sqlite3_column_type(statement, column) == SQLITE_NULL) return std::nullopt;
      return text_column(column);
    };
    const std::string stored_key = text_column(0);
    const int d = sqlite3_column_int(statement, 1);
    const auto dirs = split_dirs(text_column(2));
    const auto steps = split_steps(text_column(3));
    const auto twice_area = std::stoll(text_column(4));
    PolygonCandidate candidate;
    candidate.d = d;
    candidate.directions = dirs;
    candidate.steps = steps;
    candidate.twice_area = twice_area;
    candidate.probe_score = sqlite3_column_type(statement, 5) == SQLITE_NULL
                                ? 0.0
                                : sqlite3_column_double(statement, 5);
    const std::string status = text_column(6);
    candidate.ell = {parse_rational(text_column(7)), parse_rational(text_column(8)),
                     parse_rational(text_column(9))};
    candidate.vertices = split_points(text_column(10));
    candidate.facet_normals = split_dirs(text_column(11));
    candidate.key = stored_key;
    if (d < 3 || dirs.size() != static_cast<std::size_t>(d) ||
        steps.size() != static_cast<std::size_t>(d) ||
        candidate.vertices.size() != static_cast<std::size_t>(d) ||
        candidate.facet_normals.size() != static_cast<std::size_t>(d) ||
        candidate.key != candidate_key(candidate))
      throw std::runtime_error("inconsistent candidate geometry for key " + stored_key);
    for (int i = 0; i < d; ++i) {
      const IntPoint& p = candidate.vertices[i];
      const IntPoint& q = candidate.vertices[(i + 1) % d];
      const WideInt dx = static_cast<WideInt>(q.x) - p.x;
      const WideInt dy = static_cast<WideInt>(q.y) - p.y;
      if (dx != static_cast<WideInt>(candidate.steps[i]) * candidate.directions[i].x ||
          dy != static_cast<WideInt>(candidate.steps[i]) * candidate.directions[i].y ||
          candidate.facet_normals[i].x != -candidate.directions[i].y ||
          candidate.facet_normals[i].y != candidate.directions[i].x)
        throw std::runtime_error("inconsistent candidate geometry for key " + stored_key);
    }
    WideInt area = 0;
    for (int i = 0; i < d; ++i) {
      const auto& p = candidate.vertices[i];
      const auto& q = candidate.vertices[(i + 1) % d];
      area += static_cast<WideInt>(p.x) * q.y - static_cast<WideInt>(p.y) * q.x;
    }
    if (area != static_cast<WideInt>(candidate.twice_area))
      throw std::runtime_error("inconsistent candidate area for key " + stored_key);
    const auto stored_flags = optional_text_column(12);
    const bool has_count = sqlite3_column_type(statement, 13) != SQLITE_NULL;
    if (stored_flags.has_value() != has_count)
      throw std::runtime_error("incomplete candidate singularity metadata for key " + stored_key);
    if (stored_flags) {
      candidate.vertex_singularity_flags = split_flags(*stored_flags);
      candidate.singular_vertex_count = sqlite3_column_int(statement, 13);
      if (candidate.vertex_singularity_flags.size() != static_cast<std::size_t>(d) ||
          candidate.singular_vertex_count !=
              std::accumulate(candidate.vertex_singularity_flags.begin(),
                              candidate.vertex_singularity_flags.end(), 0))
        throw std::runtime_error("inconsistent candidate singularity metadata for key " + stored_key);
    } else if (status == "verified_unstable") {
      throw std::runtime_error("verified candidate lacks singularity metadata for key " + stored_key);
    }
    result.push_back(std::move(candidate));
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

std::uint64_t SearchDatabase::count_candidates(int dimension) const {
  if (dimension < -1 || dimension == 0)
    throw std::invalid_argument("candidate dimension must be positive or -1");
  sqlite3_stmt* statement = nullptr;
  std::string sql = "SELECT count(*) FROM candidates";
  if (dimension >= 1) sql += " WHERE d=" + std::to_string(dimension);
  sql += ";";
  sqlite3_prepare_v2(impl_->db, sql.c_str(), -1, &statement, nullptr);
  sqlite3_step(statement);
  const auto result = static_cast<std::uint64_t>(sqlite3_column_int64(statement, 0));
  sqlite3_finalize(statement);
  return result;
}

std::uint64_t SearchDatabase::count_status(const std::string& status, int dimension) const {
  if (dimension < -1 || dimension == 0)
    throw std::invalid_argument("candidate dimension must be positive or -1");
  sqlite3_stmt* statement = nullptr;
  std::string sql = "SELECT count(*) FROM candidates WHERE status=" + sql_quote(status);
  if (dimension >= 1) sql += " AND d=" + std::to_string(dimension);
  sql += ";";
  sqlite3_prepare_v2(impl_->db, sql.c_str(), -1, &statement, nullptr);
  sqlite3_step(statement);
  const auto result = static_cast<std::uint64_t>(sqlite3_column_int64(statement, 0));
  sqlite3_finalize(statement);
  return result;
}

std::uint64_t SearchDatabase::count_tested(int dimension) const {
  if (dimension < 1) throw std::invalid_argument("tested count requires a dimension");
  sqlite3_stmt* statement = nullptr;
  const std::string sql =
      "SELECT count(DISTINCT a.candidate_key) FROM attempts AS a "
      "JOIN candidates AS c ON c.key=a.candidate_key WHERE c.d=" +
      std::to_string(dimension) + ";";
  if (sqlite3_prepare_v2(impl_->db, sql.c_str(), -1, &statement, nullptr) != SQLITE_OK)
    throw std::runtime_error("cannot prepare tested count query");
  sqlite3_step(statement);
  const auto result = static_cast<std::uint64_t>(sqlite3_column_int64(statement, 0));
  sqlite3_finalize(statement);
  return result;
}

std::uint64_t SearchDatabase::count_verified(int dimension) const {
  if (dimension < 1) throw std::invalid_argument("verified count requires a dimension");
  sqlite3_stmt* statement = nullptr;
  const std::string sql =
      "SELECT count(*) FROM candidates AS c WHERE c.d=" + std::to_string(dimension) +
      " AND (EXISTS (SELECT 1 FROM candidate_validations AS v "
      "WHERE v.candidate_key=c.key AND v.status='verified_unstable') "
      "OR (c.status='verified_unstable' AND EXISTS (SELECT 1 FROM attempts AS a "
      "WHERE a.candidate_key=c.key AND a.status='verified_unstable' "
      "AND a.exact_a IS NOT NULL AND a.exact_b IS NOT NULL "
      "AND a.exact_c IS NOT NULL AND a.exact_value IS NOT NULL "
      "AND length(trim(a.exact_a)) > 0 AND length(trim(a.exact_b)) > 0 "
      "AND length(trim(a.exact_c)) > 0 AND length(trim(a.exact_value)) > 0)));";
  if (sqlite3_prepare_v2(impl_->db, sql.c_str(), -1, &statement, nullptr) != SQLITE_OK)
    throw std::runtime_error("cannot prepare verified count query");
  sqlite3_step(statement);
  const auto result = static_cast<std::uint64_t>(sqlite3_column_int64(statement, 0));
  sqlite3_finalize(statement);
  return result;
}

std::vector<VerifiedCandidateSummary> SearchDatabase::top_verified(
    int dimension, std::size_t limit) const {
  if (dimension < 1) throw std::invalid_argument("verified ranking requires a dimension");
  std::vector<VerifiedCandidateSummary> result;
  if (limit == 0) return result;
  sqlite3_stmt* statement = nullptr;
  const std::string sql =
      "SELECT c.key,c.twice_area,c.q_squared_exact,c.q_squared_value FROM candidates AS c WHERE c.d=" +
      std::to_string(dimension) +
      " AND (EXISTS (SELECT 1 FROM candidate_validations AS v "
      "WHERE v.candidate_key=c.key AND v.status='verified_unstable') "
      "OR (c.status='verified_unstable' AND EXISTS (SELECT 1 FROM attempts AS a "
      "WHERE a.candidate_key=c.key AND a.status='verified_unstable' "
      "AND a.exact_a IS NOT NULL AND a.exact_b IS NOT NULL "
      "AND a.exact_c IS NOT NULL AND a.exact_value IS NOT NULL "
      "AND length(trim(a.exact_a)) > 0 AND length(trim(a.exact_b)) > 0 "
      "AND length(trim(a.exact_c)) > 0 AND length(trim(a.exact_value)) > 0))) "
      "ORDER BY CAST(c.twice_area AS INTEGER), c.key LIMIT " +
      std::to_string(limit) + ";";
  if (sqlite3_prepare_v2(impl_->db, sql.c_str(), -1, &statement, nullptr) != SQLITE_OK)
    throw std::runtime_error("cannot prepare verified ranking query");
  while (sqlite3_step(statement) == SQLITE_ROW) {
    const auto* key_text = sqlite3_column_text(statement, 0);
    if (!key_text) throw std::runtime_error("NULL verified candidate key");
    const auto* area_text = sqlite3_column_text(statement, 1);
    if (!area_text) throw std::runtime_error("NULL verified candidate area");
    const auto* q_text = sqlite3_column_text(statement, 2);
    result.push_back({std::string(reinterpret_cast<const char*>(key_text)),
                      std::stoll(reinterpret_cast<const char*>(area_text)),
                      q_text ? std::string(reinterpret_cast<const char*>(q_text)) : std::string(),
                      sqlite3_column_type(statement, 3) == SQLITE_NULL
                          ? 0.0 : sqlite3_column_double(statement, 3)});
  }
  sqlite3_finalize(statement);
  return result;
}

std::optional<std::int64_t> SearchDatabase::min_verified_twice_area(int dimension) const {
  const auto top = top_verified(dimension, 1);
  if (top.empty()) return std::nullopt;
  return top.front().twice_area;
}

void SearchDatabase::write_report(const std::filesystem::path& path,
                                  const SearchSummary& summary, int dimension) const {
  if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
  std::ofstream output(path);
  if (!output) throw std::runtime_error("cannot write search report: " + path.string());
  if (!summary.have_verified || summary.best_verified_key.empty()) {
    output << "没找到\n";
    return;
  }

  PolygonCandidate selected;
  bool found = false;
  for (const auto& candidate : load_candidates(dimension)) {
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
         << "vertex_singularity_flags=" << join_flags(selected.vertex_singularity_flags) << '\n'
         << "singular_vertex_count=" << selected.singular_vertex_count << '\n'
         << "twice_area=" << selected.twice_area << '\n'
         << "area=" << selected.twice_area << "/2\n"
         << "q_squared_exact=";
  sqlite3_stmt* candidate_q_statement = nullptr;
  const std::string candidate_q_query =
      "SELECT q_squared_exact,q_squared_value FROM candidates WHERE key=" +
      sql_quote(selected.key) + ";";
  if (sqlite3_prepare_v2(impl_->db, candidate_q_query.c_str(), -1,
                         &candidate_q_statement, nullptr) == SQLITE_OK &&
      sqlite3_step(candidate_q_statement) == SQLITE_ROW) {
    const auto* q_text = sqlite3_column_text(candidate_q_statement, 0);
    output << (q_text ? reinterpret_cast<const char*>(q_text) : "");
    output << '\n' << "q_squared_value="
           << (sqlite3_column_type(candidate_q_statement, 1) == SQLITE_NULL
                   ? 0.0 : sqlite3_column_double(candidate_q_statement, 1));
  } else {
    output << "\nq_squared_value=";
  }
  sqlite3_finalize(candidate_q_statement);
  output
         << "boundary_length_dsigma=" << rat_string(boundary.length) << '\n'
         << "boundary_ix=" << rat_string(boundary.ix) << '\n'
         << "boundary_iy=" << rat_string(boundary.iy) << '\n'
         << "ell_P=" << rat_string(selected.ell[0]) << " + ("
         << rat_string(selected.ell[1]) << ")x + ("
         << rat_string(selected.ell[2]) << ")y\n"
         << "first_verified_key=" << summary.first_verified_key << '\n'
         << "best_twice_area=" << summary.best_twice_area << '\n'
         << "database_candidates=" << count_candidates(dimension) << '\n'
         << "database_verified=" << count_status("verified_unstable", dimension) << '\n'
         << "database_unverified=" << count_status("unverified", dimension) << '\n';

  const std::string sql =
      "SELECT profile,status,value,witness_ux,witness_uy,witness_t,"
      "exact_a,exact_b,exact_c,exact_value,q_squared_exact,q_squared_value FROM attempts WHERE candidate_key=" +
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
      if (sqlite3_column_type(statement, 10) != SQLITE_NULL) {
        output << "attempt_q_squared_exact=" << text_column(10) << '\n'
               << "attempt_q_squared_value=" << sqlite3_column_double(statement, 11) << '\n';
      }
    }
  }
  sqlite3_finalize(statement);
}

SearchSummary run_search(const AreaSearchOptions& options) {
  if (options.d < 3 || options.initial_N < 1 || options.initial_M < 1 ||
      options.time_limit_seconds <= 0 || options.shell_seconds <= 0)
    throw std::invalid_argument("invalid search options");
  SearchDatabase database(options.database);
  database.ensure_generator_revision("incremental-convex-v1");
  const std::string validation_profile = current_validation_profile(options);
  const std::string rng_state_name = dimension_state_name(options.d, "rng");
  // Shell is intentionally reset to 0 on every run. Random search does not
  // exhaust a shell, so resuming from a previously saved high shell would skip
  // the small-volume region. The database still deduplicates candidates and
  // resumes their validation stages.
  std::mt19937_64 rng(options.seed);
  restore_rng(rng, database.load_state(rng_state_name));
  std::unordered_set<std::string> queued;
  std::vector<PolygonCandidate> area_frontier;
  std::vector<PolygonCandidate> score_frontier;
  const auto previous_minimum = database.min_verified_twice_area(options.d);
  std::unordered_set<std::string> new_tested_keys;
  std::unordered_set<std::string> new_verified_keys;
  SearchSummary summary;
  const auto started = std::chrono::steady_clock::now();
  const auto deadline = started + std::chrono::duration<double>(options.time_limit_seconds);
  // Always start from shell=0. The database itself records which candidates
  // have already been tested; random search never exhausts a shell, so there
  // is no benefit in resuming from a previously saved high shell.
  std::uint64_t shell = 0;
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
  const auto note_tested = [&](const PolygonCandidate& candidate) {
    new_tested_keys.insert(candidate.key);
  };
  const auto note_verified = [&](const PolygonCandidate& candidate) {
    new_verified_keys.insert(candidate.key);
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
    summary.unverified = database.count_status("unverified", options.d);
    summary.total_tested = database.count_tested(options.d);
    summary.new_tested = new_tested_keys.size();
    summary.total_verified_unstable = database.count_verified(options.d);
    summary.new_verified_unstable = new_verified_keys.size();
    summary.verified = summary.total_verified_unstable;
    summary.top_verified = database.top_verified(options.d, 5);
    const auto current_minimum = database.min_verified_twice_area(options.d);
    summary.smaller_volume_found =
        current_minimum.has_value() &&
        (!previous_minimum.has_value() || current_minimum.value() < previous_minimum.value());
    if (!summary.top_verified.empty()) {
      summary.have_verified = true;
      summary.best_twice_area = summary.top_verified.front().twice_area;
      summary.best_verified_key = summary.top_verified.front().key;
    }
    database.write_report(options.output_directory / "k_stability_search_result.txt",
                          summary, options.d);
    return summary;
  };

  const auto run_probe = [&](PolygonCandidate& candidate, bool resume_existing) {
    if (std::chrono::steady_clock::now() >= deadline) {
      stop_requested = true;
      return false;
    }
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
    note_tested(candidate);
    DetectorOutcome outcome = detect_candidate(candidate, profile,
                                                options.certify_max_denominator);
    outcome.profile = stage_profile(DetectorTier::probe);
    database.save_attempt(candidate, outcome);
    ++summary.probes;
    candidate.probe_score = outcome.numerical.witness.value;
    database.update_probe_score(candidate.key, candidate.probe_score);
    if (outcome.verified_unstable) {
      note_verified(candidate);
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
    if (std::chrono::steady_clock::now() >= deadline) {
      stop_requested = true;
      return;
    }
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
    note_tested(candidate);
    DetectorOutcome outcome = detect_candidate(candidate, profile,
                                                options.certify_max_denominator);
    outcome.profile = profile_key;
    database.save_attempt(candidate, outcome);
    ++summary.finals;
    if (outcome.verified_unstable) {
      note_verified(candidate);
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
    if (std::chrono::steady_clock::now() >= deadline) {
      stop_requested = true;
      return;
    }
    generated.insert(candidate.key);
    if (count_as_generated) ++summary.generated;
    if (options.smooth_only && !candidate_is_smooth(candidate)) {
      ++summary.skipped;
      return;
    }
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

  for (auto candidate : database.load_candidates(options.d)) {
    prepare_loaded_candidate(std::move(candidate), false);
  }
  if (options.stop_on_first && summary.have_verified) stop_requested = true;

  while (!stop_requested && std::chrono::steady_clock::now() < deadline) {
    const auto n_bound = shell_value(options.initial_N, shell, true);
    const auto m_bound = shell_value(options.initial_M, shell, false);
    if (!n_bound || !m_bound) {
      if (options.verbose)
        std::cerr << "search stopped: shell bounds exceed int range\n";
      break;
    }
    if (options.verbose)
      std::cerr << "search shell=" << shell << " N=" << *n_bound
                << " M=" << *m_bound << '\n';
    const auto shell_deadline = std::min(deadline, std::chrono::steady_clock::now() +
        std::chrono::duration<double>(options.shell_seconds));
    if (std::chrono::steady_clock::now() >= deadline) {
      stop_requested = true;
      break;
    }
    const auto directions = primitive_directions(*n_bound);
    auto sample_step = [&](std::int64_t upper_bound) {
      const std::int64_t small = std::min<std::int64_t>(upper_bound, 2);
      if (std::uniform_int_distribution<int>(0, 99)(rng) < 90)
        return std::uniform_int_distribution<std::int64_t>(1, small)(rng);
      return std::uniform_int_distribution<std::int64_t>(1, upper_bound)(rng);
    };
    while (!stop_requested && std::chrono::steady_clock::now() < shell_deadline &&
           std::chrono::steady_clock::now() < deadline) {
      for (int batch = 0; batch < options.beam_width && !stop_requested &&
           std::chrono::steady_clock::now() < shell_deadline; ++batch) {
        std::vector<Direction> dirs{{1, 0}};
        std::vector<std::int64_t> steps{sample_step(*m_bound)};
        std::vector<IntPoint> vertices{{0, 0}, {steps.front(), 0}};
        for (int i = 1; i < options.d - 1; ++i) {
          struct FeasibleDirection {
            Direction direction;
            std::int64_t max_step;
          };
          std::vector<FeasibleDirection> feasible;
          const Direction previous = dirs.back();
          for (const Direction p : directions) {
            if (cross(previous, p) <= 0 || cross_to_origin(p, vertices.back()) <= 0)
              continue;
            if (options.smooth_only && abs_wide(cross(previous, p)) != 1)
              continue;
            std::int64_t max_step = *m_bound;
            if (p.y < 0) {
              const WideInt numerator = static_cast<WideInt>(vertices.back().y) - 1;
              const WideInt denominator = -static_cast<WideInt>(p.y);
              max_step = static_cast<std::int64_t>(numerator / denominator);
              max_step = std::min(max_step, static_cast<std::int64_t>(*m_bound));
            }
            if (max_step >= 1) feasible.push_back({p, max_step});
          }
          if (feasible.empty()) break;
          std::vector<double> direction_weights;
          direction_weights.reserve(feasible.size());
          for (const auto& option : feasible) {
            const auto coordinate_norm = std::max(std::llabs(option.direction.x),
                                                  std::llabs(option.direction.y));
            const double denominator = 1.0 + static_cast<double>(coordinate_norm);
            direction_weights.push_back(1.0 / (denominator * denominator));
          }
          const auto pick = std::discrete_distribution<std::size_t>(
              direction_weights.begin(), direction_weights.end())(rng);
          const auto selected = feasible[pick];
          const std::int64_t step = sample_step(selected.max_step);
          const WideInt x = static_cast<WideInt>(vertices.back().x) +
                            static_cast<WideInt>(step) * selected.direction.x;
          const WideInt y = static_cast<WideInt>(vertices.back().y) +
                            static_cast<WideInt>(step) * selected.direction.y;
          if (x < std::numeric_limits<std::int64_t>::min() ||
              x > std::numeric_limits<std::int64_t>::max() ||
              y < std::numeric_limits<std::int64_t>::min() ||
              y > std::numeric_limits<std::int64_t>::max() || y <= 0)
            break;
          dirs.push_back(selected.direction);
          steps.push_back(step);
          vertices.push_back({static_cast<std::int64_t>(x),
                              static_cast<std::int64_t>(y)});
        }
        if (dirs.size() != static_cast<std::size_t>(options.d - 1) ||
            steps.size() != static_cast<std::size_t>(options.d - 1)) {
          ++summary.rejected;
          continue;
        }
        PolygonCandidate candidate;
        if (build_candidate(options.d, dirs, steps, candidate)) {
          if (options.smooth_only && !candidate_is_smooth(candidate)) {
            ++summary.rejected;
            continue;
          }
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
      for (int take = 0; take < 8 && !stop_requested &&
           std::chrono::steady_clock::now() < deadline; ++take) {
        const bool use_score = take == 7;
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
          note_tested(*candidate);
          DetectorOutcome outcome = detect_candidate(*candidate, profile,
                                                     options.certify_max_denominator);
          outcome.profile = confirm_key;
          database.save_attempt(*candidate, outcome);
          ++summary.confirms;
          if (outcome.verified_unstable) {
            note_verified(*candidate);
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
      database.save_state(rng_state_name, serialize_rng(rng));
    }
    ++shell;
  }
  return finish();
}

}  // namespace kstab
