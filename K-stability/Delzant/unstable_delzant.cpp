#include "unstable_delzant.hpp"

#include <sqlite3.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace delzant {
namespace {

using Wide = __int128_t;

Wide cross(Direction a, Direction b) {
  return static_cast<Wide>(a.x) * b.y - static_cast<Wide>(a.y) * b.x;
}

std::int64_t gcd_abs(std::int64_t a, std::int64_t b) {
  const auto aa = a < 0 ? static_cast<std::uint64_t>(-(static_cast<Wide>(a)))
                        : static_cast<std::uint64_t>(a);
  const auto bb = b < 0 ? static_cast<std::uint64_t>(-(static_cast<Wide>(b)))
                        : static_cast<std::uint64_t>(b);
  return static_cast<std::int64_t>(std::gcd(aa, bb));
}

std::string candidate_key(const Candidate& candidate) {
  return canonical_candidate_key(candidate.fan, candidate.edge_lengths);
}

std::int64_t linf_diameter(const std::vector<kstab::IntPoint>& points) {
  Wide diameter = 0;
  for (std::size_t i = 0; i < points.size(); ++i) {
    for (std::size_t j = i + 1; j < points.size(); ++j) {
      Wide dx = static_cast<Wide>(points[i].x) - points[j].x;
      Wide dy = static_cast<Wide>(points[i].y) - points[j].y;
      if (dx < 0) dx = -dx;
      if (dy < 0) dy = -dy;
      diameter = std::max(diameter, std::max(dx, dy));
    }
  }
  if (diameter > std::numeric_limits<std::int64_t>::max())
    return std::numeric_limits<std::int64_t>::max();
  return static_cast<std::int64_t>(diameter);
}

constexpr const char* kDetectorProfile =
    "delzant-staged-v1|probe=32x16:norefine|confirm=128x64:refine|"
    "final=360x256:refine|negative-cutoffs=-1e-8,-1e-7|audit=fnv1a-5pct";

std::uint64_t stable_hash(const std::string& text) {
  std::uint64_t hash = 14695981039346656037ULL;
  for (const unsigned char byte : text) {
    hash ^= byte;
    hash *= 1099511628211ULL;
  }
  return hash;
}

bool audit_candidate(const std::string& key) {
  return should_audit_candidate(key);
}

bool within_bound(const std::vector<Direction>& fan, int bound) {
  return std::all_of(fan.begin(), fan.end(), [bound](Direction p) {
    return std::llabs(p.x) <= bound && std::llabs(p.y) <= bound;
  });
}

bool insert_random_blowups(std::vector<Direction>& fan, int d, int N,
                           std::mt19937_64& rng) {
  while (static_cast<int>(fan.size()) < d) {
    if (fan.size() > 10000) return false;
    std::uniform_int_distribution<std::size_t> pick(0, fan.size() - 1);
    const std::size_t i = pick(rng);
    const Direction a = fan[i];
    const Direction b = fan[(i + 1) % fan.size()];
    const Wide x = static_cast<Wide>(a.x) + b.x;
    const Wide y = static_cast<Wide>(a.y) + b.y;
    if (x < -N || x > N || y < -N || y > N) return false;
    fan.insert(fan.begin() + static_cast<std::ptrdiff_t>(i + 1),
               {static_cast<std::int64_t>(x), static_cast<std::int64_t>(y)});
  }
  return is_smooth_complete_fan(fan);
}

std::string points_text(const std::vector<kstab::IntPoint>& points) {
  std::ostringstream out;
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (i) out << ';';
    out << points[i].x << ':' << points[i].y;
  }
  return out.str();
}

std::string steps_text(const std::vector<std::int64_t>& steps) {
  std::ostringstream out;
  for (std::size_t i = 0; i < steps.size(); ++i) {
    if (i) out << ',';
    out << steps[i];
  }
  return out.str();
}

std::string fan_text(const std::vector<Direction>& fan) {
  std::ostringstream out;
  for (std::size_t i = 0; i < fan.size(); ++i) {
    if (i) out << ';';
    out << fan[i].x << ':' << fan[i].y;
  }
  return out.str();
}

std::string rat_text(const kstab::Rational& r) {
  std::ostringstream out;
  out << r;
  return out.str();
}

class Database {
 public:
  explicit Database(const std::filesystem::path& path) {
    if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());
    if (sqlite3_open(path.string().c_str(), &db_) != SQLITE_OK)
      throw std::runtime_error("cannot open Delzant database: " +
                               std::string(sqlite3_errmsg(db_)));
    exec("PRAGMA journal_mode=WAL;");
    exec("CREATE TABLE IF NOT EXISTS candidates("
         "key TEXT PRIMARY KEY,d INTEGER,fan TEXT,steps TEXT,vertices TEXT,"
         "twice_area INTEGER,ell0 TEXT,ell1 TEXT,ell2 TEXT,status TEXT,"
         "witness_a TEXT,witness_b TEXT,witness_c TEXT,df_value TEXT,"
         "detector_profile TEXT NOT NULL DEFAULT '',stage TEXT NOT NULL DEFAULT 'legacy_complete',"
         "probe_normalized REAL,confirm_normalized REAL,final_normalized REAL,"
         "probe_evaluations INTEGER NOT NULL DEFAULT 0,"
         "confirm_evaluations INTEGER NOT NULL DEFAULT 0,"
         "final_evaluations INTEGER NOT NULL DEFAULT 0);");
    ensure_column("detector_profile", "TEXT NOT NULL DEFAULT ''");
    ensure_column("stage", "TEXT NOT NULL DEFAULT 'legacy_complete'");
    ensure_column("probe_normalized", "REAL");
    ensure_column("confirm_normalized", "REAL");
    ensure_column("final_normalized", "REAL");
    ensure_column("probe_evaluations", "INTEGER NOT NULL DEFAULT 0");
    ensure_column("confirm_evaluations", "INTEGER NOT NULL DEFAULT 0");
    ensure_column("final_evaluations", "INTEGER NOT NULL DEFAULT 0");
    exec("CREATE TABLE IF NOT EXISTS state(name TEXT PRIMARY KEY,value TEXT);");
  }
  ~Database() { if (db_) sqlite3_close(db_); }
  Database(const Database&) = delete;
  Database& operator=(const Database&) = delete;

  void exec(const std::string& sql) {
    char* error = nullptr;
    if (sqlite3_exec(db_, sql.c_str(), nullptr, nullptr, &error) != SQLITE_OK) {
      const std::string message = error ? error : sqlite3_errmsg(db_);
      sqlite3_free(error);
      throw std::runtime_error("Delzant database error: " + message);
    }
  }

  void save_candidate(const Candidate& c, const std::string& status,
                      const kstab::CertifyResult* cert,
                      const std::string& profile, const std::string& stage,
                      const std::optional<double>& probe_value,
                      const std::optional<double>& confirm_value,
                      const std::optional<double>& final_value,
                      std::uint64_t probe_evaluations,
                      std::uint64_t confirm_evaluations,
                      std::uint64_t final_evaluations) {
    sqlite3_stmt* stmt = nullptr;
    const char* sql = "INSERT INTO candidates(key,d,fan,steps,vertices,twice_area,"
                      "ell0,ell1,ell2,status,witness_a,witness_b,witness_c,df_value,"
                      "detector_profile,stage,probe_normalized,confirm_normalized,"
                      "final_normalized,probe_evaluations,confirm_evaluations,"
                      "final_evaluations) "
                      "VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) "
                      "ON CONFLICT(key) DO UPDATE SET status=excluded.status,"
                      "witness_a=excluded.witness_a,witness_b=excluded.witness_b,"
                      "witness_c=excluded.witness_c,df_value=excluded.df_value,"
                      "detector_profile=excluded.detector_profile,stage=excluded.stage,"
                      "probe_normalized=excluded.probe_normalized,"
                      "confirm_normalized=excluded.confirm_normalized,"
                      "final_normalized=excluded.final_normalized,"
                      "probe_evaluations=excluded.probe_evaluations,"
                      "confirm_evaluations=excluded.confirm_evaluations,"
                      "final_evaluations=excluded.final_evaluations";
    if (sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr) != SQLITE_OK)
      throw std::runtime_error(sqlite3_errmsg(db_));
    const std::string key = c.key;
    const std::string fan = fan_text(c.fan);
    const std::string steps = steps_text(c.edge_lengths);
    const std::string vertices = points_text(c.vertices);
    const std::string e0 = rat_text(c.ell[0]);
    const std::string e1 = rat_text(c.ell[1]);
    const std::string e2 = rat_text(c.ell[2]);
    sqlite3_bind_text(stmt, 1, key.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_int(stmt, 2, static_cast<int>(c.fan.size()));
    sqlite3_bind_text(stmt, 3, fan.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 4, steps.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 5, vertices.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_int64(stmt, 6, c.twice_area);
    sqlite3_bind_text(stmt, 7, e0.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 8, e1.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 9, e2.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 10, status.c_str(), -1, SQLITE_TRANSIENT);
    if (cert && cert->certified) {
      for (int i = 0; i < 3; ++i) {
        const std::string coeff = rat_text(cert->coefficients[i]);
        sqlite3_bind_text(stmt, 11 + i, coeff.c_str(), -1, SQLITE_TRANSIENT);
      }
      const std::string value = rat_text(cert->value);
      sqlite3_bind_text(stmt, 14, value.c_str(), -1, SQLITE_TRANSIENT);
    } else {
      for (int i = 11; i <= 14; ++i) sqlite3_bind_null(stmt, i);
    }
    sqlite3_bind_text(stmt, 15, profile.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 16, stage.c_str(), -1, SQLITE_TRANSIENT);
    auto bind_optional = [&](int index, const std::optional<double>& value) {
      if (value) sqlite3_bind_double(stmt, index, *value);
      else sqlite3_bind_null(stmt, index);
    };
    bind_optional(17, probe_value);
    bind_optional(18, confirm_value);
    bind_optional(19, final_value);
    sqlite3_bind_int64(stmt, 20, static_cast<sqlite3_int64>(probe_evaluations));
    sqlite3_bind_int64(stmt, 21, static_cast<sqlite3_int64>(confirm_evaluations));
    sqlite3_bind_int64(stmt, 22, static_cast<sqlite3_int64>(final_evaluations));
    if (sqlite3_step(stmt) != SQLITE_DONE) {
      const std::string error = sqlite3_errmsg(db_);
      sqlite3_finalize(stmt);
      throw std::runtime_error(error);
    }
    sqlite3_finalize(stmt);
  }

  std::unordered_set<std::string> keys(int d, const std::string& profile) const {
    std::unordered_set<std::string> result;
    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db_, "SELECT key FROM candidates WHERE d=? AND "
                               "(detector_profile=? OR (detector_profile='' AND "
                               "status IN ('unverified','verified_unstable')))", -1,
                           &stmt, nullptr) != SQLITE_OK)
      throw std::runtime_error(sqlite3_errmsg(db_));
    sqlite3_bind_int(stmt, 1, d);
    sqlite3_bind_text(stmt, 2, profile.c_str(), -1, SQLITE_TRANSIENT);
    while (sqlite3_step(stmt) == SQLITE_ROW) {
      const auto* value = sqlite3_column_text(stmt, 0);
      if (value) result.emplace(reinterpret_cast<const char*>(value));
    }
    sqlite3_finalize(stmt);
    return result;
  }

  std::string state(const std::string& name) const {
    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db_, "SELECT value FROM state WHERE name=?", -1,
                           &stmt, nullptr) != SQLITE_OK)
      throw std::runtime_error(sqlite3_errmsg(db_));
    sqlite3_bind_text(stmt, 1, name.c_str(), -1, SQLITE_TRANSIENT);
    std::string result;
    if (sqlite3_step(stmt) == SQLITE_ROW) {
      const auto* value = sqlite3_column_text(stmt, 0);
      if (value) result = reinterpret_cast<const char*>(value);
    }
    sqlite3_finalize(stmt);
    return result;
  }

  void save_state(const std::string& name, const std::string& value) {
    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db_, "INSERT INTO state(name,value) VALUES(?,?) "
                                "ON CONFLICT(name) DO UPDATE SET value=excluded.value",
                           -1, &stmt, nullptr) != SQLITE_OK)
      throw std::runtime_error(sqlite3_errmsg(db_));
    sqlite3_bind_text(stmt, 1, name.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 2, value.c_str(), -1, SQLITE_TRANSIENT);
    if (sqlite3_step(stmt) != SQLITE_DONE) {
      const std::string error = sqlite3_errmsg(db_);
      sqlite3_finalize(stmt);
      throw std::runtime_error(error);
    }
    sqlite3_finalize(stmt);
  }

  bool first_verified(int d, std::int64_t max_diameter, Candidate& candidate,
                      kstab::CertifyResult& certificate) const {
    sqlite3_stmt* stmt = nullptr;
    const char* sql = "SELECT fan,steps,vertices,twice_area,ell0,ell1,ell2,"
                      "witness_a,witness_b,witness_c,df_value,key "
                      "FROM candidates WHERE d=? AND status='verified_unstable' "
                      "ORDER BY rowid";
    if (sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr) != SQLITE_OK)
      throw std::runtime_error(sqlite3_errmsg(db_));
    sqlite3_bind_int(stmt, 1, d);
    while (sqlite3_step(stmt) == SQLITE_ROW) {
    candidate = Candidate{};
    certificate = kstab::CertifyResult{};
    auto column = [&](int index) {
      const auto* value = sqlite3_column_text(stmt, index);
      return value ? std::string(reinterpret_cast<const char*>(value)) : std::string();
    };
    auto parse_rational = [](const std::string& value) {
      std::istringstream in(value);
      kstab::Rational result;
      in >> result;
      if (!in) throw std::runtime_error("invalid rational in Delzant database");
      return result;
    };
    std::istringstream fan_in(column(0));
    std::string item;
    while (std::getline(fan_in, item, ';')) {
      const auto colon = item.find(':');
      if (colon == std::string::npos) throw std::runtime_error("invalid stored fan");
      candidate.fan.push_back({std::stoll(item.substr(0, colon)),
                               std::stoll(item.substr(colon + 1))});
    }
    std::istringstream steps_in(column(1));
    while (std::getline(steps_in, item, ',')) candidate.edge_lengths.push_back(std::stoll(item));
    std::istringstream vertices_in(column(2));
    while (std::getline(vertices_in, item, ';')) {
      const auto colon = item.find(':');
      if (colon == std::string::npos) throw std::runtime_error("invalid stored vertices");
      candidate.vertices.push_back({std::stoll(item.substr(0, colon)),
                                    std::stoll(item.substr(colon + 1))});
    }
    candidate.twice_area = sqlite3_column_int64(stmt, 3);
    for (int i = 0; i < 3; ++i) candidate.ell[i] = parse_rational(column(4 + i));
    certificate.certified = true;
    for (int i = 0; i < 3; ++i) certificate.coefficients[i] = parse_rational(column(7 + i));
    certificate.value = parse_rational(column(10));
    candidate.key = column(11);
    if (linf_diameter(candidate.vertices) > max_diameter) {
      continue;
    }
    sqlite3_finalize(stmt);
    return true;
    }
    sqlite3_finalize(stmt);
    return false;
  }

 private:
  void ensure_column(const std::string& name, const std::string& declaration) {
    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db_, "PRAGMA table_info(candidates)", -1,
                           &stmt, nullptr) != SQLITE_OK)
      throw std::runtime_error(sqlite3_errmsg(db_));
    bool found = false;
    while (sqlite3_step(stmt) == SQLITE_ROW) {
      const auto* column = sqlite3_column_text(stmt, 1);
      if (column && name == reinterpret_cast<const char*>(column)) {
        found = true;
        break;
      }
    }
    sqlite3_finalize(stmt);
    if (!found) exec("ALTER TABLE candidates ADD COLUMN " + name + " " + declaration);
  }

  sqlite3* db_ = nullptr;
};

std::string serialize_rng(const std::mt19937_64& rng) {
  std::ostringstream out;
  out << rng;
  return out.str();
}

void restore_rng(std::mt19937_64& rng, const std::string& state) {
  if (!state.empty()) {
    std::istringstream in(state);
    in >> rng;
    if (!in) throw std::runtime_error("invalid RNG state in Delzant database");
  }
}

bool write_result(const SearchOptions& options, const Candidate& candidate,
                  const kstab::CertifyResult& cert) {
  std::filesystem::create_directories(options.output_directory);
  const auto stem = options.output_directory /
                    ("unstable_d" + std::to_string(options.d));
  const auto polygon_path = stem.string() + ".polygon";
  std::ofstream polygon(polygon_path);
  if (!polygon) return false;
  polygon << "# relatively K-unstable smooth lattice polygon\n"
          << "# twice_area=" << candidate.twice_area << "\n"
          << "# linf_diameter=" << linf_diameter(candidate.vertices) << "\n"
          << "# ell_P(x,y)=" << candidate.ell[0] << " + (" << candidate.ell[1]
          << ")x + (" << candidate.ell[2] << ")y\n"
          << "# certified_witness=max{" << cert.coefficients[0] << "x + "
          << cert.coefficients[1] << "y + " << cert.coefficients[2] << ",0}\n"
          << "# certified_M_l=" << cert.value << "\n";
  for (const auto& p : candidate.vertices) polygon << p.x << ' ' << p.y << '\n';
  polygon.close();
  if (!polygon) return false;
  const double a = kstab::rational_to_double(cert.coefficients[0]);
  const double b = kstab::rational_to_double(cert.coefficients[1]);
  const double c = kstab::rational_to_double(cert.coefficients[2]);
  const double norm = std::hypot(a, b);
  kstab::Witness witness;
  witness.ux = a / norm;
  witness.uy = b / norm;
  witness.t = -c / norm;
  witness.value = kstab::rational_to_double(cert.value);
  if (!kstab::write_svg(stem.string() + ".svg", candidate.vertices,
                        candidate.ell, &witness)) return false;
  std::ofstream report(stem.string() + ".txt");
  if (!report) return false;
  report << "status=verified_unstable\nvertices=" << candidate.vertices.size()
         << "\ntwice_area=" << candidate.twice_area << "\nfan="
         << fan_text(candidate.fan) << "\nedge_lengths="
         << steps_text(candidate.edge_lengths) << "\nell_P=" << candidate.ell[0]
         << " + (" << candidate.ell[1] << ")x + (" << candidate.ell[2]
         << ")y\nwitness=" << cert.coefficients[0] << " "
         << cert.coefficients[1] << " " << cert.coefficients[2]
         << "\ncertified_M_l=" << cert.value << '\n';
  return static_cast<bool>(report);
}

}  // namespace

std::string canonical_candidate_key(
    const std::vector<Direction>& fan,
    const std::vector<std::int64_t>& edge_lengths) {
  if (fan.empty() || fan.size() != edge_lengths.size())
    throw std::invalid_argument("fan and edge-length dimensions do not match");
  std::int64_t common = 0;
  for (const auto length : edge_lengths) {
    if (length <= 0) throw std::invalid_argument("edge lengths must be positive");
    common = std::gcd(common, length);
  }
  std::string best;
  const std::size_t d = fan.size();
  for (std::size_t offset = 0; offset < d; ++offset) {
    std::ostringstream out;
    for (std::size_t j = 0; j < d; ++j) {
      const std::size_t i = (offset + j) % d;
      if (j) out << ';';
      out << fan[i].x << ':' << fan[i].y << ':' << edge_lengths[i] / common;
    }
    const std::string current = out.str();
    if (best.empty() || current < best) best = current;
  }
  return "d" + std::to_string(d) + "|" + best;
}

bool should_audit_candidate(const std::string& canonical_key) {
  return stable_hash(canonical_key) % 20 == 0;
}

bool is_smooth_complete_fan(const std::vector<Direction>& fan) {
  if (fan.size() < 3) return false;
  long double total_turn = 0.0L;
  for (std::size_t i = 0; i < fan.size(); ++i) {
    const Direction a = fan[i];
    const Direction b = fan[(i + 1) % fan.size()];
    if ((a.x == 0 && a.y == 0) || gcd_abs(a.x, a.y) != 1) return false;
    if (cross(a, b) != 1) return false;
    const long double dot = static_cast<long double>(a.x) * b.x +
                            static_cast<long double>(a.y) * b.y;
    total_turn += std::atan2(1.0L, dot);
  }
  constexpr long double two_pi = 6.283185307179586476925286766559005768L;
  return std::fabs(total_turn - two_pi) < 1e-12L;
}

bool is_smooth_lattice_polygon(const Candidate& candidate) {
  const std::size_t d = candidate.fan.size();
  if (d < 3 || candidate.edge_lengths.size() != d ||
      candidate.vertices.size() != d || !is_smooth_complete_fan(candidate.fan))
    return false;
  Wide sx = 0;
  Wide sy = 0;
  for (std::size_t i = 0; i < d; ++i) {
    if (candidate.edge_lengths[i] <= 0) return false;
    sx += static_cast<Wide>(candidate.edge_lengths[i]) * candidate.fan[i].x;
    sy += static_cast<Wide>(candidate.edge_lengths[i]) * candidate.fan[i].y;
    const auto& p = candidate.vertices[i];
    const auto& q = candidate.vertices[(i + 1) % d];
    if (static_cast<Wide>(q.x) - p.x !=
            static_cast<Wide>(candidate.edge_lengths[i]) * candidate.fan[i].x ||
        static_cast<Wide>(q.y) - p.y !=
            static_cast<Wide>(candidate.edge_lengths[i]) * candidate.fan[i].y)
      return false;
  }
  if (sx != 0 || sy != 0) return false;
  for (std::size_t i = 0; i < d; ++i)
    if (cross(candidate.fan[i], candidate.fan[(i + 1) % d]) != 1)
      return false;
  return true;
}

bool generate_fan(int d, int N, std::mt19937_64& rng,
                  std::vector<Direction>& fan) {
  if (d < 3 || N < 1) return false;
  if (d == 3) {
    fan = {{1, 0}, {0, 1}, {-1, -1}};
    return within_bound(fan, N);
  }
  std::uniform_int_distribution<int> choose_base(0, N);
  // Include P^2 and Hirzebruch fans F_n; repeated corner blowups generate
  // smooth complete fans while preserving det(v_i,v_{i+1})=1.
  for (int trial = 0; trial < 128; ++trial) {
    const bool use_p2 = d > 4 && choose_base(rng) % 2 == 0;
    if (use_p2) fan = {{1, 0}, {0, 1}, {-1, -1}};
    else {
      const int n = choose_base(rng);
      fan = {{1, 0}, {0, 1}, {-1, n}, {0, -1}};
    }
    if (static_cast<int>(fan.size()) > d || !within_bound(fan, N)) continue;
    if (insert_random_blowups(fan, d, N, rng)) return true;
  }
  fan.clear();
  return false;
}

bool build_candidate(const std::vector<Direction>& fan, int M,
                     std::mt19937_64& rng, Candidate& candidate,
                     std::int64_t max_diameter, bool* exceeded_max_diameter) {
  if (exceeded_max_diameter) *exceeded_max_diameter = false;
  const std::size_t d = fan.size();
  if (M < 1 || max_diameter < 1 || d < 3 || !is_smooth_complete_fan(fan)) return false;
  std::vector<std::int64_t> lengths(d, 0);
  std::uniform_int_distribution<int> choose_step(1, M);
  Wide sx = 0;
  Wide sy = 0;
  for (std::size_t i = 0; i + 2 < d; ++i) {
    lengths[i] = choose_step(rng);
    sx += static_cast<Wide>(lengths[i]) * fan[i].x;
    sy += static_cast<Wide>(lengths[i]) * fan[i].y;
  }
  const Direction a = fan[d - 2];
  const Direction b = fan[d - 1];
  const Wide det = cross(a, b);
  if (det != 1) return false;
  const Wide last_a = -(sx * b.y - sy * b.x);
  const Wide last_b = -(static_cast<Wide>(a.x) * sy -
                        static_cast<Wide>(a.y) * sx);
  if (last_a <= 0 || last_b <= 0 ||
      last_a > std::numeric_limits<std::int64_t>::max() ||
      last_b > std::numeric_limits<std::int64_t>::max()) return false;
  lengths[d - 2] = static_cast<std::int64_t>(last_a);
  lengths[d - 1] = static_cast<std::int64_t>(last_b);
  std::int64_t common = 0;
  for (auto length : lengths) common = std::gcd(common, length);
  if (common > 1) for (auto& length : lengths) length /= common;

  // Rotate the fan and its lengths to a canonical representation.
  std::size_t best_offset = 0;
  std::string best;
  for (std::size_t offset = 0; offset < d; ++offset) {
    std::ostringstream key;
    for (std::size_t j = 0; j < d; ++j) {
      const auto i = (offset + j) % d;
      if (j) key << ';';
      key << fan[i].x << ':' << fan[i].y << ':' << lengths[i];
    }
    if (best.empty() || key.str() < best) {
      best = key.str();
      best_offset = offset;
    }
  }
  Candidate result;
  result.fan.reserve(d);
  result.edge_lengths.reserve(d);
  for (std::size_t j = 0; j < d; ++j) {
    const std::size_t i = (best_offset + j) % d;
    result.fan.push_back(fan[i]);
    result.edge_lengths.push_back(lengths[i]);
  }
  result.vertices.push_back({0, 0});
  Wide area = 0;
  for (std::size_t i = 0; i < d; ++i) {
    const auto& p = result.vertices[i];
    const Wide x = static_cast<Wide>(p.x) +
                   static_cast<Wide>(result.edge_lengths[i]) * result.fan[i].x;
    const Wide y = static_cast<Wide>(p.y) +
                   static_cast<Wide>(result.edge_lengths[i]) * result.fan[i].y;
    if (i + 1 < d) {
      if (x < std::numeric_limits<std::int64_t>::min() ||
          x > std::numeric_limits<std::int64_t>::max() ||
          y < std::numeric_limits<std::int64_t>::min() ||
          y > std::numeric_limits<std::int64_t>::max()) return false;
      result.vertices.push_back({static_cast<std::int64_t>(x),
                                 static_cast<std::int64_t>(y)});
    } else if (x != 0 || y != 0) return false;
  }
  if (linf_diameter(result.vertices) > max_diameter) {
    if (exceeded_max_diameter) *exceeded_max_diameter = true;
    return false;
  }
  for (std::size_t i = 0; i < d; ++i)
    area += static_cast<Wide>(result.vertices[i].x) *
                result.vertices[(i + 1) % d].y -
            static_cast<Wide>(result.vertices[i].y) *
                result.vertices[(i + 1) % d].x;
  if (area <= 0 || area > std::numeric_limits<std::int64_t>::max()) return false;
  result.twice_area = static_cast<std::int64_t>(area);
  result.key = candidate_key(result);
  if (!is_smooth_lattice_polygon(result)) return false;
  result.ell = kstab::compute_ell_p(result.vertices);
  candidate = std::move(result);
  return true;
}

SearchSummary run_search(const SearchOptions& options) {
  const auto search_start = std::chrono::steady_clock::now();
  Database database(options.database);
  SearchSummary summary;
  if (options.max_diameter < 1)
    throw std::invalid_argument("max_diameter must be positive");
  if (database.first_verified(options.d, options.max_diameter,
                              summary.candidate, summary.certificate)) {
    summary.found = true;
    if (!write_result(options, summary.candidate, summary.certificate))
      throw std::runtime_error("failed to write existing certified Delzant result files");
    summary.total_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - search_start).count();
    return summary;
  }
  const std::string profile = std::string(kDetectorProfile) +
      "|certify=" + std::to_string(options.certify_max_denominator) +
      "|final=" + std::to_string(options.theta_steps) + "x" +
      std::to_string(options.t_steps);
  std::mt19937_64 rng(options.seed);
  restore_rng(rng, database.state("d" + std::to_string(options.d) + "|rng"));
  const auto known = database.keys(options.d, profile);
  std::unordered_set<std::string> seen = known;
  const auto deadline = std::chrono::steady_clock::now() +
                        std::chrono::duration<double>(options.time_limit_seconds);
  kstab::SearchOptions probe_options;
  probe_options.theta_steps = 32;
  probe_options.t_steps = 16;
  probe_options.refine = false;
  kstab::SearchOptions confirm_options;
  confirm_options.theta_steps = 128;
  confirm_options.t_steps = 64;
  confirm_options.refine = true;
  kstab::SearchOptions final_options;
  final_options.theta_steps = options.theta_steps;
  final_options.t_steps = options.t_steps;
  final_options.refine = true;

  while (std::chrono::steady_clock::now() < deadline) {
    const auto generation_start = std::chrono::steady_clock::now();
    std::vector<Direction> fan;
    if (!generate_fan(options.d, options.N, rng, fan)) {
      summary.generation_seconds += std::chrono::duration<double>(
          std::chrono::steady_clock::now() - generation_start).count();
      continue;
    }
    Candidate candidate;
    bool oversized = false;
    if (!build_candidate(fan, options.M, rng, candidate,
                         options.max_diameter, &oversized)) {
      summary.generation_seconds += std::chrono::duration<double>(
          std::chrono::steady_clock::now() - generation_start).count();
      if (oversized) ++summary.oversized;
      continue;
    }
    summary.generation_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - generation_start).count();
    ++summary.generated;
    if (!seen.insert(candidate.key).second) {
      ++summary.duplicates;
      continue;
    }
    ++summary.tested;
    kstab::CertifyResult certificate;
    std::optional<double> probe_value;
    std::optional<double> confirm_value;
    std::optional<double> final_value;
    std::uint64_t probe_evaluations = 0;
    std::uint64_t confirm_evaluations = 0;
    std::uint64_t final_evaluations = 0;
    std::string stage = "probe";
    const bool audit = audit_candidate(candidate.key);

    auto stage_start = std::chrono::steady_clock::now();
    const auto probe = kstab::search_witness(candidate.vertices, candidate.ell,
                                             probe_options);
    const double probe_elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - stage_start).count();
    summary.probe_seconds += probe_elapsed;
    ++summary.probed;
    probe_value = probe.witness.normalized;
    probe_evaluations = static_cast<std::uint64_t>(std::max(0L, probe.evaluations));
    summary.df_evaluations += probe_evaluations;
    stage = "probe_complete";

    if (*probe_value <= -1e-8 || audit) {
      stage_start = std::chrono::steady_clock::now();
      const auto confirm = kstab::search_witness(candidate.vertices, candidate.ell,
                                                  confirm_options);
      summary.confirm_seconds += std::chrono::duration<double>(
          std::chrono::steady_clock::now() - stage_start).count();
      ++summary.confirmed;
      confirm_value = confirm.witness.normalized;
      confirm_evaluations = static_cast<std::uint64_t>(std::max(0L, confirm.evaluations));
      summary.df_evaluations += confirm_evaluations;
      stage = "confirm_complete";
      if (*confirm_value < 0.0) {
        certificate = kstab::certify_witness(candidate.vertices, candidate.ell,
                                              confirm.witness,
                                              options.certify_max_denominator);
      }
      const bool need_final = audit || *confirm_value <= -1e-7 ||
                              (*confirm_value < 0.0 && !certificate.certified);
      if (need_final) {
        stage_start = std::chrono::steady_clock::now();
        const auto final = kstab::search_witness(candidate.vertices, candidate.ell,
                                                  final_options);
        summary.final_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - stage_start).count();
        ++summary.finalized;
        final_value = final.witness.normalized;
        final_evaluations = static_cast<std::uint64_t>(std::max(0L, final.evaluations));
        summary.df_evaluations += final_evaluations;
        stage = "final_complete";
        if (*final_value < 0.0) {
          const auto final_certificate = kstab::certify_witness(
              candidate.vertices, candidate.ell, final.witness,
              options.certify_max_denominator);
          if (final_certificate.certified) certificate = final_certificate;
        }
      }
    }

    if (certificate.certified) ++summary.certified_candidates;
    database.save_candidate(candidate,
                            certificate.certified ? "verified_unstable" : "unverified",
                            certificate.certified ? &certificate : nullptr,
                            profile, stage, probe_value, confirm_value, final_value,
                            probe_evaluations, confirm_evaluations, final_evaluations);
    if (certificate.certified) {
      summary.found = true;
      summary.candidate = std::move(candidate);
      summary.certificate = certificate;
      break;
    }
    database.save_state("d" + std::to_string(options.d) + "|rng",
                        serialize_rng(rng));
  }
  database.save_state("d" + std::to_string(options.d) + "|rng",
                      serialize_rng(rng));
  if (summary.found && !write_result(options, summary.candidate, summary.certificate))
    throw std::runtime_error("failed to write certified Delzant result files");
  summary.total_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - search_start).count();
  return summary;
}

}  // namespace delzant
