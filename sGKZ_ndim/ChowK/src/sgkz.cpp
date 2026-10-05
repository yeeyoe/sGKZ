#include "chowk/sgkz.hpp"

#include <Eigen/QR>
#include <CGAL/Delaunay_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Regular_triangulation.h>
#include <CGAL/NewKernel_d/Cartesian_base.h>
#include <CGAL/NewKernel_d/Wrapper/Cartesian_wrap.h>
#include <CGAL/NewKernel_d/Kernel_d_interface.h>
#include <CGAL/NewKernel_d/Types/Weighted_point.h>
#include <CGAL/QP_functions.h>
#include <CGAL/QP_models.h>
#include <gmpxx.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>

namespace chowk_sgkz {
namespace {
using ExactDelaunayKernel = CGAL::Cartesian_d<Rational>;
using ExactDelaunay = CGAL::Delaunay_d<ExactDelaunayKernel>;
using NumericKernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using NumericRT = CGAL::Regular_triangulation<NumericKernel>;
using ExactBase = CGAL::Cartesian_base_d<Rational, CGAL::Dynamic_dimension_tag>;
using ExactWrap = CGAL::Cartesian_wrap<ExactBase>;
using ExactKernel = CGAL::Kernel_d_interface<ExactWrap>;
using ExactRT = CGAL::Regular_triangulation<ExactKernel>;

void require(bool ok, const std::string& message) { if (!ok) throw std::invalid_argument(message); }

Integer determinant(std::vector<std::vector<Integer>> a) {
  const std::size_t n = a.size();
  if (n == 0) return Integer(1);
  for (const auto& row : a) require(row.size() == n, "matrix is not square");
  Integer previous = 1;
  int sign = 1;
  for (std::size_t col = 0; col < n; ++col) {
    std::size_t pivot = col;
    while (pivot < n && a[pivot][col] == 0) ++pivot;
    if (pivot == n) return Integer(0);
    if (pivot != col) { std::swap(a[pivot], a[col]); sign = -sign; }
    const Integer p = a[col][col];
    if (col + 1 == n) return sign < 0 ? -p : p;
    for (std::size_t row = col + 1; row < n; ++row) {
      const Integer f = a[row][col];
      for (std::size_t j = col + 1; j < n; ++j) a[row][j] = (a[row][j] * p - f * a[col][j]) / previous;
    }
    previous = p;
  }
  return Integer(0);
}

Integer abs_integer(Integer x) { return x < 0 ? -x : x; }

Integer simplex_det(const std::vector<IntVector>& vertices) {
  const int d = static_cast<int>(vertices.front().size());
  std::vector<std::vector<Integer>> m(d, std::vector<Integer>(d));
  for (int i = 0; i < d; ++i) for (int j = 0; j < d; ++j)
    m[i][j] = Integer(std::to_string(vertices[i + 1][j])) - Integer(std::to_string(vertices[0][j]));
  return abs_integer(determinant(std::move(m)));
}

std::string key(const IntVector& p) { std::ostringstream s; for (auto x : p) s << x << ','; return s.str(); }

std::vector<IntVector> read_rows(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open input file: " + path.string());
  std::vector<IntVector> rows; std::string line; int d = -1;
  while (std::getline(in, line)) {
    const auto hash = line.find('#'); if (hash != std::string::npos) line.resize(hash);
    for (char& c : line) if (c == ',') c = ' ';
    std::istringstream ss(line); IntVector row; std::int64_t x;
    while (ss >> x) row.push_back(x);
    if (row.empty()) continue;
    if (d < 0) d = static_cast<int>(row.size());
    if (static_cast<int>(row.size()) != d) throw std::invalid_argument("input rows have different dimensions");
    rows.push_back(std::move(row));
  }
  if (rows.empty()) throw std::invalid_argument("input has no points");
  return rows;
}

struct TriangulationData { std::vector<IntVector> hull; std::vector<std::vector<std::size_t>> cells; Integer volume = 0; };

TriangulationData triangulate(const std::vector<IntVector>& points) {
  const int d = static_cast<int>(points.front().size());
  require(d >= 2, "dimension must be at least 2");
  ExactDelaunay dt(d);
  for (const auto& p : points) {
    std::vector<Rational> q; q.reserve(p.size()); for (auto x : p) q.emplace_back(x);
    dt.insert(ExactDelaunayKernel::Point_d(d, q.begin(), q.end()));
  }
  require(dt.current_dimension() == d, "points are not full dimensional");
  TriangulationData result;
  std::unordered_map<std::string, std::size_t> index;
  for (std::size_t i = 0; i < points.size(); ++i) index[key(points[i])] = i;
  std::vector<std::size_t> hull_indices;
  for (auto vh : dt.all_vertices(ExactDelaunay::FURTHEST)) {
    IntVector p(d); auto q = dt.associated_point(vh); for (int j = 0; j < d; ++j) {
      const Rational c = q[j]; require(c.denominator() == 1, "nonintegral hull coordinate");
      p[j] = ([](const Integer& z){ std::ostringstream ss; ss << z; return std::stoll(ss.str()); })(c.numerator());
    }
    result.hull.push_back(p); hull_indices.push_back(index.at(key(p)));
  }
  for (auto sh : dt.all_simplices()) {
    std::vector<std::size_t> cell; cell.reserve(d + 1); std::vector<IntVector> verts;
    for (int i = 0; i <= d; ++i) {
      auto q = dt.point_of_simplex(sh, i); IntVector p(d);
      for (int j = 0; j < d; ++j) { const Rational c = q[j]; require(c.denominator() == 1, "nonintegral simplex coordinate"); p[j] = ([](const Integer& z){ std::ostringstream ss; ss << z; return std::stoll(ss.str()); })(c.numerator()); }
      cell.push_back(index.at(key(p))); verts.push_back(std::move(p));
    }
    result.volume += simplex_det(verts); result.cells.push_back(std::move(cell));
  }
  require(!result.cells.empty() && result.volume > 0, "cannot triangulate full-dimensional hull");
  return result;
}

bool solve_linear(std::vector<std::vector<Rational>> a, std::vector<Rational> b, std::vector<Rational>& x) {
  const std::size_t n = a.size(); if (b.size() != n) return false;
  for (const auto& r : a) if (r.size() != n) return false;
  for (std::size_t c = 0; c < n; ++c) {
    std::size_t p = c; while (p < n && a[p][c] == 0) ++p; if (p == n) return false;
    if (p != c) { std::swap(a[p], a[c]); std::swap(b[p], b[c]); }
    Rational q = a[c][c]; for (std::size_t j = c; j < n; ++j) a[c][j] /= q; b[c] /= q;
    for (std::size_t r = 0; r < n; ++r) if (r != c && a[r][c] != 0) {
      Rational f = a[r][c]; for (std::size_t j = c; j < n; ++j) a[r][j] -= f * a[c][j]; b[r] -= f * b[c];
    }
  }
  x = std::move(b); return true;
}

bool inside_simplex(const IntVector& p, const std::vector<IntVector>& v) {
  const int d = static_cast<int>(p.size());
  std::vector<std::vector<Rational>> m(d + 1, std::vector<Rational>(d + 1));
  std::vector<Rational> rhs(d + 1); for (int r = 0; r < d; ++r) { rhs[r] = p[r]; for (int c = 0; c <= d; ++c) m[r][c] = v[c][r]; }
  rhs[d] = 1; for (int c = 0; c <= d; ++c) m[d][c] = 1;
  std::vector<Rational> lambda; if (!solve_linear(m, rhs, lambda)) return false;
  for (const auto& x : lambda) if (x < 0) return false;
  return true;
}

std::vector<IntVector> enumerate_lattice(const std::vector<IntVector>& vertices, std::int64_t k, std::size_t max_points) {
  const int d = static_cast<int>(vertices.front().size());
  auto base = triangulate(vertices); std::vector<std::vector<IntVector>> cells;
  for (const auto& cell : base.cells) { std::vector<IntVector> c; for (auto i : cell) { IntVector p = vertices[i]; for (auto& x : p) x *= k; c.push_back(std::move(p)); } cells.push_back(std::move(c)); }
  IntVector lo(d, std::numeric_limits<std::int64_t>::max()), hi(d, std::numeric_limits<std::int64_t>::min());
  for (const auto& p : vertices) for (int j = 0; j < d; ++j) { lo[j] = std::min(lo[j], p[j] * k); hi[j] = std::max(hi[j], p[j] * k); }
  std::vector<IntVector> result; IntVector current(d); std::function<void(int)> rec = [&](int axis) {
    if (axis == d) { for (const auto& cell : cells) if (inside_simplex(current, cell)) { result.push_back(current); if (max_points && result.size() > max_points) throw std::runtime_error("lattice point count exceeds --max-points"); return; } return; }
    for (std::int64_t x = lo[axis]; x <= hi[axis]; ++x) { current[axis] = x; rec(axis + 1); if (x == std::numeric_limits<std::int64_t>::max()) break; }
  }; rec(0); return result;
}

GkzVector build_vector(const PointConfiguration& c, const std::vector<std::vector<std::size_t>>& cells, std::size_t visible, std::size_t hidden, bool keep) {
  GkzVector out; out.numerators.assign(c.size(), Integer(0)); out.cells = cells.size(); out.visible_vertices = visible; out.hidden_vertices = hidden;
  if (keep) out.triangulation = std::make_shared<const std::vector<std::vector<std::size_t>>>(cells);
  for (const auto& cell : cells) { std::vector<IntVector> v; for (auto i : cell) v.push_back(c.points()[i]); const Integer det = simplex_det(v); for (auto i : cell) out.numerators[i] += det; }
  const Integer denom = Integer(c.dimension() + 1) * c.determinant_volume();
  out.values = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(c.size()));
  for (std::size_t i = 0; i < c.size(); ++i) out.values[static_cast<Eigen::Index>(i)] = CGAL::to_double(Rational(out.numerators[i], denom));
  return out;
}

GkzVector oracle_2d(const PointConfiguration& c,
                    const std::vector<Rational>& heights, bool keep) {
  const auto& points = c.points();
  require(heights.size() == points.size(), "height size mismatch");
  std::set<std::vector<std::size_t>> lower_faces;
  for (std::size_t i = 0; i < points.size(); ++i) {
    for (std::size_t j = i + 1; j < points.size(); ++j) {
      for (std::size_t k = j + 1; k < points.size(); ++k) {
        std::vector<std::vector<Rational>> a = {
            {Rational(points[i][0]), Rational(points[i][1]), Rational(1)},
            {Rational(points[j][0]), Rational(points[j][1]), Rational(1)},
            {Rational(points[k][0]), Rational(points[k][1]), Rational(1)}};
        std::vector<Rational> rhs = {heights[i], heights[j], heights[k]};
        std::vector<Rational> plane;
        if (!solve_linear(std::move(a), std::move(rhs), plane)) continue;
        bool supporting = true;
        std::vector<std::size_t> face;
        for (std::size_t q = 0; q < points.size(); ++q) {
          const Rational value = plane[0] * Rational(points[q][0]) +
                                 plane[1] * Rational(points[q][1]) + plane[2];
          if (value > heights[q]) {
            supporting = false;
            break;
          }
          if (value == heights[q]) face.push_back(q);
        }
        if (supporting) lower_faces.insert(std::move(face));
      }
    }
  }

  std::vector<std::vector<std::size_t>> cells;
  for (const auto& face : lower_faces) {
    if (face.size() < 3) continue;
    std::vector<IntVector> face_points;
    face_points.reserve(face.size());
    for (auto index : face) face_points.push_back(points[index]);
    bool full_dimensional = false;
    for (std::size_t i = 0; i < face_points.size() && !full_dimensional; ++i)
      for (std::size_t j = i + 1; j < face_points.size() && !full_dimensional; ++j)
        for (std::size_t k = j + 1; k < face_points.size(); ++k)
          if (simplex_det({face_points[i], face_points[j], face_points[k]}) != 0) {
            full_dimensional = true;
            break;
          }
    if (!full_dimensional) continue;
    const auto triangulation = triangulate(face_points);
    for (const auto& local_cell : triangulation.cells) {
      std::vector<std::size_t> cell;
      cell.reserve(3);
      for (auto local_index : local_cell) cell.push_back(face[local_index]);
      cells.push_back(std::move(cell));
    }
  }

  Integer covered = 0;
  for (const auto& cell : cells) {
    std::vector<IntVector> vertices;
    for (auto index : cell) vertices.push_back(points[index]);
    covered += simplex_det(vertices);
  }
  if (covered != c.determinant_volume()) {
    std::ostringstream message;
    message << "2D lower subdivision does not cover the convex hull (volume "
            << covered << ", expected " << c.determinant_volume()
            << ", lower faces " << lower_faces.size() << ", triangles "
            << cells.size() << ')';
    throw std::runtime_error(message.str());
  }
  std::set<std::size_t> used;
  for (const auto& cell : cells) used.insert(cell.begin(), cell.end());
  return build_vector(c, cells, used.size(), c.size() - used.size(), keep);
}

template<class Kernel, class RT>
GkzVector oracle_impl(const PointConfiguration& c, const std::vector<Rational>& heights, bool keep) {
  if (c.dimension() == 2) return oracle_2d(c, heights, keep);
  const int d = c.dimension(); RT rt(d); std::vector<typename RT::Weighted_point> weighted;
  weighted.reserve(c.size());
  for (std::size_t i = 0; i < c.size(); ++i) {
    std::vector<typename Kernel::FT> coords; for (auto x : c.points()[i]) coords.emplace_back(x);
    typename Kernel::Point_d p(coords.begin(), coords.end());
    typename Kernel::FT squared_norm = 0;
    for (auto x : c.points()[i]) {
      typename Kernel::FT coordinate;
      if constexpr (std::is_same_v<typename Kernel::FT, Rational>)
        coordinate = Rational(std::to_string(x));
      else
        coordinate = static_cast<double>(x);
      squared_norm += coordinate * coordinate;
    }
    typename Kernel::FT h;
    if constexpr (std::is_same_v<typename Kernel::FT, Rational>)
      h = heights[i];
    else
      h = CGAL::to_double(heights[i]);
    weighted.emplace_back(p, squared_norm - h);
  }
  rt.insert(weighted.begin(), weighted.end());
  require(rt.current_dimension() == d, "regular triangulation has wrong dimension");
  std::unordered_map<std::string, std::size_t> lookup; for (std::size_t i = 0; i < c.size(); ++i) lookup[key(c.points()[i])] = i;
  std::vector<std::vector<std::size_t>> cells;
  for (auto it = rt.finite_full_cells_begin(); it != rt.finite_full_cells_end(); ++it) {
    std::vector<std::size_t> cell; for (int j = 0; j <= d; ++j) { auto p = it->vertex(j)->point().point(); IntVector q(d); for (int k = 0; k < d; ++k) { auto z = p[k]; if constexpr (std::is_same_v<typename Kernel::FT, Rational>) { require(z.denominator() == 1, "regular cell coordinate is not integral"); std::ostringstream ss; ss << z.numerator(); q[k] = std::stoll(ss.str()); } else q[k] = static_cast<std::int64_t>(std::llround(static_cast<double>(z))); } cell.push_back(lookup.at(key(q))); } cells.push_back(std::move(cell));
  }
  Integer volume = 0; for (const auto& cell : cells) { std::vector<IntVector> v; for (auto i : cell) v.push_back(c.points()[i]); volume += simplex_det(v); }
  if (volume != c.determinant_volume()) {
    std::ostringstream msg;
    msg << "regular triangulation does not cover the convex hull (volume " << volume
        << ", expected " << c.determinant_volume() << ", cells " << cells.size() << "):";
    for (const auto& cell : cells) {
      msg << " [";
      for (auto index : cell) msg << index << ',';
      msg << "]";
    }
    throw std::runtime_error(msg.str());
  }
  return build_vector(c, cells, rt.number_of_vertices(), rt.number_of_hidden_vertices(), keep);
}

Eigen::MatrixXd gram(const std::vector<GkzVector>& a) { Eigen::Index n = static_cast<Eigen::Index>(a.size()); Eigen::MatrixXd g(n,n); for (Eigen::Index i=0;i<n;++i) for (Eigen::Index j=0;j<=i;++j) g(i,j)=g(j,i)=a[i].values.dot(a[j].values); return g; }
Eigen::VectorXd combine(const std::vector<GkzVector>& a, const Eigen::VectorXd& x) { Eigen::VectorXd r=Eigen::VectorXd::Zero(a.front().values.size()); for (std::size_t i=0;i<a.size();++i) r += x[static_cast<Eigen::Index>(i)]*a[i].values; return r; }

Eigen::VectorXd active_qp(const Eigen::MatrixXd& G, const Eigen::VectorXd& start, const SolverOptions& o) {
  const Eigen::Index m=G.rows(); Eigen::VectorXd x=start; if (x.size()!=m) x=Eigen::VectorXd::Constant(m,1.0/m); x=x.cwiseMax(0); x/=x.sum();
  std::vector<Eigen::Index> work; for (Eigen::Index i=0;i<m;++i) if (x[i]>o.prune_tolerance) work.push_back(i); if (work.empty()) work.push_back(0);
  for (int step=0;step<o.max_correction_steps;++step) {
    Eigen::Index r=static_cast<Eigen::Index>(work.size()); Eigen::MatrixXd K=Eigen::MatrixXd::Zero(r+1,r+1); Eigen::VectorXd b=Eigen::VectorXd::Zero(r+1); b[r]=1;
    for(Eigen::Index i=0;i<r;++i){for(Eigen::Index j=0;j<r;++j)K(i,j)=G(work[i],work[j]);K(i,r)=-1;K(r,i)=1;}
    Eigen::VectorXd sol=K.completeOrthogonalDecomposition().solve(b); Eigen::VectorXd y=sol.head(r);
    if (y.minCoeff() < -o.correction_tolerance) { double theta=1; for(Eigen::Index i=0;i<r;++i) if(y[i]<0) theta=std::min(theta,x[work[i]]/(x[work[i]]-y[i])); for(Eigen::Index i=0;i<r;++i)x[work[i]]=(1-theta)*x[work[i]]+theta*y[i]; work.erase(std::remove_if(work.begin(),work.end(),[&](Eigen::Index i){return x[i]<=o.prune_tolerance;}),work.end()); if(work.empty())work.push_back(0); continue; }
    x.setZero(); for(Eigen::Index i=0;i<r;++i)x[work[i]]=std::max(0.0,y[i]); x/=x.sum(); Eigen::VectorXd grad=G*x; double mult=x.dot(grad); Eigen::Index enter=-1; double best=std::numeric_limits<double>::infinity(); for(Eigen::Index i=0;i<m;++i) if(std::find(work.begin(),work.end(),i)==work.end()&&grad[i]<best){best=grad[i];enter=i;} if(enter<0||best>=mult-o.correction_tolerance*std::max(1.0,std::abs(mult))) return x; work.push_back(enter);
  }
  throw std::runtime_error("active-set QP exceeded correction limit");
}

ExactCertificate certify(const PointConfiguration& c, const std::vector<GkzVector>& active) {
  ExactCertificate out; if(active.empty()){out.message="empty active set";return out;} const std::size_t m=active.size(), n=c.size(); using Program=CGAL::Quadratic_program<Integer>; Program p(CGAL::EQUAL,true,Integer(0),false,Integer(0)); p.set_b(0,Integer(1));
  for(std::size_t i=0;i<m;++i){p.set_a(static_cast<int>(i),0,Integer(1)); for(std::size_t j=0;j<=i;++j){Integer s=0;for(std::size_t k=0;k<n;++k)s+=active[i].numerators[k]*active[j].numerators[k];p.set_d(static_cast<int>(i),static_cast<int>(j),2*s);}}
  auto sol=CGAL::solve_quadratic_program(p,Integer()); if(!sol.is_optimal()){out.message="exact active QP failed";return out;} auto it=sol.variable_values_begin(); out.coefficients.reserve(m); Rational sum=0; for(std::size_t i=0;i<m;++i,++it){out.coefficients.emplace_back(it->numerator(),it->denominator()); if(out.coefficients.back()<0){out.message="negative exact coefficient";return out;}sum+=out.coefficients.back();} if(sum!=1){out.message="exact coefficients do not sum to one";return out;}
  const Rational den(Integer(c.dimension()+1)*c.determinant_volume()); out.sigma.assign(n,Rational(0)); for(std::size_t i=0;i<m;++i)for(std::size_t k=0;k<n;++k)out.sigma[k]+=out.coefficients[i]*Rational(active[i].numerators[k]) / Rational(den); out.norm_squared=0; for(auto x:out.sigma)out.norm_squared+=x*x;
  std::vector<Rational> heights=out.sigma; RegularTriangulationOracle oracle(c); auto w=oracle.minimize_exact(heights,false); Rational dot=0; for(std::size_t k=0;k<n;++k)dot+=out.sigma[k]*Rational(w.numerators[k]) / Rational(den); if(out.norm_squared==dot){out.certified=true;out.message="exact global certificate";} else {out.has_witness=true;out.witness=std::move(w);out.message="exact oracle found a better GKZ vector";} return out;
}

} // namespace

std::string rational_string(const Rational& v){std::ostringstream s;if(v.denominator()==1)s<<v.numerator();else s<<v;return s.str();}
std::string integer_string(const Integer& v){std::ostringstream s; s<<v; return s.str();}

PointConfiguration::PointConfiguration(std::vector<IntVector> points,std::vector<IntVector> hull,std::vector<std::vector<std::size_t>> cells,Integer volume,std::int64_t level,Integer base):dimension_(static_cast<int>(points.front().size())),points_(std::move(points)),hull_vertices_(std::move(hull)),hull_simplices_(std::move(cells)),determinant_volume_(std::move(volume)),level_(level),base_determinant_volume_(std::move(base)){}
PointConfiguration PointConfiguration::from_points(std::vector<IntVector> points){std::set<IntVector> u(points.begin(),points.end());if(u.size()!=points.size())throw std::invalid_argument("duplicate point");auto t=triangulate(points);return PointConfiguration(std::move(points),std::move(t.hull),std::move(t.cells),std::move(t.volume),0,Integer(0));}
PointConfiguration PointConfiguration::from_points_file(const std::filesystem::path& p){return from_points(read_rows(p));}
PointConfiguration PointConfiguration::from_lattice_polytope(std::vector<IntVector> vertices,std::int64_t k,std::size_t max_points){require(k>0,"k must be positive");std::set<IntVector> u(vertices.begin(),vertices.end());if(u.size()!=vertices.size())throw std::invalid_argument("duplicate polytope vertex");auto base=triangulate(vertices);auto points=enumerate_lattice(vertices,k,max_points);if(points.empty())throw std::runtime_error("polytope contains no lattice points");auto t=triangulate(points);Integer scaled=base.volume;Integer kk(std::to_string(k));for(int i=0;i<static_cast<int>(vertices.front().size());++i)scaled*=kk;return PointConfiguration(std::move(points),std::move(t.hull),std::move(t.cells),std::move(t.volume),k,std::move(base.volume));}
PointConfiguration PointConfiguration::from_polytope_file(const std::filesystem::path& p,std::int64_t k,std::size_t max_points){return from_lattice_polytope(read_rows(p),k,max_points);}
Rational PointConfiguration::volume()const{ Rational v(determinant_volume_); for(int i=2;i<=dimension_;++i) v/=i; return v; }
Rational PointConfiguration::base_volume()const{ Rational v(base_determinant_volume_); for(int i=2;i<=dimension_;++i) v/=i; return v; }
std::vector<Rational> PointConfiguration::centroid()const{std::vector<Rational> c(dimension_,0);Integer total=0;for(const auto& cell:hull_simplices_){std::vector<IntVector> v;for(auto i:cell)v.push_back(points_[i]);Integer d=simplex_det(v);total+=d;for(int j=0;j<dimension_;++j){Integer sum=0;for(const auto& p:v)sum+=Integer(std::to_string(p[j]));c[j]+=Rational(d*sum,Integer(dimension_+1));}}for(auto& x:c)x/=Rational(total);return c;}

AffineFunction compute_ell(const PointConfiguration& c){const int d=c.dimension();const Rational scale=c.is_polytope_level()?Rational(c.level()):Rational(1);std::vector<std::vector<Rational>> m(d+1,std::vector<Rational>(d+1));for(const auto& p:c.points()){std::vector<Rational>x(d);for(int j=0;j<d;++j)x[j]=Rational(p[j])/scale;m[0][0]+=1;for(int j=0;j<d;++j){m[0][j+1]+=x[j];m[j+1][0]+=x[j];for(int k=0;k<d;++k)m[j+1][k+1]+=x[j]*x[k];}}auto cent=c.centroid();for(auto&x:cent)x/=scale;std::vector<Rational> rhs(d+1);rhs[0]=1;for(int j=0;j<d;++j)rhs[j+1]=cent[j];std::vector<Rational> sol;if(!solve_linear(m,rhs,sol))throw std::runtime_error("ell_A moment system is singular");AffineFunction out;out.exact_coefficients=sol;out.coefficients.reserve(sol.size());for(auto&x:sol)out.coefficients.push_back(CGAL::to_double(x));out.values=Eigen::VectorXd::Zero(c.size());for(std::size_t i=0;i<c.size();++i){Rational v=sol[0];for(int j=0;j<d;++j)v+=sol[j+1]*Rational(c.points()[i][j])/scale;out.exact_values.push_back(v);out.values[static_cast<Eigen::Index>(i)]=CGAL::to_double(v);}return out;}

RegularTriangulationOracle::RegularTriangulationOracle(const PointConfiguration& c):configuration_(c){}
GkzVector RegularTriangulationOracle::minimize(const Eigen::VectorXd& h,bool keep)const{if(h.size()!=static_cast<Eigen::Index>(configuration_.size()))throw std::invalid_argument("height size mismatch");std::vector<Rational> q(h.size());for(Eigen::Index i=0;i<h.size();++i)q[static_cast<std::size_t>(i)]=Rational(h[i]);return oracle_impl<NumericKernel,NumericRT>(configuration_,q,keep);}
GkzVector RegularTriangulationOracle::minimize_exact(const std::vector<Rational>& h,bool keep)const{if(h.size()!=configuration_.size())throw std::invalid_argument("height size mismatch");return oracle_impl<ExactKernel,ExactRT>(configuration_,h,keep);}
GkzVector RegularTriangulationOracle::minimize_exact_integer(const std::vector<Integer>& h,bool keep)const{std::vector<Rational> q;for(const auto&x:h)q.emplace_back(x);return minimize_exact(q,keep);}

SolverResult ShortestGkzSolver::solve(const PointConfiguration& c)const{RegularTriangulationOracle oracle(c);SolverResult r;std::vector<GkzVector> active;Eigen::VectorXd coeff;Eigen::VectorXd h=Eigen::VectorXd::Zero(c.size());active.push_back(oracle.minimize(h,false));coeff=Eigen::VectorXd::Ones(1);for(int iter=0;iter<=options_.max_iterations;++iter){Eigen::MatrixXd G=gram(active);coeff=active_qp(G,coeff,options_);Eigen::VectorXd x=combine(active,coeff);auto v=oracle.minimize(x,false);double norm=x.squaredNorm();double gap=norm-x.dot(v.values);if(gap<0&&gap>-1e-10)gap=0;r.iterations=iter;r.sigma=x;r.norm_squared=norm;r.gap=gap;r.l2_error_bound=std::sqrt(std::max(0.0,2*gap));if(gap<=std::max(options_.absolute_tolerance,options_.tolerance*std::max(1.0,norm))){r.converged=true;break;}bool seen=false;for(const auto&a:active)if(a.same_numerators(v)){seen=true;break;}if(seen){r.converged=false;break;}active.push_back(std::move(v));if(options_.verbose)std::cerr<<"iteration="<<iter<<" active="<<active.size()<<" gap="<<std::setprecision(17)<<gap<<'\n';}r.active_vectors=active;r.coefficients=coeff;if(options_.exact_certification&&(options_.exact_max_active==0||static_cast<int>(active.size())<=options_.exact_max_active)){for(int pass=0;pass<options_.max_iterations&&!r.exact.certified;++pass){r.exact=certify(c,active);if(r.exact.certified){r.sigma.resize(c.size());for(std::size_t i=0;i<c.size();++i)r.sigma[static_cast<Eigen::Index>(i)]=CGAL::to_double(r.exact.sigma[i]);r.norm_squared=CGAL::to_double(r.exact.norm_squared);r.gap=0;r.l2_error_bound=0;break;}if(!r.exact.has_witness)break;bool seen=false;for(const auto&a:active)if(a.same_numerators(r.exact.witness)){seen=true;break;}if(seen)break;active.push_back(std::move(r.exact.witness));r.active_vectors=active;Eigen::MatrixXd gg=gram(active);coeff=active_qp(gg,coeff,options_);r.coefficients=coeff;r.sigma=combine(active,coeff);r.norm_squared=r.sigma.squaredNorm();r.gap=0;}}else r.exact.message="exact certification skipped by active-set limit";return r;}

void write_result_csv(const std::filesystem::path& path,const PointConfiguration& c,const SolverResult& r,bool include_ell){std::ofstream out(path);if(!out)throw std::runtime_error("cannot open output file: "+path.string());AffineFunction ell;if(include_ell)ell=compute_ell(c);for(int j=0;j<c.dimension();++j){if(j)out<<',';out<<"x"<<j+1;}out<<",sigma,sigma_exact";if(include_ell)out<<",ell_A,ell_A_exact";out<<'\n'<<std::setprecision(17);for(std::size_t i=0;i<c.size();++i){for(int j=0;j<c.dimension();++j){if(j)out<<',';out<<c.points()[i][j];}out<<','<<r.sigma[static_cast<Eigen::Index>(i)]<<',';if(r.exact.certified)out<<rational_string(r.exact.sigma[i]);if(include_ell)out<<','<<ell.values[static_cast<Eigen::Index>(i)]<<','<<rational_string(ell.exact_values[i]);out<<'\n';}}

} // namespace sgkz
