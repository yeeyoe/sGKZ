#include "delaunay_adapter.hpp"

#include <CGAL/Cartesian_d.h>
#include <CGAL/Delaunay_d.h>

#include <algorithm>

namespace kstab_highdim::detail {
namespace {

// CGAL 6.2's exact dynamic-dimensional triangulation interface is the
// deprecated Delaunay_d API. Keeping it here isolates that compatibility
// dependency from the integration and search code.
using ExactKernel = CGAL::Cartesian_d<Rational>;
using ExactDelaunay = CGAL::Delaunay_d<ExactKernel>;

Vector point_vector(const ExactKernel::Point_d& point) {
  Vector result(point.dimension());
  for (int i = 0; i < point.dimension(); ++i) result[i] = point[i];
  return result;
}

ExactKernel::Point_d make_point(const Vector& point) {
  return ExactKernel::Point_d(static_cast<int>(point.size()), point.begin(), point.end());
}

}  // namespace

ExactTriangulation triangulate_exact(const std::vector<Vector>& points,
                                     int dimension) {
  ExactTriangulation result;
  if (points.size() < static_cast<std::size_t>(dimension + 1)) return result;
  if (dimension == 0) {
    Simplex simplex;
    simplex.vertices.push_back(points.front());
    result.cells.push_back(std::move(simplex));
    return result;
  }
  if (dimension == 1) {
    std::vector<Vector> sorted = points;
    std::sort(sorted.begin(), sorted.end(), [](const Vector& lhs, const Vector& rhs) {
      return lhs.front() < rhs.front();
    });
    for (std::size_t i = 1; i < sorted.size(); ++i) {
      if (sorted[i - 1] == sorted[i]) continue;
      result.cells.push_back({{sorted[i - 1], sorted[i]}});
    }
    if (!sorted.empty()) {
      result.boundary_facets.push_back({sorted.front()});
      if (sorted.size() > 1) result.boundary_facets.push_back({sorted.back()});
    }
    return result;
  }
  ExactDelaunay delaunay(dimension);
  for (const auto& point : points) delaunay.insert(make_point(point));
  for (const auto vertex : delaunay.all_vertices(ExactDelaunay::FURTHEST)) {
    result.hull_points.push_back(point_vector(delaunay.associated_point(vertex)));
  }
  for (const auto simplex : delaunay.all_simplices()) {
    Simplex cell;
    cell.vertices.reserve(dimension + 1);
    for (int i = 0; i <= dimension; ++i) {
      cell.vertices.push_back(point_vector(delaunay.point_of_simplex(simplex, i)));
    }
    result.cells.push_back(std::move(cell));
    for (int omitted = 0; omitted <= dimension; ++omitted) {
      if (delaunay.opposite_simplex(simplex, omitted) != ExactDelaunay::Simplex_handle()) continue;
      std::vector<Vector> facet;
      for (int i = 0; i <= dimension; ++i) {
        if (i != omitted) facet.push_back(point_vector(delaunay.point_of_simplex(simplex, i)));
      }
      result.boundary_facets.push_back(std::move(facet));
    }
  }
  return result;
}

}  // namespace kstab_highdim::detail
