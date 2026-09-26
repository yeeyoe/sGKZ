#pragma once

#include "k_stability_highdim.hpp"

namespace kstab_highdim::detail {

struct ExactTriangulation {
  std::vector<Simplex> cells;
  std::vector<std::vector<Vector>> boundary_facets;
  std::vector<Vector> hull_points;
};

ExactTriangulation triangulate_exact(const std::vector<Vector>& points,
                                     int dimension);

}  // namespace kstab_highdim::detail
