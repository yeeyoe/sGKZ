#pragma once

#include "k_stability.hpp"

#include <optional>
#include <string>
#include <utility>

namespace kstab {

// A rational convex polygon whose boundary measure is given edge by edge.
// edge_measures[i] is the total dσ measure on vertices_ccw[i] -> i+1.
struct MeasuredPolygon {
  std::vector<QPoint> vertices_ccw;
  std::vector<Rational> edge_measures;
};

MeasuredPolygon make_measured_polygon(
    const std::vector<IntPoint>& vertices_ccw,
    const std::vector<bool>& null_measure_edges = {});
BoundaryMoments boundary_moments(const MeasuredPolygon& polygon);
std::array<Rational, 3> compute_ell_p(const MeasuredPolygon& polygon);
Rational df_simple_exact(const MeasuredPolygon& polygon,
                         const std::array<Rational, 3>& ell,
                         const Rational& a, const Rational& b,
                         const Rational& c);

struct MeasuredSearchResult {
  bool unstable = false;
  Witness witness;
  double normalization = 1.0;
  long evaluations = 0;
};

MeasuredSearchResult search_witness(const MeasuredPolygon& polygon,
                                    const std::array<Rational, 3>& ell,
                                    const SearchOptions& options);

// Split a parent polygon with a chord from edge_i(s) to edge_j(t), where
// edge_i(s) = p_i + s(p_{i+1} - p_i).  The two endpoints must differ and
// 0 <= s,t <= 1.  Every new chord has zero boundary measure.
struct ChordSplit {
  QPoint first;
  QPoint second;
  MeasuredPolygon first_piece;
  MeasuredPolygon second_piece;
};

ChordSplit split_by_chord(const MeasuredPolygon& parent, std::size_t edge_i,
                          const Rational& s, std::size_t edge_j,
                          const Rational& t);

struct DecomposeOptions {
  SearchOptions witness_search;
  int root_starts = 7;
  int root_iterations = 60;
  double root_tolerance = 1e-10;
  std::int64_t rational_max_denominator = 1048576;
};

struct DecompositionCandidate {
  std::size_t first_edge = 0;
  std::size_t second_edge = 0;
  double first_parameter = 0.0;
  double second_parameter = 0.0;
  double continuity_residual = 0.0;
  bool concave = false;
  bool certified_rational = false;
  QPoint first_endpoint;
  QPoint second_endpoint;
  std::optional<Rational> first_parameter_rational;
  std::optional<Rational> second_parameter_rational;
  std::array<Rational, 3> first_ell{};
  std::array<Rational, 3> second_ell{};
  MeasuredPolygon first_piece;
  MeasuredPolygon second_piece;
  MeasuredSearchResult first_search;
  MeasuredSearchResult second_search;
};

std::vector<DecompositionCandidate> find_chord_decompositions(
    const MeasuredPolygon& parent, const DecomposeOptions& options);

// Contract for the later user-supplied multi-chord tree mode.  It deliberately
// describes only combinatorics and boundary placement; this first release
// solves the one-chord (zero-interior-node) instance only.
enum class DecompositionNodeKind { boundary_leaf, interior };

struct DecompositionTreeNode {
  DecompositionNodeKind kind = DecompositionNodeKind::boundary_leaf;
  // Boundary leaves use exactly one of these fields.  boundary_edge means a
  // free parameter along that parent edge; boundary_vertex is a fixed leaf.
  std::optional<std::size_t> boundary_edge;
  std::optional<std::size_t> boundary_vertex;
  // An interior node must list its three incident tree-edge indices in cyclic
  // order.  Boundary leaves leave this empty.
  std::vector<std::size_t> cyclic_edges;
};

struct DecompositionTreeTopology {
  std::vector<DecompositionTreeNode> nodes;
  std::vector<std::pair<std::size_t, std::size_t>> edges;
};

struct TreeTopologySummary {
  std::size_t interior_nodes = 0;
  std::size_t boundary_leaves = 0;
  std::size_t free_variables = 0;
  std::size_t continuity_equations = 0;
};

// Validates a connected acyclic graph with degree-3 interior nodes and
// degree-1 boundary leaves.  parent_edge_count is the number of original
// polygon edges.  For no fixed leaves, both count fields equal 3I+2.
TreeTopologySummary validate_tree_topology(
    const DecompositionTreeTopology& topology, std::size_t parent_edge_count);

std::string rational_string(const Rational& value);

}  // namespace kstab
