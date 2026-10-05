"""Finite-dimensional convex P1 minimization on polygon triangulations.

The module implements the numerical experiment described in
``affiness_research_ideas.md``.  It intentionally keeps geometry and the
quadratic program separate so that meshes can be tested independently.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import scipy.sparse as sp

try:
    import osqp
except ImportError as exc:  # pragma: no cover - exercised in CLI environments
    osqp = None
    _OSQP_IMPORT_ERROR = exc


Array = np.ndarray


def _cross2(a: Array, b: Array) -> float:
    """Scalar two-dimensional cross product (NumPy 2.x removed np.cross for 2D)."""
    return float(a[0] * b[1] - a[1] * b[0])


@dataclass(frozen=True)
class PolygonSample:
    name: str
    vertices: Array
    metadata: dict[str, str]


@dataclass
class Mesh:
    nodes: Array
    triangles: Array
    boundary_edges: Array
    boundary_parent: Array
    level: int
    polygon: PolygonSample


@dataclass
class LevelResult:
    values_u: Array
    psi: float
    theta_values: Array
    triangle_gradients: Array
    status: str
    iterations: int
    primal_residual: float
    dual_residual: float
    solve_time: float
    max_constraint_violation: float
    affine_fit_residual: float
    gradient_jump_max: float
    gradient_jump_mean: float
    gradient_jump_quantiles: tuple[float, float, float]
    gradient_jump_count: int
    gradient_jump_length: float
    l2_difference: float | None


def load_polygon(path: str | Path) -> PolygonSample:
    """Load an integer-vertex polygon and its optional comment metadata."""
    path = Path(path)
    metadata: dict[str, str] = {}
    points: list[tuple[float, float]] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            text = line[1:].strip()
            if "=" in text:
                key, value = text.split("=", 1)
                metadata[key.strip()] = value.strip()
            continue
        x, y = line.split()[:2]
        points.append((float(x), float(y)))
    if len(points) < 3:
        raise ValueError(f"{path}: expected at least three vertices")
    vertices = np.asarray(points, dtype=float)
    area2 = float(np.sum(vertices[:, 0] * np.roll(vertices[:, 1], -1)
                         - vertices[:, 1] * np.roll(vertices[:, 0], -1)))
    if area2 <= 0:
        raise ValueError(f"{path}: vertices must be counter-clockwise with positive area")
    edge_lengths = np.linalg.norm(np.roll(vertices, -1, axis=0) - vertices, axis=1)
    if np.any(edge_lengths <= 0):
        raise ValueError(f"{path}: polygon has a degenerate edge")
    edge_cross = []
    for i in range(len(vertices)):
        a, b, c = vertices[i], vertices[(i + 1) % len(vertices)], vertices[(i + 2) % len(vertices)]
        edge_cross.append(_cross2(b - a, c - b))
    if min(edge_cross) <= 0:
        raise ValueError(f"{path}: polygon is not strictly convex")
    return PolygonSample(path.stem, vertices, metadata)


def _key(point: Sequence[float]) -> tuple[int, int]:
    """Use a dyadic integer key for exact midpoint deduplication."""
    return (int(round(float(point[0]) * 2**40)), int(round(float(point[1]) * 2**40)))


def _make_mesh(nodes: list[Array], triangles: list[tuple[int, int, int]],
               polygon: PolygonSample, level: int) -> Mesh:
    node_array = np.asarray(nodes, dtype=float)
    tri_array = np.asarray(triangles, dtype=np.int64)
    edge_to_tri: dict[tuple[int, int], list[int]] = {}
    for ti, tri in enumerate(tri_array):
        for a, b in ((int(tri[0]), int(tri[1])), (int(tri[1]), int(tri[2])), (int(tri[2]), int(tri[0]))):
            edge_to_tri.setdefault(tuple(sorted((a, b))), []).append(ti)
    boundary: list[tuple[int, int]] = []
    parent: list[int] = []
    for edge, owners in edge_to_tri.items():
        if len(owners) != 1:
            continue
        a, b = edge
        mid = (node_array[a] + node_array[b]) / 2
        parent_edge = -1
        for i, (p, q) in enumerate(zip(polygon.vertices, np.roll(polygon.vertices, -1, axis=0))):
            if abs(_cross2(q - p, node_array[a] - p)) < 1e-8 and abs(_cross2(q - p, node_array[b] - p)) < 1e-8:
                lo = min(np.dot(node_array[a] - p, q - p), np.dot(node_array[b] - p, q - p))
                hi = max(np.dot(node_array[a] - p, q - p), np.dot(node_array[b] - p, q - p))
                norm2 = float(np.dot(q - p, q - p))
                if lo >= -1e-8 and hi <= norm2 + 1e-8:
                    parent_edge = i
                    break
        if parent_edge < 0:
            raise ValueError("could not identify parent polygon edge for boundary edge")
        boundary.append((a, b))
        parent.append(parent_edge)
    return Mesh(node_array, tri_array, np.asarray(boundary, dtype=np.int64),
                np.asarray(parent, dtype=np.int64), level, polygon)


def build_mesh(sample: PolygonSample, level: int) -> Mesh:
    """Build a fan triangulation followed by ``level`` red refinements."""
    if level < 0:
        raise ValueError("level must be non-negative")
    vertices = [p.copy() for p in sample.vertices]
    triangles = [(0, i, i + 1) for i in range(1, len(vertices) - 1)]
    mesh = _make_mesh(vertices, triangles, sample, 0)
    for lev in range(1, level + 1):
        nodes = [p.copy() for p in mesh.nodes]
        node_map = {_key(p): i for i, p in enumerate(nodes)}

        def midpoint(a: int, b: int) -> int:
            point = (mesh.nodes[a] + mesh.nodes[b]) / 2
            k = _key(point)
            if k not in node_map:
                node_map[k] = len(nodes)
                nodes.append(point)
            return node_map[k]

        new_triangles: list[tuple[int, int, int]] = []
        for a, b, c in mesh.triangles:
            a, b, c = int(a), int(b), int(c)
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            new_triangles.extend([(a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)])
        mesh = _make_mesh(nodes, new_triangles, sample, lev)
    return mesh


def _triangle_area_gradients(mesh: Mesh, values: Array) -> tuple[Array, Array]:
    gradients = np.empty((len(mesh.triangles), 2), dtype=float)
    areas = np.empty(len(mesh.triangles), dtype=float)
    for i, tri in enumerate(mesh.triangles):
        x = mesh.nodes[tri]
        mat = np.column_stack((x, np.ones(3)))
        coeff = np.linalg.solve(mat, values[tri])
        gradients[i] = coeff[:2]
        areas[i] = abs(np.linalg.det(np.column_stack((x[1] - x[0], x[2] - x[0])))) / 2
    return areas, gradients


def _triangle_areas(mesh: Mesh) -> Array:
    areas = np.empty(len(mesh.triangles), dtype=float)
    for i, tri in enumerate(mesh.triangles):
        x = mesh.nodes[tri]
        areas[i] = abs(_cross2(x[1] - x[0], x[2] - x[0])) / 2
    return areas


def _mass_norm(mesh: Mesh, values: Array) -> float:
    """Return the L2(P) norm of a P1 nodal vector without assembling M."""
    local = np.array([[2., 1., 1.], [1., 2., 1.], [1., 1., 2.]])
    total = 0.0
    for tri, area in zip(mesh.triangles, _triangle_areas(mesh)):
        v = values[tri]
        total += float(area / 12 * (v @ local @ v))
    return math.sqrt(max(0.0, total))


def assemble_qp(mesh: Mesh):
    """Assemble ``P, q, A, lower, upper`` for OSQP."""
    n = len(mesh.nodes)
    rows, cols, data = [], [], []
    q = np.zeros(n, dtype=float)
    areas = _triangle_areas(mesh)
    local = np.array([[2., 1., 1.], [1., 2., 1.], [1., 1., 2.]])
    for tri, area in zip(mesh.triangles, areas):
        block = area / 12 * local
        for i in range(3):
            for j in range(3):
                rows.append(int(tri[i])); cols.append(int(tri[j])); data.append(block[i, j])
    P = sp.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsc()
    for (a, b), parent in zip(mesh.boundary_edges, mesh.boundary_parent):
        p, r = mesh.polygon.vertices[parent], mesh.polygon.vertices[(parent + 1) % len(mesh.polygon.vertices)]
        edge = mesh.nodes[b] - mesh.nodes[a]
        density = math.gcd(abs(round(r[0] - p[0])), abs(round(r[1] - p[1]))) / np.linalg.norm(r - p)
        q[a] += density * np.linalg.norm(edge) / 2
        q[b] += density * np.linalg.norm(edge) / 2

    edge_to_tri: dict[tuple[int, int], list[int]] = {}
    for ti, tri in enumerate(mesh.triangles):
        for a, b in ((int(tri[0]), int(tri[1])), (int(tri[1]), int(tri[2])), (int(tri[2]), int(tri[0]))):
            edge_to_tri.setdefault(tuple(sorted((a, b))), []).append(ti)
    crow, ccol, cdata = [], [], []
    row_count = 0
    for edge, owners in edge_to_tri.items():
        if len(owners) != 2:
            continue
        a, b = edge
        t0, t1 = owners
        tri0, tri1 = mesh.triangles[t0], mesh.triangles[t1]
        c = next(int(x) for x in tri0 if int(x) not in edge)
        d = next(int(x) for x in tri1 if int(x) not in edge)
        # affine extension from (a,b,c) evaluated at d, represented linearly in nodal values
        xy = mesh.nodes[[a, b, c]]
        bary = np.linalg.solve(np.vstack((xy.T, np.ones(3))), np.r_[mesh.nodes[d], 1.])
        # Do not materialize a length-n row for every inner edge: red
        # refinement produces O(n) such edges, each with only four entries.
        row = row_count
        for j, value in zip((a, b, c, d), (-bary[0], -bary[1], -bary[2], 1.)):
            if value:
                crow.append(row); ccol.append(j); cdata.append(value)
        row_count += 1
    A = sp.coo_matrix((cdata, (crow, ccol)), shape=(row_count, n)).tocsc()
    return P, q, A, np.zeros(row_count), np.full(row_count, np.inf)


def _affine_residual(mesh: Mesh, values: Array) -> float:
    X = np.column_stack((mesh.nodes, np.ones(len(mesh.nodes))))
    coef, *_ = np.linalg.lstsq(X, values, rcond=None)
    return float(np.sqrt(np.mean((X @ coef - values) ** 2)))


def solve_level(mesh: Mesh, previous_solution: Array | None = None) -> LevelResult:
    if osqp is None:
        raise RuntimeError("OSQP is required; install dependencies with `python3 -m pip install -r requirements.txt`") from _OSQP_IMPORT_ERROR
    P, q, A, lower, upper = assemble_qp(mesh)
    n = len(mesh.nodes)
    solver = osqp.OSQP()
    solver.setup(P=P, q=q, A=A, l=lower, u=upper, eps_abs=1e-8, eps_rel=1e-8,
                 max_iter=200000, adaptive_rho=True, polishing=True, verbose=False)
    if previous_solution is not None and len(previous_solution) == n:
        solver.warm_start(x=previous_solution)
    started = time.perf_counter()
    result = solver.solve()
    elapsed = time.perf_counter() - started
    if result.x is None:
        raise RuntimeError(f"OSQP failed: {result.info.status}")
    values = np.asarray(result.x, dtype=float)
    if not np.all(np.isfinite(values)):
        raise RuntimeError(f"OSQP failed: {result.info.status}")
    psi = float(0.5 * values @ (P @ values) + q @ values)
    areas, gradients = _triangle_area_gradients(mesh, values)
    theta = -values
    jumps, lengths = [], []
    edge_to_tri: dict[tuple[int, int], list[int]] = {}
    for ti, tri in enumerate(mesh.triangles):
        for a, b in ((int(tri[0]), int(tri[1])), (int(tri[1]), int(tri[2])), (int(tri[2]), int(tri[0]))):
            edge_to_tri.setdefault(tuple(sorted((a, b))), []).append(ti)
    for (a, b), owners in edge_to_tri.items():
        if len(owners) == 2:
            jump = float(np.linalg.norm(gradients[owners[0]] - gradients[owners[1]]))
            jumps.append(jump); lengths.append(float(np.linalg.norm(mesh.nodes[b] - mesh.nodes[a])))
    jump_arr = np.asarray(jumps) if jumps else np.zeros(1)
    threshold = max(1e-9, 1e-6 * float(np.max(np.linalg.norm(gradients, axis=1))))
    active = jump_arr > threshold
    violation = np.maximum(0, -(A @ values)) if A.shape[0] else np.zeros(1)
    info = result.info
    return LevelResult(values, psi, theta, gradients, str(info.status), int(info.iter),
                       float(info.prim_res), float(info.dual_res), elapsed,
                       float(np.max(violation)), _affine_residual(mesh, values),
                       float(np.max(jump_arr)), float(np.mean(jump_arr)),
                       tuple(float(x) for x in np.quantile(jump_arr, [0.5, .9, .99])),
                       int(np.sum(active)), float(np.sum(np.asarray(lengths)[active])) if lengths else 0., None)


def run_experiment(input_dir: str | Path, output_dir: str | Path, max_level: int = 8,
                   samples: Iterable[str] | None = None) -> list[dict]:
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    names = set(samples) if samples is not None else None
    paths = sorted(p for p in input_dir.iterdir() if p.is_file() and not p.suffix and p.name.startswith("d"))
    rows: list[dict] = []
    for path in paths:
        if names is not None and path.stem not in names and path.name not in names:
            continue
        sample = load_polygon(path)
        previous = None
        previous_mesh = None
        for level in range(max_level + 1):
            mesh = build_mesh(sample, level)
            prolong = None
            if previous is not None and previous_mesh is not None:
                prolong = _prolongate(previous_mesh, mesh, previous)
            result = solve_level(mesh, prolong)
            if previous is not None and previous_mesh is not None:
                diff = result.values_u - prolong
                result.l2_difference = _mass_norm(mesh, diff)
            stem = f"{sample.name}_level{level}"
            np.savez_compressed(output_dir / f"{stem}.npz", nodes=mesh.nodes, triangles=mesh.triangles,
                                values_u=result.values_u, values_theta=result.theta_values,
                                triangle_gradients=result.triangle_gradients)
            metadata = {"sample": sample.name, "level": level, "status": result.status,
                        "iterations": result.iterations, "primal_residual": result.primal_residual,
                        "dual_residual": result.dual_residual, "solve_time": result.solve_time,
                        "psi": result.psi, "max_constraint_violation": result.max_constraint_violation,
                        "affine_fit_residual": result.affine_fit_residual,
                        "gradient_jump_max": result.gradient_jump_max,
                        "gradient_jump_mean": result.gradient_jump_mean,
                        "gradient_jump_quantiles": result.gradient_jump_quantiles,
                        "gradient_jump_count": result.gradient_jump_count,
                        "gradient_jump_length": result.gradient_jump_length,
                        "l2_difference": result.l2_difference,
                        "node_count": len(mesh.nodes), "triangle_count": len(mesh.triangles),
                        "interior_edge_count": sum(len(v) == 2 for v in _edge_owners(mesh).values())}
            (output_dir / f"{stem}.json").write_text(json.dumps(metadata, indent=2))
            rows.append(metadata)
            previous, previous_mesh = result.values_u, mesh
    if rows:
        fields = list(rows[0])
        with (output_dir / "summary.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader(); writer.writerows(rows)
    return rows


def _edge_owners(mesh: Mesh) -> dict[tuple[int, int], list[int]]:
    owners: dict[tuple[int, int], list[int]] = {}
    for ti, tri in enumerate(mesh.triangles):
        for a, b in ((int(tri[0]), int(tri[1])), (int(tri[1]), int(tri[2])), (int(tri[2]), int(tri[0]))):
            owners.setdefault(tuple(sorted((a, b))), []).append(ti)
    return owners


def _prolongate(old_mesh: Mesh, new_mesh: Mesh, old_values: Array) -> Array:
    """Exact P1 prolongation from a red-refined mesh to its child mesh."""
    if len(old_values) != len(old_mesh.nodes):
        raise ValueError("old_values has the wrong length for old_mesh")
    values = np.empty(len(new_mesh.nodes), dtype=float)
    old_node_map = {_key(point): i for i, point in enumerate(old_mesh.nodes)}
    midpoint_map: dict[tuple[int, int], tuple[int, int]] = {}
    for a, b in _edge_owners(old_mesh):
        midpoint_map[_key((old_mesh.nodes[a] + old_mesh.nodes[b]) / 2)] = (a, b)
    for i, point in enumerate(new_mesh.nodes):
        key = _key(point)
        if key in old_node_map:
            values[i] = old_values[old_node_map[key]]
            continue
        # Every new node in standard red refinement is the midpoint of one
        # old edge.  Looking up the edge by its dyadic midpoint is exact.
        edge = midpoint_map.get(key)
        if edge is None:
            raise ValueError("new mesh is not a red refinement of old_mesh")
        a, b = edge
        values[i] = 0.5 * (old_values[a] + old_values[b])
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default="unstable_polytope")
    parser.add_argument("--output-dir", default="results")
    parser.add_argument("--max-level", type=int, default=8)
    parser.add_argument("--sample", action="append", dest="samples")
    args = parser.parse_args()
    run_experiment(args.input_dir, args.output_dir, args.max_level, args.samples)


if __name__ == "__main__":
    main()
