"""Locate a finite-element crease and solve its ell-continuity equations exactly.

The ``convex_p1.py`` run is used only as a geometric guide.  Once a connected
set of large gradient-jump edges has identified the two parent boundary edges,
this program constructs the two sub-polygons symbolically and solves

    ell_left(U) = ell_right(U),
    ell_left(V) = ell_right(V)

over exact rational parameters and isolates the remaining algebraic roots by
exact rational intervals.  No Newton iteration or rational reconstruction is
used for the endpoint equations.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict, deque
from dataclasses import dataclass
from math import gcd
from pathlib import Path
from typing import Iterable

import numpy as np
import sympy as sp

from algebraic_unstable_scan import scan_polygon


@dataclass(frozen=True)
class Polygon:
    vertices: tuple[tuple[int, int], ...]

    @property
    def edge_count(self) -> int:
        return len(self.vertices)

    def edge_measure(self, edge: int) -> int:
        p = self.vertices[edge]
        q = self.vertices[(edge + 1) % self.edge_count]
        return gcd(abs(q[0] - p[0]), abs(q[1] - p[1]))


@dataclass(frozen=True)
class CreaseGuide:
    edges: tuple[tuple[int, int], ...]
    threshold: float
    line_point: tuple[float, float]
    line_direction: tuple[float, float]
    boundary_edge_pairs: tuple[tuple[int, int], ...]
    approximate_endpoints: tuple[tuple[float, float], tuple[float, float]]


def load_polygon(path: Path) -> Polygon:
    points: list[tuple[int, int]] = []
    for raw in path.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        fields = line.replace(",", " ").split()
        if len(fields) < 2:
            continue
        x, y = (int(fields[0]), int(fields[1]))
        points.append((x, y))
    if len(points) < 3:
        raise ValueError(f"{path}: expected at least three vertices")
    area2 = sum(x * Y - y * X for (x, y), (X, Y) in
                zip(points, points[1:] + points[:1]))
    if area2 <= 0:
        raise ValueError(f"{path}: vertices must be counter-clockwise")
    return Polygon(tuple(points))


def _vertex_singularity_flags(path: Path, count: int) -> list[bool]:
    """Read optional vertex_singularity_flags metadata from a polygon file."""
    flags = None
    for raw in path.read_text().splitlines():
        match = re.match(r"\s*#\s*vertex_singularity_flags\s*=\s*([^#]+)", raw)
        if match:
            fields = [x.strip() for x in match.group(1).split(",")]
            if len(fields) == count and all(x in {"0", "1"} for x in fields):
                flags = [x == "1" for x in fields]
                break
    return flags if flags is not None else [False] * count


def _write_crease_svg(path: Path, polygon: Polygon,
                      candidates: list[dict], singular_flags: list[bool]) -> None:
    """Render accepted exact crease candidates in the style of sample SVGs."""
    xs = [p[0] for p in polygon.vertices]
    ys = [p[1] for p in polygon.vertices]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    span = max(xmax - xmin, ymax - ymin, 1.0)
    pad = 0.10 * span
    plot_left, plot_right = 90.0, 1070.0
    plot_top, plot_bottom = 100.0, 670.0
    xlo, xhi = xmin - pad, xmax + pad
    ylo, yhi = ymin - pad, ymax + pad
    sx = (plot_right - plot_left) / (xhi - xlo)
    sy = (plot_bottom - plot_top) / (yhi - ylo)
    scale = min(sx, sy)
    xcenter = (xlo + xhi) / 2
    ycenter = (ylo + yhi) / 2
    xmid = (plot_left + plot_right) / 2
    ymid = (plot_top + plot_bottom) / 2

    def point(x: float, y: float) -> tuple[float, float]:
        return (xmid + (x - xcenter) * scale,
                ymid - (y - ycenter) * scale)

    projected = [point(float(x), float(y)) for x, y in polygon.vertices]
    polygon_points = " ".join(f"{x:.2f},{y:.2f}" for x, y in projected)
    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1120 820" width="1120" height="820">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="560.00" y="30.00" fill="#111827" font-family="Arial, sans-serif" font-size="20px" font-weight="bold" text-anchor="middle">Accepted crease candidates</text>',
        f'<polygon points="{polygon_points}" fill="#eff6ff" stroke="#111827" stroke-width="2"/>',
    ]
    for candidate in candidates:
        edges = candidate["edges"]
        s = float(candidate["s_approximation"])
        t = float(candidate["t_approximation"])
        p = polygon.vertices[edges[0]]
        q = polygon.vertices[(edges[0] + 1) % polygon.edge_count]
        u = (p[0] + s * (q[0] - p[0]), p[1] + s * (q[1] - p[1]))
        p = polygon.vertices[edges[1]]
        q = polygon.vertices[(edges[1] + 1) % polygon.edge_count]
        v = (p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1]))
        a, b = point(*u), point(*v)
        lines.append(f'<line x1="{a[0]:.2f}" y1="{a[1]:.2f}" x2="{b[0]:.2f}" y2="{b[1]:.2f}" stroke="#dc2626" stroke-width="2"/>')
    for index, ((x, y), singular) in enumerate(zip(polygon.vertices, singular_flags), 1):
        px, py = point(float(x), float(y))
        fill = "#a855f7" if singular else "#ffffff"
        lines.append(f'<circle cx="{px:.2f}" cy="{py:.2f}" r="2" fill="{fill}" stroke="#7c3aed" stroke-width="2"/>')
        lines.append(f'<text x="{px + 7:.2f}" y="{py - 7:.2f}" fill="#7c3aed" font-family="Arial, sans-serif" font-size="11px" font-weight="bold">v{index}=({x},{y})</text>')
    lines.extend([
        '<line x1="90" y1="728" x2="122" y2="728" stroke="#dc2626" stroke-width="2"/>',
        '<text x="130" y="733" fill="#dc2626" font-family="Arial, sans-serif" font-size="13px">accepted crease</text>',
        '<circle cx="350" cy="728" r="6" fill="#a855f7" stroke="#7c3aed" stroke-width="2"/>',
        '<text x="365" y="733" fill="#7c3aed" font-family="Arial, sans-serif" font-size="13px">singular vertex</text>',
        '<circle cx="530" cy="728" r="6" fill="#ffffff" stroke="#7c3aed" stroke-width="2"/>',
        '<text x="545" y="733" fill="#7c3aed" font-family="Arial, sans-serif" font-size="13px">smooth vertex</text>',
    ])
    if candidates:
        for i, candidate in enumerate(candidates):
            lines.append(f'<text x="90" y="{762 + 18 * i}" fill="#dc2626" font-family="Arial, sans-serif" font-size="12px">candidate {i + 1}: edges {candidate["edges"]}, s={candidate["s_approximation"]:.12g}, t={candidate["t_approximation"]:.12g}</text>')
    else:
        lines.append('<text x="90" y="762" fill="#6b7280" font-family="Arial, sans-serif" font-size="12px">No accepted crease candidate</text>')
    lines.append('</svg>')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _owners(triangles: np.ndarray) -> dict[tuple[int, int], list[int]]:
    owners: dict[tuple[int, int], list[int]] = defaultdict(list)
    for ti, tri in enumerate(triangles):
        for a, b in ((int(tri[0]), int(tri[1])),
                     (int(tri[1]), int(tri[2])),
                     (int(tri[2]), int(tri[0]))):
            owners[tuple(sorted((a, b)))].append(ti)
    return owners


def _components(edges: Iterable[tuple[int, int]]) -> list[list[tuple[int, int]]]:
    edges = list(edges)
    incident: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for edge in edges:
        incident[edge[0]].append(edge)
        incident[edge[1]].append(edge)
    unseen = set(edges)
    result: list[list[tuple[int, int]]] = []
    while unseen:
        first = next(iter(unseen))
        queue = deque([first])
        unseen.remove(first)
        component = []
        while queue:
            edge = queue.popleft()
            component.append(edge)
            for vertex in edge:
                for neighbor in incident[vertex]:
                    if neighbor in unseen:
                        unseen.remove(neighbor)
                        queue.append(neighbor)
        result.append(component)
    return sorted(result, key=len, reverse=True)


def _line_intersection(p: np.ndarray, q: np.ndarray, normal: np.ndarray,
                       constant: float) -> tuple[float, np.ndarray] | None:
    direction = q - p
    denominator = float(np.dot(normal, direction))
    if abs(denominator) < 1e-14:
        return None
    parameter = -float(np.dot(normal, p) + constant) / denominator
    if parameter < -1e-7 or parameter > 1.0 + 1e-7:
        return None
    return parameter, p + parameter * direction


def detect_crease(npz_path: Path, polygon: Polygon, threshold: float | None,
                  component_index: int = 0) -> CreaseGuide:
    with np.load(npz_path) as data:
        nodes = np.asarray(data["nodes"], dtype=float)
        triangles = np.asarray(data["triangles"], dtype=np.int64)
        gradients = np.asarray(data["triangle_gradients"], dtype=float)
    if nodes.ndim != 2 or nodes.shape[1] != 2:
        raise ValueError("NPZ nodes must have shape (N,2)")
    owners = _owners(triangles)
    scale = float(np.max(np.linalg.norm(gradients, axis=1)))
    actual_threshold = (max(1e-9, 1e-3 * scale)
                        if threshold is None else float(threshold))
    marked: list[tuple[int, int]] = []
    for edge, adjacent in owners.items():
        if len(adjacent) == 2:
            jump = float(np.linalg.norm(gradients[adjacent[0]] -
                                        gradients[adjacent[1]]))
            if jump >= actual_threshold:
                marked.append(edge)
    components = _components(marked)
    if not components or component_index >= len(components):
        raise ValueError("no crease component found at the selected threshold")
    component = components[component_index]
    point_indices = sorted({v for edge in component for v in edge})
    points = nodes[point_indices]
    center = points.mean(axis=0)
    _, _, right = np.linalg.svd(points - center, full_matrices=False)
    direction = np.asarray(right[0], dtype=float)
    direction /= np.linalg.norm(direction)
    normal = np.asarray((-direction[1], direction[0]))
    constant = -float(np.dot(normal, center))

    intersections: list[tuple[int, np.ndarray]] = []
    for edge in range(polygon.edge_count):
        p = np.asarray(polygon.vertices[edge], dtype=float)
        q = np.asarray(polygon.vertices[(edge + 1) % polygon.edge_count], dtype=float)
        hit = _line_intersection(p, q, normal, constant)
        if hit is not None:
            intersections.append((edge, hit[1]))
    if len(intersections) < 2:
        raise ValueError(
            f"crease line intersects {len(intersections)} boundary edges; "
            "increase the mesh level or adjust --threshold")
    intersections.sort(key=lambda item: float(np.dot(item[1] - center, direction)))
    # If the fitted line passes through a polygon vertex, both incident
    # boundary edges report the same geometric endpoint.  Keep all incident
    # edge choices: the exact solve will decide which pair gives a valid
    # interior chord.  This is essential for creases ending near vertices.
    endpoint_groups: list[list[tuple[int, np.ndarray]]] = []
    for item in intersections:
        if not endpoint_groups or np.linalg.norm(item[1] - endpoint_groups[-1][0][1]) > 1e-7:
            endpoint_groups.append([item])
        else:
            endpoint_groups[-1].append(item)
    if len(endpoint_groups) != 2:
        raise ValueError(
            f"crease line has {len(endpoint_groups)} boundary intersections; "
            "increase the mesh level or adjust --threshold")
    boundary_edge_pairs = tuple(
        (left[0], right[0])
        for left in endpoint_groups[0]
        for right in endpoint_groups[1]
        if left[0] != right[0])
    return CreaseGuide(tuple(component), actual_threshold,
                       tuple(float(x) for x in center),
                       tuple(float(x) for x in direction),
                       boundary_edge_pairs,
                       (tuple(float(x) for x in endpoint_groups[0][0][1]),
                        tuple(float(x) for x in endpoint_groups[1][0][1])))


def _det3(matrix: list[list[sp.Expr]]) -> sp.Expr:
    return (matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
            - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
            + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]))


def _moments(vertices: list[tuple[sp.Expr, sp.Expr]]) -> list[sp.Expr]:
    area = ix = iy = ixx = ixy = iyy = sp.Integer(0)
    for (x, y), (X, Y) in zip(vertices, vertices[1:] + vertices[:1]):
        cross = x * Y - y * X
        area += cross
        ix += (x + X) * cross
        iy += (y + Y) * cross
        ixx += (x * x + x * X + X * X) * cross
        ixy += (2 * x * y + x * Y + X * y + 2 * X * Y) * cross
        iyy += (y * y + y * Y + Y * Y) * cross
    return [area / 2, ix / 6, iy / 6, ixx / 12, ixy / 24, iyy / 12]


def _ell_numerator(vertices: list[tuple[sp.Expr, sp.Expr]],
                   measures: list[sp.Expr]) -> tuple[list[sp.Expr], sp.Expr]:
    area, ix, iy, ixx, ixy, iyy = _moments(vertices)
    length = boundary_x = boundary_y = sp.Integer(0)
    for (x, y), (X, Y), measure in zip(vertices, vertices[1:] + vertices[:1], measures):
        length += measure
        boundary_x += measure * (x + X) / 2
        boundary_y += measure * (y + Y) / 2
    matrix = [[area, ix, iy], [ix, ixx, ixy], [iy, ixy, iyy]]
    rhs = [length, boundary_x, boundary_y]
    denominator = _det3(matrix)
    numerators = []
    for column in range(3):
        replaced = [row[:] for row in matrix]
        for row in range(3):
            replaced[row][column] = rhs[row]
        numerators.append(_det3(replaced))
    return numerators, denominator


def _symbolic_split(polygon: Polygon, first_edge: int, second_edge: int,
                    s: sp.Symbol, t: sp.Symbol):
    if first_edge == second_edge:
        raise ValueError("chord endpoints must lie on different edges")
    if first_edge > second_edge:
        first_edge, second_edge = second_edge, first_edge
        s, t = t, s
    vertices = [(sp.Rational(x), sp.Rational(y)) for x, y in polygon.vertices]
    n = len(vertices)
    u = tuple(vertices[first_edge][k] + s *
              (vertices[(first_edge + 1) % n][k] - vertices[first_edge][k])
              for k in range(2))
    v = tuple(vertices[second_edge][k] + t *
              (vertices[(second_edge + 1) % n][k] - vertices[second_edge][k])
              for k in range(2))
    measures = [sp.Integer(polygon.edge_measure(i)) for i in range(n)]
    left_vertices = [u] + vertices[first_edge + 1:second_edge + 1] + [v]
    left_measures = ([measures[first_edge] * (1 - s)] +
                     measures[first_edge + 1:second_edge] +
                     [measures[second_edge] * t, sp.Integer(0)])
    right_vertices = [v] + vertices[second_edge + 1:] + vertices[:first_edge + 1] + [u]
    right_measures = ([measures[second_edge] * (1 - t)] +
                      measures[second_edge + 1:] + measures[:first_edge] +
                      [measures[first_edge] * s, sp.Integer(0)])
    return u, v, left_vertices, left_measures, right_vertices, right_measures


def continuity_polynomials(polygon: Polygon, first_edge: int, second_edge: int):
    s, t = sp.symbols("s t")
    u, v, left_vertices, left_measures, right_vertices, right_measures = \
        _symbolic_split(polygon, first_edge, second_edge, s, t)
    left_num, left_den = _ell_numerator(left_vertices, left_measures)
    right_num, right_den = _ell_numerator(right_vertices, right_measures)

    def value_difference(point):
        left_value = left_num[0] + left_num[1] * point[0] + left_num[2] * point[1]
        right_value = right_num[0] + right_num[1] * point[0] + right_num[2] * point[1]
        return sp.Poly(sp.expand(left_value * right_den - right_value * left_den), s, t)

    first = value_difference(u)
    second = value_difference(v)
    common = sp.gcd(first, second)
    if not common.is_zero and common.total_degree() > 0:
        # A shared factor makes the raw resultant contain a whole component.
        # Remove it before the zero-dimensional resultant computation, but
        # retain the factor in the diagnostic output.  The denominator check
        # rejects reconstructed roots on singular branches; a non-singular
        # common component is explicitly flagged for further analysis.
        first = first.exquo(common)
        second = second.exquo(common)
    return s, t, first, second, common, left_den, right_den, u, v


def _triangular_relation(first: sp.Poly, second: sp.Poly,
                         s: sp.Symbol, t: sp.Symbol):
    """Return ``A(s), B(s)`` for a generic relation ``A(s)t+B(s)=0``.

    The last degree-one subresultant is a certified relation at every
    projected root where its leading coefficient does not vanish.  Keeping
    this relation lets us pair the two univariate resultant roots without
    introducing a floating-point root finder.
    """
    subresultants = sp.subresultants(first.as_expr(), second.as_expr(), t)
    for expression in reversed(subresultants):
        candidate = sp.Poly(expression, t)
        if candidate.degree() != 1:
            continue
        coefficient = sp.Poly(candidate.coeff_monomial(t), s, domain="QQ")
        constant = sp.Poly(candidate.coeff_monomial(1), s, domain="QQ")
        # Clear denominators *jointly*.  Making the two coefficients
        # primitive independently would change their relative scale and
        # destroy the subresultant relation.
        denominators = [int(sp.denom(value))
                        for value in coefficient.all_coeffs() + constant.all_coeffs()]
        scale = sp.ilcm(*denominators) if denominators else 1
        coefficient = sp.Poly(sp.expand(coefficient.as_expr() * scale), s,
                              domain="ZZ")
        constant = sp.Poly(sp.expand(constant.as_expr() * scale), s,
                           domain="ZZ")
        contents = [abs(int(value)) for value in
                    coefficient.all_coeffs() + constant.all_coeffs() if value]
        common_content = 0
        for value in contents:
            common_content = gcd(common_content, value)
        if common_content > 1:
            coefficient = coefficient.quo_ground(common_content)
            constant = constant.quo_ground(common_content)
        return coefficient, constant
    return None, None


def _rational_roots(poly: sp.Poly, variable: sp.Symbol) -> list[sp.Rational]:
    return sorted({root for root, multiplicity in sp.polys.polytools.ground_roots(poly).items()
                   if root.is_Rational and 0 <= root <= 1})


def _primitive_integer_poly(poly: sp.Poly) -> sp.Poly:
    _, integer_poly = poly.clear_denoms(convert=True)
    _, primitive = integer_poly.primitive()
    return primitive


def _minimal_factor(poly: sp.Poly, interval: list[str], variable: sp.Symbol) -> dict:
    """Identify the irreducible resultant factor containing an isolated root."""
    lo, hi = (sp.Rational(x) for x in interval)
    factors = sp.factor_list(poly.as_expr())[1]
    matches = []
    for factor, multiplicity in factors:
        candidate = sp.Poly(factor, variable, domain="ZZ")
        roots = candidate.intervals(eps=sp.Rational(1, 10**18))
        if any(a <= hi and b >= lo for (a, b), _ in roots):
            matches.append((candidate, multiplicity))
    if not matches:
        return {"polynomial": str(poly.as_expr()), "degree": int(poly.degree()),
                "multiplicity": 1, "irreducible": False}
    candidate, multiplicity = matches[0]
    return {"polynomial": str(candidate.as_expr()),
            "degree": int(candidate.degree()),
            "multiplicity": int(multiplicity), "irreducible": True}


def _isolating_intervals(poly: sp.Poly, variable: sp.Symbol,
                         eps: sp.Rational = sp.Rational(1, 10**12)) -> list[dict]:
    result = []
    for interval, multiplicity in poly.intervals(eps=eps):
        if len(interval) != 2 or not all(value.is_Rational for value in interval):
            continue
        left, right = interval
        if right < 0 or left > 1:
            continue
        result.append({"interval": [str(left), str(right)],
                       "multiplicity": int(multiplicity)})
    return result


def _interval_poly(poly: sp.Poly, interval: tuple[sp.Rational, sp.Rational]):
    """Evaluate a univariate polynomial on a rational interval.

    Horner evaluation with rational endpoints is deliberately conservative;
    unlike a decimal evaluation it remains a proof that the value lies in
    the returned interval.
    """
    lo, hi = interval
    lower = upper = sp.Rational(0)
    for coefficient in poly.all_coeffs():
        products = (lower * lo, lower * hi, upper * lo, upper * hi)
        coefficient_values = (coefficient * lo, coefficient * hi)
        lower = min(products) + min(coefficient_values)
        upper = max(products) + max(coefficient_values)
    return lower, upper


def _interval_poly2(poly: sp.Poly,
                    s_interval: tuple[sp.Rational, sp.Rational],
                    t_interval: tuple[sp.Rational, sp.Rational]):
    """Conservative exact interval evaluation for a bivariate polynomial."""
    slo, shi = s_interval
    tlo, thi = t_interval

    def multiply(left, right):
        values = (left[0] * right[0], left[0] * right[1],
                  left[1] * right[0], left[1] * right[1])
        return min(values), max(values)

    result = (sp.Rational(0), sp.Rational(0))
    for (s_degree, t_degree), coefficient in poly.terms():
        s_power = (sp.Rational(1), sp.Rational(1))
        t_power = (sp.Rational(1), sp.Rational(1))
        for _ in range(s_degree):
            s_power = multiply(s_power, (slo, shi))
        for _ in range(t_degree):
            t_power = multiply(t_power, (tlo, thi))
        term = multiply(s_power, t_power)
        term = (term[0] * coefficient, term[1] * coefficient)
        if term[0] > term[1]:
            term = (term[1], term[0])
        result = (result[0] + term[0], result[1] + term[1])
    return result


def _paired_root_boxes(data: dict) -> list[dict]:
    """Pair resultant roots using the exact subresultant relation.

    A box is emitted only when ``A(s)`` and both moment denominators are
    provably nonzero on the rational box.  The displayed box therefore
    represents one admissible real algebraic solution, rather than an
    arbitrary Cartesian product of two one-dimensional root lists.
    """
    relation_a = data.get("relation_a")
    relation_b = data.get("relation_b")
    if relation_a is None or relation_b is None:
        return []
    s_poly = data["resultant_s"]
    t_poly = data["resultant_t"]
    # Tight intervals are used internally because the subresultant
    # coefficients can have degree near 100.  They are still exact rational
    # intervals in the JSON output.
    epsilon = sp.Rational(1, 10**18)
    s_intervals = s_poly.intervals(eps=epsilon)
    t_intervals = t_poly.intervals(eps=epsilon)
    boxes = []
    for (s_interval, multiplicity) in s_intervals:
        slo, shi = s_interval
        if shi < 0 or slo > 1:
            continue
        slo, shi = max(slo, sp.Rational(0)), min(shi, sp.Rational(1))
        a_lo, a_hi = _interval_poly(relation_a, (slo, shi))
        if a_lo <= 0 <= a_hi:
            # The linear subresultant is not a separating relation on this
            # interval (typically an endpoint or a singular branch).
            continue
        b_lo, b_hi = _interval_poly(relation_b, (slo, shi))
        values = (-b_lo / a_lo, -b_lo / a_hi,
                  -b_hi / a_lo, -b_hi / a_hi)
        t_lo, t_hi = max(sp.Rational(0), min(values)), min(sp.Rational(1), max(values))
        if t_lo > t_hi:
            continue
        for (candidate_interval, t_multiplicity) in t_intervals:
            candidate_lo, candidate_hi = candidate_interval
            if candidate_hi < t_lo or candidate_lo > t_hi:
                continue
            box_t = (max(candidate_lo, t_lo), min(candidate_hi, t_hi))
            if box_t[0] > box_t[1]:
                continue
            left_interval = _interval_poly2(
                sp.Poly(data["left_den"], data["s"], data["t"]),
                (slo, shi), box_t)
            right_interval = _interval_poly2(
                sp.Poly(data["right_den"], data["s"], data["t"]),
                (slo, shi), box_t)
            if left_interval[0] <= 0 <= left_interval[1] or \
                    right_interval[0] <= 0 <= right_interval[1]:
                continue
            boxes.append({
                "s_interval": [str(slo), str(shi)],
                "t_interval": [str(box_t[0]), str(box_t[1])],
                "s_approximation": float((slo + shi) / 2),
                "t_approximation": float((box_t[0] + box_t[1]) / 2),
                "s_multiplicity": int(multiplicity),
                "t_multiplicity": int(t_multiplicity),
                "minimal_polynomial_s": _minimal_factor(
                    s_poly, [str(slo), str(shi)], data["s"]),
                "minimal_polynomial_t": _minimal_factor(
                    t_poly, [str(box_t[0]), str(box_t[1])], data["t"]),
            })
    return boxes


def _numeric_piece(polygon: Polygon, first_edge: int, second_edge: int,
                   s_value: float, t_value: float):
    s, t = sp.symbols("s t")
    _, _, left_vertices, left_measures, right_vertices, right_measures = \
        _symbolic_split(polygon, first_edge, second_edge, s, t)
    substitution = {s: sp.Float(s_value, 30), t: sp.Float(t_value, 30)}
    left = np.asarray([[float(sp.N(x.subs(substitution), 25)),
                        float(sp.N(y.subs(substitution), 25))]
                       for x, y in left_vertices], dtype=float)
    right = np.asarray([[float(sp.N(x.subs(substitution), 25)),
                         float(sp.N(y.subs(substitution), 25))]
                        for x, y in right_vertices], dtype=float)
    lm = [float(sp.N(value.subs(substitution), 25)) for value in left_measures]
    rm = [float(sp.N(value.subs(substitution), 25)) for value in right_measures]
    return left, lm, right, rm


def _concavity_scan(left: np.ndarray, left_measures: list[float],
                    right: np.ndarray, right_measures: list[float]) -> dict:
    try:
        from algebraic_unstable_scan import _ell
        left_ell = _ell(left, left_measures)
        right_ell = _ell(right, right_measures)
    except ValueError as exc:
        return {"status": "inconclusive", "reason": str(exc)}
    difference = left_ell - right_ell
    left_values = np.asarray([difference[0] + difference[1] * x + difference[2] * y
                              for x, y in left])
    right_values = np.asarray([difference[0] + difference[1] * x + difference[2] * y
                               for x, y in right])
    scale = max(1.0, float(np.max(np.abs(left_values))),
                float(np.max(np.abs(right_values))))
    tolerance = 1e-8 * scale
    left_max = float(np.max(left_values))
    right_min = float(np.min(right_values))
    if left_max <= tolerance and right_min >= -tolerance:
        status = "pass"
    elif left_max > 10 * tolerance or right_min < -10 * tolerance:
        status = "fail"
    else:
        status = "inconclusive"
    return {"status": status, "left_max_difference": left_max,
            "right_min_difference": right_min, "tolerance": tolerance,
            "left_ell": [float(x) for x in left_ell],
            "right_ell": [float(x) for x in right_ell]}


def _scan_root_box(polygon: Polygon, edges: tuple[int, int], box: dict,
                   theta_steps: int, t_steps: int, precision: int,
    refine: bool) -> dict:
    # A Cartesian isolating box that is still broad compared with the
    # requested numerical precision cannot certify which side of a crease a
    # test line sees.  Keep the exact box in the report and defer the scan.
    widths = [float(sp.Rational(box[key][1]) - sp.Rational(box[key][0]))
              for key in ("s_interval", "t_interval")]
    # The downstream geometric scan is evaluated in IEEE double precision;
    # requiring a box narrower than 10^(-precision/2) would reject perfectly
    # useful exact isolating boxes (for example 10^-19 boxes at precision 40).
    # Treat boxes below the requested numerical resolution as point samples.
    # Only boxes wider than roughly the square root of machine precision, or
    # wider than the user-requested decimal resolution, are inconclusive.
    resolution = max(1e-12, 10.0 ** (-max(8, min(precision, 24))))
    if max(widths) > resolution:
        return {
            "edges": list(edges), "s_interval": box["s_interval"],
            "t_interval": box["t_interval"],
            "s_approximation": box["s_approximation"],
            "t_approximation": box["t_approximation"],
            "minimal_polynomial_s": box.get("minimal_polynomial_s"),
            "minimal_polynomial_t": box.get("minimal_polynomial_t"),
            "concavity": {"status": "inconclusive", "reason":
                           "isolating box wider than scan precision",
                           "box_widths": widths},
            "left_piece_scan": {"relative_df_scan_status": "inconclusive",
                                 "reason": "isolating box wider than scan precision"},
            "right_piece_scan": {"relative_df_scan_status": "inconclusive",
                                  "reason": "isolating box wider than scan precision"},
            "accepted": False, "rejection_reason": "root_box_too_wide",
        }
    s_value = box["s_approximation"]
    t_value = box["t_approximation"]
    left, left_measures, right, right_measures = _numeric_piece(
        polygon, edges[0], edges[1], s_value, t_value)
    concavity = _concavity_scan(left, left_measures, right, right_measures)
    left_scan = scan_polygon(left, left_measures, theta_steps, t_steps,
                             precision, refine)
    right_scan = scan_polygon(right, right_measures, theta_steps, t_steps,
                              precision, refine)
    statuses = {left_scan.status, right_scan.status}
    if concavity["status"] != "pass":
        accepted = False
        reason = "concavity_" + concavity["status"]
    elif "unstable" in statuses:
        accepted = False
        reason = "piece_unstable"
    elif "inconclusive" in statuses:
        accepted = False
        reason = "piece_scan_inconclusive"
    else:
        accepted = True
        reason = None
    return {
        "edges": list(edges),
        "s_interval": box["s_interval"], "t_interval": box["t_interval"],
        "s_approximation": s_value, "t_approximation": t_value,
        "minimal_polynomial_s": box.get("minimal_polynomial_s"),
        "minimal_polynomial_t": box.get("minimal_polynomial_t"),
        "concavity": concavity,
        "left_piece_scan": left_scan.as_dict(),
        "right_piece_scan": right_scan.as_dict(),
        "accepted": accepted, "rejection_reason": reason,
    }


def _system_data(polygon: Polygon, first_edge: int, second_edge: int) -> dict:
    if first_edge > second_edge:
        first_edge, second_edge = second_edge, first_edge
    s, t, first, second, common, left_den, right_den, u, v = \
        continuity_polynomials(polygon, first_edge, second_edge)
    first_integer = _primitive_integer_poly(first)
    second_integer = _primitive_integer_poly(second)
    resultant_s = _primitive_integer_poly(
        sp.Poly(sp.resultant(first_integer.as_expr(), second_integer.as_expr(), t), s))
    resultant_t = _primitive_integer_poly(
        sp.Poly(sp.resultant(first_integer.as_expr(), second_integer.as_expr(), s), t))
    relation_a, relation_b = _triangular_relation(first_integer, second_integer,
                                                   s, t)
    denominator_product = sp.Poly(sp.expand(left_den * right_den), s, t)
    singular_common = (not common.is_zero and common.total_degree() > 0 and
                       sp.gcd(common, denominator_product).total_degree() ==
                       common.total_degree())
    return {
        "first_edge": first_edge,
        "second_edge": second_edge,
        "s": s,
        "t": t,
        "first": first,
        "second": second,
        "first_integer": first_integer,
        "second_integer": second_integer,
        "resultant_s": resultant_s,
        "resultant_t": resultant_t,
        "common": common,
        "common_denominator_gcd": sp.gcd(common, denominator_product),
        "common_factor_is_singular": singular_common,
        "common_factor_requires_review": (common.total_degree() > 0 and
                                           not singular_common),
        "left_den": left_den,
        "right_den": right_den,
        "relation_a": relation_a,
        "relation_b": relation_b,
    }


def solve_exact(polygon: Polygon, first_edge: int, second_edge: int,
                data: dict | None = None) -> list[dict]:
    data = data or _system_data(polygon, first_edge, second_edge)
    s = data["s"]
    t = data["t"]
    first = data["first_integer"]
    second = data["second_integer"]
    left_den = data["left_den"]
    right_den = data["right_den"]
    resultant = data["resultant_s"]
    candidates: list[dict] = []
    for s_root in _rational_roots(resultant, s):
        p1 = sp.Poly(first.as_expr().subs(s, s_root), t)
        p2 = sp.Poly(second.as_expr().subs(s, s_root), t)
        common_t = sp.gcd(p1, p2)
        for t_root in _rational_roots(common_t, t):
            if not (0 <= t_root <= 1):
                continue
            substitutions = {s: s_root, t: t_root}
            # Clearing the ell denominators can introduce roots at which a
            # moment matrix is singular.  Such points are not polygons with
            # a defined ell function and must be rejected explicitly.
            if sp.simplify(left_den.subs(substitutions)) == 0 or \
                    sp.simplify(right_den.subs(substitutions)) == 0:
                continue
            candidates.append({"s": s_root, "t": t_root})
    return candidates


def system_summary(data: dict) -> dict:
    paired_boxes = _paired_root_boxes(data)
    relation = None
    if data.get("relation_a") is not None:
        relation = {
            "equation": str(data["relation_a"].as_expr()) + "*t + " +
                        str(data["relation_b"].as_expr()),
            "valid_when": str(data["relation_a"].as_expr()) + " != 0",
        }
    return {
        "equation_s": str(data["first_integer"].as_expr()),
        "equation_t": str(data["second_integer"].as_expr()),
        "common_factor_removed": str(data["common"].as_expr()),
        "common_denominator_gcd": str(data["common_denominator_gcd"].as_expr()),
        "common_factor_is_singular": data["common_factor_is_singular"],
        "common_factor_requires_review": data["common_factor_requires_review"],
        "resultant_s": str(data["resultant_s"].as_expr()),
        "resultant_t": str(data["resultant_t"].as_expr()),
        "equation_s_degree": int(data["first_integer"].total_degree()),
        "equation_t_degree": int(data["second_integer"].total_degree()),
        "resultant_s_degree": int(data["resultant_s"].degree()),
        "resultant_t_degree": int(data["resultant_t"].degree()),
        "equation_s_terms": len(data["first_integer"].terms()),
        "equation_t_terms": len(data["second_integer"].terms()),
        "s_root_isolating_intervals": _isolating_intervals(data["resultant_s"], data["s"]),
        "t_root_isolating_intervals": _isolating_intervals(data["resultant_t"], data["t"]),
        "triangular_relation_t": relation,
        "paired_root_boxes": paired_boxes,
    }


def _edge_parameter(point: tuple[float, float], polygon: Polygon,
                    edge: int) -> float:
    p = np.asarray(polygon.vertices[edge], dtype=float)
    q = np.asarray(polygon.vertices[(edge + 1) % polygon.edge_count], dtype=float)
    d = q - p
    return float(np.dot(np.asarray(point) - p, d) / np.dot(d, d))


def _exact_endpoint(polygon: Polygon, edge: int, parameter: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    p = polygon.vertices[edge]
    q = polygon.vertices[(edge + 1) % polygon.edge_count]
    return (sp.Rational(p[0]) + parameter * (q[0] - p[0]),
            sp.Rational(p[1]) + parameter * (q[1] - p[1]))


def _affine_interval(polygon: Polygon, edge: int,
                     interval: list[str]) -> list[list[str]]:
    """Map a parameter isolation interval to coordinate intervals."""
    lo, hi = (sp.Rational(value) for value in interval)
    p = polygon.vertices[edge]
    q = polygon.vertices[(edge + 1) % polygon.edge_count]
    result = []
    for coordinate in range(2):
        values = (sp.Rational(p[coordinate]) + lo * (q[coordinate] - p[coordinate]),
                  sp.Rational(p[coordinate]) + hi * (q[coordinate] - p[coordinate]))
        result.append([str(min(values)), str(max(values))])
    return result


def _interval_midpoint(interval: list[str]) -> float:
    """Return a decimal midpoint for a rational isolating interval."""
    lo, hi = (sp.Rational(value) for value in interval)
    return float((lo + hi) / 2)


def _format_float(value: float) -> str:
    return f"{value:.15g}"


def _json_value(value):
    if isinstance(value, sp.Basic):
        return str(value)
    if isinstance(value, tuple):
        return [_json_value(x) for x in value]
    return value


def _write_system_details(path: Path, system: dict) -> None:
    """Write the unabridged symbolic data for one candidate edge pair."""
    lines = [
        f"Exact crease decomposition details for edges {system['edges']}",
        "=" * 72,
        "",
        "Continuity equation at first endpoint:",
        system["equation_s"],
        "",
        "Continuity equation at second endpoint:",
        system["equation_t"],
        "",
        "Resultant in s:",
        system["resultant_s"],
        "",
        "Resultant in t:",
        system["resultant_t"],
        "",
        "Common factor removed:",
        system["common_factor_removed"],
        "",
        "Common denominator gcd:",
        system["common_denominator_gcd"],
        "",
    ]
    relation = system.get("triangular_relation_t")
    if relation:
        lines.extend([
            "Subresultant pairing relation:",
            relation["equation"],
            f"Valid when: {relation['valid_when']}",
            "",
        ])
    lines.extend([
        "s root isolating intervals in [0,1]:",
        json.dumps(system["s_root_isolating_intervals"], ensure_ascii=False, indent=2),
        "",
        "t root isolating intervals in [0,1]:",
        json.dumps(system["t_root_isolating_intervals"], ensure_ascii=False, indent=2),
        "",
    ])
    for box in system.get("paired_root_boxes", []):
        lines.extend(["Minimal polynomials for paired root:",
                      json.dumps({"s": box.get("minimal_polynomial_s"),
                                  "t": box.get("minimal_polynomial_t")},
                                 ensure_ascii=False, indent=2), ""])
    path.write_text("\n".join(lines) + "\n")


def _markdown_report(result: dict, detail_files: dict[tuple[int, int], str]) -> str:
    """Render the machine-readable result as a compact human-readable report."""
    lines = [
        "# Exact crease decomposition",
        "",
        f"- Crease edges: `{result['crease_edge_count']}`",
        f"- Threshold: `{result['crease_threshold']}`",
        f"- Approximate endpoints: `{result['approximate_endpoints']}`",
        f"- Candidate boundary edge pairs: `{result['boundary_edge_pairs']}`",
        f"- SVG diagram: `{result.get('svg_output', '')}` ({result.get('svg_candidate_count', 0)} accepted crease(s))",
        "",
        "## Selected algebraic branch",
        "",
    ]
    selected = result["selected_algebraic_root_boxes"]
    if not selected:
        lines.append("No branch passed the finite-element parameter window.")
    else:
        lines.extend([
            "| Edges | `s` approx. | `t` approx. | `s` isolating interval | `t` isolating interval | Endpoint approximations | Endpoint coordinate intervals |",
            "|---|---:|---:|---|---|---|---|",
        ])
        for box in selected:
            endpoints = "<br>".join(
                f"({coords[0]}, {coords[1]})"
                for coords in box["endpoint_coordinate_intervals"])
            endpoint_approximations = "<br>".join(
                f"({_format_float(coords[0])}, {_format_float(coords[1])})"
                for coords in box["endpoint_approximate_coordinates"])
            lines.append(
                f"| `{box['edges']}` | `{_format_float(box['s_approximation'])}` | "
                f"`{_format_float(box['t_approximation'])}` | "
                f"`{box['s_interval']}` | `{box['t_interval']}` | "
                f"{endpoint_approximations} | {endpoints} |")
    lines.extend(["", "## Candidate systems", ""])
    for system in result["exact_systems"]:
        lines.extend([
            f"### Edges `{system['edges']}`",
            "",
            f"- Continuity equation degrees: `{system['equation_s_degree']}` and "
            f"`{system['equation_t_degree']}` "
            f"({system['equation_s_terms']} and {system['equation_t_terms']} terms)",
            f"- Resultant degrees: `s={system['resultant_s_degree']}`, "
            f"`t={system['resultant_t_degree']}`",
            f"- Common factor removed: `{system['common_factor_removed']}`",
            f"- Common factor is singular: `{system['common_factor_is_singular']}`",
            f"- Common factor requires review: `{system['common_factor_requires_review']}`",
            f"- Paired root boxes: `{len(system['paired_root_boxes'])}`",
            "",
        ])
        detail_file = detail_files[tuple(system["edges"])]
        lines.extend([
            f"- Full equations, resultants and root lists: [`{detail_file}`]({detail_file})",
            "",
        ])
        if system["paired_root_boxes"]:
            lines.append("| `s` approx. | `t` approx. | min poly degrees | `s` interval | `t` interval | Multiplicities |")
            lines.append("|---:|---:|---|---|---|---|")
            for box in system["paired_root_boxes"]:
                lines.append(
                    f"| `{_format_float(box['s_approximation'])}` | "
                    f"`{_format_float(box['t_approximation'])}` | "
                    f"{box.get('minimal_polynomial_s', {}).get('degree', '?')}, "
                    f"{box.get('minimal_polynomial_t', {}).get('degree', '?')} | "
                    f"`{box['s_interval']}` | `{box['t_interval']}` | "
                    f"`{box['s_multiplicity']}, {box['t_multiplicity']}` |")
            lines.append("")
    lines.extend(["## Candidate screening", ""])
    for candidate in result.get("candidate_screening", []):
        lines.extend([f"### Edges `{candidate['edges']}`; "
                      f"s={_format_float(candidate['s_approximation'])}, "
                      f"t={_format_float(candidate['t_approximation'])}", "",
                      f"- Accepted: `{candidate['accepted']}`",
                      f"- Rejection reason: `{candidate.get('rejection_reason')}`",
                      f"- Concavity: `{candidate['concavity']['status']}`",
                      f"- Left scan: `{candidate['left_piece_scan']['relative_df_scan_status']}`",
                      f"- Right scan: `{candidate['right_piece_scan']['relative_df_scan_status']}`", ""])
    lines.extend([
        "## Machine-readable output",
        "",
        "The complete structured result is stored in the accompanying `.json` file.",
        "",
    ])
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polygon", type=Path, required=True)
    parser.add_argument("--npz", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=None)
    parser.add_argument("--component", type=int, default=0)
    parser.add_argument("--branch-window", type=float, default=0.5,
                        help="maximum endpoint-parameter distance used only to select "
                             "the crease branch (default: 0.5)")
    parser.add_argument("--edge-pair", action="append", default=[], metavar="I,J",
                        help="also solve this explicit boundary edge pair")
    parser.add_argument("--scan-theta-steps", type=int, default=180)
    parser.add_argument("--scan-t-steps", type=int, default=128)
    parser.add_argument("--scan-precision", type=int, default=40)
    parser.add_argument("--scan-refine", action="store_true")
    parser.add_argument("--keep-inconclusive", action="store_true")
    parser.add_argument("--show-crease-edges", action="store_true",
                        help="include all marked finite-element edges in the JSON output")
    parser.add_argument("--output", type=Path,
                        help="write the JSON result to this file (default: next to --npz, "
                             "with suffix _exact.json)")
    parser.add_argument("--svg-output", type=Path,
                        help="write accepted crease diagram (default: <sample>_crease.svg next to --npz)")
    args = parser.parse_args()

    if args.output is None:
        sample_stem = re.sub(r"_level\d+$", "", args.npz.stem)
        args.output = args.npz.with_name(f"{sample_stem}_exact.json")
    else:
        sample_stem = re.sub(r"_exact$", "", args.output.stem)
    if args.svg_output is None:
        args.svg_output = args.npz.with_name(f"{sample_stem}_crease.svg")

    polygon = load_polygon(args.polygon)
    guide = detect_crease(args.npz, polygon, args.threshold, args.component)
    explicit_pairs = []
    for value in args.edge_pair:
        fields = value.replace(" ", "").split(",")
        if len(fields) != 2:
            raise ValueError(f"invalid --edge-pair {value!r}; expected I,J")
        explicit_pairs.append((int(fields[0]), int(fields[1])))
    boundary_pairs = tuple(dict.fromkeys(guide.boundary_edge_pairs + tuple(explicit_pairs)))
    all_solutions = []
    selected_root_boxes = []
    systems = []
    system_cache: dict[tuple[int, int], dict] = {}
    for first_edge, second_edge in boundary_pairs:
        left_edge, right_edge = sorted((first_edge, second_edge))
        pair = (left_edge, right_edge)
        data = system_cache.get(pair)
        if data is None:
            data = _system_data(polygon, left_edge, right_edge)
            system_cache[pair] = data
            systems.append({"edges": [left_edge, right_edge],
                            **system_summary(data)})
            for box in systems[-1]["paired_root_boxes"]:
                s_mid = (float(sp.Rational(box["s_interval"][0])) +
                         float(sp.Rational(box["s_interval"][1]))) / 2
                t_mid = (float(sp.Rational(box["t_interval"][0])) +
                         float(sp.Rational(box["t_interval"][1]))) / 2
                approx_parameters = (
                    _edge_parameter(guide.approximate_endpoints[0], polygon, first_edge),
                    _edge_parameter(guide.approximate_endpoints[1], polygon, second_edge),
                )
                parameters = (s_mid, t_mid)
                if first_edge > second_edge:
                    parameters = parameters[::-1]
                if all(abs(a - b) <= args.branch_window
                       for a, b in zip(parameters, approx_parameters)):
                    selected_root_boxes.append({
                        "edges": [left_edge, right_edge],
                        "s_interval": box["s_interval"],
                        "t_interval": box["t_interval"],
                        "s_approximation": s_mid,
                        "t_approximation": t_mid,
                        "endpoint_coordinate_intervals": [
                            _affine_interval(polygon, left_edge, box["s_interval"]),
                            _affine_interval(polygon, right_edge, box["t_interval"]),
                        ],
                        "endpoint_approximate_coordinates": [
                            [_interval_midpoint(interval) for interval in coordinates]
                            for coordinates in [
                                _affine_interval(polygon, left_edge, box["s_interval"]),
                                _affine_interval(polygon, right_edge, box["t_interval"]),
                            ]
                        ],
                        "approximation_distance": [
                            abs(parameters[0] - approx_parameters[0]),
                            abs(parameters[1] - approx_parameters[1]),
                        ],
                    })
        raw_solutions = solve_exact(polygon, left_edge, right_edge, data)
        # The crease fit selects the relevant branch when the exact system
        # has additional roots elsewhere on the boundary.  The final
        # coordinates remain exact; this window is used only for branch
        # selection from the finite-element guide.
        approx_parameters = (
            _edge_parameter(guide.approximate_endpoints[0], polygon, first_edge),
            _edge_parameter(guide.approximate_endpoints[1], polygon, second_edge),
        )
        for solution in raw_solutions:
            parameters = (float(solution["s"]), float(solution["t"]))
            if first_edge > second_edge:
                parameters = parameters[::-1]
            if any(abs(a - b) > args.branch_window
                   for a, b in zip(parameters, approx_parameters)):
                continue
            all_solutions.append({
                "edges": [left_edge, right_edge],
                "exact_endpoints": [
                    _exact_endpoint(polygon, left_edge, solution["s"]),
                    _exact_endpoint(polygon, right_edge, solution["t"]),
                ],
                **solution,
            })
    output = {
        "crease_edge_count": len(guide.edges),
        "crease_threshold": guide.threshold,
        "boundary_edge_pairs": [list(pair) for pair in boundary_pairs],
        "approximate_endpoints": guide.approximate_endpoints,
        "exact_systems": systems,
        "selected_algebraic_root_boxes": selected_root_boxes,
        "exact_solutions": [{key: _json_value(value) for key, value in solution.items()}
                            for solution in all_solutions],
    }
    # Screen every exact algebraic root box.  The FE branch window only
    # selects the displayed branch; it must not suppress explicit candidates.
    screening = []
    for system in systems:
        edges = tuple(system["edges"])
        for box in system.get("paired_root_boxes", []):
            item = _scan_root_box(
                polygon, edges, box, args.scan_theta_steps, args.scan_t_steps,
                args.scan_precision, args.scan_refine)
            item["near_crease_branch"] = any(
                selected.get("edges") == list(edges) and
                selected.get("s_interval") == box.get("s_interval") and
                selected.get("t_interval") == box.get("t_interval")
                for selected in selected_root_boxes)
            screening.append(item)
    output["candidate_screening"] = screening
    output["accepted_candidate_screening"] = [
        item for item in screening
        if item["accepted"] or args.keep_inconclusive]
    # Only candidates that passed both geometric concavity and both relative-DF
    # scans are drawn.  Inconclusive candidates are retained only on request,
    # but are not represented as accepted creases.
    drawable = [item for item in screening if item["accepted"]]
    singular_flags = _vertex_singularity_flags(args.polygon, polygon.edge_count)
    args.svg_output.parent.mkdir(parents=True, exist_ok=True)
    _write_crease_svg(args.svg_output, polygon, drawable, singular_flags)
    output["svg_output"] = str(args.svg_output)
    output["svg_candidate_count"] = len(drawable)
    if args.show_crease_edges:
        output["crease_edges"] = [list(edge) for edge in guide.edges]
    encoded = json.dumps(output, indent=2)
    args.output.write_text(encoded + "\n")
    report_path = args.output.with_suffix(".md")
    detail_files = {}
    for system in systems:
        first_edge, second_edge = system["edges"]
        detail_name = f"{args.output.stem}_edges_{first_edge}_{second_edge}.txt"
        detail_path = args.output.with_name(detail_name)
        _write_system_details(detail_path, system)
        detail_files[(first_edge, second_edge)] = detail_name
    report_path.write_text(_markdown_report(output, detail_files) + "\n")
    print(f"JSON:   {args.output}")
    print(f"Report: {report_path}")
    print(f"SVG:    {args.svg_output}")


if __name__ == "__main__":
    main()
