"""Numerical relative-DF scan for polygons with algebraic vertices.

The geometry is evaluated at the midpoint of certified rational coordinate
intervals.  The scanner is deliberately a screening tool: it reports a
negative simple-convex relative DF witness, or that no such witness was found
on the requested grid.  Boundary measures are supplied explicitly so that a
new crease edge can have zero measure while inherited lattice edges retain
their gcd measure.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd, hypot, pi, cos, sin
from typing import Iterable, Sequence

import numpy as np


@dataclass(frozen=True)
class ScanResult:
    status: str
    minimum_relative_df: float | None
    minimum_normalized_relative_df: float | None
    witness: dict | None
    evaluations: int
    scan_precision: int
    reason: str | None = None

    def as_dict(self) -> dict:
        result = {
            "relative_df_scan_status": self.status,
            "minimum_relative_df": self.minimum_relative_df,
            "minimum_normalized_relative_df": self.minimum_normalized_relative_df,
            "witness": self.witness,
            "evaluations": self.evaluations,
            "scan_precision": self.scan_precision,
            "measure_policy": "inherited lattice edges=gcd(|dx|,|dy|); crease edge=0",
        }
        if self.reason:
            result["reason"] = self.reason
        return result


def _moments(vertices: np.ndarray) -> tuple[float, float, float, float, float, float]:
    area = ix = iy = ixx = ixy = iyy = 0.0
    for p, q in zip(vertices, np.roll(vertices, -1, axis=0)):
        x, y = p
        X, Y = q
        cross = x * Y - y * X
        area += cross
        ix += (x + X) * cross
        iy += (y + Y) * cross
        ixx += (x * x + x * X + X * X) * cross
        ixy += (2 * x * y + x * Y + X * y + 2 * X * Y) * cross
        iyy += (y * y + y * Y + Y * Y) * cross
    return area / 2, ix / 6, iy / 6, ixx / 12, ixy / 24, iyy / 12


def _ell(vertices: np.ndarray, measures: Sequence[float]) -> np.ndarray:
    area, ix, iy, ixx, ixy, iyy = _moments(vertices)
    length = bx = by = 0.0
    for p, q, measure in zip(vertices, np.roll(vertices, -1, axis=0), measures):
        length += measure
        bx += measure * (p[0] + q[0]) / 2
        by += measure * (p[1] + q[1]) / 2
    matrix = np.array([[area, ix, iy], [ix, ixx, ixy], [iy, ixy, iyy]], dtype=float)
    rhs = np.array([length, bx, by], dtype=float)
    try:
        return np.linalg.solve(matrix, rhs)
    except np.linalg.LinAlgError as exc:
        raise ValueError("singular moment matrix") from exc


def _clip(vertices: np.ndarray, a: float, b: float, c: float) -> np.ndarray:
    result: list[np.ndarray] = []
    if len(vertices) == 0:
        return np.empty((0, 2))
    values = [a * p[0] + b * p[1] + c for p in vertices]
    for i, p in enumerate(vertices):
        q = vertices[(i + 1) % len(vertices)]
        sp, sq = values[i], values[(i + 1) % len(vertices)]
        pin, qin = sp >= -1e-13, sq >= -1e-13
        if pin:
            result.append(p)
        if pin != qin and abs(sp - sq) > 1e-15:
            tau = sp / (sp - sq)
            result.append(p + tau * (q - p))
    if not result:
        return np.empty((0, 2))
    return np.asarray(result, dtype=float)


def _positive_boundary_integral(vertices: np.ndarray, measures: Sequence[float],
                                a: float, b: float, c: float) -> float:
    total = 0.0
    for p, q, measure in zip(vertices, np.roll(vertices, -1, axis=0), measures):
        sp = a * p[0] + b * p[1] + c
        sq = a * q[0] + b * q[1] + c
        if sp <= 0 and sq <= 0:
            continue
        lo, hi = 0.0, 1.0
        if sp < 0:
            lo = sp / (sp - sq)
        elif sq < 0:
            hi = sp / (sp - sq)
        g0 = sp + lo * (sq - sp)
        g1 = sp + hi * (sq - sp)
        total += measure * (hi - lo) * (g0 + g1) / 2
    return total


def _positive_area_integral(vertices: np.ndarray, ell: np.ndarray,
                            a: float, b: float, c: float) -> float:
    clipped = _clip(vertices, a, b, c)
    if len(clipped) < 3:
        return 0.0
    area, ix, iy, ixx, ixy, iyy = _moments(clipped)
    # Integral of (a x+b y+c) ell(x,y) over the clipped polygon.
    e0, e1, e2 = ell
    return (a * e0 * ix + b * e0 * iy + c * e0 * area +
            a * e1 * ixx + (a * e2 + b * e1) * ixy +
            b * e2 * iyy + c * e1 * ix + c * e2 * iy)


def relative_df(vertices: np.ndarray, measures: Sequence[float], ell: np.ndarray,
                a: float, b: float, c: float) -> float:
    return _positive_boundary_integral(vertices, measures, a, b, c) - \
        _positive_area_integral(vertices, ell, a, b, c)


def scan_polygon(vertices: Sequence[Sequence[float]], measures: Sequence[float],
                 theta_steps: int = 180, t_steps: int = 128,
                 precision: int = 40, refine: bool = False) -> ScanResult:
    if precision < 16:
        return ScanResult("inconclusive", None, None, None, 0, precision,
                          "scan precision must be at least 16 digits")
    vertices = np.asarray(vertices, dtype=float)
    measures = list(float(x) for x in measures)
    if len(vertices) < 3 or len(measures) != len(vertices):
        return ScanResult("inconclusive", None, None, None, 0, precision,
                          "invalid polygon or measure list")
    scale = max(1.0, float(np.max(np.abs(vertices))))
    # Double precision is used for the geometric quadrature.  Very wide
    # algebraic boxes are therefore reported as inconclusive rather than
    # being mistaken for a reliable stability decision.
    if not np.all(np.isfinite(vertices)) or scale > 1e12:
        return ScanResult("inconclusive", None, None, None, 0, precision,
                          "non-finite or excessively large algebraic coordinates")
    try:
        ell = _ell(vertices, measures)
    except ValueError as exc:
        return ScanResult("inconclusive", None, None, None, 0, precision, str(exc))
    boundary = sum(measures)
    diameter = max(hypot(*(p - q)) for p in vertices for q in vertices)
    sup_ell = max(abs(ell[0] + ell[1] * p[0] + ell[2] * p[1]) for p in vertices)
    normalization = max(boundary * max(sup_ell, 1e-300) * diameter, 1e-300)
    best = float("inf")
    best_data = None
    evaluations = 0
    for k in range(max(1, theta_steps)):
        theta = 2 * pi * k / max(1, theta_steps)
        a, b = cos(theta), sin(theta)
        projections = vertices @ np.array([a, b])
        lo, hi = float(np.min(projections)), float(np.max(projections))
        for j in range(max(1, t_steps) + 1):
            t = lo + (hi - lo) * j / max(1, t_steps)
            value = relative_df(vertices, measures, ell, a, b, -t)
            evaluations += 1
            if value < best:
                best = value
                best_data = {"a": a, "b": b, "c": -t, "theta": theta,
                             "offset": t, "relative_df": value}
    if refine and best_data is not None:
        # Small deterministic local search around the best grid point.
        for dtheta in (-2, -1, 1, 2):
            theta = best_data["theta"] + 2 * pi * dtheta / max(1, theta_steps)
            a, b = cos(theta), sin(theta)
            projections = vertices @ np.array([a, b])
            lo, hi = float(np.min(projections)), float(np.max(projections))
            for j in range(9):
                t = lo + (hi - lo) * j / 8
                value = relative_df(vertices, measures, ell, a, b, -t)
                evaluations += 1
                if value < best:
                    best = value
                    best_data = {"a": a, "b": b, "c": -t,
                                 "theta": theta, "offset": t,
                                 "relative_df": value}
    normalized = best / normalization if best_data is not None else None
    status = "unstable" if normalized is not None and normalized < -1e-7 \
        else "no_counterexample_found"
    return ScanResult(status, best if best_data is not None else None, normalized,
                      best_data, evaluations, precision)


def lattice_measure(p: Sequence[float], q: Sequence[float]) -> float:
    """Return gcd edge measure for integer inherited vertices."""
    return float(gcd(abs(round(q[0] - p[0])), abs(round(q[1] - p[1]))))
