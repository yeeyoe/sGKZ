#!/usr/bin/env python3
"""Render one stored K-stability candidate as a self-contained SVG.

The script intentionally uses only Python's standard library.  It reads the
candidate geometry and an exact witness from the search SQLite database,
draws the polygon and the two affine zero lines, and annotates every vertex.
"""

from __future__ import annotations

import argparse
import html
import math
import sqlite3
from fractions import Fraction
from pathlib import Path


def parse_pairs(text: str):
    result = []
    for item in text.split(";"):
        if not item:
            continue
        try:
            x, y = item.split(":", 1)
            result.append((int(x), int(y)))
        except ValueError as error:
            raise ValueError(f"invalid coordinate list: {text!r}") from error
    return result


def parse_flags(text: str | None, directions):
    if text:
        flags = [int(item) for item in text.split(",") if item != ""]
        if len(flags) == len(directions) and all(flag in (0, 1) for flag in flags):
            return flags
    flags = []
    for index, direction in enumerate(directions):
        previous = directions[(index - 1) % len(directions)]
        determinant = previous[0] * direction[1] - previous[1] * direction[0]
        flags.append(0 if abs(determinant) == 1 else 1)
    return flags


def fraction(text: str) -> Fraction:
    return Fraction(text.strip())


def rational_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def line_value(line, point):
    a, b, c = line
    return a * point[0] + b * point[1] + c


def line_box_intersections(line, bounds):
    """Return the two intersections of a line with an axis-aligned box."""
    a, b, c = (float(value) for value in line)
    min_x, max_x, min_y, max_y = bounds
    points = []
    epsilon = 1e-12

    def add(x, y):
        if (min_x - epsilon <= x <= max_x + epsilon and
                min_y - epsilon <= y <= max_y + epsilon):
            point = (max(min_x, min(max_x, x)), max(min_y, min(max_y, y)))
            if all(math.hypot(point[0] - old[0], point[1] - old[1]) > 1e-9
                   for old in points):
                points.append(point)

    if abs(b) > epsilon:
        add(min_x, -(a * min_x + c) / b)
        add(max_x, -(a * max_x + c) / b)
    if abs(a) > epsilon:
        add(-(b * min_y + c) / a, min_y)
        add(-(b * max_y + c) / a, max_y)
    if len(points) < 2:
        return []
    return points[:2]


def witness_from_row(row):
    exact = row[0:3]
    if all(value not in (None, "") for value in exact):
        return tuple(fraction(value) for value in exact), "exact"
    numeric = row[3:6]
    if all(value is not None for value in numeric):
        ux, uy, t = (float(value) for value in numeric)
        return (ux, uy, -t), "numeric"
    return None, None


def fetch_candidate(database: Path, key: str):
    connection = sqlite3.connect(database)
    connection.row_factory = sqlite3.Row
    try:
        candidate = connection.execute(
            "SELECT key,d,directions,steps,normals,vertices,twice_area,ell0,ell1,ell2,status,"
            "vertex_singularity_flags,singular_vertex_count "
            "FROM candidates WHERE key=?", (key,)
        ).fetchone()
        if candidate is None:
            raise ValueError(f"candidate key not found: {key}")
        attempts = connection.execute(
            "SELECT exact_a,exact_b,exact_c,witness_ux,witness_uy,witness_t "
            "FROM attempts WHERE candidate_key=? AND status='verified_unstable' "
            "ORDER BY rowid DESC", (key,)
        ).fetchall()
        witness = None
        witness_kind = None
        for attempt in attempts:
            witness, witness_kind = witness_from_row(tuple(attempt))
            if witness is not None:
                break
        if witness is None:
            raise ValueError("no verified exact or numerical witness is stored for this key")
        return candidate, witness, witness_kind
    finally:
        connection.close()


def svg_text(x, y, text, color="#111827", size=14, anchor="start", weight="normal"):
    escaped = html.escape(str(text))
    return (f'<text x="{x:.2f}" y="{y:.2f}" fill="{color}" '
            f'font-family="Arial, sans-serif" font-size="{size}px" '
            f'font-weight="{weight}" text-anchor="{anchor}" '
            'paint-order="stroke" stroke="white" stroke-width="4" '
            f'stroke-opacity="0.9">{escaped}</text>\n')


def render(candidate, witness, witness_kind, output: Path):
    key = candidate["key"]
    vertices = parse_pairs(candidate["vertices"])
    directions = parse_pairs(candidate["directions"])
    if len(vertices) < 3 or len(vertices) != candidate["d"]:
        raise ValueError("candidate has an invalid vertex list")
    if len(directions) != len(vertices):
        raise ValueError("candidate has inconsistent direction and vertex counts")
    flags = parse_flags(candidate["vertex_singularity_flags"], directions)
    ell = tuple(fraction(candidate[name]) for name in ("ell0", "ell1", "ell2"))

    xs = [point[0] for point in vertices]
    ys = [point[1] for point in vertices]
    data_width = max(max(xs) - min(xs), 1)
    data_height = max(max(ys) - min(ys), 1)
    data_pad = max(data_width, data_height) * 0.18
    bounds = (min(xs) - data_pad, max(xs) + data_pad,
              min(ys) - data_pad, max(ys) + data_pad)

    width, height = 1120, 820
    left, right, top, bottom = 90, 50, 70, 130
    plot_width = width - left - right
    plot_height = height - top - bottom
    min_x, max_x, min_y, max_y = bounds
    scale = min(plot_width / (max_x - min_x), plot_height / (max_y - min_y))

    def screen(point):
        return (left + (point[0] - min_x) * scale,
                height - bottom - (point[1] - min_y) * scale)

    polygon_points = " ".join(f"{screen(point)[0]:.2f},{screen(point)[1]:.2f}"
                              for point in vertices)
    ell_line = (ell[1], ell[2], ell[0])
    ell_segment = line_box_intersections(ell_line, bounds)
    witness_segment = line_box_intersections(witness, bounds)
    if not ell_segment:
        raise ValueError("ell_P=0 does not intersect the plotting box")
    if not witness_segment:
        raise ValueError("witness crease does not intersect the plotting box")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as stream:
        stream.write(f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}">
<rect width="100%" height="100%" fill="#ffffff"/>
''')
        stream.write(svg_text(width / 2, 30, f"K-stability candidate d={candidate['d']}",
                              size=20, anchor="middle", weight="bold"))
        stream.write(svg_text(width / 2, 52, key, color="#4b5563", size=10,
                              anchor="middle"))
        stream.write('<polygon points="' + polygon_points + '" fill="#eff6ff" '
                     'stroke="#111827" stroke-width="2"/>\n')

        ell_a, ell_b, ell_c = (float(value) for value in ell_line)
        ex1, ey1 = screen(ell_segment[0])
        ex2, ey2 = screen(ell_segment[1])
        stream.write(f'<line x1="{ex1:.2f}" y1="{ey1:.2f}" x2="{ex2:.2f}" y2="{ey2:.2f}" '
                     'stroke="#2563eb" stroke-width="2" stroke-dasharray="8 5"/>\n')

        wx1, wy1 = screen(witness_segment[0])
        wx2, wy2 = screen(witness_segment[1])
        stream.write(f'<line x1="{wx1:.2f}" y1="{wy1:.2f}" x2="{wx2:.2f}" y2="{wy2:.2f}" '
                     'stroke="#dc2626" stroke-width="3"/>\n')

        for index, point in enumerate(vertices):
            x, y = screen(point)
            singular = flags[index] == 1
            color = "#7c3aed" if singular else "#111827"
            fill = "#a855f7" if singular else "#ffffff"
            radius = 6 if singular else 4
            stream.write(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{radius}" '
                         f'fill="{fill}" stroke="{color}" stroke-width="2"/>\n')
            dx = 10 if x < width * 0.72 else -10
            anchor = "start" if dx > 0 else "end"
            stream.write(svg_text(x + dx, y - 10, f"v{index + 1}=({point[0]},{point[1]})",
                                  color=color, size=13, anchor=anchor,
                                  weight="bold" if singular else "normal"))

        legend_y = height - 92
        stream.write('<line x1="90" y1="' + str(legend_y) +
                     '" x2="122" y2="' + str(legend_y) +
                     '" stroke="#2563eb" stroke-width="2" stroke-dasharray="8 5"/>\n')
        stream.write(svg_text(130, legend_y + 5, "ell_P=0", size=13))
        stream.write('<line x1="250" y1="' + str(legend_y) +
                     '" x2="282" y2="' + str(legend_y) +
                     '" stroke="#dc2626" stroke-width="3"/>\n')
        stream.write(svg_text(290, legend_y + 5, "witness crease", color="#dc2626", size=13))
        stream.write('<circle cx="490" cy="' + str(legend_y) +
                     '" r="6" fill="#a855f7" stroke="#7c3aed" stroke-width="2"/>\n')
        stream.write(svg_text(505, legend_y + 5, "singular vertex", color="#7c3aed", size=13))
        stream.write(svg_text(90, height - 58,
                              f"ell_P(x,y) = {rational_text(ell[0])} + "
                              f"({rational_text(ell[1])})x + ({rational_text(ell[2])})y",
                              size=12))
        stream.write(svg_text(90, height - 34,
                              "witness: " + "a x + b y + c = 0; " + witness_kind,
                              color="#dc2626", size=12))
        stream.write('</svg>\n')


def line_text(line, kind):
    if kind == "exact":
        return (f"({rational_text(line[0])}) x + ({rational_text(line[1])}) y + "
                f"({rational_text(line[2])})")
    return (f"({float(line[0]):.17g}) x + ({float(line[1]):.17g}) y + "
            f"({float(line[2]):.17g})")


def write_vertex_file(candidate, witness, witness_kind, output: Path):
    flags = parse_flags(candidate["vertex_singularity_flags"],
                        parse_pairs(candidate["directions"]))
    ell = tuple(fraction(candidate[name]) for name in ("ell0", "ell1", "ell2"))
    ell_formula = (f"{rational_text(ell[0])} + ({rational_text(ell[1])}) x + "
                   f"({rational_text(ell[2])}) y")
    vertices = parse_pairs(candidate["vertices"])
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as stream:
        stream.write(f"# candidate_key={candidate['key']}\n")
        stream.write(f"# d={candidate['d']}\n")
        stream.write(f"# twice_area={candidate['twice_area']}\n")
        stream.write(f"# status={candidate['status']}\n")
        stream.write(f"# directions={candidate['directions']}\n")
        stream.write(f"# steps={candidate['steps']}\n")
        stream.write(f"# facet_normals={candidate['normals']}\n")
        stream.write(f"# ell_P(x,y)={ell_formula}\n")
        stream.write(f"# vertex_singularity_flags={','.join(map(str, flags))}\n")
        stream.write(f"# singular_vertex_count={sum(flags)}\n")
        stream.write(f"# witness_kind={witness_kind}\n")
        stream.write(f"# witness_line={line_text(witness, witness_kind)} = 0\n")
        for x, y in vertices:
            stream.write(f"{x} {y}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Render a stored candidate, ell_P=0, and its witness crease as SVG."
    )
    parser.add_argument("key", nargs="?", help="candidate key")
    parser.add_argument("--key", dest="key_option", help="candidate key")
    parser.add_argument(
        "--database",
        type=Path,
        default=Path(__file__).with_name("k_stability_search.sqlite"),
        help="search SQLite database (default: K-stability/k_stability_search.sqlite)",
    )
    parser.add_argument("--output", type=Path, help="output SVG path (required without --save)")
    parser.add_argument(
        "--save", action="store_true",
        help="save SVG and the commented vertex file under the project examples directory",
    )
    arguments = parser.parse_args()
    key = arguments.key or arguments.key_option
    if not key or (arguments.key and arguments.key_option):
        parser.error("provide the candidate key once, either positionally or with --key")
    if arguments.save and arguments.output is not None:
        parser.error("--save and --output are mutually exclusive")
    if not arguments.save and arguments.output is None:
        parser.error("--output is required unless --save is used")

    candidate, witness, witness_kind = fetch_candidate(arguments.database, key)
    if arguments.save:
        examples = Path(__file__).resolve().parent.parent / "examples"
        basename = f"d{candidate['d']}_a{candidate['twice_area']}"
        output = examples / f"{basename}.svg"
        vertex_output = examples / basename
    else:
        output = arguments.output
        vertex_output = None
    if output.exists():
        raise ValueError(f"output file already exists: {output}")
    if vertex_output is not None and vertex_output.exists():
        raise ValueError(f"vertex file already exists: {vertex_output}")
    render(candidate, witness, witness_kind, output)
    if vertex_output is not None:
        write_vertex_file(candidate, witness, witness_kind, vertex_output)
        print(f"vertices_file={vertex_output}")
    print(f"svg={output}")
    print(f"candidate_key={key}")
    print(f"singular_vertex_count={sum(parse_flags(candidate['vertex_singularity_flags'], parse_pairs(candidate['directions'])))}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, sqlite3.Error, ValueError) as error:
        raise SystemExit(f"error: {error}")
