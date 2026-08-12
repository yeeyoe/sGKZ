#!/usr/bin/env python3
"""Render the CSV files emitted by shortest_gkz --plot-prefix.

The C++ program remains independent of a plotting library. Install the small
Python plotting stack once with ``python3 -m pip install numpy matplotlib`` if
static PNG export is also desired. Interactive HTML export only needs Python's
standard library and a browser with access to the Plotly CDN.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import tempfile
from pathlib import Path


def read_surface(path: Path):
    rows = []
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            rows.append(row)
    if not rows:
        raise ValueError(f"No rows in {path}")
    x = [float(row["x"]) for row in rows]
    y = [float(row["y"]) for row in rows]
    sigma_vee = [float(row["sigma_vee"]) for row in rows]
    psi = [row["psi"] for row in rows]
    psi_values = (
        None if all(value == "" for value in psi)
        else [float(value) for value in psi]
    )
    return x, y, sigma_vee, psi_values


def read_triangles(path: Path):
    with path.open(newline="") as stream:
        return [(int(row["i"]), int(row["j"]), int(row["l"]))
                for row in csv.DictReader(stream)]


def read_subdivision(path: Path):
    cells = {}
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            cells.setdefault(int(row["cell"]), []).append(
                (float(row["x"]), float(row["y"]))
            )
    return [cells[index] for index in sorted(cells)]


def point_key(x, y):
    """Use a stable key for coordinates shared by the two CSV files."""
    return (round(x, 12), round(y, 12))


def subdivision_edges(cells, x, y, z):
    """Return lifted boundary lines for every subdivision cell."""
    heights = {
        point_key(point_x, point_y): value
        for point_x, point_y, value in zip(x, y, z)
    }
    edges = []
    for cell_index, cell in enumerate(cells):
        if len(cell) < 3:
            raise ValueError(
                f"Subdivision cell {cell_index} has fewer than 3 vertices")
        closed = cell + [cell[0]]
        line_z = []
        for point_x, point_y in closed:
            key = point_key(point_x, point_y)
            if key not in heights:
                raise ValueError(
                    f"Subdivision vertex ({point_x}, {point_y}) is absent "
                    "from surface")
            line_z.append(heights[key])
        edges.append({
            "x": [point[0] for point in closed],
            "y": [point[1] for point in closed],
            "z": line_z,
            "cell": cell_index,
        })
    return edges


def data_spans(x, y, z):
    """Return positive coordinate spans for a stable 3D display ratio."""
    return tuple(max(max(values) - min(values), 1e-12)
                 for values in (x, y, z))


def draw_interactive(prefix: Path, x, y, z, triangles, edge_lines,
                     suffix: str, title: str, z_scale: float,
                     show_color: bool):
    """Write a rotatable Plotly Mesh3d view without a Python Plotly dependency."""
    span_x, span_y, span_z = data_spans(x, y, z)
    aspect_ratio = json.dumps(
        {"x": span_x, "y": span_y, "z": span_z * z_scale},
        separators=(",", ":"),
    )
    payload = {
        "x": x,
        "y": y,
        "z": z,
        "i": [face[0] for face in triangles],
        "j": [face[1] for face in triangles],
        "k": [face[2] for face in triangles],
        "edges": edge_lines,
    }
    data = json.dumps(payload, separators=(",", ":"))
    title_json = json.dumps(title)
    edge_color = ("edgeColors[index % edgeColors.length]" if show_color
                  else json.dumps("#991b1b"))
    edge_script = f"""
    const edgeColors = ["#111827", "#b91c1c", "#047857", "#7c3aed",
                        "#c2410c", "#0e7490", "#4338ca", "#4d7c0f"];
    surface.edges.forEach((edge, index) => traces.push({{
      type: "scatter3d",
      mode: "lines",
      x: edge.x, y: edge.y, z: edge.z,
      line: {{color: {edge_color}, width: 6}},
      hoverinfo: "skip",
      showlegend: false
    }}));
""" if edge_lines else ""
    surface_style = ("      intensity: surface.z,\n"
                     "      colorscale: \"Viridis\",\n"
                     if show_color else
                     "      color: \"#9ecae1\",\n"
                     "      opacity: 0.95,\n")
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{title}</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    html, body, #plot {{ width: 100%; height: 100%; margin: 0; }}
    body {{ overflow: hidden; }}
  </style>
</head>
<body>
  <div id="plot" aria-label="Interactive 3D surface"></div>
  <script>
    const surface = {data};
    const surfaceTrace = {{
      type: "mesh3d",
      x: surface.x, y: surface.y, z: surface.z,
      i: surface.i, j: surface.j, k: surface.k,
{surface_style}
      flatshading: true,
      hovertemplate: "x=%{{x}}<br>y=%{{y}}<br>value=%{{z}}<extra></extra>"
    }};
    const traces = [surfaceTrace];
{edge_script}
    const layout = {{
      title: {title_json},
      scene: {{
        xaxis: {{title: "x"}},
        yaxis: {{title: "y"}},
        zaxis: {{title: "value"}},
        aspectmode: "manual",
        aspectratio: {aspect_ratio}
      }},
      margin: {{l: 0, r: 0, b: 0, t: 48}}
    }};
    Plotly.newPlot("plot", traces, layout, {{responsive: true, displaylogo: false}});
  </script>
</body>
</html>
"""
    output = prefix.parent / f"{prefix.name}_{suffix}.html"
    output.write_text(html, encoding="utf-8")
    return output


def draw_static(prefix: Path, x, y, z, triangles, edge_lines, suffix: str,
                title: str, z_scale: float, show: bool,
                show_color: bool):
    """Optionally export the same surface as a fixed PNG."""
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401, registers projection

    figure = plt.figure(figsize=(8, 6))
    axes = figure.add_subplot(111, projection="3d")
    span_x, span_y, span_z = data_spans(x, y, z)
    surface_style = {"cmap": "viridis"} if show_color else {"color": "#9ecae1"}
    axes.plot_trisurf(x, y, z, triangles=triangles, **surface_style,
                      linewidth=0.25, edgecolor="0.35", antialiased=True)
    edge_colors = (["#111827", "#b91c1c", "#047857", "#7c3aed",
                    "#c2410c", "#0e7490", "#4338ca", "#4d7c0f"]
                   if show_color else ["#991b1b"])
    for index, edge in enumerate(edge_lines):
        axes.plot(edge["x"], edge["y"], edge["z"],
                  color=edge_colors[index % len(edge_colors)], linewidth=2.2)
    # Keep all data coordinates and labels unchanged; only the rendered box is
    # made taller in z when requested.
    axes.set_box_aspect((span_x, span_y, span_z * z_scale))
    axes.set_xlabel("x")
    axes.set_ylabel("y")
    axes.set_zlabel("value")
    axes.zaxis.set_major_locator(MaxNLocator(nbins=5))
    axes.set_title(title)
    figure.tight_layout()
    output = prefix.parent / f"{prefix.name}_{suffix}.png"
    figure.savefig(output, dpi=180)
    if not show:
        plt.close(figure)
    return output


def draw_subdivision_svg(prefix: Path, cells, title: str):
    """Write a dependency-free 2D image of the subdivision."""
    width, height, margin = 1000, 760, 64
    points = [point for cell in cells for point in cell]
    min_x = min(point[0] for point in points)
    max_x = max(point[0] for point in points)
    min_y = min(point[1] for point in points)
    max_y = max(point[1] for point in points)
    span_x = max(max_x - min_x, 1e-12)
    span_y = max(max_y - min_y, 1e-12)
    scale = min((width - 2 * margin) / span_x,
                (height - 2 * margin) / span_y)

    def project(point):
        x, y = point
        return (margin + (x - min_x) * scale,
                height - margin - (y - min_y) * scale)

    palette = ["#2563eb", "#16a34a", "#d97706", "#9333ea", "#dc2626",
               "#0891b2", "#4f46e5", "#65a30d"]
    shapes = [
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="0 0 {width} {height}" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="64" y="34" font-family="sans-serif" font-size="24">'
        f'{title}</text>',
    ]
    for index, cell in enumerate(cells):
        polygon = " ".join(f"{x:.8f},{y:.8f}" for x, y in map(project, cell))
        shapes.append(
            f'<polygon points="{polygon}" fill="{palette[index % len(palette)]}" '
            'fill-opacity="0.28" stroke="#202020" stroke-width="1.4"/>')
    shapes.append('</svg>')
    output = prefix.parent / f"{prefix.name}_subdivision.svg"
    output.write_text("\n".join(shapes), encoding="utf-8")
    return output


def draw_subdivision_static(prefix: Path, cells, title: str, show: bool):
    """Optionally export the 2D subdivision as a PNG."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    figure, axes = plt.subplots(figsize=(8, 6))
    palette = ["#2563eb", "#16a34a", "#d97706", "#9333ea", "#dc2626",
               "#0891b2", "#4f46e5", "#65a30d"]
    for index, cell in enumerate(cells):
        axes.add_patch(Polygon(cell, closed=True,
                               facecolor=palette[index % len(palette)],
                               edgecolor="#202020", alpha=0.28,
                               linewidth=1.2))
    points = [point for cell in cells for point in cell]
    axes.scatter([point[0] for point in points], [point[1] for point in points],
                 s=8, color="#202020", zorder=3)
    axes.set_aspect("equal", adjustable="box")
    axes.set_xlabel("x")
    axes.set_ylabel("y")
    axes.set_title(title)
    figure.tight_layout()
    output = prefix.parent / f"{prefix.name}_subdivision.png"
    figure.savefig(output, dpi=180)
    if not show:
        plt.close(figure)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prefix", type=Path,
                        help="Prefix passed to shortest_gkz --plot-prefix")
    parser.add_argument("--sigma-vee", action="store_true",
                        help="Generate the sigma_A^vee view; off by default in polygon mode")
    parser.add_argument("--no-sigma-vee", action="store_true",
                        help="Do not generate the sigma_A^vee 3D view")
    parser.add_argument("--no-psi", action="store_true",
                        help="Do not generate the psi_k 3D view")
    parser.add_argument("--no-subdivision", action="store_true",
                        help="Do not generate the S(sigma_A) 2D image")
    parser.add_argument("--no-edges", action="store_true",
                        help="Do not overlay lifted subdivision edges")
    parser.add_argument("--z-scale", type=float, default=1.0,
                        help="Stretch the displayed z direction")
    parser.add_argument("--color", action="store_true",
                        help="Use colored 3D surfaces and lifted edges")
    parser.add_argument("--png", action="store_true",
                        help="Also export fixed-angle PNG files")
    parser.add_argument("--show", action="store_true",
                        help="Display the optional PNG figures after saving them")
    args = parser.parse_args()
    if not math.isfinite(args.z_scale) or args.z_scale <= 0:
        parser.error("--z-scale must be a finite positive number")

    x, y, sigma_vee, psi = read_surface(
        args.prefix.parent / f"{args.prefix.name}_surface.csv")
    if args.sigma_vee and args.no_sigma_vee:
        parser.error("--sigma-vee and --no-sigma-vee cannot be combined")
    show_sigma_vee = (psi is None or args.sigma_vee) and not args.no_sigma_vee
    needs_triangles = show_sigma_vee or (psi is not None and not args.no_psi)
    triangles = (read_triangles(
        args.prefix.parent / f"{args.prefix.name}_triangles.csv")
                 if needs_triangles else [])
    cells = None
    if not args.no_subdivision or not args.no_edges:
        cells = read_subdivision(
            args.prefix.parent / f"{args.prefix.name}_subdivision.csv")
    sigma_edges = (subdivision_edges(cells, x, y, sigma_vee)
                   if show_sigma_vee and cells is not None and not args.no_edges
                   else [])
    psi_edges = (subdivision_edges(cells, x, y, psi)
                 if cells is not None and psi is not None and not args.no_edges
                 else [])
    outputs = []
    z_suffix = ("" if args.z_scale == 1.0
                else f" (visual z x{args.z_scale:g})")
    stable_projection_plot = args.prefix.name.endswith("_stable")
    psi_label = ("stable projection" if stable_projection_plot
                 else "last QP solution")
    if show_sigma_vee:
        sigma_title = ("sigma_k^vee" if psi is not None else "sigma_A^vee") + z_suffix
        outputs.append(draw_interactive(args.prefix, x, y, sigma_vee, triangles,
                                        sigma_edges, "sigma_vee", sigma_title,
                                        args.z_scale, args.color))
    if psi is not None and not args.no_psi:
        outputs.append(draw_interactive(
            args.prefix, x, y, psi, triangles, psi_edges, "psi_k",
            f"psi_k ({psi_label})" + z_suffix, args.z_scale, args.color))
    if not args.no_subdivision:
        subdivision_title = "S(sigma_k)" if psi is not None else "S(sigma_A)"
        subdivision_title += f" ({psi_label})"
        outputs.append(draw_subdivision_svg(args.prefix, cells, subdivision_title))

    if args.png or args.show:
        try:
            cache_root = Path(tempfile.gettempdir()) / "shortest-gkz-matplotlib"
            cache_root.mkdir(parents=True, exist_ok=True)
            os.environ.setdefault("MPLCONFIGDIR", str(cache_root / "matplotlib"))
            os.environ.setdefault("XDG_CACHE_HOME", str(cache_root / "xdg"))
            import matplotlib
            if not args.show:
                matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError as error:
            parser.error("matplotlib is required for --png; install numpy and matplotlib")
            raise error
        if show_sigma_vee:
            sigma_title = (r"$\sigma_k^\vee$" if psi is not None
                           else r"$\sigma_A^\vee$") + z_suffix
            outputs.append(draw_static(args.prefix, x, y, sigma_vee, triangles,
                                       sigma_edges, "sigma_vee", sigma_title,
                                       args.z_scale, args.show, args.color))
        if psi is not None and not args.no_psi:
            outputs.append(draw_static(
                args.prefix, x, y, psi, triangles, psi_edges, "psi_k",
                rf"$\psi_k$ ({psi_label})" + z_suffix,
                args.z_scale, args.show, args.color))
        if not args.no_subdivision:
            subdivision_title = (r"$S(\sigma_k)$" if psi is not None
                                 else r"$S(\sigma_A)$")
            subdivision_title += f" ({psi_label})"
            outputs.append(draw_subdivision_static(args.prefix, cells,
                                                   subdivision_title, args.show))
    else:
        plt = None
    for output in outputs:
        print(output)
    if args.show and plt is not None:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
