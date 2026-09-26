#!/usr/bin/env python3
"""Plot active-set size and gap from a shortest_gkz verbose log."""

from __future__ import annotations

import argparse
import html
import json
import math
import os
import re
import tempfile
from pathlib import Path


FLOAT_PATTERN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
ITERATION_PATTERN = re.compile(
    rf"^iteration=(?P<iteration>\d+)\s+"
    rf"active=(?P<active>\d+)\s+"
    rf"norm2=(?P<norm2>{FLOAT_PATTERN})\s+"
    rf"gap=(?P<gap>{FLOAT_PATTERN})$"
)
PROJECTION_STABLE_PATTERN = re.compile(
    r"^projection\s+event=stable\s+iteration=(?P<iteration>\d+)\b"
)
PROJECTION_RANK_PATTERN = re.compile(
    r"^projection\s+event=rank\s+iteration=(?P<iteration>\d+)\s+"
    r"rank=(?P<rank>\d+)\b"
)
PARAMETERS_PATTERN = re.compile(
    rf"^parameters\s+tolerance=(?P<tolerance>{FLOAT_PATTERN})\s+"
    rf"absolute_tolerance=(?P<absolute_tolerance>{FLOAT_PATTERN})$"
)


def read_solver_tolerances(path: Path):
    """Read tolerances recorded by a verbose shortest_gkz run, if present."""
    with path.open(encoding="utf-8") as stream:
        for raw_line in stream:
            match = PARAMETERS_PATTERN.fullmatch(raw_line.strip())
            if match is None:
                continue
            tolerance = float(match.group("tolerance"))
            absolute_tolerance = float(match.group("absolute_tolerance"))
            if (not math.isfinite(tolerance) or tolerance < 0 or
                    not math.isfinite(absolute_tolerance) or
                    absolute_tolerance < 0):
                raise ValueError(f"Invalid tolerances recorded in {path}")
            return tolerance, absolute_tolerance
    return None


def read_history(path: Path):
    iterations = []
    active_sizes = []
    norm2_values = []
    gaps = []
    stable_projections = []
    rank_iterations = []
    ranks = []
    with path.open(encoding="utf-8") as stream:
        for line_number, raw_line in enumerate(stream, start=1):
            line = raw_line.strip()
            stable_match = PROJECTION_STABLE_PATTERN.match(line)
            if stable_match is not None:
                stable_projections.append(int(stable_match.group("iteration")))
                continue
            rank_match = PROJECTION_RANK_PATTERN.match(line)
            if rank_match is not None:
                rank_iterations.append(int(rank_match.group("iteration")))
                ranks.append(int(rank_match.group("rank")))
                continue
            if not line or not line.startswith("iteration="):
                continue
            match = ITERATION_PATTERN.fullmatch(line)
            if match is None:
                raise ValueError(
                    f"Malformed iteration record at {path}:{line_number}"
                )
            iteration = int(match.group("iteration"))
            active = int(match.group("active"))
            norm2 = float(match.group("norm2"))
            gap = float(match.group("gap"))
            if iterations and iteration <= iterations[-1]:
                raise ValueError(
                    f"Iteration numbers must be strictly increasing at "
                    f"{path}:{line_number}; the log may contain multiple runs"
                )
            if active <= 0:
                raise ValueError(f"Active-set size must be positive at {path}:{line_number}")
            if not math.isfinite(norm2) or norm2 < 0:
                raise ValueError(f"norm2 must be finite and nonnegative at {path}:{line_number}")
            if not math.isfinite(gap) or gap < 0:
                raise ValueError(f"gap must be finite and nonnegative at {path}:{line_number}")
            iterations.append(iteration)
            active_sizes.append(active)
            norm2_values.append(norm2)
            gaps.append(gap)
    if not iterations:
        raise ValueError(f"No iteration records found in {path}")
    return (iterations, active_sizes, norm2_values, gaps, stable_projections,
            rank_iterations, ranks)


def default_prefix(log_path: Path) -> Path:
    base = log_path.with_suffix("")
    if base.name.endswith(".iterations"):
        return base
    return base.with_name(base.name + "_iterations")


def draw_interactive(prefix: Path, title: str, iterations, active_sizes,
                     gaps, thresholds, stable_projections, rank_iterations,
                     ranks) -> Path:
    payload = json.dumps(
        {
            "iteration": iterations,
            "active": active_sizes,
            "gap": gaps,
            "threshold": thresholds,
            "stable_projection": stable_projections,
            "rank_iteration": rank_iterations,
            "rank": ranks,
        },
        separators=(",", ":"),
    )
    title_json = json.dumps(title)
    document = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(title)}</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    html, body, #plot {{ width: 100%; height: 100%; margin: 0; }}
    body {{ overflow: hidden; }}
  </style>
</head>
<body>
  <div id="plot" aria-label="Iteration history"></div>
  <script>
    const history = {payload};
    const common = {{type: "scatter", mode: "lines"}};
    const traces = [
      {{...common, name: "active set size", x: history.iteration,
        y: history.active, xaxis: "x", yaxis: "y",
        line: {{color: "#2563eb", width: 2}}}},
      {{...common, name: "affine rank", x: history.rank_iteration,
        y: history.rank, xaxis: "x", yaxis: "y",
        line: {{color: "#dc2626", width: 2}}}},
      {{...common, name: "gap", x: history.iteration,
        y: history.gap, xaxis: "x2", yaxis: "y2",
        line: {{color: "#b91c1c", width: 2}}}},
      {{...common, name: "stopping threshold", x: history.iteration,
        y: history.threshold, xaxis: "x2", yaxis: "y2",
        line: {{color: "#2563eb", width: 2}},
        hovertemplate: "threshold=%{{y:.3e}}<extra></extra>"}}
    ];
    const projectionShapes = history.stable_projection.map((iteration) => ({{
      type: "line", xref: "x2", yref: "paper", x0: iteration,
      x1: iteration, y0: 0, y1: 1, line: {{color: "#059669", width: 1,
      dash: "dot"}}
    }}));
    const axisStyle = {{showgrid: true, gridcolor: "#e5e7eb",
                       zeroline: false, automargin: true}};
    const layout = {{
      title: {title_json},
      showlegend: false,
      hovermode: "x",
      plot_bgcolor: "#ffffff",
      paper_bgcolor: "#ffffff",
      xaxis: {{...axisStyle, domain: [0, 1], anchor: "y",
               matches: "x2", showticklabels: false}},
      yaxis: {{...axisStyle, domain: [0.55, 1.00], title: "active set size",
               rangemode: "tozero", tickformat: "d"}},
      xaxis2: {{...axisStyle, domain: [0, 1], anchor: "y2",
                title: "iteration"}},
      yaxis2: {{...axisStyle, domain: [0.00, 0.45], title: "gap",
                type: "log", tickformat: ".0e"}},
      shapes: projectionShapes,
      margin: {{l: 92, r: 34, b: 70, t: 64}}
    }};
    Plotly.newPlot("plot", traces, layout,
                   {{responsive: true, displaylogo: false}});
  </script>
</body>
</html>
"""
    output = Path(str(prefix) + ".html")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(document, encoding="utf-8")
    return output


def draw_static(prefix: Path, title: str, iterations, active_sizes,
                gaps, thresholds, stable_projections, rank_iterations, ranks,
                show: bool) -> Path:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    figure, axes = plt.subplots(2, 1, sharex=True, figsize=(10, 7))
    series = (
        (active_sizes, "active set size", "#2563eb"),
        (gaps, "gap", "#b91c1c"),
    )
    for axes_item, (values, label, color) in zip(axes, series):
        axes_item.plot(iterations, values, color=color, linewidth=1.8)
        axes_item.set_ylabel(label)
        axes_item.grid(True, color="#e5e7eb", linewidth=0.8)
        axes_item.margins(x=0)
        for iteration in stable_projections:
            axes_item.axvline(iteration, color="#059669", linestyle=":",
                              linewidth=1.0)
    if rank_iterations:
        axes[0].plot(rank_iterations, ranks, color="#dc2626", linewidth=1.8,
                     label="affine rank")
        axes[0].legend(loc="upper left", framealpha=0.9)
    axes[0].yaxis.set_major_locator(MaxNLocator(integer=True))
    if any(value <= 0 for value in gaps):
        raise ValueError("log scale requires every recorded gap to be positive")
    axes[1].set_yscale("log")
    axes[1].yaxis.set_major_formatter(
        FuncFormatter(lambda value, _: f"{value:.0e}")
    )
    if all(value > 0 for value in thresholds):
        axes[1].plot(iterations, thresholds, color="#2563eb", linewidth=1.8,
                     label="stopping threshold")
        axes[1].legend(loc="upper right", framealpha=0.9)
    axes[1].set_xlabel("iteration")
    figure.suptitle(title)
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    output = Path(str(prefix) + ".png")
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    if not show:
        plt.close(figure)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path,
                        help="stderr log produced by shortest_gkz --verbose")
    parser.add_argument("--output-prefix", type=Path,
                        help="output path without .html/.png")
    parser.add_argument("--log-gap", action="store_true",
                        help=argparse.SUPPRESS)
    parser.add_argument("--tolerance", type=float, default=None,
                        help="relative gap tolerance (default: read from log, otherwise 1e-11)")
    parser.add_argument("--absolute-tolerance", type=float, default=None,
                        help="absolute gap tolerance (default: read from log, otherwise 1e-14)")
    parser.add_argument("--png", action="store_true",
                        help="also export a static PNG")
    parser.add_argument("--show", action="store_true",
                        help="display the optional PNG after saving it")
    args = parser.parse_args()

    try:
        recorded_tolerances = read_solver_tolerances(args.log)
        history = read_history(args.log)
    except (OSError, ValueError) as error:
        parser.error(str(error))

    recorded_tolerance, recorded_absolute = recorded_tolerances or (1e-11, 1e-14)
    tolerance = recorded_tolerance if args.tolerance is None else args.tolerance
    absolute_tolerance = (recorded_absolute if args.absolute_tolerance is None
                          else args.absolute_tolerance)
    if (not math.isfinite(tolerance) or tolerance < 0 or
            not math.isfinite(absolute_tolerance) or absolute_tolerance < 0):
        parser.error("tolerances must be finite and nonnegative")

    (iterations, active_sizes, norm2_values, gaps, stable_projections,
     rank_iterations, ranks) = history
    if any(value <= 0 for value in gaps):
        parser.error("log scale requires every recorded gap to be positive")

    # Per-iteration stopping threshold used by shortest_gkz:
    # absolute_tolerance + tolerance * norm2.
    thresholds = [
        absolute_tolerance + tolerance * max(1e-300, norm2)
        for norm2 in norm2_values
    ]

    prefix = args.output_prefix or default_prefix(args.log)
    title = f"Iteration history: {args.log.stem}"
    outputs = [draw_interactive(prefix, title, iterations, active_sizes,
                                gaps, thresholds, stable_projections,
                                rank_iterations, ranks)]

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
            parser.error("matplotlib is required for --png")
            raise error
        outputs.append(draw_static(prefix, title, iterations, active_sizes,
                                   gaps, thresholds, stable_projections,
                                   rank_iterations, ranks,
                                   args.show))
    else:
        plt = None

    for output in outputs:
        print(output)
    if args.show and plt is not None:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
