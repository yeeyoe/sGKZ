#!/usr/bin/env python3
"""Plot active-set size, norm2, and gap from a shortest_gkz verbose log."""

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


def read_history(path: Path):
    iterations = []
    active_sizes = []
    norm2_values = []
    gaps = []
    with path.open(encoding="utf-8") as stream:
        for line_number, raw_line in enumerate(stream, start=1):
            line = raw_line.strip()
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
    return iterations, active_sizes, norm2_values, gaps


def default_prefix(log_path: Path) -> Path:
    base = log_path.with_suffix("")
    if base.name.endswith(".iterations"):
        return base
    return base.with_name(base.name + "_iterations")


def draw_interactive(prefix: Path, title: str, iterations, active_sizes,
                     norm2_values, gaps, log_gap: bool) -> Path:
    payload = json.dumps(
        {
            "iteration": iterations,
            "active": active_sizes,
            "norm2": norm2_values,
            "gap": gaps,
        },
        separators=(",", ":"),
    )
    title_json = json.dumps(title)
    gap_axis_type = json.dumps("log" if log_gap else "linear")
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
      {{...common, name: "norm2", x: history.iteration,
        y: history.norm2, xaxis: "x2", yaxis: "y2",
        line: {{color: "#047857", width: 2}}}},
      {{...common, name: "gap", x: history.iteration,
        y: history.gap, xaxis: "x3", yaxis: "y3",
        line: {{color: "#b91c1c", width: 2}}}}
    ];
    const axisStyle = {{showgrid: true, gridcolor: "#e5e7eb",
                       zeroline: false, automargin: true}};
    const layout = {{
      title: {title_json},
      showlegend: false,
      hovermode: "x",
      plot_bgcolor: "#ffffff",
      paper_bgcolor: "#ffffff",
      xaxis: {{...axisStyle, domain: [0, 1], anchor: "y",
               matches: "x3", showticklabels: false}},
      yaxis: {{...axisStyle, domain: [0.70, 1.00], title: "active set size",
               rangemode: "tozero", tickformat: "d"}},
      xaxis2: {{...axisStyle, domain: [0, 1], anchor: "y2",
                matches: "x3", showticklabels: false}},
      yaxis2: {{...axisStyle, domain: [0.35, 0.65], title: "norm2"}},
      xaxis3: {{...axisStyle, domain: [0, 1], anchor: "y3",
                title: "iteration"}},
      yaxis3: {{...axisStyle, domain: [0.00, 0.30], title: "gap",
                type: {gap_axis_type}}},
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
                norm2_values, gaps, log_gap: bool, show: bool) -> Path:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    figure, axes = plt.subplots(3, 1, sharex=True, figsize=(10, 9))
    series = (
        (active_sizes, "active set size", "#2563eb"),
        (norm2_values, "norm2", "#047857"),
        (gaps, "gap", "#b91c1c"),
    )
    for axes_item, (values, label, color) in zip(axes, series):
        axes_item.plot(iterations, values, color=color, linewidth=1.8)
        axes_item.set_ylabel(label)
        axes_item.grid(True, color="#e5e7eb", linewidth=0.8)
        axes_item.margins(x=0)
    axes[0].yaxis.set_major_locator(MaxNLocator(integer=True))
    if log_gap:
        if any(value <= 0 for value in gaps):
            raise ValueError("--log-gap requires every recorded gap to be positive")
        axes[2].set_yscale("log")
    axes[2].set_xlabel("iteration")
    figure.suptitle(title)
    figure.tight_layout(rect=(0, 0, 1, 0.97))
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
                        help="use a logarithmic y axis for gap")
    parser.add_argument("--png", action="store_true",
                        help="also export a static PNG")
    parser.add_argument("--show", action="store_true",
                        help="display the optional PNG after saving it")
    args = parser.parse_args()

    try:
        history = read_history(args.log)
        if args.log_gap and any(value <= 0 for value in history[3]):
            parser.error("--log-gap requires every recorded gap to be positive")
    except (OSError, ValueError) as error:
        parser.error(str(error))

    prefix = args.output_prefix or default_prefix(args.log)
    title = f"Iteration history: {args.log.stem}"
    outputs = [draw_interactive(prefix, title, *history, args.log_gap)]

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
        outputs.append(draw_static(prefix, title, *history,
                                   args.log_gap, args.show))
    else:
        plt = None

    for output in outputs:
        print(output)
    if args.show and plt is not None:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
