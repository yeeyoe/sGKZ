#!/usr/bin/env python3
"""Plot the area--Q_P(g)^2 distribution from the search database.

The database stores the candidate-level instability indicator as
``q_squared_exact``/``q_squared_value``.  Ranking is performed with exact
fractions using Q_P(g)^2 / V_P; floating point values are only used by the
interactive chart.
"""

from __future__ import annotations

import argparse
import html
import json
import sqlite3
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path


@dataclass(frozen=True)
class CandidatePoint:
    key: str
    area: Fraction
    q_squared: Fraction

    @property
    def ratio(self) -> Fraction:
        return self.q_squared / self.area


def parse_fraction(value: str, label: str) -> Fraction:
    try:
        result = Fraction(value.strip())
    except (AttributeError, ValueError, ZeroDivisionError) as error:
        raise ValueError(f"invalid {label}: {value!r}") from error
    return result


def load_points(database: Path, max_area: Fraction) -> list[CandidatePoint]:
    if not database.exists():
        raise ValueError(f"database not found: {database}")
    connection = sqlite3.connect(database)
    connection.row_factory = sqlite3.Row
    try:
        rows = connection.execute(
            "SELECT key, twice_area, q_squared_exact "
            "FROM candidates "
            "WHERE status='verified_unstable' "
            "AND q_squared_value IS NOT NULL "
            "AND q_squared_exact IS NOT NULL"
        ).fetchall()
    except sqlite3.Error as error:
        raise ValueError(f"could not query database: {error}") from error
    finally:
        connection.close()

    points: list[CandidatePoint] = []
    for row in rows:
        twice_area = parse_fraction(row["twice_area"], "twice_area")
        area = twice_area / 2
        if area <= 0:
            continue
        if area > max_area:
            continue
        q_squared = parse_fraction(row["q_squared_exact"], "q_squared_exact")
        if q_squared < 0:
            raise ValueError(f"negative q_squared_exact for candidate {row['key']}")
        points.append(CandidatePoint(row["key"], area, q_squared))
    return points


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def filename_bound(value: Fraction) -> str:
    text = fraction_text(value).replace("/", "-")
    return text.replace("-", "m") if text.startswith("m") else text


def render_html(points: list[CandidatePoint], max_area: Fraction) -> str:
    data = [
        {
            "key": point.key,
            "area": float(point.area),
            "qSquared": float(point.q_squared),
        }
        for point in points
    ]
    encoded = json.dumps(data, ensure_ascii=True, separators=(",", ":")).replace("<", "\\u003c")
    bound = html.escape(fraction_text(max_area))
    count = f"{len(points):,}"
    return f'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Area vs Q_P(g)^2 (V &lt;= {bound})</title>
<style>
  :root {{ color-scheme: light dark; --fg: light-dark(#172033, #e7eaf0); --muted: light-dark(#526078, #aab4c5); --line: light-dark(#cbd3df, #465268); --series: light-dark(#2f73c8, #7fb5ff); --bg: light-dark(#ffffff, #171a20); }}
  html, body {{ margin: 0; padding: 16px; background: var(--bg); color: var(--fg); font-family: system-ui, sans-serif; }}
  h1 {{ font-size: 18px; font-weight: 500; margin: 0 0 4px; }}
  p {{ color: var(--muted); margin: 0 0 12px; }}
  .plot {{ width: 100%; }}
  svg {{ display: block; width: 100%; height: auto; overflow: visible; }}
  .axis path, .axis line, .frame {{ stroke: var(--line); fill: none; }}
  .axis text, .axis-title {{ fill: var(--fg); font-size: 12px; }}
  .grid line {{ stroke: var(--line); opacity: .32; }}
  .point {{ fill: var(--series); opacity: .38; }}
  .hover {{ fill: var(--fg); stroke: var(--bg); stroke-width: 1.5; }}
  .tooltip {{ position: absolute; display: none; max-width: min(460px, calc(100% - 24px)); padding: 8px 10px; background: var(--bg); color: var(--fg); border: 1px solid var(--line); pointer-events: none; overflow-wrap: anywhere; }}
</style>
</head>
<body>
<h1>Area vs. Q<sub>P</sub>(g)<sup>2</sup></h1>
<p>{count} verified unstable candidates with V<sub>P</sub> &le; {bound}</p>
<div class="plot"><svg id="chart" role="img" aria-label="Area versus squared instability indicator"></svg></div>
<div id="tooltip" class="tooltip" role="tooltip"></div>
<script src="https://cdn.jsdelivr.net/npm/d3@7.9.0/dist/d3.min.js"></script>
<script>
(() => {{
  const data = {encoded};
  const svg = d3.select('#chart');
  const tooltip = document.getElementById('tooltip');
  const plot = document.querySelector('.plot');
  let width = 0;
  function draw() {{
    width = Math.max(320, plot.clientWidth || 736);
    const compact = width < 430;
    const height = Math.max(360, Math.min(560, width * 0.68));
    const margin = {{top: 12, right: 16, bottom: compact ? 68 : 62, left: compact ? 72 : 84}};
    const innerWidth = width - margin.left - margin.right;
    const innerHeight = height - margin.top - margin.bottom;
    svg.attr('viewBox', `0 0 ${{width}} ${{height}}`).selectAll('*').remove();
    const chart = svg.append('g').attr('transform', `translate(${{margin.left}},${{margin.top}})`);
    const x = d3.scaleLog().domain(d3.extent(data, d => d.area)).nice().range([2, innerWidth - 2]);
    const y = d3.scaleSymlog().constant(1e-4).domain([0, d3.max(data, d => d.qSquared)]).nice().range([innerHeight - 2, 2]);
    const xd = x.domain(), yd = y.domain();
    const xv = compact ? d3.range(4).map(i => 10 ** (Math.log10(xd[0]) + i * (Math.log10(xd[1]) - Math.log10(xd[0])) / 3)) : null;
    const yv = compact ? d3.range(5).map(i => yd[0] + i * (yd[1] - yd[0]) / 4) : null;
    const xa = () => (xv ? d3.axisBottom(x).tickValues(xv) : d3.axisBottom(x).ticks(7));
    const ya = () => (yv ? d3.axisLeft(y).tickValues(yv) : d3.axisLeft(y).ticks(7));
    chart.append('g').attr('class', 'grid').attr('transform', `translate(0,${{innerHeight}})`).call(xa().tickSize(-innerHeight).tickFormat(''));
    chart.append('g').attr('class', 'grid').call(ya().tickSize(-innerWidth).tickFormat(''));
    chart.append('rect').attr('class', 'frame').attr('x', 0).attr('y', 0).attr('width', innerWidth).attr('height', innerHeight);
    chart.append('g').attr('class', 'axis').attr('transform', `translate(0,${{innerHeight}})`).call(xa().tickFormat(d3.format('.3~s')));
    chart.append('g').attr('class', 'axis').call(ya().tickFormat(d3.format('.3~s')));
    chart.append('text').attr('class', 'axis-title').attr('data-axis', 'x').attr('x', innerWidth / 2).attr('y', innerHeight + margin.bottom - 16).attr('text-anchor', 'middle').text('Area V_P (log scale)');
    chart.append('text').attr('class', 'axis-title').attr('data-axis', 'y').attr('transform', 'rotate(-90)').attr('x', -innerHeight / 2).attr('y', -margin.left + 20).attr('text-anchor', 'middle').text('Q_P(g)^2 (symlog; zero included)');
    const projected = data.map(d => ({{d, px: x(d.area), py: y(d.qSquared)}}));
    chart.append('g').selectAll('circle').data(projected).join('circle').attr('class', 'point').attr('cx', d => d.px).attr('cy', d => d.py).attr('r', 2.3);
    const hover = chart.append('circle').attr('class', 'hover').attr('r', 5).attr('visibility', 'hidden');
    const tree = d3.quadtree().x(d => d.px).y(d => d.py).addAll(projected);
    const overlay = chart.append('rect').attr('data-chart-hit', 'true').attr('x', 0).attr('y', 0).attr('width', innerWidth).attr('height', innerHeight).attr('fill', 'transparent');
    overlay.on('pointermove', event => {{
      const [px, py] = d3.pointer(event, overlay.node());
      const hit = tree.find(px, py, compact ? 28 : 36);
      if (!hit) {{ hover.attr('visibility', 'hidden'); tooltip.style.display = 'none'; return; }}
      hover.attr('cx', hit.px).attr('cy', hit.py).attr('visibility', 'visible');
      tooltip.textContent = `${{hit.d.key}} | Area: ${{d3.format(',.6~g')(hit.d.area)}} | Q_P(g)^2: ${{d3.format(',.6~g')(hit.d.qSquared)}}`;
      tooltip.style.display = 'block';
      const [tx, ty] = d3.pointer(event, document.body);
      tooltip.style.left = `${{Math.min(tx + 12, document.body.clientWidth - tooltip.offsetWidth - 8)}}px`;
      tooltip.style.top = `${{Math.max(8, ty + 12)}}px`;
    }}).on('pointerleave', () => {{ hover.attr('visibility', 'hidden'); tooltip.style.display = 'none'; }});
  }}
  draw();
  new ResizeObserver(draw).observe(plot);
}})();
</script>
</body>
</html>
'''


def print_ranking(points: list[CandidatePoint], top_n: int, max_area: Fraction) -> None:
    print(f"filtered_candidates={len(points)}")
    print(f"max_area={fraction_text(max_area)}")
    print("rank\tkey\tarea\tq_squared_over_area")
    for index, point in enumerate(sorted(points, key=lambda item: (-item.ratio, item.key))[:top_n], 1):
        print(f"{index}\t{point.key}\t{fraction_text(point.area)}\t{float(point.ratio):.4f}")
    if not points:
        print("no candidates satisfy the area and Q-value filters")


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot area versus Q_P(g)^2 from the K-stability search database.")
    parser.add_argument("--max-area", required=True, help="keep candidates with V_P <= this value")
    parser.add_argument("--top-n", type=int, default=10, help="number of ranked keys to print (default: 10)")
    parser.add_argument("--database", type=Path, default=Path(__file__).with_name("k_stability_search.sqlite"))
    parser.add_argument("--output", type=Path, help="HTML output path")
    args = parser.parse_args()
    try:
        max_area = parse_fraction(args.max_area, "max-area")
        if max_area <= 0:
            raise ValueError("max-area must be positive")
        if args.top_n < 1:
            raise ValueError("top-n must be at least 1")
        points = load_points(args.database, max_area)
        if not points:
            print_ranking(points, args.top_n, max_area)
            return 1
        output = args.output or (Path(__file__).resolve().parent / "Q_volume_joint_graph" /
                                 f"area-q-squared-v{filename_bound(max_area)}.html")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(render_html(points, max_area), encoding="utf-8")
        print_ranking(points, args.top_n, max_area)
        print(f"html={output}")
        return 0
    except (OSError, sqlite3.Error, ValueError) as error:
        print(f"error: {error}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
