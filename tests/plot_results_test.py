#!/usr/bin/env python3
"""Small standard-library regression test for plot_results.py."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path


def write_fixture(prefix: Path, include_psi: bool = False) -> None:
    psi_rows = (",10\n", ",11\n", ",14\n", ",12\n")
    rows = [
        "x,y,sigma,sigma_vee,psi\n",
        "0,0,0,0" + (psi_rows[0] if include_psi else ",\n"),
        "1,0,1,1" + (psi_rows[1] if include_psi else ",\n"),
        "1,1,4,4" + (psi_rows[2] if include_psi else ",\n"),
        "0,1,2,2" + (psi_rows[3] if include_psi else ",\n"),
    ]
    prefix.with_name(prefix.name + "_surface.csv").write_text(
        "".join(rows),
        encoding="utf-8",
    )
    prefix.with_name(prefix.name + "_triangles.csv").write_text(
        "i,j,l\n0,1,2\n0,2,3\n", encoding="utf-8"
    )
    prefix.with_name(prefix.name + "_subdivision.csv").write_text(
        "cell,vertex,x,y\n"
        "0,0,0,0\n0,1,1,0\n0,2,1,1\n"
        "1,0,0,0\n1,1,1,1\n1,2,0,1\n",
        encoding="utf-8",
    )


def run(script: Path, prefix: Path, *options: str) -> str:
    command = [sys.executable, str(script), str(prefix), *options]
    return subprocess.check_output(command, text=True)


def output(prefix: Path, suffix: str) -> Path:
    return prefix.with_name(prefix.name + suffix)


def main() -> int:
    script = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix="plot-results-test-") as name:
        root = Path(name)

        scaled = root / "scaled"
        write_fixture(scaled)
        run(script, scaled, "--z-scale", "10")
        html = output(scaled, "_sigma_vee.html").read_text(encoding="utf-8")
        assert 'zaxis: {title: "value"}' in html
        assert 'surface.edges.forEach' in html
        assert html.count('type: "scatter3d"') == 1
        # z remains [0, 1, 4, 2]; only the scene aspect ratio is stretched.
        assert '"z":[0.0,1.0,4.0,2.0]' in html
        assert 'aspectmode: "manual"' in html
        assert 'aspectratio: {"x":1.0,"y":1.0,"z":40.0}' in html
        assert 'color: "#9ecae1"' in html
        assert 'color: "#991b1b"' in html
        assert 'colorscale: "Viridis"' not in html
        assert output(scaled, "_subdivision.svg").exists()

        polygon = root / "polygon"
        write_fixture(polygon, include_psi=True)
        run(script, polygon, "--z-scale", "10")
        assert not output(polygon, "_sigma_vee.html").exists()
        assert output(polygon, "_psi_k.html").exists()
        psi_html = output(polygon, "_psi_k.html").read_text(
            encoding="utf-8")
        assert "psi_k (last QP solution)" in psi_html
        subdivision = output(polygon, "_subdivision.svg").read_text(
            encoding="utf-8")
        assert "S(sigma_k) (last QP solution)" in subdivision
        run(script, polygon, "--z-scale", "10", "--sigma-vee")
        assert output(polygon, "_sigma_vee.html").exists()

        stable = root / "polygon_stable"
        write_fixture(stable, include_psi=True)
        run(script, stable)
        stable_psi = output(stable, "_psi_k.html").read_text(
            encoding="utf-8")
        assert "psi_k (stable projection)" in stable_psi
        stable_subdivision = output(stable, "_subdivision.svg").read_text(
            encoding="utf-8")
        assert "S(sigma_k) (stable projection)" in stable_subdivision

        colored = root / "colored"
        write_fixture(colored)
        run(script, colored, "--color")
        colored_html = output(colored, "_sigma_vee.html").read_text(
            encoding="utf-8"
        )
        assert 'intensity: surface.z' in colored_html
        assert 'colorscale: "Viridis"' in colored_html
        assert 'color: "#9ecae1"' not in colored_html

        no_edges = root / "no_edges"
        write_fixture(no_edges)
        run(script, no_edges, "--z-scale", "10", "--no-edges")
        no_edges_html = output(no_edges, "_sigma_vee.html").read_text(
            encoding="utf-8"
        )
        assert 'surface.edges.forEach' not in no_edges_html
        assert 'type: "scatter3d"' not in no_edges_html

        no_subdivision = root / "no_subdivision"
        write_fixture(no_subdivision)
        run(script, no_subdivision, "--no-subdivision")
        assert not output(no_subdivision, "_subdivision.svg").exists()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
