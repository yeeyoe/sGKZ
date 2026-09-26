#!/usr/bin/env python3
"""Regression tests for K-stability/plot_area_q.py."""

from __future__ import annotations

import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path


def main() -> int:
    script = Path(__file__).resolve().parents[1] / "K-stability" / "plot_area_q.py"
    with tempfile.TemporaryDirectory(prefix="plot-area-q-test-") as name:
        root = Path(name)
        database = root / "fixture.sqlite"
        output = root / "area.html"
        connection = sqlite3.connect(database)
        connection.execute(
            "CREATE TABLE candidates("
            "key TEXT PRIMARY KEY, twice_area TEXT, status TEXT, "
            "q_squared_exact TEXT, q_squared_value REAL)"
        )
        connection.executemany(
            "INSERT INTO candidates VALUES(?,?,?,?,?)",
            [
                ("boundary", "20", "verified_unstable", "3", 3.0),
                ("best", "10", "verified_unstable", "3", 3.0),
                ("aaa-tie", "10", "verified_unstable", "3", 3.0),
                ("exact-order", "6", "verified_unstable", "2", 2.0),
                ("over-limit", "22", "verified_unstable", "100", 100.0),
                ("no-q", "4", "verified_unstable", None, None),
                ("pending", "2", "pending", "1000", 1000.0),
            ],
        )
        connection.commit()
        connection.close()

        command = [
            sys.executable, str(script), "--max-area", "10", "--top-n", "3",
            "--database", str(database), "--output", str(output),
        ]
        result = subprocess.run(command, text=True, capture_output=True, check=True)
        lines = result.stdout.splitlines()
        assert "filtered_candidates=4" in lines
        assert "rank\tkey\tarea\tq_squared_over_area" in lines
        assert not any("q_squared_exact" in line for line in lines)
        ranking = [line for line in lines if line.startswith(("1\t", "2\t", "3\t"))]
        assert all(len(line.split("\t")[-1].split(".")[-1]) == 4 for line in ranking)
        assert ranking[0].split("\t")[1] == "exact-order"
        assert ranking[1].split("\t")[1] == "aaa-tie"
        assert ranking[2].split("\t")[1] == "best"
        assert output.exists()
        html = output.read_text(encoding="utf-8")
        assert "4 verified unstable candidates" in html
        assert "Area V_P (log scale)" in html
        assert "Q_P(g)^2" in html
        assert "boundary" in html
        assert "over-limit" not in html

        empty = subprocess.run(
            [sys.executable, str(script), "--max-area", "1", "--database", str(database)],
            text=True, capture_output=True,
        )
        assert empty.returncode != 0
        assert "no candidates satisfy" in empty.stdout

        invalid = subprocess.run(
            [sys.executable, str(script), "--max-area", "0", "--database", str(database)],
            text=True, capture_output=True,
        )
        assert invalid.returncode != 0
        assert "must be positive" in invalid.stdout
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
