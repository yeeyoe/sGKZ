#!/usr/bin/env python3
"""Regression test for the Wang-Zhou polygon generator."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path


def main() -> int:
    generator = Path(sys.argv[1]).resolve()
    shortest_gkz = Path(sys.argv[2]).resolve()
    sys.path.insert(0, str(generator.parent))
    import generate_wang_zhou

    assert generate_wang_zhou.default_output(5).name == "Wang_Zhou_a5"
    expected = [
        "0 21",
        "1 21",
        "2 20",
        "3 18",
        "4 15",
        "5 11",
        "7 1",
        "7 -1",
        "0 -1",
    ]

    with tempfile.TemporaryDirectory(prefix="wang-zhou-test-") as name:
        output = Path(name) / "Wang_Zhou_a1"
        subprocess.run(
            [sys.executable, str(generator), "1", "--output", str(output)],
            check=True,
            capture_output=True,
            text=True,
        )
        lines = [
            line
            for line in output.read_text(encoding="ascii").splitlines()
            if line and not line.startswith("#")
        ]
        assert lines == expected

        result = subprocess.run(
            [
                str(shortest_gkz),
                "--polygon",
                str(output),
                "--k",
                "1",
                "--max-iterations",
                "0",
                "--no-exact",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode in (0, 2), result.stderr
        assert "points=" in result.stdout
        assert "base_twice_area=" in result.stdout

        invalid = subprocess.run(
            [sys.executable, str(generator), "0", "--output", str(output)],
            capture_output=True,
            text=True,
        )
        assert invalid.returncode != 0
        assert "a must be positive" in invalid.stderr

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
