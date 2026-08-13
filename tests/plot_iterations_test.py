#!/usr/bin/env python3
"""Standard-library regression test for plot_iterations.py."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path


def main() -> int:
    script = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix="plot-iterations-test-") as name:
        root = Path(name)
        log = root / "sample.log"
        log.write_text(
            "parameters tolerance=2e-3 absolute_tolerance=4e-5\n"
            "iteration=0 active=1 norm2=0.5 gap=0.25\n"
            "iteration=1 active=2 norm2=0.4 gap=1e-2\n"
            "projection event=rank iteration=1 rank=2 observations=3 stall=0\n"
            "projection event=stable iteration=1 rank=2 observations=3 "
            "p_norm2=0.3 vertex_new=true\n"
            "iteration=2 active=2 norm2=0.39 gap=1e-5\n",
            encoding="utf-8",
        )
        subprocess.run(
            [sys.executable, str(script), str(log)],
            check=True,
            capture_output=True,
            text=True,
        )
        output = root / "sample_iterations.html"
        document = output.read_text(encoding="utf-8")
        assert '"iteration":[0,1,2]' in document
        assert '"active":[1,2,2]' in document
        assert '"gap":[0.25,0.01,1e-05]' in document
        assert '"threshold":[' in document
        assert '0.00104' in document
        assert '"stable_projection":[1]' in document
        assert '"rank_iteration":[1]' in document
        assert '"rank":[2]' in document
        assert 'name: "affine rank"' in document
        assert '#dc2626' in document
        assert "projectionShapes" in document
        assert 'title: "active set size"' in document
        assert 'title: "gap"' in document
        assert 'title: "norm2"' not in document
        assert 'type: "log"' in document
        assert 'tickformat: ".0e"' in document
        assert '"stopping threshold"' in document

        bad_log = root / "multiple.log"
        bad_log.write_text(
            "iteration=0 active=1 norm2=0.5 gap=0.25\n"
            "iteration=0 active=1 norm2=0.5 gap=0.25\n",
            encoding="utf-8",
        )
        invalid = subprocess.run(
            [sys.executable, str(script), str(bad_log)],
            capture_output=True,
            text=True,
        )
        assert invalid.returncode != 0
        assert "strictly increasing" in invalid.stderr
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
