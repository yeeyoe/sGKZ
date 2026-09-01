#!/usr/bin/env python3
"""Generate the Wang-Zhou lattice polygon for a positive integer a."""

from __future__ import annotations

import argparse
from pathlib import Path


MAX_A = (1 << 63) - 1 - 20
EXAMPLES_DIR = Path(__file__).resolve().parent / "examples"


def positive_int(value: str) -> int:
    try:
        result = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("a must be an integer") from error
    if result <= 0:
        raise argparse.ArgumentTypeError("a must be positive")
    if result > MAX_A:
        raise argparse.ArgumentTypeError(
            f"a must be at most {MAX_A} so every coordinate fits in int64"
        )
    return result


def wang_zhou_vertices(a: int) -> list[tuple[int, int]]:
    return [
        (0, a + 20),
        (1, a + 20),
        (2, a + 19),
        (3, a + 17),
        (4, a + 14),
        (5, a + 10),
        (7, a),
        (7, -a),
        (0, -a),
    ]


def default_output(a: int) -> Path:
    return EXAMPLES_DIR / f"Wang_Zhou_a{a}"


def write_polygon(a: int, output: Path) -> None:
    lines = [
        f"# Wang-Zhou polygon generated with a={a}.",
        "# Vertices are listed clockwise.",
    ]
    lines.extend(f"{x} {y}" for x, y in wang_zhou_vertices(a))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="ascii")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "a",
        nargs="?",
        type=positive_int,
        help="positive integer parameter; prompted for when omitted",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="output polygon path (default: examples/Wang_Zhou_a{a})",
    )
    args = parser.parse_args()

    a = args.a
    if a is None:
        try:
            a = positive_int(input("a = ").strip())
        except (EOFError, argparse.ArgumentTypeError) as error:
            parser.error(str(error) if str(error) else "no value supplied for a")

    output = args.output if args.output is not None else default_output(a)
    write_polygon(a, output)
    print(output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
