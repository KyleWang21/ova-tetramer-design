#!/usr/bin/env python3
"""Build a compact factorial library from AF3-contacting donor mutations."""

from __future__ import annotations

import argparse
import itertools
import pathlib


def sequence(path: pathlib.Path) -> str:
    return "".join(line.strip() for line in path.read_text().splitlines() if not line.startswith(">"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", type=pathlib.Path, required=True)
    ap.add_argument("--base", type=pathlib.Path, action="append", required=True)
    ap.add_argument("--donor", type=pathlib.Path, required=True)
    ap.add_argument("--positions", default="93,94,96,155,156,192")
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--max-mutations", type=int, default=24)
    args = ap.parse_args()
    ref, donor = sequence(args.reference), sequence(args.donor)
    positions = [int(item) for item in args.positions.split(",")]
    records: list[tuple[str, str]] = []
    seen: set[str] = set()
    for base_index, base_path in enumerate(args.base, 1):
        base = sequence(base_path)
        for bits in itertools.product((0, 1), repeat=len(positions)):
            chars = list(base)
            used = []
            for position, enabled in zip(positions, bits):
                if enabled:
                    chars[position - 1] = donor[position - 1]
                    used.append(position)
            candidate = "".join(chars)
            nmut = sum(a != b for a, b in zip(ref, candidate))
            if nmut >= 25 or nmut > args.max_mutations or candidate in seen:
                continue
            if candidate[257:265] != ref[257:265]:
                continue
            seen.add(candidate)
            records.append((f"AF3REC-B{base_index}-P{'.'.join(map(str, used)) or 'none'} nmut={nmut}", candidate))
    args.out.write_text("".join(f">{header}\n{seq}\n" for header, seq in records))
    print(f"wrote {len(records)} unique recombinations to {args.out}")


if __name__ == "__main__":
    main()
