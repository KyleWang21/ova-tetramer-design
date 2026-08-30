#!/usr/bin/env python3
"""Combine the first N records from one or more FASTA files."""

from __future__ import annotations

import argparse
import pathlib


def records(path: pathlib.Path) -> list[tuple[str, str]]:
    output = []
    header = None
    sequence: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                output.append((header, "".join(sequence)))
            header, sequence = line, []
        elif header is not None:
            sequence.append(line.strip())
    if header is not None:
        output.append((header, "".join(sequence)))
    return output


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", action="append", type=pathlib.Path, required=True)
    ap.add_argument("--top-per-input", type=int, default=1)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    args = ap.parse_args()
    combined = []
    for path in args.input:
        combined.extend(records(path)[: args.top_per_input])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as handle:
        for header, sequence in combined:
            handle.write(f"{header}\n{sequence}\n")
    print(f"wrote {len(combined)} records to {args.out}")


if __name__ == "__main__":
    main()
