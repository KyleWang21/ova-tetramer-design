#!/usr/bin/env python3
"""Copy the first FASTA records with stable numbered headers."""

from __future__ import annotations

import argparse
import pathlib


def records(path: pathlib.Path) -> list[tuple[str, str]]:
    output: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                output.append((header, "".join(chunks)))
            header, chunks = line[1:], []
        elif header is not None:
            chunks.append(line.strip())
    if header is not None:
        output.append((header, "".join(chunks)))
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--top", type=int, required=True)
    parser.add_argument("--prefix", default="OVA-IPTM80-MS")
    args = parser.parse_args()
    selected = records(args.input)[: args.top]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as handle:
        for index, (old_header, sequence) in enumerate(selected, 1):
            handle.write(f">{args.prefix}{index:02d} source={old_header.split()[0]}\n{sequence}\n")
    print(f"wrote {len(selected)} records to {args.out}")


if __name__ == "__main__":
    main()
