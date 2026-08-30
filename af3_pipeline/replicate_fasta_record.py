#!/usr/bin/env python3
"""Write repeated copies of one FASTA record with unique seed annotations."""

from __future__ import annotations

import argparse
import pathlib

from combine_fasta_top import records


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=pathlib.Path, required=True)
    ap.add_argument("--record", type=int, required=True, help="one-based FASTA record")
    ap.add_argument("--copies", type=int, required=True)
    ap.add_argument("--seed-start", type=int, default=1)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    args = ap.parse_args()
    selected = records(args.input)[args.record - 1]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as handle:
        for offset in range(args.copies):
            handle.write(
                f"{selected[0]}|af3_seed={args.seed_start + offset}\n{selected[1]}\n"
            )
    print(f"replicated record {args.record} x{args.copies} to {args.out}")


if __name__ == "__main__":
    main()
