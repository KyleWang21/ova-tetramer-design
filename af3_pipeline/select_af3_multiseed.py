#!/usr/bin/env python3
"""Select complete single-seed AF3 candidates for multi-seed follow-up."""

from __future__ import annotations

import argparse
import csv
import pathlib


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--minimum", type=float, default=0.78)
    ap.add_argument("--exclude-fasta", type=pathlib.Path, action="append", default=[])
    args = ap.parse_args()
    excluded: set[str] = set()
    for path in args.exclude_fasta:
        sequence: list[str] = []
        for line in path.read_text().splitlines() + [">"]:
            if line.startswith(">"):
                if sequence:
                    excluded.add("".join(sequence))
                sequence = []
            else:
                sequence.append(line.strip())
    rows = list(csv.DictReader(args.table.open(), delimiter="\t"))
    rows = [
        row for row in rows
        if float(row["mean_iptm"]) >= args.minimum
        and int(row["n_mutations"]) < 25
        and row["sequence"] not in excluded
    ]
    rows.sort(key=lambda row: (-float(row["mean_iptm"]), -float(row["min_iptm"])))
    rows = rows[: args.top]
    if not rows:
        raise SystemExit("no candidates pass selection threshold")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as handle:
        for rank, row in enumerate(rows, 1):
            handle.write(
                f">OVA-IPTM80-MS{rank:02d} source={row['system']} "
                f"nmut={row['n_mutations']} seed1_mean={float(row['mean_iptm']):.3f} "
                f"seed1_min={float(row['min_iptm']):.3f}\n{row['sequence']}\n"
            )
    print(f"selected {len(rows)} candidates")
    for row in rows:
        print(row["system"], row["n_mutations"], f"mean={float(row['mean_iptm']):.3f}")


if __name__ == "__main__":
    main()
