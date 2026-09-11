#!/usr/bin/env python3
"""Attach full sequences to AF3 summary rows using a candidate table."""

from __future__ import annotations

import argparse
import csv
import pathlib


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=pathlib.Path, required=True)
    parser.add_argument("--sequences", type=pathlib.Path, required=True)
    parser.add_argument("--summary-key", default="candidate")
    parser.add_argument("--sequence-key", default="candidate")
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    sequences = {
        row[args.sequence_key]: row["sequence"] for row in read_tsv(args.sequences)
    }
    rows = read_tsv(args.summary)
    for row in rows:
        key = row[args.summary_key]
        if key not in sequences:
            raise SystemExit(f"missing sequence for {key}")
        if len(sequences[key]) != 386:
            raise SystemExit(f"unexpected sequence length for {key}: {len(sequences[key])}")
        row["sequence"] = sequences[key]
    if not rows:
        raise SystemExit("empty summary")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows to {args.out}")


if __name__ == "__main__":
    main()
