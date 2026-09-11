#!/usr/bin/env python3
"""Assign stable non-overlapping ky task-number blocks to staged v1.0 screens."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--start", type=int, default=1400)
    parser.add_argument("--block", type=int, default=32)
    parser.add_argument(
        "--reserved", action="append", default=[], metavar="START-END",
        help="task-number interval reserved by another workflow; may be repeated",
    )
    args = parser.parse_args()
    reserved: list[tuple[int, int]] = []
    for interval in args.reserved:
        left, right = (int(value) for value in interval.split("-", 1))
        if left > right:
            raise ValueError(f"invalid reserved interval: {interval}")
        reserved.append((left, right))

    def intersects_reserved(base: int) -> bool:
        end = base + args.block - 1
        return any(base <= right and end >= left for left, right in reserved)

    screens = sorted(path.parent for path in args.root.glob("*/STAGING_READY"))
    used = set()
    for screen in screens:
        path = screen / "v1_job_base.tsv"
        if not path.exists():
            continue
        with path.open() as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        base = int(row["job_base"])
        if intersects_reserved(base):
            raise ValueError(f"{path} overlaps a reserved task-number interval")
        used.add(base)
    next_base = args.start
    for screen in screens:
        path = screen / "v1_job_base.tsv"
        if path.exists():
            continue
        while next_base in used or intersects_reserved(next_base):
            next_base += args.block
        row = {
            "screen": screen.name,
            "job_base": next_base,
            "stoichiometry_jobs": f"{next_base}-{next_base + 7}",
            "protenix_jobs": f"{next_base + 8}-{next_base + 15}",
            "opendde_jobs": f"{next_base + 16}-{next_base + 23}",
        }
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, list(row), delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerow(row)
        print(f"allocated {screen.name}: {next_base}-{next_base + args.block - 1}")
        used.add(next_base)
        next_base += args.block


if __name__ == "__main__":
    main()
