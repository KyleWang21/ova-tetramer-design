#!/usr/bin/env python3
"""Merge recovered tied-AF2 trajectory libraries by exact sequence.

For a sequence observed in more than one run, retain the observation with the
largest weakest-state AF2 ipTM.  These scores remain generation-only evidence;
the merged library is intended for a new template-free AF3 screen.
"""

from __future__ import annotations

import argparse
import csv
import math
import pathlib


def number(value: object, default: float) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result if math.isfinite(result) else default


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", type=pathlib.Path, action="append", required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    unique: dict[str, dict[str, str]] = {}
    for path in args.library:
        for row in read_tsv(path):
            sequence = row.get("sequence", "")
            if not sequence:
                continue
            candidate = dict(row)
            candidate["recovered_library"] = str(path)
            old = unique.get(sequence)
            candidate_key = (
                number(candidate.get("af2_iptm"), -math.inf),
                -number(candidate.get("generation_average_rank"), math.inf),
            )
            old_key = (
                number(old.get("af2_iptm"), -math.inf),
                -number(old.get("generation_average_rank"), math.inf),
            ) if old else (-math.inf, -math.inf)
            if old is None or candidate_key > old_key:
                unique[sequence] = candidate

    rows = sorted(
        unique.values(),
        key=lambda row: (
            -number(row.get("af2_iptm"), -math.inf),
            number(row.get("generation_average_rank"), math.inf),
            row["sequence"],
        ),
    )
    if not rows:
        raise SystemExit("no recovered trajectory candidates")
    args.out.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    table = args.out / "merged_recovered_trajectory_library.tsv"
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    fasta = args.out / "merged_recovered_trajectory_library.fasta"
    with fasta.open("w") as handle:
        for index, row in enumerate(rows, 1):
            handle.write(
                f">RECOVERED-TRAJ-{index:04d} nmut={row.get('n_mutations', '')} "
                f"weak_af2={row.get('af2_iptm', '')}\n{row['sequence']}\n"
            )
    print(f"merged_unique={len(rows)} table={table}")


if __name__ == "__main__":
    main()
