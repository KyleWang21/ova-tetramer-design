#!/usr/bin/env python3
"""Collect diverse high-ipTM discrete sequences from AF2 dimer design trajectories."""

from __future__ import annotations

import argparse
import csv
import pathlib


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--top", type=int, default=4)
    args = ap.parse_args()
    rows = []
    for path in sorted(args.root.glob("volc_seed*/trajectory.tsv")):
        for row in csv.DictReader(path.open(), delimiter="\t"):
            row["source"] = str(path)
            rows.append(row)
    native_cys = {12, 31, 74, 121, 368, 383}
    rows = [
        row for row in rows
        if len(row["sequence"]) == 386
        and {index for index, aa in enumerate(row["sequence"], 1) if aa == "C"} == native_cys
    ]
    rows.sort(key=lambda row: (float(row["iptm"]), float(row["plddt"]), -int(row["n_mutations"])), reverse=True)
    selected = []
    for row in rows:
        if row["sequence"] in {old["sequence"] for old in selected}: continue
        if all(sum(a != b for a, b in zip(row["sequence"], old["sequence"])) >= 2 for old in selected):
            selected.append(row)
        if len(selected) == args.top: break
    args.out.mkdir(parents=True, exist_ok=True)
    fields = ["rank", *list(selected[0])]
    with (args.out / "af2_trajectory_shortlist.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t"); writer.writeheader()
        for rank, row in enumerate(selected, 1): writer.writerow({"rank": rank, **row})
    with (args.out / "af2_trajectory_shortlist.fasta").open("w") as handle:
        for rank, row in enumerate(selected, 1):
            handle.write(
                f">OVA-NC-DIMER-AF2-{rank:02d}|seed={row['seed']}|stage={row['stage']}|"
                f"iter={row['iteration']}|iptm={float(row['iptm']):.4f}|mut={row['mutations']}\n"
                f"{row['sequence']}\n"
            )
    for rank, row in enumerate(selected, 1):
        print(rank, f"seed={row['seed']}", row["stage"], f"iter={row['iteration']}",
              f"ipTM={float(row['iptm']):.3f}", f"pLDDT={float(row['plddt']):.3f}",
              f"mut={row['n_mutations']}")


if __name__ == "__main__":
    main()
