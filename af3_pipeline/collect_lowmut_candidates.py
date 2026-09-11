#!/usr/bin/env python3
"""Merge two-state AF2 design outputs into a deduplicated AF3 shortlist."""

from __future__ import annotations

import argparse
import csv
import pathlib


def one_row(path: pathlib.Path) -> dict[str, str]:
    with path.open() as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def f(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", type=pathlib.Path, required=True)
    parser.add_argument("--max-mutations", type=int, default=25)
    parser.add_argument("--top", type=int, default=16)
    parser.add_argument(
        "--warm-quota", type=int, default=0,
        help="reserve this many slots for current-C4-initialized interface-pruning runs",
    )
    args = parser.parse_args()

    rows: list[dict[str, object]] = []
    stage2_paths = sorted(args.experiment.glob("run_*/stage2_second_interface/af2diff_candidates.tsv"))
    stage2_paths += sorted(args.experiment.glob("run_*/stage2_second_interface_warm/af2diff_candidates.tsv"))
    for stage2_path in stage2_paths:
        run = stage2_path.parents[1].name
        branch = "warm" if stage2_path.parent.name.endswith("_warm") else "cold"
        stage1_path = stage2_path.parents[1] / "stage1_first_interface" / "af2diff_candidates.tsv"
        if not stage1_path.exists():
            continue
        first, second = one_row(stage1_path), one_row(stage2_path)
        nmut = int(second["n_mutations"])
        if nmut > args.max_mutations:
            continue
        rows.append({
            "candidate": f"OVA-LM-C4-{run.upper()}-{branch.upper()}",
            "run": run,
            "branch": branch,
            "n_mutations": nmut,
            "mutations": second["mutations"],
            "first_iptm": f(first["iptm"]),
            "second_iptm": f(second["iptm"]),
            "first_interface_pae": f(first["interface_pae"]),
            "second_interface_pae": f(second["interface_pae"]),
            "first_target_con": f(first["first_target_con"]),
            "second_target_con": f(second["second_target_con"]),
            "first_plddt": f(first["plddt"]),
            "second_plddt": f(second["plddt"]),
            "stage1_loss": f(first["loss"]),
            "stage2_loss": f(second["loss"]),
            "sequence": second["sequence"],
        })

    # One representative per exact sequence.  Put <=20-mutation designs first,
    # then favor the weakest of the two AF2 interface scores and fold confidence.
    unique: dict[str, dict[str, object]] = {}
    for row in rows:
        seq = str(row["sequence"])
        rank = (
            int(row["n_mutations"]) > 20,
            -min(float(row["first_iptm"]), float(row["second_iptm"])),
            -min(float(row["first_plddt"]), float(row["second_plddt"])),
            int(row["n_mutations"]),
            float(row["stage1_loss"]) + float(row["stage2_loss"]),
        )
        old = unique.get(seq)
        if old is None or rank < old["_rank"]:
            row["_rank"] = rank
            unique[seq] = row
    ranked = sorted(unique.values(), key=lambda row: row["_rank"])
    if args.warm_quota:
        selected = [row for row in ranked if row["branch"] == "warm"][: args.warm_quota]
        selected_sequences = {str(row["sequence"]) for row in selected}
        selected += [
            row for row in ranked if str(row["sequence"]) not in selected_sequences
        ][: max(0, args.top - len(selected))]
    else:
        selected = ranked[: args.top]
    if not selected:
        raise SystemExit("no completed candidates at or below the mutation limit")

    table = args.experiment / "combined_candidates.tsv"
    fasta = args.experiment / "shortlist.fasta"
    fields = [key for key in selected[0] if key != "_rank"]
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows({key: row[key] for key in fields} for row in selected)
    with fasta.open("w") as handle:
        for row in selected:
            handle.write(
                f">{row['candidate']} mutations={row['n_mutations']} "
                f"first_ipTM={float(row['first_iptm']):.3f} "
                f"second_ipTM={float(row['second_iptm']):.3f}\n{row['sequence']}\n"
            )
    print(f"selected {len(selected)} unique candidates from {len(rows)} completed runs")
    print(table)
    print(fasta)


if __name__ == "__main__":
    main()
