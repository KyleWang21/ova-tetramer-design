#!/usr/bin/env python3
"""Continue the best fully discrete sequence across four AF3 backbone contexts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import shutil
from pathlib import Path


def read_fasta(path: Path, wanted: str) -> str:
    records: dict[str, str] = {}
    name = None
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = ""
        elif name:
            records[name] += line
    return records[wanted]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--max-mutations", type=int, default=24)
    parser.add_argument("--min-surface-fraction", type=float, default=0.8)
    parser.add_argument("--run-glob", default="run_xm3_*")
    parser.add_argument("--anchor-prefix", default="CB1")
    parser.add_argument("--source-label", default="E243-BEST-DISCRETE")
    args = parser.parse_args()

    frozen = args.out / "SOURCE_SELECTION_FROZEN"
    if frozen.exists():
        print(f"{args.out} already frozen")
        return
    reference = read_fasta(args.reference_fasta, args.reference_name)
    surface = set(map(int, re.findall(r"\d+", (args.source / "surface_positions.txt").read_text())))
    glyco = [
        i for i in range(len(reference) - 2)
        if reference[i] == "N" and reference[i + 1] != "P" and reference[i + 2] in "ST"
    ]
    epitope_start = reference.index("SIINFEKL")

    observations = []
    for path in sorted(args.source.glob(f"{args.run_glob}/**/trajectory.tsv")):
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                sequence = row["sequence"]
                nmut = int(row["n_mutations"])
                fraction = float(row["surface_mutation_fraction"])
                if row["stage"] != "hard" or float(row["soft"]) < 0.999:
                    continue
                if nmut > args.max_mutations or fraction < args.min_surface_fraction:
                    continue
                if len(sequence) != len(reference):
                    continue
                if any((a == "C") != (b == "C") for a, b in zip(reference, sequence)):
                    continue
                if sequence[epitope_start:epitope_start + 8] != "SIINFEKL":
                    continue
                if any(sequence[i:i + 3] != reference[i:i + 3] for i in glyco):
                    continue
                state1 = float(row["state1_iptm"])
                state2 = float(row["state2_iptm"])
                observations.append({
                    **row,
                    "source_trajectory": str(path.relative_to(args.source)),
                    "weakest_tied_af2_iptm": min(state1, state2),
                    "mean_tied_af2_iptm": (state1 + state2) / 2,
                    "minimum_plddt": min(float(row["state1_plddt"]), float(row["state2_plddt"])),
                })
    if not observations:
        raise SystemExit("no eligible fully discrete trajectory state")

    best_by_sequence = {}
    for row in observations:
        sequence = row["sequence"]
        key = (
            row["weakest_tied_af2_iptm"], row["mean_tied_af2_iptm"],
            row["minimum_plddt"], -int(row["n_mutations"]),
        )
        if sequence not in best_by_sequence or key > best_by_sequence[sequence][0]:
            best_by_sequence[sequence] = (key, row)
    best = max((item[1] for item in best_by_sequence.values()), key=lambda row: (
        row["weakest_tied_af2_iptm"], row["mean_tied_af2_iptm"],
        row["minimum_plddt"], -int(row["n_mutations"]),
    ))

    args.out.mkdir(parents=True, exist_ok=False)
    (args.out / "design_masks").mkdir()
    with (args.out / "anchors.fasta").open("w") as handle:
        for index in range(1, 5):
            handle.write(f">{args.anchor_prefix}-SEED-{index:02d}\n{best['sequence']}\n")

    with (args.source / "target_manifest.tsv").open() as handle:
        source_rows = list(csv.DictReader(handle, delimiter="\t"))
        fieldnames = list(source_rows[0].keys())
    if len(source_rows) != 4:
        raise SystemExit(f"expected four E243 target contexts, found {len(source_rows)}")
    output_rows = []
    for index, row in enumerate(source_rows, 1):
        row = dict(row)
        anchor = f"{args.anchor_prefix}-SEED-{index:02d}"
        slug = args.anchor_prefix.lower()
        mask = args.out / "design_masks" / f"{slug}-seed-{index:02d}_joint_design.txt"
        shutil.copy2(row["joint_design"], mask)
        row["anchor"] = anchor
        row["anchor_n_mutations"] = best["n_mutations"]
        try:
            row["joint_design"] = str(mask.resolve().relative_to(Path.cwd().resolve()))
        except ValueError:
            row["joint_design"] = str(mask)
        row["surrogate_name"] = args.source_label
        row["surrogate_rank"] = "1"
        if "surrogate_predicted_af3" in row:
            row["surrogate_predicted_af3"] = "nan"
        if "surrogate_predicted_protenix" in row:
            row["surrogate_predicted_protenix"] = "nan"
        output_rows.append(row)
    with (args.out / "target_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    shutil.copy2(args.source / "surface_positions.txt", args.out / "surface_positions.txt")

    audit_fields = list(best.keys()) + ["sequence_sha256"]
    with (args.out / "source_selection.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=audit_fields, delimiter="\t")
        writer.writeheader()
        writer.writerow({**best, "sequence_sha256": hashlib.sha256(best["sequence"].encode()).hexdigest()})
    frozen.touch()
    print(
        "frozen discrete continuation source",
        best["source_trajectory"], best["stage"], best["iteration"],
        f"nmut={best['n_mutations']}",
        f"tied_AF2={float(best['state1_iptm']):.3f}/{float(best['state2_iptm']):.3f}",
    )


if __name__ == "__main__":
    main()
