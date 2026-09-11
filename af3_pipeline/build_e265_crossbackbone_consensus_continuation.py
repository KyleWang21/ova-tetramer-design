#!/usr/bin/env python3
"""Compress E261 into a <=20-mutation cross-backbone consensus and continue it."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import shutil
from collections import defaultdict
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


def row_key(row: dict[str, object]) -> tuple[float, float, float, int]:
    return (
        float(row["weakest_tied_af2_iptm"]),
        float(row["mean_tied_af2_iptm"]),
        float(row["minimum_plddt"]),
        -int(row["n_mutations"]),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--max-source-mutations", type=int, default=24)
    parser.add_argument("--max-anchor-mutations", type=int, default=20)
    parser.add_argument("--min-anchor-mutations", type=int, default=12)
    parser.add_argument("--min-surface-fraction", type=float, default=0.8)
    args = parser.parse_args()

    if (args.out / "SOURCE_SELECTION_FROZEN").exists():
        print(f"{args.out} already frozen")
        return
    reference = read_fasta(args.reference_fasta, args.reference_name)
    surface = set(map(int, re.findall(r"\d+", (args.source / "surface_positions.txt").read_text())))
    glyco = [
        i for i in range(len(reference) - 2)
        if reference[i] == "N" and reference[i + 1] != "P" and reference[i + 2] in "ST"
    ]
    epitope_start = reference.index("SIINFEKL")

    by_context: dict[int, list[dict[str, object]]] = defaultdict(list)
    pattern = re.compile(r"run_cb1_(\d{2})")
    for path in sorted(args.source.glob("run_cb1_*/**/trajectory.tsv")):
        match = pattern.search(str(path))
        if not match:
            continue
        context = int(match.group(1)) // 2
        with path.open() as handle:
            for raw in csv.DictReader(handle, delimiter="\t"):
                sequence = raw["sequence"]
                nmut = int(raw["n_mutations"])
                fraction = float(raw["surface_mutation_fraction"])
                if raw["stage"] != "hard" or float(raw["soft"]) < 0.999:
                    continue
                if nmut > args.max_source_mutations or fraction < args.min_surface_fraction:
                    continue
                if len(sequence) != len(reference):
                    continue
                if any((a == "C") != (b == "C") for a, b in zip(reference, sequence)):
                    continue
                if sequence[epitope_start:epitope_start + 8] != "SIINFEKL":
                    continue
                if any(sequence[i:i + 3] != reference[i:i + 3] for i in glyco):
                    continue
                state1 = float(raw["state1_iptm"])
                state2 = float(raw["state2_iptm"])
                row: dict[str, object] = {
                    **raw,
                    "context": context,
                    "source_trajectory": str(path.relative_to(args.source)),
                    "weakest_tied_af2_iptm": min(state1, state2),
                    "mean_tied_af2_iptm": (state1 + state2) / 2,
                    "minimum_plddt": min(float(raw["state1_plddt"]), float(raw["state2_plddt"])),
                }
                by_context[context].append(row)
    missing = sorted(set(range(4)) - set(by_context))
    if missing:
        raise SystemExit(f"no eligible fully discrete E261 state for contexts {missing}")
    selected = [max(by_context[index], key=row_key) for index in range(4)]

    support: dict[tuple[int, str], dict[str, object]] = {}
    for row in selected:
        context = int(row["context"])
        sequence = str(row["sequence"])
        for position, (old, new) in enumerate(zip(reference, sequence), 1):
            if old == new:
                continue
            key = (position, new)
            record = support.setdefault(key, {
                "position": position,
                "reference_aa": old,
                "consensus_aa": new,
                "contexts": set(),
                "weakest_iptm_sum": 0.0,
                "mean_iptm_sum": 0.0,
                "surface": int(position in surface),
            })
            record["contexts"].add(context)
            record["weakest_iptm_sum"] = float(record["weakest_iptm_sum"]) + float(row["weakest_tied_af2_iptm"])
            record["mean_iptm_sum"] = float(record["mean_iptm_sum"]) + float(row["mean_tied_af2_iptm"])

    variants = list(support.values())
    for record in variants:
        record["context_count"] = len(record["contexts"])
    variants.sort(key=lambda record: (
        -int(record["context_count"]), -int(record["surface"]),
        -float(record["weakest_iptm_sum"]), -float(record["mean_iptm_sum"]),
        int(record["position"]), str(record["consensus_aa"]),
    ))

    chosen: list[dict[str, object]] = []
    used_positions: set[int] = set()
    for record in variants:
        position = int(record["position"])
        if position in used_positions:
            continue
        if int(record["context_count"]) < 2 and len(chosen) >= args.min_anchor_mutations:
            continue
        chosen.append(record)
        used_positions.add(position)
        if len(chosen) == args.max_anchor_mutations:
            break
    if len(chosen) < args.min_anchor_mutations:
        for record in variants:
            position = int(record["position"])
            if position in used_positions:
                continue
            chosen.append(record)
            used_positions.add(position)
            if len(chosen) == args.min_anchor_mutations:
                break
    if not args.min_anchor_mutations <= len(chosen) <= args.max_anchor_mutations:
        raise SystemExit(f"consensus mutation count outside target: {len(chosen)}")
    surface_fraction = sum(int(record["surface"]) for record in chosen) / len(chosen)
    if surface_fraction < args.min_surface_fraction:
        raise SystemExit(f"consensus surface fraction {surface_fraction:.3f} below requirement")

    consensus = list(reference)
    for record in chosen:
        consensus[int(record["position"]) - 1] = str(record["consensus_aa"])
    consensus_sequence = "".join(consensus)
    if any((a == "C") != (b == "C") for a, b in zip(reference, consensus_sequence)):
        raise SystemExit("consensus changes native cysteine pattern")
    if consensus_sequence[epitope_start:epitope_start + 8] != "SIINFEKL":
        raise SystemExit("consensus changes SIINFEKL")
    if any(consensus_sequence[i:i + 3] != reference[i:i + 3] for i in glyco):
        raise SystemExit("consensus changes native glycosylation motif")

    args.out.mkdir(parents=True, exist_ok=False)
    (args.out / "design_masks").mkdir()
    with (args.out / "anchors.fasta").open("w") as handle:
        for index in range(1, 5):
            handle.write(f">CONS2-SEED-{index:02d}\n{consensus_sequence}\n")

    with (args.source / "target_manifest.tsv").open() as handle:
        source_rows = list(csv.DictReader(handle, delimiter="\t"))
        fields = list(source_rows[0])
    if len(source_rows) != 4:
        raise SystemExit(f"expected four E261 target contexts, found {len(source_rows)}")
    output_rows = []
    for index, raw in enumerate(source_rows, 1):
        row = dict(raw)
        mask = args.out / "design_masks" / f"cons2-seed-{index:02d}_joint_design.txt"
        shutil.copy2(row["joint_design"], mask)
        row["anchor"] = f"CONS2-SEED-{index:02d}"
        row["anchor_n_mutations"] = str(len(chosen))
        try:
            row["joint_design"] = str(mask.resolve().relative_to(Path.cwd().resolve()))
        except ValueError:
            row["joint_design"] = str(mask)
        row["surrogate_name"] = "E261-CROSSBACKBONE-CONSENSUS"
        row["surrogate_rank"] = "1"
        if "surrogate_predicted_af3" in row:
            row["surrogate_predicted_af3"] = "nan"
        if "surrogate_predicted_protenix" in row:
            row["surrogate_predicted_protenix"] = "nan"
        output_rows.append(row)
    with (args.out / "target_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(output_rows)
    shutil.copy2(args.source / "surface_positions.txt", args.out / "surface_positions.txt")

    chosen_keys = {(int(record["position"]), str(record["consensus_aa"])) for record in chosen}
    consensus_rows = []
    for record in variants:
        contexts = sorted(record["contexts"])
        consensus_rows.append({
            **{key: value for key, value in record.items() if key != "contexts"},
            "contexts": ",".join(map(str, contexts)),
            "selected": int((int(record["position"]), str(record["consensus_aa"])) in chosen_keys),
        })
    with (args.out / "consensus_mutations.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(consensus_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(consensus_rows)
    selection_fields = list(selected[0])
    with (args.out / "context_source_selection.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, selection_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(selected)
    with (args.out / "consensus_audit.tsv").open("w", newline="") as handle:
        row = {
            "n_mutations": len(chosen),
            "surface_mutation_fraction": surface_fraction,
            "mutations": ",".join(
                f"{record['reference_aa']}{record['position']}{record['consensus_aa']}"
                for record in sorted(chosen, key=lambda item: int(item["position"]))
            ),
            "sequence_sha256": hashlib.sha256(consensus_sequence.encode()).hexdigest(),
            "sequence": consensus_sequence,
        }
        writer = csv.DictWriter(handle, list(row), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerow(row)
    (args.out / "SOURCE_SELECTION_FROZEN").touch()
    print(
        f"frozen E265 consensus nmut={len(chosen)} surface={surface_fraction:.3f}; "
        + "contexts=" + ",".join(f"{float(row['weakest_tied_af2_iptm']):.3f}" for row in selected)
    )


if __name__ == "__main__":
    main()
