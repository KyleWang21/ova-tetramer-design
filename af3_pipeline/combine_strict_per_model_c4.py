#!/usr/bin/env python3
"""Combine one strict AF3 seed-1 table with independent seeds 2--5."""

from __future__ import annotations

import argparse
import csv
import pathlib
import re
import statistics

from run_final20_failclosed_filter import glyco_sites, reference_rsasa


NATIVE_CYS = {12, 31, 74, 121, 368, 383}
EPITOPE = set(range(258, 266))


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta_sequence(path: pathlib.Path, name: str) -> str:
    active = False
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if active:
                break
            active = line[1:].startswith(name)
        elif active:
            chunks.append(line.strip())
    sequence = "".join(chunks)
    if not sequence:
        raise ValueError(f"FASTA record {name!r} not found in {path}")
    return sequence


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed1-per-model", type=pathlib.Path, required=True)
    parser.add_argument("--seed1-system", required=True)
    parser.add_argument("--multiseed-experiment", type=pathlib.Path, required=True)
    parser.add_argument("--candidate-system", required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--reference-cif", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    manifest_rows = read_tsv(args.multiseed_experiment / "manifest.tsv")
    candidates = {
        row["sequence"] for row in manifest_rows if row["name"] == args.candidate_system
    }
    if len(candidates) != 1:
        raise SystemExit(f"expected one sequence for {args.candidate_system}; found {len(candidates)}")
    sequence = candidates.pop()
    reference = fasta_sequence(args.reference_fasta, args.reference_name)
    if len(sequence) != len(reference):
        raise SystemExit("candidate/reference length mismatch")

    rows = [
        dict(row, validation_seed=1, source_table=str(args.seed1_per_model))
        for row in read_tsv(args.seed1_per_model)
        if row["system"] == args.seed1_system
    ]
    for seed in range(2, 6):
        path = args.multiseed_experiment / f"seed{seed}_strict" / "seed1_strict_per_model.tsv"
        if not path.exists():
            raise SystemExit(f"missing strict table: {path}")
        rows.extend(
            dict(row, validation_seed=seed, source_table=str(path))
            for row in read_tsv(path)
            if row["system"] == args.candidate_system
        )

    observed = set()
    for row in rows:
        match = re.search(r"seed-(\d+)_sample-(\d+)", row["model"])
        if not match:
            raise SystemExit(f"cannot parse seed/sample from {row['model']}")
        key = (int(match.group(1)), int(match.group(2)))
        if key in observed:
            raise SystemExit(f"duplicate model {key}")
        observed.add(key)
    expected = {(seed, sample) for seed in range(1, 6) for sample in range(5)}
    complete = observed == expected

    changed = [index for index, (old, new) in enumerate(zip(reference, sequence), 1) if old != new]
    rsasa = reference_rsasa(args.reference_cif, "A", reference)
    surface_fraction = (
        sum(rsasa.get(index, 0.0) >= 0.20 for index in changed) / len(changed)
        if changed else 0.0
    )
    sequence_constraints = (
        len(changed) <= 24
        and surface_fraction >= 0.80
        and {index for index, aa in enumerate(sequence, 1) if aa == "C"} == NATIVE_CYS
        and all(sequence[index - 1] == reference[index - 1] for index in EPITOPE)
        and glyco_sites(sequence) == glyco_sites(reference)
    )
    iptm = [float(row["iptm"]) for row in rows]
    strict_geometry = all(
        int(row["atomic_clash_count_2p4"]) == 0
        and int(row["severe_ca_overlap_count_2p5"]) == 0
        and int(row["valid_c4_network"]) == 1
        and int(row["native_c74_c121_bonds"]) == 4
        and int(row["interchain_ss_bonds"]) == 0
        and int(row["two_interfaces_each_design_ge3"]) == 1
        and float(row["max_chain_ca_rmsd_to_1ova"]) <= 2.0
        for row in rows
    )
    per_seed_strict = all(
        len(seed_rows := [row for row in rows if int(row["validation_seed"]) == seed]) == 5
        and statistics.fmean(float(row["iptm"]) for row in seed_rows) >= 0.80
        and min(float(row["iptm"]) for row in seed_rows) >= 0.75
        for seed in range(1, 6)
    )
    passes = bool(
        complete and len(rows) == 25 and statistics.fmean(iptm) >= 0.80
        and sequence_constraints and strict_geometry and per_seed_strict
    )
    summary = {
        "candidate": args.candidate_system,
        "seed1_system": args.seed1_system,
        "n_models": len(rows),
        "complete_5seeds_x_5samples": int(complete),
        "n_mutations": len(changed),
        "surface_mutation_fraction": surface_fraction,
        "mean_iptm": statistics.fmean(iptm),
        "median_iptm": statistics.median(iptm),
        "min_iptm": min(iptm),
        "max_iptm": max(iptm),
        "zero_atomic_clash_models": sum(int(row["atomic_clash_count_2p4"]) == 0 for row in rows),
        "zero_ca_overlap_models": sum(int(row["severe_ca_overlap_count_2p5"]) == 0 for row in rows),
        "valid_c4_network_models": sum(int(row["valid_c4_network"]) for row in rows),
        "native_ss_complete_models": sum(int(row["native_c74_c121_bonds"]) == 4 for row in rows),
        "zero_interchain_ss_models": sum(int(row["interchain_ss_bonds"]) == 0 for row in rows),
        "two_interfaces_design_ge3_models": sum(int(row["two_interfaces_each_design_ge3"]) for row in rows),
        "gate_sequence_constraints": int(sequence_constraints),
        "gate_per_seed_strict": int(per_seed_strict),
        "gate_all_model_geometry": int(strict_geometry),
        "accepted_af3": int(passes),
        "mutations": ",".join(
            f"{reference[index - 1]}{index}{sequence[index - 1]}" for index in changed
        ),
        "sequence": sequence,
    }
    args.out.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out / "combined_per_model.tsv", rows)
    write_tsv(args.out / "combined_25model_summary.tsv", [summary])
    print(
        f"{args.candidate_system} models={len(rows)} mean={summary['mean_iptm']:.4f} "
        f"min={summary['min_iptm']:.4f} zero_clash={summary['zero_atomic_clash_models']} "
        f"accepted_af3={summary['accepted_af3']}"
    )


if __name__ == "__main__":
    main()
