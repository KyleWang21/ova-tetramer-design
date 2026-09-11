#!/usr/bin/env python3
"""Enumerate AF3-clash-directed repairs using alleles from Accepted_AF3 sequences."""

from __future__ import annotations

import argparse
import csv
import itertools
import pathlib


NATIVE_CYS = {12, 31, 74, 121, 368, 383}


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta_records(path: pathlib.Path) -> list[tuple[str, str]]:
    output: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    for line in path.read_text().splitlines() + [">"]:
        if line.startswith(">"):
            if name is not None:
                output.append((name, "".join(chunks)))
            name, chunks = line[1:], []
        elif name is not None:
            chunks.append(line.strip())
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", type=pathlib.Path, required=True)
    parser.add_argument("--strict-summary", type=pathlib.Path, required=True)
    parser.add_argument("--clash-recurrence", type=pathlib.Path, required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--accepted-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--surface-positions", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--max-parents", type=int, default=20)
    parser.add_argument("--max-edits", type=int, default=3)
    parser.add_argument("--max-mutations", type=int, default=24)
    parser.add_argument("--min-surface-fraction", type=float, default=0.80)
    parser.add_argument("--name-prefix", default="E289-REPAIR")
    parser.add_argument("--exclude-glob", action="append", default=[])
    args = parser.parse_args()

    reference = next(
        sequence for name, sequence in fasta_records(args.reference_fasta)
        if name.split()[0] == args.reference_name
    )
    accepted = [sequence for _, sequence in fasta_records(args.accepted_fasta)]
    if len(reference) != 386 or not accepted or any(len(sequence) != 386 for sequence in accepted):
        raise ValueError("expected 386-aa reference and Accepted_AF3 sequences")
    surface = {
        int(value)
        for value in args.surface_positions.read_text().replace("\n", ",").split(",")
        if value.strip()
    }
    manifest = {row["name"]: row for row in read_tsv(args.experiment / "manifest.tsv")}
    project = pathlib.Path(__file__).resolve().parents[1]
    excluded_sequences: set[str] = set()
    for pattern in args.exclude_glob:
        for path in project.glob(pattern):
            try:
                excluded_sequences.update(
                    row.get("sequence", "") for row in read_tsv(path)
                    if len(row.get("sequence", "")) == 386
                )
            except (OSError, csv.Error):
                continue
    strict = [
        row for row in read_tsv(args.strict_summary)
        if int(row["valid_c4_network_models"]) == 5
        and int(row["two_interfaces_design_ge3_models"]) == 5
        and int(row["zero_ca_overlap_models"]) == 5
        and int(row["native_ss_complete_models"]) == 5
        and int(row["zero_interchain_ss_models"]) == 5
        and float(row["mean_iptm"]) >= 0.80
    ]
    strict.sort(key=lambda row: (
        -float(row["mean_iptm"]), -int(row["zero_atomic_clash_models"]), row["system"]
    ))
    strict = strict[: args.max_parents]
    clash_positions: dict[str, set[int]] = {}
    for row in read_tsv(args.clash_recurrence):
        clash_positions.setdefault(row["system"], set()).update(
            map(int, row["residue_pair"].split("-"))
        )

    records: dict[str, dict[str, object]] = {}
    for parent_rank, parent in enumerate(strict, 1):
        name = parent["system"]
        sequence = manifest[name]["sequence"]
        direct = sorted(
            position for position in clash_positions.get(name, set())
            if sequence[position - 1] != reference[position - 1]
        )
        if direct:
            implicated = direct
            repair_basis = "direct_clash_mutation"
        else:
            mutations = {
                index for index, (old, new) in enumerate(zip(reference, sequence, strict=True), 1)
                if old != new
            }
            implicated = sorted(
                position for position in mutations
                if any(abs(position - clash) <= 4 for clash in clash_positions.get(name, set()))
            )
            repair_basis = "nearby_mutation_for_native_clash"
        if not implicated:
            continue
        allele_options = {}
        for position in implicated:
            values = {reference[position - 1]}
            values.update(candidate[position - 1] for candidate in accepted)
            values.discard(sequence[position - 1])
            values.discard("C")
            if values:
                allele_options[position] = sorted(values)
        positions = sorted(allele_options)
        for edit_count in range(1, min(args.max_edits, len(positions)) + 1):
            for chosen in itertools.combinations(positions, edit_count):
                for alleles in itertools.product(*(allele_options[position] for position in chosen)):
                    variant = list(sequence)
                    edits = []
                    for position, allele in zip(chosen, alleles, strict=True):
                        edits.append(f"{variant[position - 1]}{position}{allele}")
                        variant[position - 1] = allele
                    variant_sequence = "".join(variant)
                    if variant_sequence in excluded_sequences:
                        continue
                    mutation_positions = {
                        index for index, (old, new) in enumerate(
                            zip(reference, variant_sequence, strict=True), 1
                        ) if old != new
                    }
                    nmut = len(mutation_positions)
                    surface_fraction = (
                        len(mutation_positions & surface) / nmut if nmut else 1.0
                    )
                    if nmut > args.max_mutations or surface_fraction < args.min_surface_fraction:
                        continue
                    if {index for index, aa in enumerate(variant_sequence, 1) if aa == "C"} != NATIVE_CYS:
                        continue
                    if variant_sequence[257:265] != "SIINFEKL":
                        continue
                    record = {
                        "name": "pending",
                        "parent": name,
                        "parent_rank": parent_rank,
                        "parent_mean_iptm": float(parent["mean_iptm"]),
                        "parent_zero_clash_models": int(parent["zero_atomic_clash_models"]),
                        "repair_basis": repair_basis,
                        "observed_clash_positions": ",".join(
                            map(str, sorted(clash_positions.get(name, set())))
                        ),
                        "repair_positions": ",".join(map(str, implicated)),
                        "repair_edits": ",".join(edits),
                        "n_repairs": len(edits),
                        "n_mutations": nmut,
                        "surface_mutation_fraction": surface_fraction,
                        "sequence": variant_sequence,
                    }
                    old = records.get(variant_sequence)
                    priority = (
                        float(record["parent_mean_iptm"]),
                        int(record["parent_zero_clash_models"]),
                        -int(record["n_repairs"]),
                    )
                    if old is None or priority > (
                        float(old["parent_mean_iptm"]),
                        int(old["parent_zero_clash_models"]),
                        -int(old["n_repairs"]),
                    ):
                        records[variant_sequence] = record

    rows = sorted(records.values(), key=lambda row: (
        -float(row["parent_mean_iptm"]), -int(row["parent_zero_clash_models"]),
        int(row["n_repairs"]), int(row["n_mutations"]), str(row["sequence"]),
    ))
    if not rows:
        raise RuntimeError("no clash-directed repair sequences generated")
    for index, row in enumerate(rows, 1):
        row["name"] = f"{args.name_prefix}-{index:04d}"
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    with args.out.with_suffix(".fasta").open("w") as handle:
        for row in rows:
            handle.write(
                f">{row['name']} parent={row['parent']} repairs={row['repair_edits']} "
                f"nmut={row['n_mutations']}\n{row['sequence']}\n"
            )
    print(
        f"parents={len(strict)} parents_with_clashes="
        f"{sum(bool(clash_positions.get(row['system'])) for row in strict)} "
        f"excluded_tested={len(excluded_sequences)} unique_repairs={len(rows)}"
    )


if __name__ == "__main__":
    main()
