#!/usr/bin/env python3
"""Create low-mutation variants by reverting AF3-noninterface parent changes."""

from __future__ import annotations

import argparse
import csv
import itertools
import pathlib


NATIVE_CYS = {12, 31, 74, 121, 368, 383}
EPITOPE = set(range(258, 266))


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta_sequence(path: pathlib.Path, name: str) -> str:
    active = False
    chunks: list[str] = []
    for line in path.read_text().splitlines() + [">"]:
        if line.startswith(">"):
            if active:
                return "".join(chunks)
            active = line[1:].split()[0] == name
        elif active:
            chunks.append(line.strip())
    raise ValueError(f"missing FASTA record {name}")


def glyco_sites(sequence: str) -> set[int]:
    return {
        index + 1 for index in range(len(sequence) - 2)
        if sequence[index] == "N" and sequence[index + 1] != "P" and sequence[index + 2] in "ST"
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--per-model", type=pathlib.Path, required=True)
    parser.add_argument("--system", required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--surface-positions", type=pathlib.Path, required=True)
    parser.add_argument("--exclude-tsv", type=pathlib.Path, action="append", default=[])
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    reference = fasta_sequence(args.reference_fasta, args.reference_name)
    parent_row = next(row for row in read_tsv(args.manifest) if row["name"] == args.system)
    parent = parent_row["sequence"]
    model_rows = [row for row in read_tsv(args.per_model) if row["system"] == args.system]
    if len(model_rows) != 5:
        raise ValueError(f"expected five seed1 models for {args.system}, found {len(model_rows)}")
    mutations = {index for index, (old, new) in enumerate(zip(reference, parent), 1) if old != new}
    recurrence = {
        position: sum(
            position in {int(value) for value in row["design_interface_positions"].split(",") if value}
            for row in model_rows
        )
        for position in mutations
    }
    noninterface = sorted(position for position, count in recurrence.items() if count == 0)
    if not noninterface:
        raise SystemExit("no zero-recurrence mutated positions to revert")
    surface = {
        int(value) for value in args.surface_positions.read_text().replace(",", " ").split()
    }
    excluded = set()
    for path in args.exclude_tsv:
        excluded.update(row.get("sequence", "") for row in read_tsv(path))

    rows = []
    for size in range(1, len(noninterface) + 1):
        for reverted in itertools.combinations(noninterface, size):
            sequence = list(parent)
            for position in reverted:
                sequence[position - 1] = reference[position - 1]
            sequence = "".join(sequence)
            changed = {index for index, (old, new) in enumerate(zip(reference, sequence), 1) if old != new}
            if sequence in excluded or not changed or len(changed) > 24:
                continue
            if {index for index, aa in enumerate(sequence, 1) if aa == "C"} != NATIVE_CYS:
                continue
            if changed & EPITOPE or glyco_sites(sequence) != glyco_sites(reference):
                continue
            surface_fraction = len(changed & surface) / len(changed)
            if surface_fraction < 0.80:
                continue
            rows.append({
                "name": f"{args.system}_REV{len(rows) + 1:03d}",
                "parent": args.system,
                "reverted_positions": ",".join(map(str, reverted)),
                "n_reverted": len(reverted),
                "n_mutations": len(changed),
                "surface_mutation_fraction": surface_fraction,
                "predicted_design_interface_fraction_if_topology_preserved": (
                    sum(count > 0 for position, count in recurrence.items() if position in changed) / len(changed)
                ),
                "sequence": sequence,
            })
    if not rows:
        raise SystemExit("all reversion variants were excluded")
    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "reversion_library.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    with (args.out / "reversion_library.fasta").open("w") as handle:
        for row in rows:
            handle.write(
                f">{row['name']} nmut={row['n_mutations']} reverted={row['reverted_positions']}\n"
                f"{row['sequence']}\n"
            )
    with (args.out / "parent_interface_recurrence.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["position", "parent_mutation", "interface_models_of_5"])
        for position in sorted(mutations):
            writer.writerow([
                position, f"{reference[position - 1]}{position}{parent[position - 1]}", recurrence[position]
            ])
    print(
        f"parent={args.system} mutations={len(mutations)} zero_interface={noninterface} "
        f"novel_valid_variants={len(rows)}"
    )


if __name__ == "__main__":
    main()
