#!/usr/bin/env python3
"""Enumerate the sequence bridge between two experimentally untested AF3 hits."""

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
    parser.add_argument("--system-a", required=True)
    parser.add_argument("--system-b", required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--surface-positions", type=pathlib.Path, required=True)
    parser.add_argument("--exclude-tsv", type=pathlib.Path, action="append", default=[])
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    reference = fasta_sequence(args.reference_fasta, args.reference_name)
    manifest = {row["name"]: row for row in read_tsv(args.manifest)}
    sequence_a = manifest[args.system_a]["sequence"]
    sequence_b = manifest[args.system_b]["sequence"]
    different = [index for index, (a, b) in enumerate(zip(sequence_a, sequence_b), 1) if a != b]
    if len(different) > 16:
        raise ValueError(f"refusing 2^{len(different)} crossover space")
    surface = {
        int(value) for value in args.surface_positions.read_text().replace(",", " ").split()
    }
    excluded = {sequence_a, sequence_b}
    for path in args.exclude_tsv:
        excluded.update(row.get("sequence", "") for row in read_tsv(path))
    rows = []
    for mask in itertools.product((0, 1), repeat=len(different)):
        sequence = list(sequence_a)
        chosen_b = []
        for use_b, position in zip(mask, different):
            if use_b:
                sequence[position - 1] = sequence_b[position - 1]
                chosen_b.append(position)
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
            "name": f"{args.system_a}_{args.system_b}_X{len(rows) + 1:03d}",
            "parent_a": args.system_a,
            "parent_b": args.system_b,
            "crossover_positions": ",".join(map(str, different)),
            "crossover_b_positions": ",".join(map(str, chosen_b)),
            "n_mutations": len(changed),
            "surface_mutation_fraction": surface_fraction,
            "sequence": sequence,
        })
    if not rows:
        raise SystemExit("no novel valid crossover sequences")
    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "crossover_library.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    with (args.out / "crossover_library.fasta").open("w") as handle:
        for row in rows:
            handle.write(
                f">{row['name']} nmut={row['n_mutations']} b_positions={row['crossover_b_positions']}\n"
                f"{row['sequence']}\n"
            )
    print(
        f"parents={args.system_a},{args.system_b} differing_positions={different} "
        f"novel_valid={len(rows)}"
    )


if __name__ == "__main__":
    main()
