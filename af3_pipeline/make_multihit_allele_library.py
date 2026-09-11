#!/usr/bin/env python3
"""Enumerate observed alleles across several AF3 seed-1 hit sequences."""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import pathlib


NATIVE_CYS = {12, 31, 74, 121, 368, 383}
EPITOPE = set(range(258, 266))


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta_records(path: pathlib.Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line in path.read_text().splitlines() + [">"]:
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:].split()[0] if line[1:].strip() else None, []
        elif header is not None:
            chunks.append(line.strip())
    return records


def fasta_sequence(path: pathlib.Path, name: str) -> str:
    for header, sequence in fasta_records(path):
        if header == name:
            return sequence
    raise ValueError(f"missing FASTA record {name!r}")


def glyco_sites(sequence: str) -> set[int]:
    return {
        index + 1 for index in range(len(sequence) - 2)
        if sequence[index] == "N" and sequence[index + 1] != "P" and sequence[index + 2] in "ST"
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--system", action="append", required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--surface-positions", type=pathlib.Path, required=True)
    parser.add_argument("--exclude-tsv", type=pathlib.Path, action="append", default=[])
    parser.add_argument("--exclude-fasta", type=pathlib.Path, action="append", default=[])
    parser.add_argument("--exclude-glob", action="append", default=[])
    parser.add_argument(
        "--include-reference-position", type=int, action="append", default=[],
        help=(
            "add the P3-13R residue as an allowed allele at this 1-based "
            "position; intended for AF3 clash-directed local expansion"
        ),
    )
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--max-library", type=int, default=100000)
    parser.add_argument("--name-prefix", default="MULTIHIT")
    args = parser.parse_args()

    reference = fasta_sequence(args.reference_fasta, args.reference_name)
    manifest = {
        row.get("name") or row.get("candidate", ""): row
        for row in read_tsv(args.manifest)
    }
    missing = set(args.system) - set(manifest)
    if missing:
        raise ValueError(f"systems absent from manifest: {sorted(missing)}")
    parents = [(name, manifest[name]["sequence"]) for name in args.system]
    if any(len(sequence) != len(reference) for _, sequence in parents):
        raise ValueError("parent/reference length mismatch")
    requested_reference_positions = set(args.include_reference_position)
    if any(position < 1 or position > len(reference) for position in requested_reference_positions):
        raise ValueError("include-reference-position is outside the reference sequence")
    variable = sorted(set([
        position
        for position, residues in enumerate(zip(*(sequence for _, sequence in parents)), 1)
        if len(set(residues)) > 1
    ]) | requested_reference_positions)
    alleles = {
        position: sorted(
            {sequence[position - 1] for _, sequence in parents}
            | ({reference[position - 1]} if position in requested_reference_positions else set())
        )
        for position in variable
    }
    space = math.prod(len(alleles[position]) for position in variable)
    if space > args.max_library:
        raise ValueError(f"refusing allele space of {space} sequences")
    surface = {
        int(value) for value in args.surface_positions.read_text().replace(",", " ").split()
    }
    excluded = {sequence for _, sequence in parents}
    for path in args.exclude_tsv:
        excluded.update(row.get("sequence", "") for row in read_tsv(path))
    for path in args.exclude_fasta:
        excluded.update(sequence for _, sequence in fasta_records(path))
    project = pathlib.Path(__file__).resolve().parents[1]
    for pattern in args.exclude_glob:
        for path in project.glob(pattern):
            try:
                excluded.update(
                    row.get("sequence", "") for row in read_tsv(path)
                    if len(row.get("sequence", "")) == len(reference)
                )
            except (OSError, csv.Error):
                continue
    reference_glyco = glyco_sites(reference)

    rows: list[dict[str, object]] = []
    rejection: dict[str, int] = {}
    for choice in itertools.product(*(alleles[position] for position in variable)):
        sequence = list(parents[0][1])
        for position, residue in zip(variable, choice):
            sequence[position - 1] = residue
        sequence = "".join(sequence)
        changed = {
            position for position, (old, new) in enumerate(zip(reference, sequence), 1)
            if old != new
        }
        reason = ""
        if sequence in excluded:
            reason = "excluded"
        elif not changed or len(changed) > 24:
            reason = "mutation_budget"
        elif {position for position, aa in enumerate(sequence, 1) if aa == "C"} != NATIVE_CYS:
            reason = "cysteine"
        elif changed & EPITOPE:
            reason = "SIINFEKL"
        elif glyco_sites(sequence) != reference_glyco:
            reason = "glycosylation"
        else:
            surface_fraction = len(changed & surface) / len(changed)
            if surface_fraction < 0.80:
                reason = "surface_fraction"
        if reason:
            rejection[reason] = rejection.get(reason, 0) + 1
            continue
        rows.append({
            "name": f"{args.name_prefix}-{len(rows) + 1:04d}",
            "parents": ",".join(args.system),
            "variable_positions": ",".join(map(str, variable)),
            "chosen_alleles": ",".join(
                f"{position}{sequence[position - 1]}" for position in variable
            ),
            "n_mutations": len(changed),
            "surface_mutation_fraction": len(changed & surface) / len(changed),
            "sequence": sequence,
        })
    if not rows:
        raise SystemExit("no novel valid multi-hit allele sequence")
    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "allele_library.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    with (args.out / "allele_library.fasta").open("w") as handle:
        for row in rows:
            handle.write(
                f">{row['name']} nmut={row['n_mutations']} alleles={row['chosen_alleles']}\n"
                f"{row['sequence']}\n"
            )
    with (args.out / "rejection_counts.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["reason", "count"]); writer.writerows(sorted(rejection.items()))
    print(
        f"parents={len(parents)} variable_positions={variable} allele_space={space} "
        f"novel_valid={len(rows)}"
    )


if __name__ == "__main__":
    main()
