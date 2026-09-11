#!/usr/bin/env python3
"""Prepare a compact allele matrix around E298's recurrent R340-E341 clash."""

from __future__ import annotations

import argparse
import csv
import pathlib


NATIVE_CYS = {12, 31, 74, 121, 368, 383}


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta_sequence(path: pathlib.Path, wanted: str) -> str:
    name = None
    chunks: list[str] = []
    for line in path.read_text().splitlines() + [">"]:
        if line.startswith(">"):
            if name == wanted:
                return "".join(chunks)
            name, chunks = line[1:].split()[0], []
        elif name is not None:
            chunks.append(line.strip())
    raise ValueError(f"missing FASTA record {wanted}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=pathlib.Path, required=True)
    parser.add_argument("--system", default="ova_body_c4_02_ms25_r02")
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--surface-positions", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--exclude-glob", action="append", default=[])
    args = parser.parse_args()

    reference = fasta_sequence(args.reference_fasta, args.reference_name)
    source_rows = read_tsv(args.source / "manifest.tsv")
    base = next(row["sequence"] for row in source_rows if row["name"] == args.system)
    surface = {
        int(value)
        for value in args.surface_positions.read_text().replace("\n", ",").split(",")
        if value.strip()
    }
    project = pathlib.Path(__file__).resolve().parents[1]
    excluded: set[str] = set()
    for pattern in args.exclude_glob:
        for path in project.glob(pattern):
            try:
                excluded.update(
                    row.get("sequence", "") for row in read_tsv(path)
                    if len(row.get("sequence", "")) == 386
                )
            except (OSError, csv.Error):
                continue

    # H is observed in an Accepted_AF3 sequence; N is observed in the E284
    # 0.826 AF3 parent.  S345 is the P3 reversion and D347 is the consensus
    # Accepted_AF3 allele.  The matrix separates clash-local and basin-return
    # effects while staying at 12--14 total mutations.
    edit_sets = [
        ((341, "H"),),
        ((341, "N"),),
        ((345, "S"),),
        ((347, "D"),),
        ((341, "H"), (345, "S")),
        ((341, "N"), (345, "S")),
        ((341, "H"), (347, "D")),
        ((341, "N"), (347, "D")),
        ((345, "S"), (347, "D")),
    ]
    rows: list[dict[str, object]] = []
    seen: set[str] = set()
    for edits in edit_sets:
        variant = list(base)
        labels = []
        for position, allele in edits:
            labels.append(f"{variant[position - 1]}{position}{allele}")
            variant[position - 1] = allele
        sequence = "".join(variant)
        if sequence in excluded or sequence in seen:
            continue
        seen.add(sequence)
        mutation_positions = {
            index for index, (old, new) in enumerate(zip(reference, sequence, strict=True), 1)
            if old != new
        }
        nmut = len(mutation_positions)
        surface_fraction = len(mutation_positions & surface) / nmut
        if nmut > 24 or surface_fraction < 0.80:
            continue
        if {index for index, aa in enumerate(sequence, 1) if aa == "C"} != NATIVE_CYS:
            continue
        if sequence[257:265] != "SIINFEKL":
            continue
        rows.append({
            "name": f"E301-REPAIR-{len(rows) + 1:04d}",
            "parent": args.system,
            "parent_models": 20,
            "parent_mean_iptm": 0.820,
            "parent_zero_clash_models": 19,
            "observed_clash": "R340-E341;seed3_sample3",
            "repair_edits": ",".join(labels),
            "n_repairs": len(edits),
            "n_mutations": nmut,
            "surface_mutation_fraction": surface_fraction,
            "sequence": sequence,
        })
    if not rows:
        raise RuntimeError("all E301 repairs were previously tested or failed constraints")
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
    print(f"novel_repairs={len(rows)} excluded_tested={len(excluded)}")


if __name__ == "__main__":
    main()
