#!/usr/bin/env python3
"""Add the AF3-state-discriminating P94C/Y213C pair to double-ring designs."""

from __future__ import annotations

import argparse
import csv
import pathlib


def mutate(sequence: str, changes: dict[int, str]) -> str:
    chars = list(sequence)
    for position, amino_acid in changes.items():
        chars[position - 1] = amino_acid
    return "".join(chars)


def read_fasta(path: pathlib.Path) -> list[tuple[str, str]]:
    records = []
    header = None
    sequence = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(sequence)))
            header, sequence = line[1:], []
        else:
            sequence.append(line.strip())
    if header is not None:
        records.append((header, "".join(sequence)))
    return records


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--parent-manifest", type=pathlib.Path, required=True)
    ap.add_argument("--ai-fasta", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--ai-top", type=int, default=3)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    parent = next(
        row["sequence"] for row in csv.DictReader(args.parent_manifest.open(), delimiter="\t")
        if row["name"] == "ova_body_c4_04"
    )
    records = [("parent_double_ring", parent), *read_fasta(args.ai_fasta)[:args.ai_top]]
    output = []
    for index, (source, sequence) in enumerate(records, 1):
        if len(sequence) != 386:
            raise SystemExit(f"{source}: expected 386 aa")
        triple = mutate(sequence, {94: "C", 213: "C"})
        required = {94, 96, 155, 190, 196, 213}
        if any(triple[position - 1] != "C" for position in required):
            raise SystemExit(f"{source}: failed to retain all six designed cysteines")
        output.append((f"OVA-BODY-SS3-{index:02d}|source={source}|add=P94C,Y213C", triple))
    with (args.out / "triple_ring_candidates.fasta").open("w") as handle:
        for header, sequence in output:
            handle.write(f">{header}\n{sequence}\n")
    with (args.out / "triple_ring_candidates.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["candidate", "source", "length", "sequence"])
        for index, (header, sequence) in enumerate(output, 1):
            writer.writerow([f"OVA-BODY-SS3-{index:02d}", header, len(sequence), sequence])
    print(f"wrote {len(output)} triple-ring candidates")


if __name__ == "__main__":
    main()
