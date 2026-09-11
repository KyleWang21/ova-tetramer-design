#!/usr/bin/env python3
"""Audit sequence mutations against the two intended C4 interface masks."""

from __future__ import annotations

import argparse
import csv
import pathlib


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


def positions(path: pathlib.Path) -> set[int]:
    return {int(value) for value in path.read_text().replace(",", " ").split()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--interface1-mask", type=pathlib.Path, required=True)
    parser.add_argument("--interface2-mask", type=pathlib.Path, required=True)
    parser.add_argument("--interface1-name", default="AC")
    parser.add_argument("--interface2-name", default="AD")
    parser.add_argument("--hotspots", default="99,155,157,182,184,334,336,338")
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    reference = fasta_sequence(args.reference_fasta, args.reference_name)
    mask1 = positions(args.interface1_mask)
    mask2 = positions(args.interface2_mask)
    hotspots = {int(value) for value in args.hotspots.split(",") if value}
    rows = []
    with args.manifest.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sequence = row["sequence"]
            changed = {
                index for index, (old, new) in enumerate(zip(reference, sequence), 1)
                if old != new
            }
            in1 = changed & mask1
            in2 = changed & mask2
            in_hotspots = changed & hotspots
            rows.append({
                "candidate": row.get("name") or row.get("candidate"),
                "n_mutations": len(changed),
                f"{args.interface1_name}_designed_mutation_count": len(in1),
                f"{args.interface1_name}_designed_mutation_positions": ",".join(map(str, sorted(in1))),
                f"{args.interface2_name}_designed_mutation_count": len(in2),
                f"{args.interface2_name}_designed_mutation_positions": ",".join(map(str, sorted(in2))),
                "both_design_masks_ge3": int(len(in1) >= 3 and len(in2) >= 3),
                "original_hotspot_mutation_count": len(in_hotspots),
                "original_hotspot_mutation_positions": ",".join(map(str, sorted(in_hotspots))),
                "sequence": sequence,
            })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    print(
        f"candidates={len(rows)} both_masks_ge3={sum(int(row['both_design_masks_ge3']) for row in rows)} "
        f"AC_min={min(int(row[f'{args.interface1_name}_designed_mutation_count']) for row in rows)} "
        f"AD_min={min(int(row[f'{args.interface2_name}_designed_mutation_count']) for row in rows)}"
    )


if __name__ == "__main__":
    main()
