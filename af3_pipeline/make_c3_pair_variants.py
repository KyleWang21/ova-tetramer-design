#!/usr/bin/env python3
"""Create targeted variants that disrupt C3-exclusive AF3 contact pairs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


VARIANTS = [
    (),
    ((210, "E"),),
    ((342, "E"),),
    ((210, "E"), (342, "E")),
    ((95, "D"), (210, "E"), (342, "E")),
    ((98, "H"), (210, "E"), (342, "E")),
    ((98, "R"), (210, "E"), (342, "E")),
    ((183, "A"), (210, "E"), (342, "E")),
    ((95, "D"), (98, "H"), (210, "E"), (342, "E")),
    ((95, "E"), (98, "R"), (210, "D"), (342, "K")),
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--system", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(args.manifest.open(), delimiter="\t"))
    source = next(row for row in rows if row["name"] == args.system)
    native = source["sequence"]
    positions = [95, 98, 183, 210, 342]
    output = []
    for index, changes in enumerate(VARIANTS, 1):
        sequence = list(native)
        for position, amino_acid in changes:
            sequence[position - 1] = amino_acid
        sequence = "".join(sequence)
        mutations = [
            f"{native[position - 1]}{position}{sequence[position - 1]}"
            for position in positions if sequence[position - 1] != native[position - 1]
        ]
        output.append({
            "variant": index,
            "mutations": ",".join(mutations) or "base",
            "n_mutations": len(mutations),
            "interface_sequence": "".join(sequence[position - 1] for position in positions),
            "full_sequence": sequence,
        })
    fields = list(output[0])
    with (args.out / "c3_pair_variants.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output)
    with (args.out / "c3_pair_variants.fasta").open("w") as handle:
        for row in output:
            handle.write(f">C3PAIR-{row['variant']:02d}|mut={row['mutations']}\n{row['full_sequence']}\n")
    print(f"wrote {len(output)} variants from {args.system}")


if __name__ == "__main__":
    main()
