#!/usr/bin/env python3
"""Build pure-OVA candidates from AF3 success/failure disulfide discriminator pairs."""

from __future__ import annotations

import argparse
import csv
import pathlib


LOCKS = [
    ("S152C_P208C", {152: "C", 208: "C"}),
    ("S99C_Q195C", {99: "C", 195: "C"}),
    ("P94C_Q195C", {94: "C", 195: "C"}),
    ("S99C_Q195C_S152C_P208C", {99: "C", 195: "C", 152: "C", 208: "C"}),
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--parent-manifest", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    parent = next(
        row["sequence"] for row in csv.DictReader(args.parent_manifest.open(), delimiter="\t")
        if row["name"] == "ova_body_c4_04"
    )
    records = []
    for index, (name, changes) in enumerate(LOCKS, 1):
        sequence = list(parent)
        for position, amino_acid in changes.items(): sequence[position - 1] = amino_acid
        sequence = "".join(sequence)
        required = {96, 155, 190, 196, *changes}
        assert len(sequence) == 386 and all(sequence[position - 1] == "C" for position in required)
        records.append((f"OVA-BODY-DLOCK-{index:02d}|{name}", sequence))
    with (args.out / "disulfide_lock_scan.fasta").open("w") as handle:
        for header, sequence in records: handle.write(f">{header}\n{sequence}\n")
    with (args.out / "disulfide_lock_scan.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t"); writer.writerow(["candidate", "sequence"])
        writer.writerows(records)
    print(f"wrote {len(records)} discriminator-lock candidates")


if __name__ == "__main__":
    main()
