#!/usr/bin/env python3
"""Intersect homomer-position-normalized AF3 and Protenix ProLIF fingerprints."""

from __future__ import annotations

import argparse
import csv
import pathlib
from collections import defaultdict


NONVDW = {"Hydrophobic", "HydrogenBond", "SaltBridge", "PiStacking", "CationPi"}


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def normalized_key(row: dict[str, str]) -> tuple[int, int, str]:
    positions = sorted((int(row["position_a"]), int(row["position_b"])))
    return positions[0], positions[1], row["category"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--unique", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()
    by_candidate_system: dict[tuple[str, str], set[tuple[int, int, str]]] = defaultdict(set)
    for row in read_tsv(args.unique):
        if row["category"] not in NONVDW:
            continue
        if row["source"].endswith("__AF3"):
            system = "AF3"
        elif row["source"].endswith("__Protenix"):
            system = "Protenix-v2"
        else:
            continue
        by_candidate_system[(row["candidate"], system)].add(normalized_key(row))

    candidates = sorted({candidate for candidate, _ in by_candidate_system})
    rows = []
    for candidate in candidates:
        af3 = by_candidate_system[(candidate, "AF3")]
        protenix = by_candidate_system[(candidate, "Protenix-v2")]
        shared = af3 & protenix
        union = af3 | protenix
        rows.append({
            "candidate": candidate,
            "af3_nonvdw_position_type_count": len(af3),
            "protenix_nonvdw_position_type_count": len(protenix),
            "shared_nonvdw_position_type_count": len(shared),
            "nonvdw_position_type_jaccard": len(shared) / len(union) if union else 0.0,
            "shared_nonvdw_position_types": ";".join(
                f"{left}-{right}:{category}" for left, right, category in sorted(shared)
            ),
        })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    print(f"wrote normalized cross-model fingerprints for {len(rows)} candidates")


if __name__ == "__main__":
    main()
