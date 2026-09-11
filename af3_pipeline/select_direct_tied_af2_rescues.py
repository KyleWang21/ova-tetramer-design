#!/usr/bin/env python3
"""Select trajectory candidates by measured tied-AF2 evidence, not AF3 surrogate rank."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def finite_float(value: str) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if number == number else None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranking", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--top", type=int, default=8)
    parser.add_argument("--max-mutations", type=int, default=24)
    parser.add_argument("--min-surface-fraction", type=float, default=0.8)
    args = parser.parse_args()

    with args.ranking.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    eligible = []
    for row in rows:
        weakest = finite_float(row.get("af2_iptm", ""))
        plddt = finite_float(row.get("af2_plddt", ""))
        pae = finite_float(row.get("af2_interface_pae", ""))
        surface = finite_float(row.get("surface_mutation_fraction", ""))
        if row.get("method") != "joint_tied_af2" or None in (weakest, plddt, pae, surface):
            continue
        if int(row["n_mutations"]) > args.max_mutations or surface < args.min_surface_fraction:
            continue
        eligible.append(row)
    eligible.sort(
        key=lambda row: (
            -float(row["af2_iptm"]),
            -float(row["af2_plddt"]),
            float(row["af2_interface_pae"]),
            int(row["n_mutations"]),
            row["name"],
        )
    )
    selected = eligible[: args.top]
    if len(selected) != args.top:
        raise SystemExit(f"needed {args.top} eligible tied-AF2 rescues, found {len(selected)}")

    args.out.mkdir(parents=True, exist_ok=True)
    out_tsv = args.out / "direct_tied_af2_rescues.tsv"
    with out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=selected[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(selected)
    with (args.out / "direct_tied_af2_rescues.fasta").open("w") as handle:
        for rank, row in enumerate(selected, 1):
            handle.write(
                f">E259-DIRECT-AF2-{rank:02d} source={row['name']} "
                f"nmut={row['n_mutations']} weakest_tied_af2_iptm={float(row['af2_iptm']):.4f}\n"
                f"{row['sequence']}\n"
            )
    print(f"selected {len(selected)} direct tied-AF2 rescue candidates")
    for rank, row in enumerate(selected, 1):
        print(rank, row["name"], row["n_mutations"], row["af2_iptm"], row["source"])


if __name__ == "__main__":
    main()
