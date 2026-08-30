#!/usr/bin/env python3
"""Extract symmetry-recurrent heavy-atom interface positions from a homooligomer PDB."""

from __future__ import annotations

import argparse
import csv
import pathlib

import numpy as np
from Bio.PDB import PDBParser


def residue_xyz(residue) -> np.ndarray:  # noqa: ANN001
    return np.stack([
        np.asarray(atom.coord, dtype=np.float32) for atom in residue
        if atom.element != "H" and not atom.name.startswith("H")
    ])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--cutoff", type=float, default=6.0)
    ap.add_argument("--recurrence", type=float, default=0.75)
    ap.add_argument("--chain-length", type=int, default=386)
    ap.add_argument("--exclude", default="1,12,31,74,121,257-264,368,383")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    excluded: set[int] = set()
    for item in args.exclude.split(","):
        if "-" in item:
            lo, hi = map(int, item.split("-", 1)); excluded.update(range(lo, hi + 1))
        elif item:
            excluded.add(int(item))

    model = next(PDBParser(QUIET=True).get_structure("oligomer", str(args.pdb)).get_models())
    chains = [
        chain for chain in model
        if sum(1 for residue in chain if residue.id[0] == " " and int(residue.id[1]) <= args.chain_length)
        == args.chain_length
    ]
    if len(chains) < 2:
        raise SystemExit("expected at least two full-length protein chains")
    residues = {
        str(chain.id): {
            int(residue.id[1]): residue_xyz(residue) for residue in chain
            if residue.id[0] == " " and int(residue.id[1]) <= args.chain_length
        }
        for chain in chains
    }
    rows = []
    for position in range(1, args.chain_length + 1):
        per_chain = []
        partner_counts = []
        for chain in chains:
            source = residues[str(chain.id)][position]
            best = float("inf"); partners = 0
            for partner in chains:
                if partner is chain:
                    continue
                target = np.concatenate(list(residues[str(partner.id)].values()), axis=0)
                distance = float(np.sqrt(np.square(source[:, None] - target[None]).sum(-1)).min())
                best = min(best, distance)
                partners += int(distance < args.cutoff)
            per_chain.append(best); partner_counts.append(partners)
        recurrence = sum(distance < args.cutoff for distance in per_chain) / len(chains)
        selected = recurrence >= args.recurrence and position not in excluded
        rows.append({
            "position": position,
            "mean_min_heavy_distance": float(np.mean(per_chain)),
            "max_min_heavy_distance": float(max(per_chain)),
            "chain_recurrence": recurrence,
            "mean_partner_count": float(np.mean(partner_counts)),
            "excluded": int(position in excluded),
            "selected": int(selected),
        })
    fields = list(rows[0])
    with (args.out / "interface_positions.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t"); writer.writeheader(); writer.writerows(rows)
    selected = [int(row["position"]) for row in rows if row["selected"]]
    (args.out / "design_positions.txt").write_text(",".join(map(str, selected)) + "\n")
    print(f"chains={len(chains)} selected={len(selected)} positions={','.join(map(str, selected))}")


if __name__ == "__main__":
    main()
