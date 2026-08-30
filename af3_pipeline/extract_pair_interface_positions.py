#!/usr/bin/env python3
"""Extract homomer positions recurrently contacting specified chain pairs."""

from __future__ import annotations

import argparse
import csv
import pathlib

import numpy as np
from Bio.PDB import PDBParser


def parse_ranges(text: str) -> set[int]:
    output: set[int] = set()
    for item in text.split(","):
        if not item:
            continue
        if "-" in item:
            lo, hi = map(int, item.split("-", 1)); output.update(range(lo, hi + 1))
        else:
            output.add(int(item))
    return output


def xyz(residue) -> np.ndarray:  # noqa: ANN001
    return np.stack([
        np.asarray(atom.coord, dtype=float) for atom in residue
        if atom.element != "H" and not atom.name.startswith("H")
    ])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", type=pathlib.Path, required=True)
    ap.add_argument("--pairs", required=True, help="comma-separated pairs, e.g. A-B,C-D")
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--recurrence", type=float, default=0.75)
    ap.add_argument("--exclude", default="1,12,31,74,121,257-264,368,383")
    args = ap.parse_args()
    excluded = parse_ranges(args.exclude)
    model = next(PDBParser(QUIET=True).get_structure("oligomer", str(args.pdb)).get_models())
    chains = {str(chain.id): chain for chain in model}
    pairs = [tuple(item.split("-", 1)) for item in args.pairs.split(",")]
    if any(a not in chains or b not in chains for a, b in pairs):
        raise SystemExit("requested pair chain is missing")
    residues = {
        chain_id: {int(r.id[1]): xyz(r) for r in chain if r.id[0] == " "}
        for chain_id, chain in chains.items()
    }
    rows = []
    for position in range(1, 387):
        distances = []
        for left, right in pairs:
            right_all = np.concatenate(list(residues[right].values()))
            left_all = np.concatenate(list(residues[left].values()))
            distances.append(float(np.linalg.norm(
                residues[left][position][:, None] - right_all[None], axis=-1
            ).min()))
            distances.append(float(np.linalg.norm(
                residues[right][position][:, None] - left_all[None], axis=-1
            ).min()))
        recurrence = sum(value < args.cutoff for value in distances) / len(distances)
        selected = recurrence >= args.recurrence and position not in excluded
        rows.append({
            "position": position,
            "mean_min_heavy_distance": float(np.mean(distances)),
            "max_min_heavy_distance": float(max(distances)),
            "recurrence": recurrence,
            "excluded": int(position in excluded),
            "selected": int(selected),
        })
    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "pair_interface_positions.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    selected = [row["position"] for row in rows if row["selected"]]
    (args.out / "design_positions.txt").write_text(",".join(map(str, selected)) + "\n")
    print(f"pairs={args.pairs} selected={len(selected)} positions={','.join(map(str, selected))}")


if __name__ == "__main__":
    main()
