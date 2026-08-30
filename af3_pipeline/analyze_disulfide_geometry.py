#!/usr/bin/env python3
"""Measure candidate directed C4 disulfide geometries across passing AF3 models."""

from __future__ import annotations

import argparse
import csv
import itertools
import pathlib

import numpy as np
from Bio.PDB import MMCIFParser


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", type=pathlib.Path, action="append", required=True)
    ap.add_argument("--system", action="append", required=True)
    ap.add_argument("--pairs", default="94-213,96-190,155-196")
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--iptm-min", type=float, default=0.50)
    ap.add_argument("--weak-min", type=float, default=0.45)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    systems = set(args.system)
    score_rows = []
    for path in args.scores:
        score_rows.extend(
            row for row in csv.DictReader(path.open(), delimiter="\t")
            if row["system"] in systems
            and float(row["iptm"]) >= args.iptm_min
            and float(row["min_incident_iptm"]) >= args.weak_min
            and float(row["has_clash"]) == 0
        )
    pairs = [tuple(int(value) for value in item.split("-", 1)) for item in args.pairs.split(",")]
    parser = MMCIFParser(QUIET=True)
    rows = []
    for score in score_rows:
        cif = pathlib.Path(score["cif_path"])
        model = next(parser.get_structure(cif.stem, str(cif)).get_models())
        chains = [str(chain.id) for chain in model]
        cb = {}
        for chain in model:
            for residue in chain:
                if residue.id[0] == " " and "CA" in residue:
                    atom = residue["CB"] if "CB" in residue else residue["CA"]
                    cb[(str(chain.id), int(residue.id[1]))] = np.asarray(atom.coord)
        anchor = chains[0]
        cycles = [(anchor, *order) for order in itertools.permutations(chains[1:])]
        for left, right in pairs:
            choices = []
            for cycle in cycles:
                distances = [
                    float(np.linalg.norm(cb[(cycle[index], left)] - cb[(cycle[(index + 1) % 4], right)]))
                    for index in range(4)
                ]
                choices.append((max(distances), float(np.mean(distances)), cycle, distances))
            maximum, mean, cycle, distances = min(choices)
            rows.append({
                "system": score["system"],
                "model": score["model"],
                "pair": f"{left}-{right}",
                "best_cycle": "-".join(cycle),
                "mean_cb_distance_a": mean,
                "max_cb_distance_a": maximum,
                "std_cb_distance_a": float(np.std(distances)),
                "all_four_within_6a": int(maximum <= 6.0),
                "cif_path": str(cif),
            })
    with (args.out / "disulfide_pair_geometry.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    for pair in sorted({row["pair"] for row in rows}):
        selected = [row for row in rows if row["pair"] == pair]
        means = [float(row["mean_cb_distance_a"]) for row in selected]
        maxima = [float(row["max_cb_distance_a"]) for row in selected]
        print(
            pair, f"n={len(selected)}", f"mean_median={np.median(means):.3f}",
            f"max_median={np.median(maxima):.3f}",
            f"all4<=6A={sum(int(row['all_four_within_6a']) for row in selected)}/{len(selected)}",
        )


if __name__ == "__main__":
    main()
