#!/usr/bin/env python3
"""Report the exact inter-chain heavy-atom pairs below the AF3 2.4 A gate."""

from __future__ import annotations

import argparse
import csv
import itertools
import pathlib
from collections import defaultdict

import numpy as np
from scipy.spatial import cKDTree

from screen_final20_af3_structures import structure_arrays


def read_reference(path: pathlib.Path, name: str) -> str:
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


def write_tsv(path: pathlib.Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", type=pathlib.Path, required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--cutoff", type=float, default=2.4)
    args = parser.parse_args()

    reference = read_reference(args.reference_fasta, args.reference_name)
    manifest = {
        row["name"]: row
        for row in csv.DictReader((args.experiment / "manifest.tsv").open(), delimiter="\t")
    }
    detail: list[dict[str, object]] = []
    model_sets: dict[tuple[str, int, int], set[str]] = defaultdict(set)
    atom_pair_counts: dict[tuple[str, int, int], int] = defaultdict(int)
    residue_names: dict[tuple[str, int, int], set[str]] = defaultdict(set)
    mutation_flags: dict[tuple[str, int, int], bool] = {}

    for cif in sorted(args.experiment.glob("out_s*/**/seed-*/*_model.cif")):
        stem = cif.name.removesuffix("_model.cif")
        system = stem.split("_seed-", 1)[0]
        if system not in manifest:
            continue
        sequence = manifest[system]["sequence"]
        mutations = {
            index for index, (old, new) in enumerate(zip(reference, sequence, strict=True), 1)
            if old != new
        }
        _, chains = structure_arrays(cif)
        model = cif.parent.name
        for chain_a, chain_b in itertools.combinations(sorted(chains), 2):
            left, right = chains[chain_a], chains[chain_b]
            neighbors = cKDTree(left["coords"]).query_ball_tree(
                cKDTree(right["coords"]), args.cutoff
            )
            for atom_a, hits in enumerate(neighbors):
                for atom_b in hits:
                    residue_a = int(left["resids"][atom_a])
                    residue_b = int(right["resids"][atom_b])
                    distance = float(np.linalg.norm(
                        np.asarray(left["coords"])[atom_a] - np.asarray(right["coords"])[atom_b]
                    ))
                    key = (system, min(residue_a, residue_b), max(residue_a, residue_b))
                    model_sets[key].add(model)
                    atom_pair_counts[key] += 1
                    residue_names[key].update({
                        str(left["resnames"][atom_a]), str(right["resnames"][atom_b])
                    })
                    mutation_flags[key] = residue_a in mutations or residue_b in mutations
                    detail.append({
                        "system": system, "model": model,
                        "chain_a": chain_a, "residue_a": residue_a,
                        "resname_a": str(left["resnames"][atom_a]),
                        "atom_a": str(left["names"][atom_a]),
                        "chain_b": chain_b, "residue_b": residue_b,
                        "resname_b": str(right["resnames"][atom_b]),
                        "atom_b": str(right["names"][atom_b]),
                        "distance_A": distance,
                        "mutation_involved": int(mutation_flags[key]),
                        "cif": str(cif),
                    })

    detail.sort(key=lambda row: (
        str(row["system"]), str(row["model"]), float(row["distance_A"]),
        str(row["chain_a"]), int(row["residue_a"]), str(row["atom_a"]),
    ))
    summary = []
    for key in sorted(model_sets, key=lambda item: (
        item[0], -len(model_sets[item]), -atom_pair_counts[item], item[1], item[2]
    )):
        system, residue_a, residue_b = key
        summary.append({
            "system": system,
            "residue_pair": f"{residue_a}-{residue_b}",
            "residue_names": ",".join(sorted(residue_names[key])),
            "models_with_clash": len(model_sets[key]),
            "model_names": ",".join(sorted(model_sets[key])),
            "atom_pair_count": atom_pair_counts[key],
            "mutation_involved": int(mutation_flags[key]),
        })
    write_tsv(args.out / "clash_atom_pairs.tsv", [
        "system", "model", "chain_a", "residue_a", "resname_a", "atom_a",
        "chain_b", "residue_b", "resname_b", "atom_b", "distance_A",
        "mutation_involved", "cif",
    ], detail)
    write_tsv(args.out / "clash_residue_pair_recurrence.tsv", [
        "system", "residue_pair", "residue_names", "models_with_clash",
        "model_names", "atom_pair_count", "mutation_involved",
    ], summary)
    print(f"models={len({(row['system'], row['model']) for row in detail})} "
          f"atom_pairs={len(detail)} residue_pairs={len(summary)}")


if __name__ == "__main__":
    main()
