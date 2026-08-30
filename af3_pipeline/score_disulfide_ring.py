#!/usr/bin/env python3
"""Score intended inter-chain disulfide rings in AF3 pure-OVA outputs."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib

import numpy as np
from Bio.PDB import MMCIFParser


PAIR_BY_INDEX = {
    1: [(94, 213)],
    2: [(96, 190)],
    3: [(155, 196)],
    4: [(96, 190), (155, 196)],
    5: [(96, 190)],
    6: [(155, 196)],
    7: [(96, 190)],
    8: [(155, 196)],
}


def largest_component(chains: list[str], edges: list[tuple[str, str]]) -> int:
    adjacent = {chain: set() for chain in chains}
    for a, b in edges:
        adjacent[a].add(b); adjacent[b].add(a)
    best = 0
    for start in chains:
        seen = {start}; stack = [start]
        while stack:
            current = stack.pop()
            for other in adjacent[current]:
                if other not in seen:
                    seen.add(other); stack.append(other)
        best = max(best, len(seen))
    return best


def score_cif(cif: pathlib.Path, pairs: list[tuple[int, int]], cutoff: float) -> dict[str, object]:
    model = next(MMCIFParser(QUIET=True).get_structure(cif.stem, str(cif)).get_models())
    chains = [str(chain.id) for chain in model]
    sg: dict[tuple[str, int], np.ndarray] = {}
    for chain in model:
        for residue in chain:
            if residue.id[0] == " " and "SG" in residue:
                sg[(str(chain.id), int(residue.id[1]))] = np.asarray(residue["SG"].coord)
    bonds: list[tuple[str, str, int, int, float]] = []
    for left, right in pairs:
        candidates = []
        for chain_a in chains:
            for chain_b in chains:
                if chain_a == chain_b:
                    continue
                key_a, key_b = (chain_a, left), (chain_b, right)
                if key_a not in sg or key_b not in sg:
                    continue
                distance = float(np.linalg.norm(sg[key_a] - sg[key_b]))
                if distance <= cutoff:
                    candidates.append((distance, key_a, key_b))
        used: set[tuple[str, int]] = set()
        for distance, key_a, key_b in sorted(candidates):
            if key_a in used or key_b in used:
                continue
            used.update((key_a, key_b))
            bonds.append((key_a[0], key_b[0], left, right, distance))
    edges = [(bond[0], bond[1]) for bond in bonds]
    expected = len(chains) * len(pairs) if len(chains) > 1 else 0
    all_close: list[tuple[str, str, int, int, float]] = []
    sg_items = sorted(sg.items())
    for index, ((chain_a, pos_a), xyz_a) in enumerate(sg_items):
        for (chain_b, pos_b), xyz_b in sg_items[index + 1:]:
            distance = float(np.linalg.norm(xyz_a - xyz_b))
            if distance <= cutoff:
                all_close.append((chain_a, chain_b, pos_a, pos_b, distance))
    intended_sets = [{left, right} for left, right in pairs]
    native = [bond for bond in all_close
              if bond[0] == bond[1] and {bond[2], bond[3]} == {74, 121}]
    intended_close = [bond for bond in all_close
                      if bond[0] != bond[1] and {bond[2], bond[3]} in intended_sets]
    unexpected = [bond for bond in all_close if bond not in native and bond not in intended_close]
    interchain = [bond for bond in all_close if bond[0] != bond[1]]
    component = largest_component(chains, edges)
    full_ring = len(chains) > 1 and len(bonds) >= expected and component == len(chains)
    return {
        "intended_pairs": ";".join(f"{left}-{right}" for left, right in pairs),
        "expected_ss_bonds": expected,
        "observed_ss_bonds": len(bonds),
        "chain_component": component,
        "full_ss_ring": int(full_ring),
        "native_c74_c121_bonds": len(native),
        "unexpected_ss_bonds": len(unexpected),
        "chemically_clean_ring": int(full_ring and len(native) == len(chains) and not unexpected),
        "interchain_ss_bonds": len(interchain),
        "no_interchain_ss": int(not interchain),
        "noncovalent_chemically_clean": int(len(native) == len(chains) and not interchain),
        "median_ss_distance_a": float(np.median([bond[4] for bond in bonds])) if bonds else float("nan"),
        "bonds_json": json.dumps(bonds),
        "unexpected_bonds_json": json.dumps(unexpected),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--experiment", type=pathlib.Path, required=True)
    ap.add_argument("--cutoff", type=float, default=2.6)
    ap.add_argument("--pairs", help="Use the same comma-separated residue pairs for every system")
    ap.add_argument("--no-designed-pairs", action="store_true",
                    help="Audit a strictly noncovalent design: expect zero inter-chain S-S bonds")
    ap.add_argument("--pairs-map", type=pathlib.Path,
                    help="TSV with system and comma-separated pairs columns")
    args = ap.parse_args()
    pairs_map = {}
    if args.pairs_map:
        for row in csv.DictReader(args.pairs_map.open(), delimiter="\t"):
            pairs_map[row["system"]] = [
                tuple(int(value) for value in item.split("-", 1))
                for item in row["pairs"].split(",")
            ]
    rows: list[dict[str, object]] = []
    for summary in sorted(args.experiment.glob("out_s*/*/seed-*/*_summary_confidences.json")):
        system = summary.parent.parent.name
        if args.no_designed_pairs:
            pairs = []
        elif system in pairs_map:
            pairs = pairs_map[system]
        elif args.pairs:
            pairs = [tuple(int(value) for value in item.split("-", 1)) for item in args.pairs.split(",")]
        else:
            index = int(system.rsplit("_", 1)[1])
            pairs = PAIR_BY_INDEX[index]
        cif = pathlib.Path(str(summary).replace("_summary_confidences.json", "_model.cif"))
        confidence = json.loads(summary.read_text())
        rows.append({
            "system": system,
            "model": summary.parent.name,
            "iptm": confidence.get("iptm") or 0.0,
            "has_clash": confidence.get("has_clash", 0.0),
            **score_cif(cif, pairs, args.cutoff),
            "cif_path": str(cif),
        })
    if not rows:
        raise SystemExit("no AF3 seed/sample outputs found")
    with (args.experiment / "disulfide_model_scores.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    grouped = []
    for system in sorted({str(row["system"]) for row in rows}):
        selected = [row for row in rows if row["system"] == system]
        grouped.append({
            "system": system,
            "n_models": len(selected),
            "full_ss_ring_models": sum(int(row["full_ss_ring"]) for row in selected),
            "full_ss_ring_fraction": sum(int(row["full_ss_ring"]) for row in selected) / len(selected),
            "chemically_clean_models": sum(int(row["chemically_clean_ring"]) for row in selected),
            "chemically_clean_fraction": sum(int(row["chemically_clean_ring"]) for row in selected) / len(selected),
            "median_observed_ss_bonds": float(np.median([int(row["observed_ss_bonds"]) for row in selected])),
            "median_native_c74_c121_bonds": float(np.median([int(row["native_c74_c121_bonds"]) for row in selected])),
            "median_unexpected_ss_bonds": float(np.median([int(row["unexpected_ss_bonds"]) for row in selected])),
            "zero_interchain_ss_models": sum(int(row["no_interchain_ss"]) for row in selected),
            "zero_interchain_ss_fraction": sum(int(row["no_interchain_ss"]) for row in selected) / len(selected),
            "noncovalent_chemically_clean_models": sum(int(row["noncovalent_chemically_clean"]) for row in selected),
            "noncovalent_chemically_clean_fraction": sum(int(row["noncovalent_chemically_clean"]) for row in selected) / len(selected),
            "median_iptm": float(np.median([float(row["iptm"]) for row in selected])),
        })
    with (args.experiment / "disulfide_system_scores.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(grouped[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(grouped)
    for row in grouped:
        print(
            row["system"], f"ring={row['full_ss_ring_models']}/{row['n_models']}",
            f"bonds={row['median_observed_ss_bonds']}", f"ipTM={row['median_iptm']:.3f}",
        )


if __name__ == "__main__":
    main()
