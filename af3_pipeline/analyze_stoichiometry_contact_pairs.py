#!/usr/bin/env python3
"""Compare residue-pair contact recurrence between C4 and C3 AF3 ensembles."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import MMCIFParser


def read_rows(path: Path, system: str | None, iptm_min: float, weak_min: float) -> list[dict[str, str]]:
    rows = []
    for row in csv.DictReader(path.open(), delimiter="\t"):
        if system and row.get("system") != system:
            continue
        weak = float(row.get("min_incident_iptm", row.get("weakest_chain", 0)))
        if float(row["iptm"]) >= iptm_min and weak >= weak_min and int(float(row.get("has_clash", 0))) == 0:
            rows.append(row)
    if not rows:
        raise SystemExit(f"no passing rows in {path}")
    return rows


def chain_atoms(chain, end: int) -> tuple[list[int], np.ndarray]:  # noqa: ANN001
    positions, coordinates = [], []
    for residue in chain:
        if residue.id[0] != " " or int(residue.id[1]) > end or "CA" not in residue:
            continue
        atom = residue["CB"] if "CB" in residue else residue["CA"]
        positions.append(int(residue.id[1]))
        coordinates.append(np.asarray(atom.coord, dtype=np.float32))
    return positions, np.stack(coordinates)


def pair_recurrence(rows: list[dict[str, str]], end: int, cutoff: float) -> dict[tuple[int, int], float]:
    parser = MMCIFParser(QUIET=True)
    hits: defaultdict[tuple[int, int], int] = defaultdict(int)
    for row in rows:
        structure = parser.get_structure("model", row["cif_path"])
        chains = list(next(structure.get_models()).get_chains())
        atoms = {str(chain.id): chain_atoms(chain, end) for chain in chains}
        model_pairs: set[tuple[int, int]] = set()
        for index, chain_a in enumerate(chains):
            posa, xyza = atoms[str(chain_a.id)]
            for chain_b in chains[index + 1:]:
                posb, xyzb = atoms[str(chain_b.id)]
                delta = xyza[:, None, :] - xyzb[None, :, :]
                ai, bj = np.nonzero(np.einsum("ijk,ijk->ij", delta, delta) < cutoff ** 2)
                model_pairs.update(tuple(sorted((posa[int(i)], posb[int(j)]))) for i, j in zip(ai, bj))
        for pair in model_pairs:
            hits[pair] += 1
    return {pair: count / len(rows) for pair, count in hits.items()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c4-scores", type=Path, required=True)
    parser.add_argument("--c3-scores", type=Path, required=True)
    parser.add_argument("--c4-system")
    parser.add_argument("--c3-system", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--ova-end", type=int, default=386)
    parser.add_argument("--contact-cutoff", type=float, default=8.0)
    parser.add_argument("--iptm-min", type=float, default=0.50)
    parser.add_argument("--weak-min", type=float, default=0.45)
    parser.add_argument("--c3-min", type=float, default=0.70)
    parser.add_argument("--c4-max", type=float, default=0.20)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    c4_rows = read_rows(args.c4_scores, args.c4_system, args.iptm_min, args.weak_min)
    c3_rows = read_rows(args.c3_scores, args.c3_system, args.iptm_min, args.weak_min)
    c4 = pair_recurrence(c4_rows, args.ova_end, args.contact_cutoff)
    c3 = pair_recurrence(c3_rows, args.ova_end, args.contact_cutoff)
    pairs = sorted(set(c4) | set(c3))
    records = [{
        "position_i": pair[0],
        "position_j": pair[1],
        "c4_recurrence": c4.get(pair, 0.0),
        "c3_recurrence": c3.get(pair, 0.0),
        "c3_minus_c4": c3.get(pair, 0.0) - c4.get(pair, 0.0),
        "c3_exclusive": int(c3.get(pair, 0.0) >= args.c3_min and c4.get(pair, 0.0) <= args.c4_max),
    } for pair in pairs]
    records.sort(key=lambda row: (row["c3_exclusive"], row["c3_minus_c4"], row["c3_recurrence"]), reverse=True)
    with (args.out / "contact_pair_recurrence.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    exclusive = [row for row in records if row["c3_exclusive"]]
    with (args.out / "c3_exclusive_pairs.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(exclusive)

    top = exclusive[:40]
    if top:
        labels = [f"{row['position_i']}–{row['position_j']}" for row in top][::-1]
        values = [row["c3_minus_c4"] for row in top][::-1]
        fig, ax = plt.subplots(figsize=(8.4, max(4.6, 0.22 * len(top))))
        ax.barh(labels, values, color="#b23628")
        ax.set_xlabel("C3 − C4 contact recurrence")
        ax.set_ylabel("OVA residue pair")
        ax.set_title("C3-exclusive contact pairs (C3≥0.70, C4≤0.20)")
        ax.grid(axis="x", alpha=0.2)
        fig.tight_layout()
        fig.savefig(args.out / "c3_exclusive_contact_pairs.png", dpi=220)
        plt.close(fig)
    print(f"C4 models={len(c4_rows)} C3 models={len(c3_rows)} exclusive pairs={len(exclusive)}")
    for row in exclusive[:30]:
        print(f"{row['position_i']}-{row['position_j']} C3={row['c3_recurrence']:.2f} C4={row['c4_recurrence']:.2f}")


if __name__ == "__main__":
    main()
