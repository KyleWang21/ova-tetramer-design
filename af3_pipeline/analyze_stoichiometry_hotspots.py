#!/usr/bin/env python3
"""Find OVA surface positions that distinguish AF3 C4 and C3 ensembles.

Residue contact recurrence is measured per chain, not per structure, so the
three- and four-chain ensembles are directly comparable.  The resulting
positive and negative hotspots are intended as ProteinMPNN design positions;
AF3 remains the independent acceptance test.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import MMCIFParser, Model, PDBIO, Select, ShrakeRupley, Structure


AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


class OvaSelect(Select):
    def __init__(self, end: int) -> None:
        self.end = end

    def accept_residue(self, residue) -> bool:  # noqa: ANN001
        return residue.id[0] == " " and int(residue.id[1]) <= self.end


def parse_ranges(spec: str) -> set[int]:
    result: set[int] = set()
    for part in filter(None, spec.split(",")):
        if "-" in part:
            lo, hi = (int(value) for value in part.split("-", 1))
            result.update(range(lo, hi + 1))
        else:
            result.add(int(part))
    return result


def read_rows(path: Path, system: str | None, iptm_min: float, weak_min: float) -> list[dict[str, str]]:
    rows = []
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if system and row.get("system") != system:
                continue
            iptm = float(row["iptm"])
            weak = float(row.get("min_incident_iptm", row.get("weakest_chain", 0)))
            clash = int(float(row.get("has_clash", 0)))
            if iptm >= iptm_min and weak >= weak_min and clash == 0:
                rows.append(row)
    if not rows:
        raise SystemExit(f"no passing models in {path}")
    return rows


def representative_atoms(chain, end: int) -> tuple[list[int], np.ndarray, dict[int, float]]:  # noqa: ANN001
    positions, coordinates, confidence = [], [], {}
    for residue in chain:
        position = int(residue.id[1])
        if residue.id[0] != " " or position > end or "CA" not in residue:
            continue
        atom = residue["CB"] if "CB" in residue else residue["CA"]
        positions.append(position)
        coordinates.append(np.asarray(atom.coord, dtype=np.float32))
        confidence[position] = float(residue["CA"].bfactor)
    return positions, np.stack(coordinates), confidence


def ensemble_contacts(
    rows: list[dict[str, str]], end: int, cutoff: float
) -> tuple[dict[int, float], dict[int, float], dict[int, float], object]:
    parser = MMCIFParser(QUIET=True)
    chain_hits: defaultdict[int, int] = defaultdict(int)
    contact_counts: defaultdict[int, int] = defaultdict(int)
    plddt: defaultdict[int, list[float]] = defaultdict(list)
    denominator = 0
    first_chain = None
    for row in rows:
        cif = Path(row["cif_path"])
        structure = parser.get_structure(cif.stem, str(cif))
        chains = list(next(structure.get_models()).get_chains())
        denominator += len(chains)
        per_chain_hits = {str(chain.id): set() for chain in chains}
        coordinates = {}
        for chain in chains:
            if first_chain is None:
                first_chain = chain
            positions, xyz, confidence = representative_atoms(chain, end)
            coordinates[str(chain.id)] = (positions, xyz)
            for position, value in confidence.items():
                plddt[position].append(value)
        for index, chain_a in enumerate(chains):
            posa, xyza = coordinates[str(chain_a.id)]
            for chain_b in chains[index + 1:]:
                posb, xyzb = coordinates[str(chain_b.id)]
                delta = xyza[:, None, :] - xyzb[None, :, :]
                ai, bj = np.nonzero(np.einsum("ijk,ijk->ij", delta, delta) < cutoff ** 2)
                for atom_index in ai:
                    position = posa[int(atom_index)]
                    per_chain_hits[str(chain_a.id)].add(position)
                    contact_counts[position] += 1
                for atom_index in bj:
                    position = posb[int(atom_index)]
                    per_chain_hits[str(chain_b.id)].add(position)
                    contact_counts[position] += 1
        for hits in per_chain_hits.values():
            for position in hits:
                chain_hits[position] += 1
    assert first_chain is not None
    recurrence = {position: chain_hits[position] / denominator for position in range(1, end + 1)}
    mean_contacts = {position: contact_counts[position] / denominator for position in range(1, end + 1)}
    confidence = {
        position: float(np.mean(plddt[position])) if plddt[position] else 0.0
        for position in range(1, end + 1)
    }
    return recurrence, mean_contacts, confidence, first_chain


def isolated_sasa(chain, end: int) -> dict[int, float]:  # noqa: ANN001
    structure = Structure.Structure("isolated")
    model = Model.Model(0)
    copied = copy.deepcopy(chain)
    for residue in list(copied):
        if residue.id[0] != " " or int(residue.id[1]) > end:
            copied.detach_child(residue.id)
    model.add(copied)
    structure.add(model)
    ShrakeRupley(probe_radius=1.4, n_points=100).compute(structure, level="R")
    return {int(residue.id[1]): float(residue.sasa) for residue in copied}


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
    parser.add_argument("--recurrence-min", type=float, default=0.50)
    parser.add_argument("--delta-min", type=float, default=0.20)
    parser.add_argument("--sasa-min", type=float, default=20.0)
    parser.add_argument("--plddt-min", type=float, default=70.0)
    parser.add_argument("--top-positive", type=int, default=10)
    parser.add_argument("--top-negative", type=int, default=10)
    parser.add_argument("--exclude", default="14,69,77,79,94,96,155,190-193,196,213,257-264,289-295")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    c4_rows = read_rows(args.c4_scores, args.c4_system, args.iptm_min, args.weak_min)
    c3_rows = read_rows(args.c3_scores, args.c3_system, args.iptm_min, args.weak_min)
    c4_rec, c4_contacts, c4_plddt, first_chain = ensemble_contacts(
        c4_rows, args.ova_end, args.contact_cutoff
    )
    c3_rec, c3_contacts, c3_plddt, _ = ensemble_contacts(
        c3_rows, args.ova_end, args.contact_cutoff
    )
    sasa = isolated_sasa(first_chain, args.ova_end)
    excluded = parse_ranges(args.exclude)
    sequence = {
        int(residue.id[1]): AA3.get(residue.resname, "X")
        for residue in first_chain
        if residue.id[0] == " " and int(residue.id[1]) <= args.ova_end
    }

    records = []
    for position in range(1, args.ova_end + 1):
        eligible = (
            position not in excluded
            and sequence.get(position) not in {"C", "G", "P", "X"}
            and sasa.get(position, 0.0) >= args.sasa_min
            and min(c4_plddt[position], c3_plddt[position]) >= args.plddt_min
        )
        records.append({
            "position": position,
            "aa": sequence.get(position, "X"),
            "c4_contact_recurrence": c4_rec[position],
            "c3_contact_recurrence": c3_rec[position],
            "c4_minus_c3": c4_rec[position] - c3_rec[position],
            "c3_minus_c4": c3_rec[position] - c4_rec[position],
            "c4_mean_contacts": c4_contacts[position],
            "c3_mean_contacts": c3_contacts[position],
            "isolated_sasa_a2": sasa.get(position, 0.0),
            "c4_mean_plddt": c4_plddt[position],
            "c3_mean_plddt": c3_plddt[position],
            "eligible": int(eligible),
        })
    positive = sorted(
        (row for row in records if row["eligible"]
         and row["c4_contact_recurrence"] >= args.recurrence_min
         and row["c4_minus_c3"] >= args.delta_min),
        key=lambda row: (row["c4_minus_c3"], row["c4_contact_recurrence"], row["isolated_sasa_a2"]),
        reverse=True,
    )[:args.top_positive]
    negative = sorted(
        (row for row in records if row["eligible"]
         and row["c3_contact_recurrence"] >= args.recurrence_min
         and row["c3_minus_c4"] >= args.delta_min),
        key=lambda row: (row["c3_minus_c4"], row["c3_contact_recurrence"], row["isolated_sasa_a2"]),
        reverse=True,
    )[:args.top_negative]
    positive_positions = {int(row["position"]) for row in positive}
    negative_positions = {int(row["position"]) for row in negative}
    selected = sorted(positive_positions | negative_positions)
    for row in records:
        position = int(row["position"])
        row["positive_hotspot"] = int(position in positive_positions)
        row["negative_hotspot"] = int(position in negative_positions)
        row["selected"] = int(position in selected)

    fields = list(records[0])
    with (args.out / "stoichiometry_hotspots.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    (args.out / "design_positions.txt").write_text(",".join(map(str, selected)) + "\n")
    (args.out / "summary.json").write_text(json.dumps({
        "c4_models": len(c4_rows),
        "c3_models": len(c3_rows),
        "positive_positions": sorted(positive_positions),
        "negative_positions": sorted(negative_positions),
        "all_design_positions": selected,
    }, indent=2) + "\n")
    for label, rows in (("c4", c4_rows), ("c3", c3_rows)):
        best = max(rows, key=lambda row: (
            float(row["iptm"]),
            float(row.get("min_incident_iptm", row.get("weakest_chain", 0))),
        ))
        structure = MMCIFParser(QUIET=True).get_structure(label, best["cif_path"])
        writer = PDBIO()
        writer.set_structure(structure)
        writer.save(str(args.out / f"representative_{label}_ova_only.pdb"), OvaSelect(args.ova_end))

    x = np.arange(1, args.ova_end + 1)
    c4 = np.array([c4_rec[position] for position in x])
    c3 = np.array([c3_rec[position] for position in x])
    delta = c4 - c3
    fig, axes = plt.subplots(2, 1, figsize=(12, 6.2), sharex=True, constrained_layout=True)
    axes[0].plot(x, c4, color="#2369bd", lw=1.25, label="C4 (40 AF3 models)")
    axes[0].plot(x, c3, color="#d54a3a", lw=1.25, label="C3 (15 AF3 models)")
    axes[0].set_ylabel("inter-chain contact recurrence")
    axes[0].set_ylim(0, 1.05)
    axes[0].legend(frameon=False)
    axes[1].axhline(0, color="#777777", lw=0.8)
    axes[1].bar(x, delta, width=1.0, color=np.where(delta >= 0, "#2369bd", "#d54a3a"))
    axes[1].scatter(sorted(positive_positions), delta[np.array(sorted(positive_positions)) - 1],
                    s=28, color="#082f6b", label="C4-positive design")
    axes[1].scatter(sorted(negative_positions), delta[np.array(sorted(negative_positions)) - 1],
                    s=28, color="#7f160c", marker="s", label="C3-negative design")
    axes[1].set_ylabel("C4 − C3 recurrence")
    axes[1].set_xlabel("OVA residue")
    axes[1].legend(frameon=False, ncol=2)
    fig.suptitle("OVA-body stoichiometry hotspots from AF3 ensembles")
    fig.savefig(args.out / "c4_vs_c3_hotspots.png", dpi=220)
    plt.close(fig)
    print(f"C4 models={len(c4_rows)} C3 models={len(c3_rows)}")
    print("positive positions:", ",".join(map(str, sorted(positive_positions))))
    print("negative positions:", ",".join(map(str, sorted(negative_positions))))
    print("all design positions:", ",".join(map(str, selected)))


if __name__ == "__main__":
    main()
