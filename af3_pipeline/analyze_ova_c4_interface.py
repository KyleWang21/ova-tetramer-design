#!/usr/bin/env python3
"""Extract the recurrent OVA--OVA interface from high-confidence C4 models.

The input C4 models may contain a C-terminal oligomerization module.  Only OVA
residues 1..386 are considered for contacts, exposure, the cropped backbone,
and the proposed ProteinMPNN design positions.
"""

from __future__ import annotations

import argparse
import copy
import csv
import math
import pathlib
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import MMCIFParser, Model, PDBIO, Select, Structure
from Bio.PDB.SASA import ShrakeRupley


class OvaSelect(Select):
    def __init__(self, end: int) -> None:
        self.end = end

    def accept_residue(self, residue) -> bool:  # noqa: ANN001
        return int(residue.id[1]) <= self.end and residue.id[0] == " "


def representative_atoms(chain, end: int) -> tuple[list[int], np.ndarray, dict[int, float]]:  # noqa: ANN001
    positions: list[int] = []
    coords: list[np.ndarray] = []
    plddt: dict[int, float] = {}
    for residue in chain:
        pos = int(residue.id[1])
        if residue.id[0] != " " or pos > end or "CA" not in residue:
            continue
        atom = residue["CB"] if "CB" in residue else residue["CA"]
        positions.append(pos)
        coords.append(np.asarray(atom.coord, dtype=np.float32))
        plddt[pos] = float(residue["CA"].bfactor)
    return positions, np.stack(coords), plddt


def isolated_chain_sasa(chain, end: int) -> dict[int, float]:  # noqa: ANN001
    isolated = Structure.Structure("ova")
    model = Model.Model(0)
    copied = copy.deepcopy(chain)
    for residue in list(copied):
        if residue.id[0] != " " or int(residue.id[1]) > end:
            copied.detach_child(residue.id)
    model.add(copied)
    isolated.add(model)
    ShrakeRupley(probe_radius=1.4, n_points=100).compute(isolated, level="R")
    return {int(res.id[1]): float(res.sasa) for res in copied if hasattr(res, "sasa")}


def read_rows(
    paths: list[pathlib.Path], iptm_min: float, weak_min: float, systems: set[str]
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in paths:
        for row in csv.DictReader(path.open(), delimiter="\t"):
            if (
                (not systems or row["system"] in systems)
                and
                float(row["iptm"]) >= iptm_min
                and float(row["min_incident_iptm"]) >= weak_min
                and int(float(row["has_clash"])) == 0
                and int(row["n_chains"]) == 4
            ):
                rows.append(row)
    if not rows:
        raise SystemExit("no models pass the requested C4 confidence gates")
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", type=pathlib.Path, action="append", required=True)
    ap.add_argument("--system", action="append", default=[], help="Exact system name to retain")
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--ova-end", type=int, default=386)
    ap.add_argument("--contact-cutoff", type=float, default=8.0)
    ap.add_argument("--iptm-min", type=float, default=0.50)
    ap.add_argument("--weak-min", type=float, default=0.45)
    ap.add_argument("--recurrence-min", type=float, default=0.60)
    ap.add_argument("--sasa-min", type=float, default=25.0)
    ap.add_argument("--plddt-min", type=float, default=70.0)
    ap.add_argument("--top", type=int, default=24)
    ap.add_argument("--exclude", default="257-264,289-295")
    ap.add_argument("--label", default="candidate")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    excluded: set[int] = set()
    for part in args.exclude.split(","):
        lo, hi = (int(x) for x in part.split("-", 1))
        excluded.update(range(lo, hi + 1))

    rows = read_rows(args.scores, args.iptm_min, args.weak_min, set(args.system))
    parser = MMCIFParser(QUIET=True)
    contact_chain_hits: defaultdict[int, int] = defaultdict(int)
    contact_counts: defaultdict[int, int] = defaultdict(int)
    plddt_values: defaultdict[int, list[float]] = defaultdict(list)
    aa: dict[int, str] = {}
    three_to_one = {
        "ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D", "CYS":"C",
        "GLN":"Q", "GLU":"E", "GLY":"G", "HIS":"H", "ILE":"I",
        "LEU":"L", "LYS":"K", "MET":"M", "PHE":"F", "PRO":"P",
        "SER":"S", "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V",
    }
    best = max(rows, key=lambda row: (float(row["iptm"]), float(row["min_incident_iptm"])))
    first_chain = None
    per_model_interfaces: list[set[int]] = []
    for row in rows:
        cif = pathlib.Path(row["cif_path"])
        structure = parser.get_structure(cif.stem, str(cif))
        chains = list(next(structure.get_models()).get_chains())
        if len(chains) != 4:
            raise ValueError(f"{cif}: expected four chains, found {len(chains)}")
        by_chain = {}
        chain_hits = {str(chain.id): set() for chain in chains}
        model_interface: set[int] = set()
        for chain in chains:
            positions, coords, plddt = representative_atoms(chain, args.ova_end)
            by_chain[str(chain.id)] = (positions, coords)
            for pos, value in plddt.items():
                plddt_values[pos].append(value)
            if first_chain is None:
                first_chain = chain
                for residue in chain:
                    if residue.id[0] == " " and int(residue.id[1]) <= args.ova_end:
                        aa[int(residue.id[1])] = three_to_one.get(residue.resname, "X")
        for i, chain_a in enumerate(chains):
            posa, xyza = by_chain[str(chain_a.id)]
            for chain_b in chains[i + 1:]:
                posb, xyzb = by_chain[str(chain_b.id)]
                delta = xyza[:, None, :] - xyzb[None, :, :]
                mask = np.einsum("ijk,ijk->ij", delta, delta) < args.contact_cutoff ** 2
                ai, bj = np.nonzero(mask)
                for index in ai:
                    pos = posa[int(index)]
                    chain_hits[str(chain_a.id)].add(pos)
                    model_interface.add(pos)
                    contact_counts[pos] += 1
                for index in bj:
                    pos = posb[int(index)]
                    chain_hits[str(chain_b.id)].add(pos)
                    model_interface.add(pos)
                    contact_counts[pos] += 1
        for hits in chain_hits.values():
            for pos in hits:
                contact_chain_hits[pos] += 1
        per_model_interfaces.append(model_interface)

    assert first_chain is not None
    sasa = isolated_chain_sasa(first_chain, args.ova_end)
    denom = len(rows) * 4
    output: list[dict[str, object]] = []
    for pos in range(1, args.ova_end + 1):
        recurrence = contact_chain_hits[pos] / denom
        mean_contacts = contact_counts[pos] / denom
        mean_plddt = float(np.mean(plddt_values[pos])) if plddt_values[pos] else 0.0
        surface = sasa.get(pos, 0.0)
        score = recurrence * math.log1p(mean_contacts) * min(surface / 50.0, 2.0)
        eligible = (
            recurrence >= args.recurrence_min
            and surface >= args.sasa_min
            and mean_plddt >= args.plddt_min
            and pos not in excluded
            and aa.get(pos) not in {"G", "P", "X"}
        )
        output.append({
            "position": pos,
            "wt_aa": aa.get(pos, "X"),
            "contact_chain_recurrence": round(recurrence, 6),
            "mean_cross_chain_cb_contacts": round(mean_contacts, 4),
            "isolated_chain_sasa_a2": round(surface, 3),
            "mean_plddt": round(mean_plddt, 3),
            "design_score": round(score, 6),
            "eligible": int(eligible),
            "excluded_functional_region": int(pos in excluded),
        })
    ranked = sorted((row for row in output if row["eligible"]), key=lambda row: float(row["design_score"]), reverse=True)
    selected = sorted(int(row["position"]) for row in ranked[:args.top])
    for row in output:
        row["selected"] = int(int(row["position"]) in selected)
    fields = list(output[0])
    with (args.out / "interface_positions.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader(); writer.writerows(output)
    (args.out / "design_positions.txt").write_text(",".join(map(str, selected)) + "\n")

    best_cif = pathlib.Path(best["cif_path"])
    best_structure = parser.get_structure("best_c4", str(best_cif))
    io = PDBIO(); io.set_structure(best_structure)
    io.save(str(args.out / "c4_ova_only_backbone.pdb"), OvaSelect(args.ova_end))

    x = np.arange(1, args.ova_end + 1)
    recurrence = np.array([float(row["contact_chain_recurrence"]) for row in output])
    surface = np.array([float(row["isolated_chain_sasa_a2"]) for row in output])
    fig, axes = plt.subplots(2, 1, figsize=(12, 5.8), sharex=True, constrained_layout=True)
    axes[0].plot(x, recurrence, color="#1769aa", lw=1.2)
    axes[0].axhline(args.recurrence_min, color="#777777", ls="--", lw=1)
    axes[0].scatter(selected, recurrence[np.array(selected) - 1], color="#d62728", s=25, zorder=3)
    axes[0].set_ylabel("C4 contact recurrence")
    axes[0].set_ylim(0, 1.05)
    axes[1].plot(x, surface, color="#4b8f29", lw=1.1)
    axes[1].axhline(args.sasa_min, color="#777777", ls="--", lw=1)
    axes[1].scatter(selected, surface[np.array(selected) - 1], color="#d62728", s=25, zorder=3)
    axes[1].set_ylabel("Monomer SASA (Å²)")
    axes[1].set_xlabel("OVA residue")
    fig.suptitle(f"Pure-OVA interface prior from {len(rows)} passing {args.label} models")
    fig.savefig(args.out / "ova_c4_interface_recurrence.png", dpi=180)
    plt.close(fig)

    (args.out / "analysis_summary.tsv").write_text(
        "passing_models\tbest_iptm\tbest_weakest_chain\tbest_cif\tselected_positions\n"
        f"{len(rows)}\t{best['iptm']}\t{best['min_incident_iptm']}\t{best_cif}\t"
        f"{','.join(map(str, selected))}\n"
    )
    print(f"passing models: {len(rows)}")
    print(f"best backbone: ipTM={best['iptm']}, weak={best['min_incident_iptm']}, {best_cif}")
    print("selected OVA positions:", ",".join(map(str, selected)))


if __name__ == "__main__":
    main()
