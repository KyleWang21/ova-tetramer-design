#!/usr/bin/env python3
"""Rank OVA surface patches for AI interface design.

The score deliberately combines independent evidence: solvent exposure on a
well-folded AF3 chain, AF3 per-residue confidence, MSA conservation, and the
frequency with which a residue participates in the otherwise low-confidence
P3 dimer predictions.  The dimer frequency is used as a surface-location hint,
not as evidence that the dimer itself is real.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
import re
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select, ShrakeRupley


AA = set("ACDEFGHIKLMNPQRSTVWY")
AROMATIC = set("FWY")
HYDROPHOBIC = set("FILMVWY")
MAX_ASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "Q": 225.0, "E": 223.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}


class OneChain(Select):
    def __init__(self, chain_id: str):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id


def clean_a3m_sequence(seq: str) -> str:
    return "".join(c for c in seq if not c.islower() and c != "\x00")


def read_a3m(path: pathlib.Path, length: int) -> list[str]:
    seqs: list[str] = []
    current: list[str] = []
    for raw in path.read_text(errors="replace").replace("\x00", "").splitlines():
        if raw.startswith(">"): 
            if current:
                seqs.append(clean_a3m_sequence("".join(current)))
                current = []
        else:
            current.append(raw.strip())
    if current:
        seqs.append(clean_a3m_sequence("".join(current)))
    return [s[:length] for s in seqs if len(s) >= length]


def msa_stats(seqs: list[str], length: int) -> tuple[np.ndarray, np.ndarray]:
    coverage = np.zeros(length, dtype=float)
    conservation = np.zeros(length, dtype=float)
    for i in range(length):
        column = [s[i] for s in seqs if i < len(s) and s[i] in AA]
        coverage[i] = len(column) / max(1, len(seqs))
        if not column:
            continue
        counts = np.array(list(Counter(column).values()), dtype=float)
        probs = counts / counts.sum()
        entropy = -float(np.sum(probs * np.log(probs)))
        conservation[i] = 1.0 - entropy / math.log(20.0)
    return coverage, np.clip(conservation, 0.0, 1.0)


def dimer_contact_frequency(path: pathlib.Path, system: str) -> dict[int, float]:
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    sets = [set(json.loads(r["interface_residues_json"])) for r in rows if r["system"] == system]
    counts: Counter[int] = Counter()
    for residues in sets:
        counts.update(int(x) for x in residues)
    return {r: n / len(sets) for r, n in counts.items()} if sets else {}


def normalise(values: np.ndarray) -> np.ndarray:
    lo, hi = np.nanpercentile(values, [5, 95])
    if hi <= lo:
        return np.zeros_like(values)
    return np.clip((values - lo) / (hi - lo), 0.0, 1.0)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", type=pathlib.Path, required=True)
    ap.add_argument("--a3m", type=pathlib.Path, required=True)
    ap.add_argument("--model-scores", type=pathlib.Path, required=True)
    ap.add_argument("--system", default="p3_dimer_msa")
    ap.add_argument("--chain", default="A")
    ap.add_argument("--dssp", type=pathlib.Path, help="deprecated; Shrake-Rupley SASA is used")
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--n-patches", type=int, default=3)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    structure = MMCIFParser(QUIET=True).get_structure("P3", str(args.cif))
    model = next(structure.get_models())
    chain = model[args.chain]
    residues = [r for r in chain if r.id[0] == " " and "CA" in r]
    if not residues:
        raise SystemExit(f"no protein residues found in chain {args.chain}")
    seq_length = max(int(r.id[1]) for r in residues)

    target_pdb = args.out / "p3_af3_monomer_chainA.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(target_pdb), OneChain(args.chain))
    # The vendored old mkdssp binary segfaults on current AF3 PDBs.  Biopython's
    # deterministic Shrake-Rupley implementation avoids an external dependency.
    sasa_model = next(PDBParser(QUIET=True).get_structure("P3_sasa", str(target_pdb)).get_models())
    ShrakeRupley(n_points=200).compute(sasa_model, level="R")
    sasa_by_residue = {
        int(r.id[1]): float(getattr(r, "sasa", 0.0))
        for r in sasa_model[args.chain]
        if r.id[0] == " "
    }

    seqs = read_a3m(args.a3m, seq_length)
    coverage, conservation = msa_stats(seqs, seq_length)
    dimer_freq = dimer_contact_frequency(args.model_scores, args.system)

    records: list[dict[str, object]] = []
    for residue in residues:
        rid = int(residue.id[1])
        aa = residue.resname
        one = {
            "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q",
            "GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K",
            "MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W",
            "TYR":"Y","VAL":"V",
        }.get(aa, "X")
        rsa = min(1.5, sasa_by_residue.get(rid, 0.0) / MAX_ASA.get(one, 200.0))
        plddt = float(residue["CA"].bfactor)
        coord = np.asarray(residue["CB"].coord if "CB" in residue else residue["CA"].coord)
        excluded_reason = ""
        # Full-sequence numbering.  Protect signal/N-terminal region, OVA257-264
        # SIINFEKL with a margin, the native N292 glycan neighbourhood, and the
        # flexible extreme C terminus.
        if rid <= 20:
            excluded_reason = "N_terminal_or_signal"
        elif 250 <= rid <= 271:
            excluded_reason = "SIINFEKL_margin"
        elif 286 <= rid <= 298:
            excluded_reason = "N292_glycan_margin"
        elif rid >= seq_length - 10:
            excluded_reason = "C_terminal"
        elif plddt < 75:
            excluded_reason = "low_plddt"
        elif rsa < 0.15:
            excluded_reason = "buried"
        records.append({
            "residue": rid,
            "aa": one,
            "rsa": rsa,
            "plddt": plddt,
            "msa_coverage": float(coverage[rid - 1]),
            "conservation": float(conservation[rid - 1]),
            "dimer_contact_frequency": float(dimer_freq.get(rid, 0.0)),
            "x": float(coord[0]), "y": float(coord[1]), "z": float(coord[2]),
            "excluded_reason": excluded_reason,
        })

    rsa_n = normalise(np.array([float(r["rsa"]) for r in records]))
    plddt_n = normalise(np.array([float(r["plddt"]) for r in records]))
    for i, r in enumerate(records):
        anchor = 1.0 if r["aa"] in AROMATIC else (0.6 if r["aa"] in HYDROPHOBIC else 0.2)
        r["residue_score"] = (
            0.45 * float(r["dimer_contact_frequency"])
            + 0.20 * rsa_n[i]
            + 0.15 * plddt_n[i]
            + 0.10 * float(r["conservation"])
            + 0.10 * anchor
        ) if not r["excluded_reason"] else 0.0

    coords = np.array([[r["x"], r["y"], r["z"]] for r in records], dtype=float)
    patch_candidates = []
    for i, center in enumerate(records):
        if center["excluded_reason"]:
            continue
        distances = np.linalg.norm(coords - coords[i], axis=1)
        members = [
            r for r, d in zip(records, distances)
            if d <= 11.0 and not r["excluded_reason"] and float(r["rsa"]) >= 0.20
        ]
        if len(members) < 5:
            continue
        # Use the strongest six-to-ten exposed residues as BindCraft hotspots.
        hotspots = sorted(members, key=lambda x: float(x["residue_score"]), reverse=True)[:10]
        score = float(np.mean([float(r["residue_score"]) for r in hotspots]))
        recurrence = float(np.mean([float(r["dimer_contact_frequency"]) for r in hotspots]))
        patch_candidates.append({
            "center": int(center["residue"]),
            "center_xyz": coords[i],
            "score": score,
            "recurrence": recurrence,
            "members": sorted(int(r["residue"]) for r in hotspots),
        })
    patch_candidates.sort(key=lambda x: (x["score"], x["recurrence"]), reverse=True)
    selected = []
    for patch in patch_candidates:
        if all(
            np.linalg.norm(patch["center_xyz"] - p["center_xyz"]) >= 22.0
            and set(patch["members"]).isdisjoint(p["members"])
            for p in selected
        ):
            selected.append(patch)
        if len(selected) == args.n_patches:
            break
    for rank, patch in enumerate(selected, 1):
        patch["rank"] = rank
        patch["name"] = f"patch{rank}_r{patch['center']}"
        patch["bindcraft_hotspots"] = ",".join(f"A{x}" for x in patch["members"])
        patch.pop("center_xyz")

    columns = list(records[0])
    with (args.out / "residue_hotspot_scores.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t")
        writer.writeheader(); writer.writerows(records)
    (args.out / "selected_patches.json").write_text(json.dumps({
        "source_cif": str(args.cif),
        "target_pdb": str(target_pdb),
        "msa_sequences": len(seqs),
        "score_definition": "0.45*dimer_contact_frequency + 0.20*RSA + 0.15*pLDDT + 0.10*conservation + 0.10*anchor_chemistry",
        "caveat": "dimer contact recurrence locates candidate surfaces; the dimer AF3 confidence itself failed",
        "patches": selected,
    }, indent=2) + "\n")

    xyz = coords
    score = np.array([float(r["residue_score"]) for r in records])
    fig = plt.figure(figsize=(11, 8))
    axes = [fig.add_subplot(121, projection="3d"), fig.add_subplot(122, projection="3d")]
    for ax, azim in zip(axes, (35, 215)):
        ax.scatter(xyz[:,0], xyz[:,1], xyz[:,2], c=score, cmap="viridis", s=15, alpha=.55,
                   vmin=0, vmax=max(0.01, float(score.max())))
        colors = ["#e41a1c", "#ff7f00", "#984ea3"]
        for patch, color in zip(selected, colors):
            inds = [i for i,r in enumerate(records) if int(r["residue"]) in patch["members"]]
            ax.scatter(xyz[inds,0], xyz[inds,1], xyz[inds,2], color=color, s=60,
                       edgecolor="black", linewidth=.4, label=patch["name"])
        ax.view_init(elev=18, azim=azim)
        ax.set_axis_off()
    axes[0].legend(loc="upper left", frameon=False)
    fig.suptitle("P3-13R candidate interface hotspots\nAF3 surface + MSA + repeated low-confidence dimer contacts")
    fig.tight_layout()
    fig.savefig(args.out / "hotspot_map.png", dpi=220, bbox_inches="tight")
    fig.savefig(args.out / "hotspot_map.svg", bbox_inches="tight")

    print(json.dumps(selected, indent=2))
    print(f"MSA sequences: {len(seqs)}")
    print(f"target PDB: {target_pdb}")


if __name__ == "__main__":
    main()
