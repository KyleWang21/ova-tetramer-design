#!/usr/bin/env python3
"""Compare AF3 P3-13R chain folds with experimental uncleaved ovalbumin."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from Bio import Align
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1
from Bio.SVDSuperimposer import SVDSuperimposer


AA_OVERRIDES = {"SEP": "S"}  # 1OVA contains two phosphoserines.


def chain_records(cif: Path, chain_id: str):
    structure = MMCIFParser(QUIET=True).get_structure(cif.stem, str(cif))
    chain = structure[0][chain_id]
    records = []
    for residue in chain:
        if not is_aa(residue, standard=False) or "CA" not in residue:
            continue
        aa = AA_OVERRIDES.get(residue.resname, protein_letters_3to1.get(residue.resname, "X"))
        records.append((aa, residue["CA"].coord.astype(float), float(residue["CA"].bfactor)))
    return records


def paired_indices(ref_seq: str, query_seq: str):
    aligner = Align.PairwiseAligner(mode="global")
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(ref_seq, query_seq)[0]
    pairs = []
    for (r0, r1), (q0, q1) in zip(*alignment.aligned):
        assert r1 - r0 == q1 - q0
        pairs.extend(zip(range(r0, r1), range(q0, q1)))
    return pairs


def superpose(ref_records, query_records):
    ref_seq = "".join(row[0] for row in ref_records)
    query_seq = "".join(row[0] for row in query_records)
    pairs = paired_indices(ref_seq, query_seq)
    ref = np.asarray([ref_records[i][1] for i, _ in pairs])
    query = np.asarray([query_records[j][1] for _, j in pairs])
    plddt = np.asarray([query_records[j][2] for _, j in pairs])
    positions = np.asarray([j + 1 for _, j in pairs])
    sup = SVDSuperimposer()
    sup.set(ref, query)
    sup.run()
    transformed = sup.get_transformed()
    distances = np.linalg.norm(transformed - ref, axis=1)
    identity = sum(ref_records[i][0] == query_records[j][0] for i, j in pairs) / len(pairs)
    high = plddt >= 70.0
    return {
        "n": len(pairs),
        "identity": identity,
        "rmsd": float(np.sqrt(np.mean(distances ** 2))),
        "rmsd_plddt70": float(np.sqrt(np.mean(distances[high] ** 2))),
        "within_2a": float(np.mean(distances <= 2.0)),
        "within_3a": float(np.mean(distances <= 3.0)),
        "within_5a": float(np.mean(distances <= 5.0)),
        "positions": positions,
        "distances": distances,
        "ref_coords": ref,
        "query_coords": transformed,
        "plddt": plddt,
    }


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = ["model", "chain", "matched_ca", "sequence_identity", "ca_rmsd_a",
              "ca_rmsd_plddt70_a", "within_2a", "within_3a", "within_5a", "cif_path"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--reference-chain", default="A")
    ap.add_argument("--predictions", type=Path, required=True)
    ap.add_argument("--chains", default="A,B", help="comma-separated AF3 chain IDs")
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--label", default="P3-13R")
    ap.add_argument("--mutation-positions", default="14,69,77,79,191,192,193")
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    ref_records = chain_records(args.reference, args.reference_chain)
    results = []
    details = []
    for cif in sorted(args.predictions.glob("seed-*/*_model.cif")):
        for chain_id in args.chains.split(","):
            comparison = superpose(ref_records, chain_records(cif, chain_id))
            details.append((cif, chain_id, comparison))
            results.append({
                "model": cif.parent.name,
                "chain": chain_id,
                "matched_ca": comparison["n"],
                "sequence_identity": round(comparison["identity"], 5),
                "ca_rmsd_a": round(comparison["rmsd"], 4),
                "ca_rmsd_plddt70_a": round(comparison["rmsd_plddt70"], 4),
                "within_2a": round(comparison["within_2a"], 5),
                "within_3a": round(comparison["within_3a"], 5),
                "within_5a": round(comparison["within_5a"], 5),
                "cif_path": str(cif),
            })
    write_tsv(args.output_dir / "model_comparison.tsv", results)

    rmsd = [float(row["ca_rmsd_a"]) for row in results]
    rmsd_hc = [float(row["ca_rmsd_plddt70_a"]) for row in results]
    summary = {
        "reference": str(args.reference),
        "reference_chain": args.reference_chain,
        "reference_matched_residues": details[0][2]["n"],
        "n_af3_chains": len(results),
        "sequence_identity": details[0][2]["identity"],
        "median_ca_rmsd_a": statistics.median(rmsd),
        "range_ca_rmsd_a": [min(rmsd), max(rmsd)],
        "median_ca_rmsd_plddt70_a": statistics.median(rmsd_hc),
        "median_within_2a": statistics.median(float(row["within_2a"]) for row in results),
        "median_within_3a": statistics.median(float(row["within_3a"]) for row in results),
        "median_within_5a": statistics.median(float(row["within_5a"]) for row in results),
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    best_cif, best_chain, best = min(details, key=lambda row: row[2]["rmsd"])
    fig, axes = plt.subplots(2, 1, figsize=(9.2, 6.5), gridspec_kw={"height_ratios": [1, 2]})
    axes[0].hist(rmsd, bins=12, color="#356ec9", alpha=0.85)
    axes[0].axvline(summary["median_ca_rmsd_a"], color="black", ls="--", lw=1.5,
                    label=f"median {summary['median_ca_rmsd_a']:.2f} Å")
    axes[0].set_xlabel("Cα RMSD to 1OVA chain A (Å)")
    axes[0].set_ylabel("AF3 chains")
    axes[0].legend(frameon=False)
    axes[1].plot(best["positions"], best["distances"], lw=1.2, color="#356ec9")
    mutation_positions = [int(value) for value in args.mutation_positions.split(",") if value]
    for pos in mutation_positions:
        axes[1].axvline(pos, color="#d13c3c", lw=0.8, alpha=0.5)
    axes[1].axhline(2.0, color="#777777", ls=":", lw=1)
    axes[1].set_xlabel(f"{args.label} residue position")
    axes[1].set_ylabel("Best-model Cα deviation (Å)")
    axes[1].set_title(
        f"Best: {best_cif.parent.name} chain {best_chain}; red = {len(mutation_positions)} mutation sites"
    )
    axes[1].set_xlim(1, 386)
    fig.suptitle(f"{args.label} AF3 chain fold vs 1.95 Å uncleaved ovalbumin crystal structure (1OVA)")
    fig.tight_layout()
    fig.savefig(args.output_dir / "experimental_comparison.png", dpi=220)

    fig = plt.figure(figsize=(10, 4.8))
    for idx, azim in enumerate((30, 120), 1):
        ax = fig.add_subplot(1, 2, idx, projection="3d")
        ax.plot(*best["ref_coords"].T, color="#555555", lw=1.2, label="1OVA crystal")
        ax.plot(*best["query_coords"].T, color="#2b70d6", lw=1.0, alpha=0.9, label=f"AF3 {args.label}")
        ax.view_init(elev=20, azim=azim)
        ax.set_axis_off()
        if idx == 1:
            ax.legend(frameon=False, loc="upper left")
    fig.suptitle(f"Experimental overlay: {best['n']} Cα; RMSD {best['rmsd']:.2f} Å")
    fig.tight_layout()
    fig.savefig(args.output_dir / "experimental_overlay.png", dpi=220)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
