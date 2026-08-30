#!/usr/bin/env python3
"""Rank ProteinMPNN sequences for one target PDB state against decoy PDB states."""

from __future__ import annotations

import argparse
import csv
import pathlib
import statistics

from Bio.PDB import PDBIO, PDBParser, Select
from colabdesign.mpnn import mk_mpnn_model


AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


class ChainSelect(Select):
    def __init__(self, chain_id: str) -> None:
        self.chain_id = chain_id

    def accept_chain(self, chain) -> bool:  # noqa: ANN001
        return str(chain.id) == self.chain_id

    def accept_residue(self, residue) -> bool:  # noqa: ANN001
        return residue.id[0] == " " and int(residue.id[1]) <= 386


def chain_ids(pdb: pathlib.Path) -> list[str]:
    structure = PDBParser(QUIET=True).get_structure(pdb.stem, str(pdb))
    return [
        str(chain.id) for chain in next(structure.get_models()).get_chains()
        if sum(1 for residue in chain if residue.id[0] == " ") >= 386
    ]


def native_sequence(pdb: pathlib.Path) -> str:
    structure = PDBParser(QUIET=True).get_structure(pdb.stem, str(pdb))
    chain = next(chain for chain in structure.get_chains()
                 if sum(1 for residue in chain if residue.id[0] == " ") >= 386)
    return "".join(AA3[res.resname] for res in chain if res.id[0] == " " and int(res.id[1]) <= 386)


def score_state(model, pdb: pathlib.Path, sequences: list[str], positions: list[int], rm_aa: str) -> list[float]:  # noqa: ANN001
    chains = chain_ids(pdb)
    if not chains:
        raise ValueError(f"no 386-residue protein chain in {pdb}")
    spec = ",".join(f"{chain}{position}" for chain in chains for position in positions)
    model.prep_inputs(
        pdb_filename=str(pdb), chain=",".join(chains), homooligomer=len(chains) > 1,
        fix_pos=spec, inverse=True, rm_aa=rm_aa, verbose=False,
    )
    return [float(model.score(seq=sequence, key=model.key())["score"]) for sequence in sequences]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--positive-pdb", action="append", type=pathlib.Path, required=True)
    ap.add_argument("--negative-pdb", action="append", type=pathlib.Path, required=True)
    ap.add_argument("--candidates", type=pathlib.Path, required=True)
    ap.add_argument("--positions", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20260830)
    ap.add_argument("--rm-aa", default="C,P")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    positions = [int(value) for value in args.positions.read_text().strip().split(",")]
    rows = list(csv.DictReader(args.candidates.open(), delimiter="\t"))
    sequences = [row["full_sequence"] for row in rows]
    if not rows or any(len(sequence) != 386 for sequence in sequences):
        raise SystemExit("candidate table must contain 386-aa full_sequence records")

    native = native_sequence(args.positive_pdb[0])
    if len(native) != 386:
        raise SystemExit(f"target backbone chain has {len(native)} residues")
    monomer = args.out / "target_chain_a_monomer.pdb"
    structure = PDBParser(QUIET=True).get_structure("target", str(args.positive_pdb[0]))
    first_chain = next(chain for chain in structure.get_chains()
                       if sum(1 for residue in chain if residue.id[0] == " ") >= 386)
    io = PDBIO(); io.set_structure(structure); io.save(str(monomer), ChainSelect(str(first_chain.id)))

    model = mk_mpnn_model(
        model_name="v_48_020", backbone_noise=0.0, seed=args.seed,
        verbose=False, weights="original",
    )
    all_sequences = [native, *sequences]
    scores: dict[str, list[float]] = {}
    for index, pdb in enumerate(args.positive_pdb, 1):
        scores[f"positive_{index}"] = score_state(model, pdb, all_sequences, positions, args.rm_aa)
    for index, pdb in enumerate(args.negative_pdb, 1):
        scores[f"negative_{index}"] = score_state(model, pdb, all_sequences, positions, args.rm_aa)
    scores["monomer"] = score_state(model, monomer, all_sequences, positions, args.rm_aa)

    native_positive = statistics.mean(values[0] for key, values in scores.items() if key.startswith("positive_"))
    native_monomer = scores["monomer"][0]
    for row_index, row in enumerate(rows, 1):
        positive = [values[row_index] for key, values in scores.items() if key.startswith("positive_")]
        negative = [values[row_index] for key, values in scores.items() if key.startswith("negative_")]
        monomer_penalty = scores["monomer"][row_index] - native_monomer
        specificity = min(negative) - max(positive)
        row["positive_mean_score"] = statistics.mean(positive)
        row["positive_worst_score"] = max(positive)
        row["negative_best_score"] = min(negative)
        row["decoy_specificity"] = specificity
        row["positive_gain_vs_native"] = native_positive - statistics.mean(positive)
        row["monomer_penalty_vs_native"] = monomer_penalty
        row["multistate_objective"] = statistics.mean(positive) - 0.75 * specificity + 0.50 * max(0.0, monomer_penalty)
        for key, values in scores.items():
            row[f"{key}_score"] = values[row_index]
    rows.sort(key=lambda row: (float(row["multistate_objective"]), -float(row["decoy_specificity"])))
    for rank, row in enumerate(rows, 1):
        row["multistate_rank"] = rank
    fields = ["multistate_rank", *[field for field in rows[0] if field != "multistate_rank"]]
    with (args.out / "multistate_all.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t"); writer.writeheader(); writer.writerows(rows)

    selected: list[dict[str, str]] = []
    for row in rows:
        if float(row["positive_gain_vs_native"]) < 0.25 or float(row["monomer_penalty_vs_native"]) > 0.30:
            continue
        interface = row["interface_sequence"]
        if all(sum(a != b for a, b in zip(interface, old["interface_sequence"])) >= 2 for old in selected):
            selected.append(row)
        if len(selected) == args.top:
            break
    if len(selected) < args.top:
        for row in rows:
            if row not in selected:
                selected.append(row)
            if len(selected) == args.top:
                break
    with (args.out / "multistate_shortlist.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t"); writer.writeheader(); writer.writerows(selected)
    with (args.out / "multistate_shortlist.fasta").open("w") as handle:
        for index, row in enumerate(selected, 1):
            handle.write(
                f">OVA-NC-C4-S32-AI{index:02d}|rank={row['multistate_rank']}|mut={row['mutations']}|"
                f"spec={float(row['decoy_specificity']):.4f}\n{row['full_sequence']}\n"
            )
    with (args.out / "native_state_scores.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t"); writer.writerow(["state", "native_score"])
        for key, values in scores.items(): writer.writerow([key, values[0]])
    for index, row in enumerate(selected, 1):
        print(index, row["mutations"], f"target={float(row['positive_mean_score']):.3f}",
              f"spec={float(row['decoy_specificity']):.3f}",
              f"mono={float(row['monomer_penalty_vs_native']):.3f}")


if __name__ == "__main__":
    main()
