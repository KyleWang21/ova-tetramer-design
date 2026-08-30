#!/usr/bin/env python3
"""Rank OVA double-disulfide designs on successful vs failed AF3 C4 states."""

from __future__ import annotations

import argparse
import csv
import pathlib
import statistics

from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select
from colabdesign.mpnn import mk_mpnn_model


AA3 = {
    "ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D", "CYS":"C",
    "GLN":"Q", "GLU":"E", "GLY":"G", "HIS":"H", "ILE":"I",
    "LEU":"L", "LYS":"K", "MET":"M", "PHE":"F", "PRO":"P",
    "SER":"S", "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V",
}


class ProteinSelect(Select):
    def __init__(self, chains: set[str] | None = None) -> None:
        self.chains = chains

    def accept_chain(self, chain) -> bool:  # noqa: ANN001
        return self.chains is None or str(chain.id) in self.chains

    def accept_residue(self, residue) -> bool:  # noqa: ANN001
        return residue.id[0] == " " and int(residue.id[1]) <= 386


def cif_to_pdb(cif: pathlib.Path, pdb: pathlib.Path, chains: set[str] | None = None) -> None:
    structure = MMCIFParser(QUIET=True).get_structure(pdb.stem, str(cif))
    io = PDBIO(); io.set_structure(structure); io.save(str(pdb), ProteinSelect(chains))


def chain_ids(pdb: pathlib.Path) -> list[str]:
    structure = PDBParser(QUIET=True).get_structure(pdb.stem, str(pdb))
    return [str(chain.id) for chain in next(structure.get_models()).get_chains()]


def native_sequence(pdb: pathlib.Path) -> str:
    structure = PDBParser(QUIET=True).get_structure(pdb.stem, str(pdb))
    chain = next(structure.get_chains())
    return "".join(AA3[res.resname] for res in chain if res.id[0] == " ")


def score_state(model, pdb: pathlib.Path, sequences: list[str], positions: list[int], rm_aa: str) -> list[float]:  # noqa: ANN001
    chains = chain_ids(pdb)
    spec = ",".join(f"{chain}{pos}" for chain in chains for pos in positions)
    model.prep_inputs(
        pdb_filename=str(pdb), chain=",".join(chains), homooligomer=True,
        fix_pos=spec, inverse=True, rm_aa=rm_aa, verbose=False,
    )
    return [float(model.score(seq=sequence, key=model.key())["score"]) for sequence in sequences]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--positive-cif", action="append", type=pathlib.Path, required=True)
    ap.add_argument("--negative-cif", action="append", type=pathlib.Path, required=True)
    ap.add_argument("--candidates", type=pathlib.Path, required=True)
    ap.add_argument("--positions", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20260831)
    ap.add_argument("--rm-aa", default="C,P")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    backbone_dir = args.out / "backbones"; backbone_dir.mkdir(exist_ok=True)
    positions = [int(value) for value in args.positions.read_text().strip().split(",")]
    rows = list(csv.DictReader(args.candidates.open(), delimiter="\t"))
    sequences = [row["full_sequence"] for row in rows]
    if not rows or any(len(sequence) != 386 for sequence in sequences):
        raise SystemExit("candidate table must contain 386-aa full_sequence records")

    positive = []
    for index, cif in enumerate(args.positive_cif, 1):
        pdb = backbone_dir / f"positive_{index}.pdb"; cif_to_pdb(cif, pdb); positive.append(pdb)
    negative = []
    for index, cif in enumerate(args.negative_cif, 1):
        pdb = backbone_dir / f"negative_{index}.pdb"; cif_to_pdb(cif, pdb); negative.append(pdb)
    monomer = backbone_dir / "positive_chain_a_monomer.pdb"
    cif_to_pdb(args.positive_cif[0], monomer, {"A"})
    native = native_sequence(positive[0])
    if len(native) != 386:
        raise SystemExit(f"positive backbone chain has {len(native)} residues")

    model = mk_mpnn_model(
        model_name="v_48_020", backbone_noise=0.0, seed=args.seed,
        verbose=False, weights="original",
    )
    all_sequences = [native, *sequences]
    state_scores: dict[str, list[float]] = {}
    for index, pdb in enumerate(positive, 1):
        state_scores[f"positive_{index}"] = score_state(model, pdb, all_sequences, positions, args.rm_aa)
    for index, pdb in enumerate(negative, 1):
        state_scores[f"negative_{index}"] = score_state(model, pdb, all_sequences, positions, args.rm_aa)
    state_scores["monomer"] = score_state(model, monomer, all_sequences, positions, args.rm_aa)

    native_pos = statistics.mean(values[0] for key, values in state_scores.items() if key.startswith("positive_"))
    native_mono = state_scores["monomer"][0]
    for row_index, row in enumerate(rows, 1):
        pos = [values[row_index] for key, values in state_scores.items() if key.startswith("positive_")]
        neg = [values[row_index] for key, values in state_scores.items() if key.startswith("negative_")]
        monomer_score = state_scores["monomer"][row_index]
        row["positive_mean_score"] = statistics.mean(pos)
        row["positive_worst_score"] = max(pos)
        row["negative_best_score"] = min(neg)
        row["failed_state_specificity"] = min(neg) - max(pos)
        row["positive_gain_vs_native"] = native_pos - statistics.mean(pos)
        row["monomer_penalty_vs_native"] = monomer_score - native_mono
        row["lock_objective"] = (
            statistics.mean(pos)
            - 0.75 * (min(neg) - max(pos))
            + 0.50 * max(0.0, monomer_score - native_mono)
        )
        for key, values in state_scores.items():
            row[f"{key}_score"] = values[row_index]
    rows.sort(key=lambda row: (float(row["lock_objective"]), -float(row["failed_state_specificity"])))
    for rank, row in enumerate(rows, 1):
        row["lock_rank"] = rank
    fields = ["lock_rank", *[field for field in rows[0] if field != "lock_rank"]]
    with (args.out / "double_ring_multistate_all.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t"); writer.writeheader(); writer.writerows(rows)

    selected = []
    for row in rows:
        if float(row["positive_gain_vs_native"]) < 0.25:
            continue
        if float(row["monomer_penalty_vs_native"]) > 0.30:
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
    with (args.out / "double_ring_multistate_shortlist.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t"); writer.writeheader(); writer.writerows(selected)
    with (args.out / "double_ring_multistate_shortlist.fasta").open("w") as handle:
        for index, row in enumerate(selected, 1):
            handle.write(
                f">OVA-BODY-SS-AI-{index:02d}|rank={row['lock_rank']}|mut={row['mutations']}|"
                f"spec={float(row['failed_state_specificity']):.4f}\n{row['full_sequence']}\n"
            )
    with (args.out / "native_state_scores.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t"); writer.writerow(["state", "native_score"])
        for key, values in state_scores.items(): writer.writerow([key, values[0]])
    for index, row in enumerate(selected, 1):
        print(index, row["mutations"], f"pos={float(row['positive_mean_score']):.3f}",
              f"spec={float(row['failed_state_specificity']):.3f}",
              f"mono={float(row['monomer_penalty_vs_native']):.3f}")


if __name__ == "__main__":
    main()
