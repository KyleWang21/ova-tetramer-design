#!/usr/bin/env python3
"""Rank pure-OVA interface designs by C4-positive/C2-C3-dimer-negative MPNN scores."""

from __future__ import annotations

import argparse
import csv
import pathlib

from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select
from colabdesign.mpnn import mk_mpnn_model


class CropSelect(Select):
    def __init__(self, end: int, chains: set[str] | None = None) -> None:
        self.end = end
        self.chains = chains

    def accept_chain(self, chain) -> bool:  # noqa: ANN001
        return self.chains is None or str(chain.id) in self.chains

    def accept_residue(self, residue) -> bool:  # noqa: ANN001
        return residue.id[0] == " " and int(residue.id[1]) <= self.end


def cif_to_pdb(cif: pathlib.Path, pdb: pathlib.Path, end: int, chains: set[str] | None = None) -> None:
    structure = MMCIFParser(QUIET=True).get_structure(pdb.stem, str(cif))
    io = PDBIO(); io.set_structure(structure); io.save(str(pdb), CropSelect(end, chains))


def chains_in_pdb(pdb: pathlib.Path) -> list[str]:
    structure = PDBParser(QUIET=True).get_structure(pdb.stem, str(pdb))
    return [str(chain.id) for chain in next(structure.get_models()).get_chains()]


def score_state(
    model, pdb: pathlib.Path, sequences: list[str], positions: list[int]
) -> list[float]:  # noqa: ANN001
    chains = chains_in_pdb(pdb)
    design_spec = ",".join(f"{chain}{pos}" for chain in chains for pos in positions)
    model.prep_inputs(
        pdb_filename=str(pdb), chain=",".join(chains), homooligomer=True,
        fix_pos=design_spec, inverse=True, rm_aa="C,P", verbose=False,
    )
    return [float(model.score(seq=sequence, key=model.key())["score"]) for sequence in sequences]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--c4-pdb", type=pathlib.Path, required=True)
    ap.add_argument("--c2-cif", type=pathlib.Path, required=True)
    ap.add_argument("--c3-cif", type=pathlib.Path, required=True)
    ap.add_argument("--dimer-cif", type=pathlib.Path, required=True)
    ap.add_argument("--candidates", type=pathlib.Path, required=True)
    ap.add_argument("--positions", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--ova-end", type=int, default=386)
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20260830)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    positions = [int(value) for value in args.positions.read_text().strip().split(",")]
    rows = list(csv.DictReader(args.candidates.open(), delimiter="\t"))
    sequences = [row["full_sequence"] for row in rows]
    if not rows or any(len(sequence) != args.ova_end for sequence in sequences):
        raise SystemExit("candidate table must contain non-empty 386-aa full_sequence records")

    backbones = args.out / "backbones"; backbones.mkdir(exist_ok=True)
    c2_pdb = backbones / "c2_ova_only.pdb"
    c3_pdb = backbones / "c3_ova_only.pdb"
    dimer_pdb = backbones / "p3_dimer_ova_only.pdb"
    monomer_pdb = backbones / "c4_chain_a_monomer.pdb"
    cif_to_pdb(args.c2_cif, c2_pdb, args.ova_end)
    cif_to_pdb(args.c3_cif, c3_pdb, args.ova_end)
    cif_to_pdb(args.dimer_cif, dimer_pdb, args.ova_end)
    c4_structure = PDBParser(QUIET=True).get_structure("c4", str(args.c4_pdb))
    io = PDBIO(); io.set_structure(c4_structure)
    io.save(str(monomer_pdb), CropSelect(args.ova_end, {"A"}))

    native_structure = PDBParser(QUIET=True).get_structure("native", str(args.c4_pdb))
    aa3 = {
        "ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D", "CYS":"C",
        "GLN":"Q", "GLU":"E", "GLY":"G", "HIS":"H", "ILE":"I",
        "LEU":"L", "LYS":"K", "MET":"M", "PHE":"F", "PRO":"P",
        "SER":"S", "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V",
    }
    native_chain = next(native_structure.get_chains())
    native = "".join(aa3[res.resname] for res in native_chain if res.id[0] == " ")
    if len(native) != args.ova_end:
        raise SystemExit(f"native C4 backbone chain has length {len(native)}, expected {args.ova_end}")

    states = {
        "c4": args.c4_pdb,
        "c2": c2_pdb,
        "c3": c3_pdb,
        "dimer": dimer_pdb,
        "monomer": monomer_pdb,
    }
    model = mk_mpnn_model(
        model_name="v_48_020", backbone_noise=0.0, seed=args.seed,
        verbose=False, weights="original",
    )
    native_scores: dict[str, float] = {}
    for state, pdb in states.items():
        values = score_state(model, pdb, [native, *sequences], positions)
        native_scores[state] = values[0]
        for row, value in zip(rows, values[1:], strict=True):
            row[f"{state}_mpnn_score"] = value

    for row in rows:
        c4 = float(row["c4_mpnn_score"])
        deltas = [float(row[f"{state}_mpnn_score"]) - c4 for state in ("c2", "c3", "dimer")]
        row["c2_minus_c4"] = deltas[0]
        row["c3_minus_c4"] = deltas[1]
        row["dimer_minus_c4"] = deltas[2]
        row["min_negative_delta"] = min(deltas)
        row["c4_gain_vs_native"] = native_scores["c4"] - c4
        row["monomer_penalty_vs_native"] = float(row["monomer_mpnn_score"]) - native_scores["monomer"]
        row["multistate_objective"] = (
            c4
            - 0.65 * float(row["min_negative_delta"])
            + 0.50 * max(0.0, float(row["monomer_penalty_vs_native"]))
        )
    rows.sort(key=lambda row: (
        float(row["multistate_objective"]),
        -float(row["min_negative_delta"]),
        float(row["c4_mpnn_score"]),
    ))
    for rank, row in enumerate(rows, 1):
        row["multistate_rank"] = rank
    fields = ["multistate_rank", *[field for field in rows[0] if field != "multistate_rank"]]
    with (args.out / "multistate_all.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)

    selected: list[dict[str, str]] = []
    for row in rows:
        if float(row["c4_gain_vs_native"]) < 0.20:
            continue
        if float(row["monomer_penalty_vs_native"]) > 0.35:
            continue
        sequence = row["interface_sequence"]
        if all(sum(a != b for a, b in zip(sequence, old["interface_sequence"])) >= 2 for old in selected):
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
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader(); writer.writerows(selected)
    with (args.out / "multistate_shortlist.fasta").open("w") as handle:
        for index, row in enumerate(selected, 1):
            handle.write(
                f">OVA-BODY-C4-{index:02d}|msrank={row['multistate_rank']}|mut={row['mutations']}|"
                f"delta={float(row['min_negative_delta']):.4f}\n{row['full_sequence']}\n"
            )
    with (args.out / "native_state_scores.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["state", "native_mpnn_score", "backbone"])
        for state, pdb in states.items():
            writer.writerow([state, native_scores[state], pdb])
    for index, row in enumerate(selected, 1):
        print(
            index, row["mutations"],
            f"C4={float(row['c4_mpnn_score']):.4f}",
            f"minNegDelta={float(row['min_negative_delta']):.4f}",
            f"monoPenalty={float(row['monomer_penalty_vs_native']):.4f}",
        )


if __name__ == "__main__":
    main()
