#!/usr/bin/env python3
"""Tied ProteinMPNN redesign of the recurrent, pure-OVA C4 interface."""

from __future__ import annotations

import argparse
import csv
import pathlib

from Bio.PDB import PDBParser
from colabdesign.mpnn import mk_mpnn_model


AA3 = {
    "ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D", "CYS":"C",
    "GLN":"Q", "GLU":"E", "GLY":"G", "HIS":"H", "ILE":"I",
    "LEU":"L", "LYS":"K", "MET":"M", "PHE":"F", "PRO":"P",
    "SER":"S", "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V",
}


def charge(sequence: str) -> int:
    return sum(sequence.count(aa) for aa in "KR") - sum(sequence.count(aa) for aa in "DE")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", type=pathlib.Path, required=True)
    ap.add_argument("--positions", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--chains", default="A,B,C,D")
    ap.add_argument("--temperatures", default="0.05,0.10,0.20")
    ap.add_argument("--samples-per-temperature", type=int, default=128)
    ap.add_argument("--seed", type=int, default=20260829)
    ap.add_argument("--shortlist", type=int, default=48)
    ap.add_argument("--min-mutations", type=int, default=4)
    ap.add_argument("--max-mutations", type=int, default=12)
    ap.add_argument("--max-abs-charge", type=int, default=7)
    ap.add_argument("--max-hydrophobic", type=int, default=11)
    ap.add_argument("--rm-aa", default="C,P", help="Amino acids excluded from designed positions")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    positions = [int(x) for x in args.positions.read_text().strip().split(",")]
    chains = args.chains.split(",")

    structure = PDBParser(QUIET=True).get_structure("ova_c4", str(args.pdb))
    first = next(structure.get_chains())
    native = "".join(AA3[res.resname] for res in first if res.id[0] == " " and res.resname in AA3)
    if len(native) != 386:
        raise SystemExit(f"expected a 386-aa OVA chain, found {len(native)}")
    design_spec = ",".join(f"{chain}{pos}" for chain in chains for pos in positions)
    model = mk_mpnn_model(
        model_name="v_48_020", backbone_noise=0.0, seed=args.seed,
        verbose=False, weights="original",
    )
    model.prep_inputs(
        pdb_filename=str(args.pdb), chain=args.chains, homooligomer=True,
        fix_pos=design_spec, inverse=True, rm_aa=args.rm_aa, verbose=True,
    )
    native_score = float(model.score(seq=native, key=model.key())["score"])
    native_interface = "".join(native[pos - 1] for pos in positions)
    rows: list[dict[str, object]] = []
    for temperature in [float(value) for value in args.temperatures.split(",")]:
        result = model.sample_parallel(
            batch=args.samples_per_temperature, temperature=temperature, rescore=True
        )
        for sequence, score, seqid in zip(result["seq"], result["score"], result["seqid"]):
            sequence = str(sequence)
            interface = "".join(sequence[pos - 1] for pos in positions)
            mutations = [
                f"{native[pos - 1]}{pos}{sequence[pos - 1]}"
                for pos in positions if sequence[pos - 1] != native[pos - 1]
            ]
            rows.append({
                "temperature": temperature,
                "mpnn_c4_score": float(score),
                "mpnn_seqid_design_positions": float(seqid),
                "n_mutations": len(mutations),
                "mutations": ",".join(mutations),
                "interface_sequence": interface,
                "interface_net_charge": charge(interface),
                "interface_hydrophobic_count": sum(interface.count(aa) for aa in "AFILMVWY"),
                "full_sequence": sequence,
            })
    deduplicated: dict[str, dict[str, object]] = {}
    for row in rows:
        key = str(row["interface_sequence"])
        if key not in deduplicated or float(row["mpnn_c4_score"]) < float(deduplicated[key]["mpnn_c4_score"]):
            deduplicated[key] = row
    rows = sorted(deduplicated.values(), key=lambda row: float(row["mpnn_c4_score"]))
    for rank, row in enumerate(rows, 1):
        row["c4_rank"] = rank
    fields = ["c4_rank", *[field for field in rows[0] if field != "c4_rank"]]
    with (args.out / "mpnn_all.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)

    shortlist: list[dict[str, object]] = []
    for row in rows:
        n_mutations = int(row["n_mutations"])
        if not (args.min_mutations <= n_mutations <= args.max_mutations):
            continue
        if abs(int(row["interface_net_charge"])) > args.max_abs_charge:
            continue
        if int(row["interface_hydrophobic_count"]) > args.max_hydrophobic:
            continue
        interface = str(row["interface_sequence"])
        if all(
            sum(a != b for a, b in zip(interface, str(old["interface_sequence"]))) >= 3
            for old in shortlist
        ):
            shortlist.append(row)
        if len(shortlist) == args.shortlist:
            break
    with (args.out / "mpnn_shortlist.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader(); writer.writerows(shortlist)
    with (args.out / "mpnn_shortlist.fasta").open("w") as handle:
        for index, row in enumerate(shortlist, 1):
            handle.write(
                f">ova_c4body_{index:03d}|score={float(row['mpnn_c4_score']):.4f}|"
                f"mut={row['mutations']}\n{row['full_sequence']}\n"
            )
    (args.out / "native_control.tsv").write_text(
        "mpnn_c4_score\tinterface_sequence\tpositions\n"
        f"{native_score}\t{native_interface}\t{','.join(map(str, positions))}\n"
    )
    print(f"native C4 score: {native_score:.4f}; interface: {native_interface}")
    print(f"unique designs: {len(rows)}; shortlist: {len(shortlist)}")
    for index, row in enumerate(shortlist[:12], 1):
        print(index, round(float(row["mpnn_c4_score"]), 4), row["n_mutations"], row["mutations"])


if __name__ == "__main__":
    main()
