#!/usr/bin/env python3
"""ProteinMPNN redesign of a tied homotetramer module on an AF3 backbone."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib

from Bio.PDB import MMCIFParser, PDBIO
from colabdesign.mpnn import mk_mpnn_model


AA3 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
    "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
    "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}


def identity(a: str, b: str) -> float:
    return sum(x == y for x, y in zip(a, b)) / len(a)


def net_charge(seq: str) -> int:
    return sum(seq.count(x) for x in "KR") - sum(seq.count(x) for x in "DE")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--fixed-end", type=int, default=391)
    ap.add_argument("--module-start", type=int, default=392)
    ap.add_argument("--module-end", type=int, default=424)
    ap.add_argument("--chains", default="A,B,C,D")
    ap.add_argument("--temperatures", default="0.05,0.10,0.20")
    ap.add_argument("--samples-per-temperature", type=int, default=48)
    ap.add_argument("--seed", type=int, default=20260829)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    pdb = args.out / "af3_tetramer_backbone.pdb"
    structure = MMCIFParser(QUIET=True).get_structure("tetramer", str(args.cif))
    io = PDBIO(); io.set_structure(structure); io.save(str(pdb))
    chain_ids = args.chains.split(",")
    fix_pos = ",".join(f"{chain}1-{chain}{args.fixed_end}" for chain in chain_ids)

    model = mk_mpnn_model(
        model_name="v_48_020", backbone_noise=0.0, seed=args.seed,
        verbose=False, weights="original",
    )
    model.prep_inputs(
        pdb_filename=str(pdb), chain=args.chains, homooligomer=True,
        fix_pos=fix_pos, rm_aa="C", verbose=True,
    )
    first_chain = next(structure.get_chains())
    native = "".join(AA3[r.resname] for r in first_chain if r.id[0] == " " and r.resname in AA3)
    native_module = native[args.module_start - 1:args.module_end]
    native_score = float(model.score(seq=native, key=model.key())["score"])
    (args.out / "native_control.json").write_text(json.dumps({
        "full_sequence": native, "module_sequence": native_module,
        "mpnn_score": native_score,
    }, indent=2) + "\n")
    rows: list[dict[str, object]] = []
    for temperature in [float(x) for x in args.temperatures.split(",")]:
        result = model.sample_parallel(
            batch=args.samples_per_temperature, temperature=temperature, rescore=True
        )
        for seq, score, seqid in zip(result["seq"], result["score"], result["seqid"]):
            module = str(seq)[args.module_start - 1:args.module_end]
            rows.append({
                "temperature": temperature,
                "mpnn_score": float(score),
                "mpnn_seqid_design_positions": float(seqid),
                "module_identity_to_pLI": identity(module, native_module),
                "module_mutations": sum(a != b for a, b in zip(module, native_module)),
                "module_net_charge": net_charge(module),
                "module_hydrophobic_fraction": sum(module.count(x) for x in "AFILMVWY") / len(module),
                "module_sequence": module,
                "full_sequence": str(seq),
            })
    dedup = {}
    for row in rows:
        seq = str(row["module_sequence"])
        if seq not in dedup or float(row["mpnn_score"]) < float(dedup[seq]["mpnn_score"]):
            dedup[seq] = row
    rows = sorted(dedup.values(), key=lambda r: float(r["mpnn_score"]))
    for rank, row in enumerate(rows, 1):
        row["rank"] = rank
    fields = ["rank", *[k for k in rows[0] if k != "rank"]]
    with (args.out / "mpnn_all.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)

    shortlist = []
    for row in rows:
        ident = float(row["module_identity_to_pLI"])
        charge = abs(int(row["module_net_charge"]))
        if not (0.35 <= ident <= 0.90 and charge <= 12):
            continue
        module = str(row["module_sequence"])
        if all(sum(a != b for a, b in zip(module, str(old["module_sequence"]))) >= 3 for old in shortlist):
            shortlist.append(row)
        if len(shortlist) == 12:
            break
    with (args.out / "mpnn_shortlist.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader(); writer.writerows(shortlist)
    with (args.out / "mpnn_shortlist.fasta").open("w") as handle:
        for i, row in enumerate(shortlist, 1):
            handle.write(
                f">mpnn_pLI_{i:02d}|score={float(row['mpnn_score']):.4f}|"
                f"mut={row['module_mutations']}|module={row['module_sequence']}\n"
                f"{row['full_sequence']}\n"
            )
    print(f"native module: {native_module}; MPNN score: {native_score:.4f}")
    print(f"unique designs: {len(rows)}; shortlist: {len(shortlist)}")
    for i, row in enumerate(shortlist, 1):
        print(i, round(float(row["mpnn_score"]), 4), row["module_mutations"], row["module_sequence"])


if __name__ == "__main__":
    main()
