#!/usr/bin/env python3
"""Rank module sequences by ProteinMPNN preference for C4 over C3 backbones."""

from __future__ import annotations

import argparse
import csv
import pathlib

from Bio.PDB import MMCIFParser, PDBIO
from colabdesign.mpnn import mk_mpnn_model


def to_pdb(cif: pathlib.Path, pdb: pathlib.Path) -> None:
    structure=MMCIFParser(QUIET=True).get_structure("x",str(cif))
    io=PDBIO();io.set_structure(structure);io.save(str(pdb))


def score_backbone(model, pdb: pathlib.Path, chains: str, sequences: list[str], fixed_end: int) -> list[float]:
    ids=chains.split(",");fix=",".join(f"{c}1-{c}{fixed_end}" for c in ids)
    model.prep_inputs(pdb_filename=str(pdb),chain=chains,homooligomer=True,
                      fix_pos=fix,rm_aa="C")
    return [float(model.score(seq=seq)["score"]) for seq in sequences]


def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument("--c3-cif",type=pathlib.Path,required=True)
    ap.add_argument("--c4-cif",type=pathlib.Path,required=True)
    ap.add_argument("--c4-low-cif",type=pathlib.Path,
                    help="Optional competing low-confidence C4 backbone for robust multistate ranking")
    ap.add_argument("--candidates",type=pathlib.Path,required=True)
    ap.add_argument("--out",type=pathlib.Path,required=True)
    ap.add_argument("--fixed-end",type=int,default=391)
    ap.add_argument("--top",type=int,default=12)
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    rows=list(csv.DictReader(args.candidates.open(),delimiter="\t"))
    sequences=[r["full_sequence"] for r in rows]
    p3=args.out/"c3_backbone.pdb";p4=args.out/"c4_backbone.pdb"
    to_pdb(args.c3_cif,p3);to_pdb(args.c4_cif,p4)
    model=mk_mpnn_model(model_name="v_48_020",backbone_noise=0.0,seed=20260831,weights="original")
    c3=score_backbone(model,p3,"A,B,C",sequences,args.fixed_end)
    c4=score_backbone(model,p4,"A,B,C,D",sequences,args.fixed_end)
    c4_low = None
    if args.c4_low_cif:
        p4_low=args.out/"c4_low_backbone.pdb";to_pdb(args.c4_low_cif,p4_low)
        c4_low=score_backbone(model,p4_low,"A,B,C,D",sequences,args.fixed_end)
    for idx,(r,s3,s4) in enumerate(zip(rows,c3,c4)):
        r["c3_mpnn_score"]=s3;r["c4_rescored_mpnn_score"]=s4
        if c4_low is not None:
            r["c4_low_mpnn_score"]=c4_low[idx]
            r["c4_worst_mpnn_score"]=max(s4,c4_low[idx])
            r["c4_mean_mpnn_score"]=(s4+c4_low[idx])/2
            r["c4_preference_delta_c3_minus_c4_worst"]=s3-max(s4,c4_low[idx])
        r["c4_preference_delta_c3_minus_c4"]=s3-s4
    if c4_low is not None:
        rows.sort(key=lambda r:(float(r["c4_preference_delta_c3_minus_c4_worst"]),
                                -float(r["c4_worst_mpnn_score"])),reverse=True)
    else:
        rows.sort(key=lambda r:(float(r["c4_preference_delta_c3_minus_c4"]),-float(r["c4_rescored_mpnn_score"])),reverse=True)
    for i,r in enumerate(rows,1):r["specificity_rank"]=i
    fields=["specificity_rank",*[k for k in rows[0] if k!="specificity_rank"]]
    with (args.out/"specificity_all.tsv").open("w",newline="") as h:
        w=csv.DictWriter(h,fields,delimiter="\t");w.writeheader();w.writerows(rows)
    shortlist=rows[:args.top]
    with (args.out/"specificity_shortlist.tsv").open("w",newline="") as h:
        w=csv.DictWriter(h,fields,delimiter="\t");w.writeheader();w.writerows(shortlist)
    with (args.out/"specificity_shortlist.fasta").open("w") as h:
        for i,r in enumerate(shortlist,1):
            h.write(f">mpnn_c4spec_{i:02d}|delta={float(r['c4_preference_delta_c3_minus_c4']):.4f}|module={r['module_sequence']}\n{r['full_sequence']}\n")
    for i,r in enumerate(shortlist,1):
        print(i,round(float(r["c4_preference_delta_c3_minus_c4"]),4),
              round(float(r["c3_mpnn_score"]),4),round(float(r["c4_rescored_mpnn_score"]),4),r["module_sequence"])


if __name__=="__main__":main()
