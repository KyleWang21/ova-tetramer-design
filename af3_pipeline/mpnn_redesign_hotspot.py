#!/usr/bin/env python3
"""Tied ProteinMPNN point redesign of the OVA surface hotspot in a C4 model."""

from __future__ import annotations

import argparse
import csv
import pathlib

from Bio.PDB import MMCIFParser, PDBIO
from colabdesign.mpnn import mk_mpnn_model


AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
       "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
       "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--positions", default="99,155,157,182,184,334,336,338")
    ap.add_argument("--chains", default="A,B,C,D")
    ap.add_argument("--temperatures", default="0.05,0.10,0.20")
    ap.add_argument("--samples-per-temperature", type=int, default=64)
    ap.add_argument("--seed", type=int, default=20260830)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    structure = MMCIFParser(QUIET=True).get_structure("c4", str(args.cif))
    pdb = args.out / "af3_c4_backbone.pdb"
    io=PDBIO();io.set_structure(structure);io.save(str(pdb))
    chains=args.chains.split(","); positions=[int(x) for x in args.positions.split(",")]
    design_spec=",".join(f"{c}{p}" for c in chains for p in positions)
    model=mk_mpnn_model(model_name="v_48_020",backbone_noise=0.0,seed=args.seed,weights="original")
    # inverse=True means the listed hotspot positions are the only designable sites.
    model.prep_inputs(pdb_filename=str(pdb),chain=args.chains,homooligomer=True,
                      fix_pos=design_spec,inverse=True,rm_aa="C",verbose=True)
    first=next(structure.get_chains())
    native="".join(AA3[r.resname] for r in first if r.id[0]==" " and r.resname in AA3)
    native_hot="".join(native[p-1] for p in positions)
    native_score=float(model.score(seq=native)["score"])
    rows=[]
    for temp in [float(x) for x in args.temperatures.split(",")]:
        result=model.sample_parallel(batch=args.samples_per_temperature,temperature=temp,rescore=True)
        for seq,score,seqid in zip(result["seq"],result["score"],result["seqid"]):
            seq=str(seq); hot="".join(seq[p-1] for p in positions)
            muts=[f"{native[p-1]}{p}{seq[p-1]}" for p in positions if native[p-1]!=seq[p-1]]
            rows.append({"temperature":temp,"mpnn_score":float(score),"seqid":float(seqid),
                         "n_mutations":len(muts),"mutations":",".join(muts),
                         "hotspot_sequence":hot,"hotspot_hydrophobic_count":sum(hot.count(x) for x in "AFILMVWY"),
                         "full_sequence":seq})
    dedup={}
    for r in rows:
        key=r["hotspot_sequence"]
        if key not in dedup or float(r["mpnn_score"])<float(dedup[key]["mpnn_score"]):dedup[key]=r
    rows=sorted(dedup.values(),key=lambda r:float(r["mpnn_score"]))
    for i,r in enumerate(rows,1):r["rank"]=i
    fields=["rank",*[k for k in rows[0] if k!="rank"]]
    with (args.out/"mpnn_hotspot_all.tsv").open("w",newline="") as h:
        w=csv.DictWriter(h,fields,delimiter="\t");w.writeheader();w.writerows(rows)
    shortlist=[]
    for r in rows:
        if 2<=int(r["n_mutations"])<=6 and int(r["hotspot_hydrophobic_count"])<=5:
            hot=str(r["hotspot_sequence"])
            if all(sum(a!=b for a,b in zip(hot,str(x["hotspot_sequence"])))>=2 for x in shortlist):
                shortlist.append(r)
        if len(shortlist)==12:break
    with (args.out/"mpnn_hotspot_shortlist.tsv").open("w",newline="") as h:
        w=csv.DictWriter(h,fields,delimiter="\t");w.writeheader();w.writerows(shortlist)
    with (args.out/"mpnn_hotspot_shortlist.fasta").open("w") as h:
        for i,r in enumerate(shortlist,1):h.write(f">mpnn_hotspot_{i:02d}|{r['mutations']}|score={float(r['mpnn_score']):.4f}\n{r['full_sequence']}\n")
    print("positions",positions,"native",native_hot,"native_score",round(native_score,4))
    print("unique",len(rows),"shortlist",len(shortlist))
    for i,r in enumerate(shortlist,1):print(i,round(float(r["mpnn_score"]),4),r["mutations"],r["hotspot_sequence"])


if __name__=="__main__":main()
