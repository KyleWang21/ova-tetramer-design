#!/usr/bin/env python3
"""Prepare AF3 stoichiometry controls for a small designed binder."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib

from prepare_inputs import protein


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binder-fasta", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--stoichiometries", default="2,3,4,5")
    ap.add_argument("--seeds", type=int, default=3)
    args = ap.parse_args()
    if (args.out / "manifest.tsv").exists(): raise SystemExit(f"experiment exists: {args.out}")
    seq = "".join(x.strip() for x in args.binder_fasta.read_text().splitlines()
                  if not x.startswith(">"))
    args.out.mkdir(parents=True); rows=[]; alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for shard, n in enumerate(map(int, args.stoichiometries.split(","))):
        name=f"bindcraft_patch1_binder_n{n}"; ids=alphabet[:n]
        d=args.out/f"in_s{shard}"; d.mkdir()
        payload={"name":name,"modelSeeds":list(range(1,args.seeds+1)),"dialect":"alphafold3",
                 "version":2,"sequences":[protein(c,seq) for c in ids]}
        jp=d/f"{name}.json"; jp.write_text(json.dumps(payload,indent=2)+"\n")
        rows.append({"name":name,"shard":shard,"kind":"binder_stoichiometry_control",
                     "n_chains":n,"chain_ids":ids,"chain_length":",".join([str(len(seq))]*n),
                     "ova_ranges":";".join(["0-0"]*n),
                     "module_ranges":";".join([f"1-{len(seq)}"]*n),
                     "architecture":f"BindCraft-binder x{n}","seeds":args.seeds,
                     "expected_samples":args.seeds*5,"description":"Binder homo-oligomerization control",
                     "input_json":str(jp),"sequence":seq})
    with (args.out/"manifest.tsv").open("w",newline="") as h:
        w=csv.DictWriter(h,list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
    (args.out/"representative_chain.fasta").write_text(f">BindCraft_patch1_binder\n{seq}\n")
    print(f"prepared {len(rows)} binder stoichiometry controls")


if __name__ == "__main__": main()
