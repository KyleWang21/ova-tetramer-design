#!/usr/bin/env python3
"""Prepare chain-count competition from one existing homomer AF3 input."""

from __future__ import annotations

import argparse
import csv
import copy
import json
import pathlib


def first_range(spec: str) -> str:
    return spec.split(";", 1)[0]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-json", type=pathlib.Path, required=True)
    ap.add_argument("--source-manifest", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--stoichiometries", default="2,3,4,5")
    ap.add_argument("--seeds", default="1,2,3")
    args = ap.parse_args()
    if (args.out / "manifest.tsv").exists(): raise SystemExit(f"experiment exists: {args.out}")
    source = json.loads(args.source_json.read_text()); old_name = source["name"]
    source_rows = list(csv.DictReader(args.source_manifest.open(), delimiter="\t"))
    base = next(r for r in source_rows if r["name"] == old_name)
    first_protein = source["sequences"][0]; seed_values = [int(x) for x in args.seeds.split(",")]
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"; args.out.mkdir(parents=True); rows=[]
    for shard, n in enumerate(map(int, args.stoichiometries.split(","))):
        name=f"{old_name}_n{n}"; d=args.out/f"in_s{shard}"; d.mkdir()
        proteins=[]
        for cid in alphabet[:n]:
            p=copy.deepcopy(first_protein); p["protein"]["id"]=cid; proteins.append(p)
        payload={"name":name,"modelSeeds":seed_values,"dialect":"alphafold3",
                 "version":source.get("version",2),"sequences":proteins}
        jp=d/f"{name}.json"; jp.write_text(json.dumps(payload,indent=2)+"\n")
        row=dict(base); row.update({"name":name,"shard":str(shard),"n_chains":str(n),
            "candidate":base.get("candidate", old_name),
            "chain_ids":alphabet[:n],"chain_length":",".join([str(len(first_protein['protein']['sequence']))]*n),
            "ova_ranges":";".join([first_range(base["ova_ranges"])]*n),
            "module_ranges":";".join([first_range(base["module_ranges"])]*n),
            "seeds":str(len(seed_values)),"expected_samples":str(len(seed_values)*5),
            "description":base["description"]+f"; forced n={n}","input_json":str(jp)})
        if "binder_ranges" in row: row["binder_ranges"]=";".join([first_range(base["binder_ranges"])]*n)
        rows.append(row)
    with (args.out/"manifest.tsv").open("w",newline="") as h:
        w=csv.DictWriter(h,list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
    print(f"prepared stoichiometries {[r['n_chains'] for r in rows]}")


if __name__ == "__main__": main()
