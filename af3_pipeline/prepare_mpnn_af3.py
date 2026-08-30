#!/usr/bin/env python3
"""Prepare four-chain AF3 jobs for ProteinMPNN module redesigns."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib

from prepare_inputs import combine_domain_a3ms, portable, protein


def read_fasta(path: pathlib.Path) -> list[tuple[str, str]]:
    records = []; name = None; parts = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None: records.append((name, "".join(parts)))
            name = line[1:].split("|", 1)[0]; parts = []
        else: parts.append(line.strip())
    if name is not None: records.append((name, "".join(parts)))
    return records


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", type=pathlib.Path, required=True)
    ap.add_argument("--ova-msa", type=pathlib.Path, required=True)
    ap.add_argument("--module-msa", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--top", type=int, default=6)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--shards", type=int, default=6)
    ap.add_argument("--ova-end", type=int, default=386)
    ap.add_argument("--module-start", type=int, default=392)
    ap.add_argument("--module-end", type=int, default=424)
    ap.add_argument("--prefix", default="p3_pli_mpnn")
    ap.add_argument("--kind", default="mpnn_tetramer_candidate")
    ap.add_argument("--label", default="ProteinMPNN tied C4 module redesign")
    args = ap.parse_args()
    if (args.out / "manifest.tsv").exists(): raise SystemExit(f"experiment exists: {args.out}")
    records = read_fasta(args.fasta)[:args.top]
    args.out.mkdir(parents=True)
    for shard in range(args.shards): (args.out / f"in_s{shard}").mkdir()
    msas = args.out / "msas"; msas.mkdir()
    ova = args.ova_msa.read_text(); module = args.module_msa.read_text()
    rows=[]
    for index, (source_name, seq) in enumerate(records):
        if len(seq) != args.module_end: raise SystemExit(f"{source_name}: {len(seq)} aa")
        name = f"{args.prefix}_{index + 1:02d}"
        shard = index % args.shards
        a3m = msas / f"{name}.a3m"
        a3m.write_text(combine_domain_a3ms(
            ova, seq, (1,args.ova_end), module, (args.module_start,args.module_end)
        ))
        ids="ABCD"
        payload={"name":name,"modelSeeds":list(range(1,args.seeds+1)),"dialect":"alphafold3","version":2,
                 "sequences":[protein(c,seq,portable(a3m)) for c in ids]}
        jp=args.out/f"in_s{shard}"/f"{name}.json"; jp.write_text(json.dumps(payload,indent=2)+"\n")
        rows.append({"name":name,"shard":shard,"kind":args.kind,"n_chains":4,
                     "chain_ids":ids,"chain_length":",".join([str(len(seq))]*4),
                     "ova_ranges":";".join([f"1-{args.ova_end}"]*4),
                     "module_ranges":";".join([f"{args.module_start}-{args.module_end}"]*4),
                     "architecture":f"P3-GGGGS-{source_name}","seeds":args.seeds,
                     "expected_samples":args.seeds*5,"description":f"{args.label} {source_name}",
                     "input_json":str(jp),"sequence":seq})
    with (args.out/"manifest.tsv").open("w",newline="") as h:
        w=csv.DictWriter(h,list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
    with (args.out/"representative_chains.fasta").open("w") as h:
        for r in rows:h.write(f">{r['name']} | {r['architecture']}\n{r['sequence']}\n")
    print(f"prepared {len(rows)} MPNN candidates across {args.shards} shards")


if __name__ == "__main__": main()
