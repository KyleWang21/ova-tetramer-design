#!/usr/bin/env python3
"""Prepare AF3 1--6 chain stoichiometry competition for one fusion sequence."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib

from prepare_inputs import combine_domain_a3ms, portable, protein


def read_first_fasta(path: pathlib.Path) -> str:
    parts = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if parts:
                break
        else:
            parts.append(line.strip())
    return "".join(parts)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sequence-fasta", type=pathlib.Path, required=True)
    ap.add_argument("--ova-msa", type=pathlib.Path, required=True)
    ap.add_argument("--module-msa", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--ova-end", type=int, default=386)
    ap.add_argument("--module-start", type=int, default=392)
    ap.add_argument("--module-end", type=int, default=424)
    ap.add_argument("--stoichiometries", default="1,2,3,4,5,6")
    ap.add_argument("--seeds", type=int, default=3)
    args = ap.parse_args()
    if (args.out / "manifest.tsv").exists():
        raise SystemExit(f"experiment exists: {args.out}")
    seq = read_first_fasta(args.sequence_fasta)
    if len(seq) != args.module_end:
        raise SystemExit(f"expected {args.module_end} aa, found {len(seq)}")
    args.out.mkdir(parents=True)
    msa_dir = args.out / "msas"; msa_dir.mkdir()
    combined = combine_domain_a3ms(
        args.ova_msa.read_text(), seq, (1, args.ova_end),
        args.module_msa.read_text(), (args.module_start, args.module_end),
    )
    msa_path = msa_dir / "p3_pli_f05.a3m"; msa_path.write_text(combined)
    rows = []
    chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for shard, n in enumerate(int(x) for x in args.stoichiometries.split(",")):
        name = f"p3_pli_f05_n{n}"
        indir = args.out / f"in_s{shard}"; indir.mkdir()
        ids = chain_ids[:n]
        payload = {
            "name": name, "modelSeeds": list(range(1, args.seeds + 1)),
            "dialect": "alphafold3", "version": 2,
            "sequences": [protein(cid, seq, portable(msa_path)) for cid in ids],
        }
        input_json = indir / f"{name}.json"
        input_json.write_text(json.dumps(payload, indent=2) + "\n")
        rows.append({
            "name": name, "shard": shard, "kind": "stoichiometry_competition",
            "n_chains": n, "chain_ids": ids, "chain_length": ",".join([str(len(seq))] * n),
            "ova_ranges": ";".join([f"1-{args.ova_end}"] * n),
            "module_ranges": ";".join([f"{args.module_start}-{args.module_end}"] * n),
            "architecture": "P3-GGGGS-pLI", "seeds": args.seeds,
            "expected_samples": args.seeds * 5,
            "description": f"P3-GGGGS-pLI forced input stoichiometry n={n}",
            "input_json": str(input_json),
        })
    with (args.out / "manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    (args.out / "representative_chain.fasta").write_text(f">P3-GGGGS-pLI\n{seq}\n")
    print(f"prepared stoichiometries {[r['n_chains'] for r in rows]} in {args.out}")


if __name__ == "__main__":
    main()
