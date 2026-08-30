#!/usr/bin/env python3
"""Prepare an independent AF3 check of a BindCraft OVA binder."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib

from prepare_inputs import combine_domain_a3ms, portable, protein
from build_tetramer_design import make_constructs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binder-fasta", type=pathlib.Path, required=True)
    ap.add_argument("--ova-msa", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--seeds", type=int, default=3)
    args = ap.parse_args()
    if (args.out / "manifest.tsv").exists():
        raise SystemExit(f"experiment exists: {args.out}")

    constructs, _, _ = make_constructs()
    p3 = next(str(x["aa"]) for x in constructs if x["id"] == "D0-P3-13R")
    binder = "".join(
        line.strip() for line in args.binder_fasta.read_text().splitlines()
        if not line.startswith(">")
    )
    args.out.mkdir(parents=True)
    (args.out / "in_s0").mkdir()
    (args.out / "msas").mkdir()
    ova_a3m = args.out / "msas" / "p3.a3m"
    ova_a3m.write_text(combine_domain_a3ms(args.ova_msa.read_text(), p3, (1, len(p3))))

    systems = [
        {
            "name": "bindcraft_patch1_binder_monomer",
            "kind": "bindcraft_binder_fold",
            "proteins": [protein("A", binder)],
            "chain_ids": "A",
            "lengths": [len(binder)],
            "ova_ranges": [(0, 0)],
            "module_ranges": [(1, len(binder))],
            "description": "AF3 fold check of the ProteinMPNN-refined BindCraft binder",
        },
        {
            "name": "p3_bindcraft_patch1_complex",
            "kind": "bindcraft_target_complex",
            "proteins": [protein("A", p3, portable(ova_a3m)), protein("B", binder)],
            "chain_ids": "AB",
            "lengths": [len(p3), len(binder)],
            "ova_ranges": [(1, len(p3)), (0, 0)],
            "module_ranges": [(0, 0), (1, len(binder))],
            "description": "Independent full-length P3 plus BindCraft binder complex",
        },
    ]
    rows = []
    for system in systems:
        payload = {
            "name": system["name"],
            "modelSeeds": list(range(1, args.seeds + 1)),
            "dialect": "alphafold3",
            "version": 2,
            "sequences": system["proteins"],
        }
        jp = args.out / "in_s0" / f"{system['name']}.json"
        jp.write_text(json.dumps(payload, indent=2) + "\n")
        rows.append({
            "name": system["name"], "shard": 0, "kind": system["kind"],
            "n_chains": len(system["proteins"]), "chain_ids": system["chain_ids"],
            "chain_length": ",".join(map(str, system["lengths"])),
            "ova_ranges": ";".join(f"{a}-{b}" for a, b in system["ova_ranges"]),
            "module_ranges": ";".join(f"{a}-{b}" for a, b in system["module_ranges"]),
            "architecture": system["name"], "seeds": args.seeds,
            "expected_samples": args.seeds * 5, "description": system["description"],
            "input_json": str(jp), "sequence": binder,
        })
    with (args.out / "manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    (args.out / "representative_chains.fasta").write_text(
        f">P3-13R\n{p3}\n>BindCraft_patch1_mpnn1\n{binder}\n"
    )
    print(f"prepared {len(systems)} systems; binder={len(binder)} aa")


if __name__ == "__main__":
    main()
