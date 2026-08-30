#!/usr/bin/env python3
"""Prepare four-chain AF3 screens of a BindCraft binder fused to P3-13R."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib

from prepare_inputs import GCN4_PLI, combine_domain_a3ms, portable, protein
from build_tetramer_design import make_constructs


G5 = "GGGGS"
G10 = G5 * 2
G15 = G5 * 3


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binder-fasta", type=pathlib.Path, required=True)
    ap.add_argument("--ova-msa", type=pathlib.Path, required=True)
    ap.add_argument("--pli-msa", type=pathlib.Path, required=True)
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

    def describe(name: str, parts: list[tuple[str, str]]) -> dict[str, object]:
        seq = "".join(x[1] for x in parts)
        cursor = 1; ranges: dict[str, tuple[int, int]] = {}
        for label, part in parts:
            lo, hi = cursor, cursor + len(part) - 1
            if label in {"ova", "binder", "pli"}:
                ranges[label] = (lo, hi)
            cursor = hi + 1
        return {"name": name, "seq": seq, "ranges": ranges,
                "architecture": "-".join(x[0] for x in parts)}

    systems = [
        describe("p3_bc70_pli_mid_g5", [("ova", p3), ("g5", G5), ("binder", binder), ("g5", G5), ("pli", GCN4_PLI)]),
        describe("bc70_p3_pli_n_g5", [("binder", binder), ("g5", G5), ("ova", p3), ("g5", G5), ("pli", GCN4_PLI)]),
        describe("bc70_p3_pli_n_g15", [("binder", binder), ("g15", G15), ("ova", p3), ("g5", G5), ("pli", GCN4_PLI)]),
        describe("bc70_p3_n_g15", [("binder", binder), ("g15", G15), ("ova", p3)]),
        describe("p3_bc70_c_g5", [("ova", p3), ("g5", G5), ("binder", binder)]),
        describe("p3_bc70_c_g10", [("ova", p3), ("g10", G10), ("binder", binder)]),
    ]

    args.out.mkdir(parents=True)
    msas = args.out / "msas"; msas.mkdir()
    for i in range(len(systems)):
        (args.out / f"in_s{i}").mkdir()
    ova_msa = args.ova_msa.read_text(); pli_msa = args.pli_msa.read_text()
    rows = []
    for shard, system in enumerate(systems):
        seq = str(system["seq"]); ranges = dict(system["ranges"])
        a3m = msas / f"{system['name']}.a3m"
        a3m.write_text(combine_domain_a3ms(
            ova_msa, seq, ranges["ova"],
            pli_msa if "pli" in ranges else None,
            ranges.get("pli", (0, 0)),
        ))
        ids = "ABCD"
        payload = {"name": system["name"], "modelSeeds": list(range(1, args.seeds + 1)),
                   "dialect": "alphafold3", "version": 2,
                   "sequences": [protein(c, seq, portable(a3m)) for c in ids]}
        jp = args.out / f"in_s{shard}" / f"{system['name']}.json"
        jp.write_text(json.dumps(payload, indent=2) + "\n")
        primary_module = ranges.get("pli", ranges["binder"])
        rows.append({
            "name": system["name"], "shard": shard, "kind": "bindcraft_fusion_tetramer",
            "n_chains": 4, "chain_ids": ids, "chain_length": ",".join([str(len(seq))] * 4),
            "ova_ranges": ";".join([f"{ranges['ova'][0]}-{ranges['ova'][1]}"] * 4),
            "module_ranges": ";".join([f"{primary_module[0]}-{primary_module[1]}"] * 4),
            "binder_ranges": ";".join([f"{ranges['binder'][0]}-{ranges['binder'][1]}"] * 4),
            "architecture": system["architecture"], "seeds": args.seeds,
            "expected_samples": args.seeds * 5,
            "description": "Four-chain P3 fusion with AI-designed hotspot binder",
            "input_json": str(jp), "sequence": seq,
        })
    with (args.out / "manifest.tsv").open("w", newline="") as h:
        w = csv.DictWriter(h, list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    with (args.out / "representative_chains.fasta").open("w") as h:
        for row in rows: h.write(f">{row['name']} | {row['architecture']}\n{row['sequence']}\n")
    print(f"prepared {len(rows)} BindCraft fusion candidates")


if __name__ == "__main__":
    main()
