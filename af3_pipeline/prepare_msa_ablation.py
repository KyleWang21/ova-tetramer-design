#!/usr/bin/env python3
"""Prepare an AF3 control with OVA-only MSA and a gapped designed module."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def portable(path: Path) -> str:
    return str(path.resolve()).replace(
        "/vepfs-mlp2/c20250508/400083", "/root/400083", 1
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-json", type=Path, required=True)
    ap.add_argument("--source-manifest", type=Path, required=True)
    ap.add_argument("--ova-only-msa", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--seeds", default="1,2,3")
    args = ap.parse_args()
    if (args.out / "manifest.tsv").exists():
        raise SystemExit(f"experiment exists: {args.out}")

    payload = json.loads(args.source_json.read_text())
    old_name = payload["name"]
    query = payload["sequences"][0]["protein"]["sequence"]
    lines = args.ova_only_msa.read_text().splitlines()
    assert lines[0].startswith(">") and len(query) == len(lines[1])
    lines[1] = query

    msa_dir = args.out / "msas"
    input_dir = args.out / "in_s0"
    msa_dir.mkdir(parents=True)
    input_dir.mkdir()
    msa_path = msa_dir / f"{args.name}.a3m"
    msa_path.write_text("\n".join(lines) + "\n")
    payload["name"] = args.name
    payload["modelSeeds"] = [int(seed) for seed in args.seeds.split(",")]
    for entity in payload["sequences"]:
        protein = entity["protein"]
        protein["unpairedMsaPath"] = portable(msa_path)
        protein["pairedMsa"] = ""
        protein["templates"] = []
    json_path = input_dir / f"{args.name}.json"
    json_path.write_text(json.dumps(payload, indent=2) + "\n")

    rows = list(csv.DictReader(args.source_manifest.open(), delimiter="\t"))
    row = next(row for row in rows if row["name"] == old_name)
    row.update({
        "name": args.name,
        "shard": "0",
        "kind": "ova_only_msa_ablation",
        "seeds": str(len(payload["modelSeeds"])),
        "expected_samples": str(5 * len(payload["modelSeeds"])),
        "description": row["description"] + "; module MSA removed",
        "input_json": str(json_path),
    })
    with (args.out / "manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    print(f"prepared {args.name}; MSA lines={len(lines)}; query length={len(query)}")


if __name__ == "__main__":
    main()
