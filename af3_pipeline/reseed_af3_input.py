#!/usr/bin/env python3
"""Create an independent AF3 reseed experiment from one existing input."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-json", type=pathlib.Path, required=True)
    ap.add_argument("--source-manifest", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--seeds", required=True, help="comma-separated AF3 seeds")
    args = ap.parse_args()
    if (args.out / "manifest.tsv").exists(): raise SystemExit(f"experiment exists: {args.out}")
    seed_values = [int(x) for x in args.seeds.split(",")]
    payload = json.loads(args.source_json.read_text()); old_name = payload["name"]
    payload["name"] = args.name; payload["modelSeeds"] = seed_values
    d = args.out / "in_s0"; d.mkdir(parents=True)
    jp = d / f"{args.name}.json"; jp.write_text(json.dumps(payload, indent=2) + "\n")
    source_rows = list(csv.DictReader(args.source_manifest.open(), delimiter="\t"))
    row = next(r for r in source_rows if r["name"] == old_name)
    row.update({"name": args.name, "shard": "0", "seeds": str(len(seed_values)),
                "expected_samples": str(len(seed_values) * 5), "input_json": str(jp),
                "description": row["description"] + f"; independent seeds {seed_values}"})
    with (args.out / "manifest.tsv").open("w", newline="") as h:
        w = csv.DictWriter(h, list(row), delimiter="\t"); w.writeheader(); w.writerow(row)
    print(f"prepared {args.name} with seeds {seed_values}")


if __name__ == "__main__": main()
