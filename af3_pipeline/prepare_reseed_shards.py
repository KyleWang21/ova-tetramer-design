#!/usr/bin/env python3
"""Split independent AF3 seeds for one system across one-GPU shards."""

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
    ap.add_argument("--seeds", required=True)
    args = ap.parse_args()
    if (args.out / "manifest.tsv").exists():
        raise SystemExit(f"experiment exists: {args.out}")
    seeds = [int(value) for value in args.seeds.split(",")]
    if not seeds or len(seeds) > 8:
        raise SystemExit("one to eight seeds are required")
    payload = json.loads(args.source_json.read_text())
    old_name = payload["name"]
    source_rows = list(csv.DictReader(args.source_manifest.open(), delimiter="\t"))
    source = next(row for row in source_rows if row["name"] == old_name)
    rows = []
    for shard, seed in enumerate(seeds):
        input_dir = args.out / f"in_s{shard}"; input_dir.mkdir(parents=True, exist_ok=True)
        current = dict(payload)
        current["name"] = args.name
        current["modelSeeds"] = [seed]
        json_path = input_dir / f"{args.name}.json"
        json_path.write_text(json.dumps(current, indent=2) + "\n")
        row = dict(source)
        row.update({
            "name": args.name,
            "shard": str(shard),
            "seeds": str(len(seeds)),
            "expected_samples": str(len(seeds) * 5),
            "input_json": str(json_path),
            "description": source["description"] + f"; independent seed {seed}",
        })
        rows.append(row)
    with (args.out / "manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    print(f"prepared {args.name}: {len(seeds)} one-seed shards, seeds={seeds}")


if __name__ == "__main__":
    main()
