#!/usr/bin/env python3
"""Prepare independent one-seed-per-shard AF3 jobs from an existing manifest row."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-manifest", type=pathlib.Path, required=True)
    parser.add_argument("--source-name", required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--output-name", required=True)
    parser.add_argument("--seeds", default="2,3,4,5")
    parser.add_argument("--task-name", required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit(f"experiment exists: {args.out}")

    with args.source_manifest.open() as handle:
        source_rows = list(csv.DictReader(handle, delimiter="\t"))
    matches = [row for row in source_rows if row["name"] == args.source_name]
    if len(matches) != 1:
        raise SystemExit(
            f"expected one source row named {args.source_name}; found {len(matches)}"
        )
    source = matches[0]
    with pathlib.Path(source["input_json"]).open() as handle:
        payload = json.load(handle)
    seeds = [int(value) for value in args.seeds.split(",")]
    if not seeds or len(seeds) > 8 or len(set(seeds)) != len(seeds):
        raise SystemExit("--seeds must contain 1..8 unique integers")

    args.out.mkdir(parents=True)
    rows: list[dict[str, object]] = []
    for shard, seed in enumerate(seeds):
        input_dir = args.out / f"in_s{shard}"
        input_dir.mkdir()
        shard_payload = dict(payload)
        shard_payload["name"] = args.output_name
        shard_payload["modelSeeds"] = [seed]
        json_path = input_dir / f"{args.output_name}.json"
        with json_path.open("w") as handle:
            json.dump(shard_payload, handle, indent=2)
            handle.write("\n")
        row = dict(source)
        row.update({
            "name": args.output_name,
            "shard": shard,
            "seeds": len(seeds),
            "expected_samples": len(seeds) * 5,
            "description": f"{source['description']}; independent seed {seed}",
            "input_json": str(json_path),
        })
        rows.append(row)

    with (args.out / "manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with (args.out / "jobs.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["task_name", "mode", "status", "note"])
        writer.writerow([
            args.task_name,
            "local_af3_multiseed_failfast",
            "QUEUED",
            f"{args.source_name} independent seeds {','.join(map(str, seeds))}",
        ])
    print(f"prepared {args.output_name}: seeds={seeds} shards={len(seeds)}")


if __name__ == "__main__":
    main()
