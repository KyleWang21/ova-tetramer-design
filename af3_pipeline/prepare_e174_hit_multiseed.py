#!/usr/bin/env python3
"""Prepare independent AF3 seeds 2-5 for one strict E174 seed-1 hit."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import re


FIELDS = [
    "name", "shard", "kind", "n_chains", "chain_ids", "chain_length",
    "ova_ranges", "module_ranges", "architecture", "seeds", "expected_samples",
    "description", "input_json", "sequence",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=pathlib.Path, required=True)
    parser.add_argument("--source-system", required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--task-name", required=True)
    parser.add_argument("--note", required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"ky-\d{8}-\d{3}", args.task_name):
        raise SystemExit(f"invalid task name: {args.task_name}")
    rows = list(csv.DictReader((args.source / "manifest.tsv").open(), delimiter="\t"))
    matches = [row for row in rows if row["name"] == args.source_system]
    if len(matches) != 1:
        raise SystemExit(f"expected one source row for {args.source_system}, found {len(matches)}")
    source_row = matches[0]
    source_json = json.loads(pathlib.Path(source_row["input_json"]).read_text())
    sequence = source_row["sequence"]
    for entity in source_json["sequences"]:
        if entity["protein"]["sequence"] != sequence:
            raise SystemExit("source JSON sequence disagrees with manifest")

    args.out.mkdir(parents=True, exist_ok=True)
    output_rows = []
    for shard, seed in enumerate(range(2, 6)):
        input_dir = args.out / f"in_s{shard}"
        input_dir.mkdir(exist_ok=True)
        payload = json.loads(json.dumps(source_json))
        payload["name"] = args.name
        payload["modelSeeds"] = [seed]
        input_json = input_dir / f"{args.name}.json"
        input_json.write_text(json.dumps(payload, indent=2) + "\n")
        output_rows.append({
            **{field: source_row[field] for field in FIELDS if field in source_row},
            "name": args.name,
            "shard": shard,
            "seeds": 4,
            "expected_samples": 20,
            "description": f"{source_row['description']}; independent seed {seed}",
            "input_json": str(input_json),
            "sequence": sequence,
        })

    with (args.out / "manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    with (args.out / "jobs.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["task_name", "mode", "status", "note"])
        writer.writerow([args.task_name, "local_af3_multiseed", "QUEUED", args.note])
    print(f"prepared {args.name}: seeds2-5 in {args.out}")


if __name__ == "__main__":
    main()
