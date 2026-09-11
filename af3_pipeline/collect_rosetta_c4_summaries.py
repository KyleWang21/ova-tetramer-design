#!/usr/bin/env python3
"""Collect per-candidate Rosetta summaries, preserving missing jobs fail-closed."""

import argparse
import csv
import pathlib


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--root", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()
    manifest = list(csv.DictReader(args.manifest.open(), delimiter="\t"))
    rows = []
    for item in manifest:
        path = args.root / item["candidate"] / f"{item['candidate']}_{item['source_model']}_rosetta_summary.tsv"
        if path.exists():
            rows.append(next(csv.DictReader(path.open(), delimiter="\t")))
        else:
            rows.append({
                "candidate": item["candidate"], "source_model": item["source_model"],
                "nrelax": 0, "production_ensemble_complete": 0, "rosetta_pass": 0,
                "failed_rosetta_gates": "MISSING",
            })
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows({field: row.get(field, "") for field in fields} for row in rows)
    print(f"collected {len(rows)} summaries; completed={sum(int(row['nrelax']) > 0 for row in rows)}")


if __name__ == "__main__":
    main()
