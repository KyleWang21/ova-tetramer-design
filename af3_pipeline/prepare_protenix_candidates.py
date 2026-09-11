#!/usr/bin/env python3
"""Prepare counted homotetramer Protenix inputs for final AF3 candidates."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import re


def replace_query(a3m: str, sequence: str) -> str:
    lines = a3m.replace("\x00", "").splitlines()
    for index, line in enumerate(lines):
        if line.startswith(">"):
            query_index = index + 1
            while query_index < len(lines) and not lines[query_index].strip():
                query_index += 1
            if query_index == len(lines) or lines[query_index].startswith(">"):
                raise ValueError("first A3M record has no sequence")
            lines[query_index] = sequence
            return "\n".join(lines) + "\n"
    raise ValueError("A3M contains no FASTA header")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--validation", type=pathlib.Path, required=True)
    parser.add_argument("--base-a3m", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--passing-only", action="store_true")
    parser.add_argument("--top", type=int, default=0, help="prepare only the first N retained rows")
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    msa_dir = args.out / "msa"
    input_dir = args.out / "input"
    msa_dir.mkdir(exist_ok=True)
    input_dir.mkdir(exist_ok=True)
    base_a3m = args.base_a3m.read_text()

    rows = list(csv.DictReader(args.validation.open(), delimiter="\t"))
    retained = 0
    for row in rows:
        if args.passing_only and row.get("passes_final") != "1":
            continue
        source = row.get("candidate") or row.get("seed1_system") or row.get("system") or "candidate"
        tag = re.sub(r"[^a-z0-9]+", "_", source.lower()).strip("_")
        name = f"ova_iptm80_{tag}_protenix"
        sequence = row["sequence"]
        msa_path = (msa_dir / f"{name}_1024.a3m").resolve()
        msa_path.write_text(replace_query(base_a3m, sequence))
        payload = [{
            "name": name,
            "sequences": [{"proteinChain": {
                "sequence": sequence,
                "count": 4,
                "unpairedMsaPath": str(msa_path),
            }}],
        }]
        input_path = input_dir / f"{name}.json"
        input_path.write_text(json.dumps(payload, indent=2) + "\n")
        print(input_path)
        retained += 1
        if args.top and retained >= args.top:
            break


if __name__ == "__main__":
    main()
