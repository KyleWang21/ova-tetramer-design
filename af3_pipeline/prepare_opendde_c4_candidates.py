#!/usr/bin/env python3
"""Prepare generic 10-seed OpenDDE C4 inputs from an accepted-candidate table."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import re

from prepare_inputs import parse_a3m_records


def write_a3m(path: pathlib.Path, records: list[tuple[str, str]]) -> None:
    path.write_text("".join(f"{header}\n{sequence}\n" for header, sequence in records))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--final", type=pathlib.Path, required=True)
    parser.add_argument("--msa-dir", type=pathlib.Path, required=True)
    parser.add_argument("--runtime-msa-dir", required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--seeds", default="101,102,103,104,105,106,107,108,109,110")
    parser.add_argument("--shards", type=int, default=8)
    parser.add_argument("--max-msa-records", type=int, default=1024)
    args = parser.parse_args()

    manifest_path = args.out / "manifest.tsv"
    if manifest_path.exists():
        raise SystemExit(f"experiment already prepared: {args.out}")
    rows = list(csv.DictReader(args.final.open(), delimiter="\t"))
    if not rows:
        raise ValueError("candidate table is empty")
    seeds = [int(value) for value in args.seeds.split(",")]
    if len(seeds) != 10 or len(set(seeds)) != 10:
        raise ValueError("OpenDDE production screen requires 10 unique seeds")
    if args.max_msa_records < 1:
        raise ValueError("--max-msa-records must be positive")

    args.out.mkdir(parents=True)
    input_dir = args.out / "input"
    input_dir.mkdir()
    msa_dir = args.out / "msa"
    msa_dir.mkdir()
    manifest = []
    msa_manifest = []
    for index, row in enumerate(rows):
        candidate = row["candidate"]
        sequence = row["sequence"]
        token = re.sub(r"[^a-z0-9]+", "_", candidate.lower()).strip("_")
        matches = sorted(args.msa_dir.glob(f"*{token}*_1024.a3m"))
        if len(matches) != 1:
            raise ValueError(f"{candidate}: expected one MSA, found {matches}")
        source_records = parse_a3m_records(matches[0].read_text())
        query = source_records[0][1]
        if query != sequence:
            raise ValueError(f"{candidate}: MSA query differs from candidate sequence")
        capped_records = source_records[: args.max_msa_records]
        capped_name = f"{token}_opendde_{len(capped_records)}.a3m"
        capped_path = msa_dir / capped_name
        write_a3m(capped_path, capped_records)
        job = f"{token}_opendde_c4"
        runtime_msa = f"{args.runtime_msa_dir.rstrip('/')}/{capped_name}"
        payload = [{
            "name": job,
            "modelSeeds": seeds,
            "sequences": [{"proteinChain": {
                "sequence": sequence,
                "count": 4,
                "unpairedMsaPath": runtime_msa,
            }}],
        }]
        input_json = input_dir / f"{job}.json"
        input_json.write_text(json.dumps(payload, indent=2) + "\n")
        manifest.append({
            "shard": index % args.shards,
            "name": candidate,
            "kind": "accepted_af3_candidate",
            "n_seeds": len(seeds),
            "input_json": str(input_json.relative_to(args.out)),
        })
        msa_manifest.append({
            "candidate": candidate,
            "source_msa": str(matches[0]),
            "source_records": len(source_records),
            "kept_records": len(capped_records),
            "capped_msa": str(capped_path.relative_to(args.out)),
            "runtime_msa": runtime_msa,
        })

    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(manifest[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest)
    with (args.out / "msa_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, list(msa_manifest[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(msa_manifest)
    print(
        f"prepared {len(manifest)} OpenDDE C4 systems with {len(seeds)} seeds each; "
        f"MSA records capped at {args.max_msa_records}"
    )


if __name__ == "__main__":
    main()
