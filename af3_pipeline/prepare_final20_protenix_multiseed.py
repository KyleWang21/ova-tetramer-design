#!/usr/bin/env python3
"""Prepare four additional Protenix-v2 seeds for each final OVA candidate."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import shutil


def count_a3m_records(path: pathlib.Path) -> int:
    return sum(line.startswith(">") for line in path.read_text().splitlines())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--final", type=pathlib.Path, required=True)
    parser.add_argument("--input-dir", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--seeds", default="11,21,41,51")
    parser.add_argument("--shards", type=int, default=8)
    parser.add_argument(
        "--expected-candidates",
        type=int,
        default=20,
        help="fail unless this many candidates were materialized; use 0 to disable",
    )
    parser.add_argument(
        "--msa-path-prefix",
        help="absolute directory used inside generated JSONs (for a remote compute cluster)",
    )
    args = parser.parse_args()
    if (args.out / "volc_manifest.tsv").exists():
        raise SystemExit(f"experiment already prepared: {args.out}")
    candidates = [
        row["candidate"]
        for row in csv.DictReader(args.final.open(), delimiter="\t")
    ]
    inputs = sorted(args.input_dir.glob("*.json"))
    by_candidate = {}
    for candidate in candidates:
        token = candidate.lower().replace("-", "_")
        matches = [path for path in inputs if token in path.stem]
        if len(matches) != 1:
            raise ValueError(f"{candidate}: expected one input JSON, found {matches}")
        by_candidate[candidate] = str(matches[0].resolve())
    if args.expected_candidates and len(by_candidate) != args.expected_candidates:
        raise ValueError(
            f"expected {args.expected_candidates} candidates, found {len(by_candidate)}"
        )
    seeds = [int(item) for item in args.seeds.split(",")]
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "logs").mkdir(exist_ok=True)
    input_out = args.out / "input"
    input_out.mkdir()
    msa_out = args.out / "msa"
    msa_out.mkdir()
    materialized = {}
    msa_audit: dict[tuple[str, str], dict[str, object]] = {}
    for candidate, input_path in by_candidate.items():
        source = pathlib.Path(input_path)
        payload = json.loads(source.read_text())
        systems = payload if isinstance(payload, list) else [payload]
        for system in systems:
            for sequence in system["sequences"]:
                protein = sequence.get("proteinChain", sequence.get("protein"))
                if protein is None:
                    raise ValueError(f"{source}: unsupported protein-chain schema")
                msa_path = protein.pop("unpairedMsaPath", None)
                if msa_path:
                    msa_source = pathlib.Path(msa_path)
                    if not msa_source.exists():
                        msa_source = args.input_dir.parent / "msa" / msa_source.name
                    msa_target = msa_out / msa_source.name
                    shutil.copyfile(msa_source, msa_target)
                    msa_audit[(candidate, msa_target.name)] = {
                        "candidate": candidate,
                        "msa_file": str(msa_target.relative_to(args.out)),
                        "record_count": count_a3m_records(msa_target),
                    }
                    prefix = pathlib.Path(args.msa_path_prefix) if args.msa_path_prefix else msa_out.resolve()
                    protein["unpairedMsaPath"] = str(prefix / msa_target.name)
        target = input_out / source.name
        target.write_text(json.dumps(payload, indent=2) + "\n")
        materialized[candidate] = target
    rows = []
    job_index = 0
    for seed in seeds:
        for candidate in sorted(by_candidate):
            rows.append({
                "shard": job_index % args.shards,
                "candidate": candidate,
                "seed": seed,
                "input_json": str(materialized[candidate].relative_to(args.out)),
            })
            job_index += 1
    with (args.out / "volc_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    if msa_audit:
        audit_rows = list(msa_audit.values())
        with (args.out / "msa_manifest.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, list(audit_rows[0]), delimiter="\t", lineterminator="\n"
            )
            writer.writeheader(); writer.writerows(audit_rows)
    print(
        f"prepared {len(rows)} Protenix runs "
        f"({len(seeds)} seeds x {len(by_candidate)} candidates) across {args.shards} shards; "
        f"audited_msa_files={len(msa_audit)}"
    )


if __name__ == "__main__":
    main()
