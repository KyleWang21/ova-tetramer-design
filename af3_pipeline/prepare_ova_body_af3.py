#!/usr/bin/env python3
"""Prepare no-module, four-chain AF3 jobs for pure-OVA interface designs."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib

from prepare_inputs import parse_a3m_records, portable, protein


def read_fasta(path: pathlib.Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(sequence)))
            header, sequence = line[1:], []
        else:
            sequence.append(line.strip())
    if header is not None:
        records.append((header, "".join(sequence)))
    return records


def replace_query(base_a3m: str, sequence: str) -> str:
    records = parse_a3m_records(base_a3m)
    if len(records[0][1]) != len(sequence):
        raise ValueError(f"MSA query length {len(records[0][1])} != sequence length {len(sequence)}")
    output = [records[0][0], sequence]
    for header, aligned in records[1:]:
        output.extend((header, aligned))
    return "\n".join(output) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", type=pathlib.Path, required=True)
    ap.add_argument("--ova-msa", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--seed-start", type=int, default=1)
    ap.add_argument(
        "--seed-offset",
        type=int,
        default=0,
        help="add index*offset to each candidate's first model seed",
    )
    ap.add_argument("--shards", type=int, default=8)
    ap.add_argument("--n-chains", type=int, default=4, choices=range(1, 7))
    args = ap.parse_args()
    if (args.out / "manifest.tsv").exists():
        raise SystemExit(f"experiment exists: {args.out}")
    records = read_fasta(args.fasta)[:args.top]
    if not records:
        raise SystemExit("no FASTA records")
    if args.shards < 1 or args.shards > 8:
        raise SystemExit("shards must be 1..8")
    if any(len(sequence) != 386 for _, sequence in records):
        raise SystemExit("all pure-OVA candidates must be exactly 386 aa")
    # Queue/watchdog logs may intentionally pre-create the experiment directory
    # before candidate generation finishes.  The manifest guard above remains
    # the authoritative idempotency check.
    args.out.mkdir(parents=True, exist_ok=True)
    for shard in range(args.shards):
        (args.out / f"in_s{shard}").mkdir()
    msa_dir = args.out / "msas"; msa_dir.mkdir()
    base_a3m = args.ova_msa.read_text()
    rows: list[dict[str, object]] = []
    for index, (source_header, sequence) in enumerate(records):
        shard = index % args.shards
        name = f"ova_body_c{args.n_chains}_{index + 1:02d}"
        first_seed = args.seed_start + index * args.seed_offset
        model_seeds = list(range(first_seed, first_seed + args.seeds))
        msa_path = msa_dir / f"{name}.a3m"
        msa_path.write_text(replace_query(base_a3m, sequence))
        payload = {
            "name": name,
            "modelSeeds": model_seeds,
            "dialect": "alphafold3",
            "version": 2,
            "sequences": [
                protein(chain, sequence, portable(msa_path))
                for chain in "ABCDEF"[:args.n_chains]
            ],
        }
        json_path = args.out / f"in_s{shard}" / f"{name}.json"
        json_path.write_text(json.dumps(payload, indent=2) + "\n")
        rows.append({
            "name": name,
            "shard": shard,
            "kind": f"ova_body_c{args.n_chains}_candidate",
            "n_chains": args.n_chains,
            "chain_ids": "ABCDEF"[:args.n_chains],
            "chain_length": ",".join(["386"] * args.n_chains),
            "ova_ranges": ";".join(["1-386"] * args.n_chains),
            "module_ranges": ";".join(["0-0"] * args.n_chains),
            "architecture": f"pure OVA P3-13R C{args.n_chains} interface redesign; no linker/module",
            "seeds": args.seeds,
            "expected_samples": args.seeds * 5,
            "description": source_header,
            "input_json": str(json_path),
            "sequence": sequence,
        })
    with (args.out / "manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    with (args.out / "representative_chains.fasta").open("w") as handle:
        for row in rows:
            handle.write(f">{row['name']} | {row['description']}\n{row['sequence']}\n")
    print(f"prepared {len(rows)} pure-OVA C{args.n_chains} candidates across {args.shards} shards")


if __name__ == "__main__":
    main()
