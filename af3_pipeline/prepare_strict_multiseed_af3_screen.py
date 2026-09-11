#!/usr/bin/env python3
"""Build the frozen 25-model AF3 screening inputs for reseeded E158 hits."""

from __future__ import annotations

import argparse
import csv
import pathlib
import re


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing empty output: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fasta_sequence(path: pathlib.Path, name: str) -> str:
    active = False
    chunks: list[str] = []
    for line in path.read_text().splitlines() + [">"]:
        if line.startswith(">"):
            if active:
                return "".join(chunks)
            active = line[1:].split()[0] == name
        elif active:
            chunks.append(line.strip())
    raise ValueError(f"FASTA record {name!r} not found")


def mutations(reference: str, sequence: str) -> str:
    return ",".join(
        f"{old}{position}{new}"
        for position, (old, new) in enumerate(zip(reference, sequence), 1)
        if old != new
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed1-experiment", type=pathlib.Path, required=True)
    parser.add_argument(
        "--multiseed-experiment", type=pathlib.Path, action="append", required=True
    )
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    reference = fasta_sequence(args.reference_fasta, args.reference_name)
    seed1_manifest = read_tsv(args.seed1_experiment / "manifest.tsv")
    seed1_scores = read_tsv(args.seed1_experiment / "model_scores.tsv")
    seed1_by_sequence = {row["sequence"]: row for row in seed1_manifest}
    final_rows: list[dict[str, object]] = []
    model_rows: list[dict[str, object]] = []

    for rank, experiment in enumerate(args.multiseed_experiment, 1):
        manifest = read_tsv(experiment / "manifest.tsv")
        unique = {(row["name"], row["sequence"]) for row in manifest}
        if len(unique) != 1:
            raise ValueError(f"{experiment}: expected one reseeded sequence, found {len(unique)}")
        candidate, sequence = next(iter(unique))
        source = seed1_by_sequence.get(sequence)
        if source is None:
            raise ValueError(f"{candidate}: sequence absent from seed-1 manifest")
        combined = [row for row in seed1_scores if row["system"] == source["name"]]
        combined += [
            row for row in read_tsv(experiment / "model_scores.tsv")
            if row["system"] == candidate
        ]
        seeds: dict[int, int] = {}
        for row in combined:
            match = re.fullmatch(r"seed-(\d+)_sample-(\d+)", row["model"])
            if match is None:
                raise ValueError(f"{candidate}: unexpected model name {row['model']!r}")
            seed = int(match.group(1))
            seeds[seed] = seeds.get(seed, 0) + 1
            model_rows.append({
                "candidate": candidate,
                "seed": seed,
                "model": row["model"],
                "iptm": row["iptm"],
                "min_incident_iptm": row["min_incident_iptm"],
                "contact_pair_pae_median": row["contact_pair_pae_median"],
                "cif_path": row["cif_path"],
            })
        if len(combined) != 25 or seeds != {1: 5, 2: 5, 3: 5, 4: 5, 5: 5}:
            raise ValueError(f"{candidate}: expected 25 models as five per seed, found {seeds}")
        final_rows.append({
            "candidate": candidate,
            "rank": rank,
            "mutations": mutations(reference, sequence),
            "sequence": sequence,
        })

    args.out.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out / "final_candidates.tsv", final_rows)
    write_tsv(args.out / "per_model_input.tsv", model_rows)
    print(f"prepared {len(final_rows)} candidates and {len(model_rows)} models in {args.out}")


if __name__ == "__main__":
    main()
