#!/usr/bin/env python3
"""Prepare seeds 2--5 for the top strict seed-1 C4 hits.

One AF3 seed is placed in each shard.  Four consecutive shard numbers are
reserved per biological candidate, which makes capacity-limited/fail-fast
submission straightforward while keeping every Volc job to one GPU.
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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
    raise ValueError(f"FASTA record {name!r} not found in {path}")


def positions(path: pathlib.Path) -> set[int]:
    return {
        int(value)
        for value in path.read_text().replace(",", " ").split()
        if value.strip()
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=pathlib.Path, required=True)
    parser.add_argument("--strict-summary", type=pathlib.Path)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--surface-positions", type=pathlib.Path, required=True)
    parser.add_argument("--top", type=int, default=4)
    parser.add_argument(
        "--skip", type=int, default=0,
        help="skip this many higher-ranked strict seed1 hits before selecting --top",
    )
    args = parser.parse_args()

    strict_path = args.strict_summary or (
        args.source / "seed1_strict_final" / "seed1_strict_summary.tsv"
    )
    source_manifest_path = args.source / "manifest.tsv"
    strict = [
        row for row in read_tsv(strict_path)
        if int(row["seed1_strict_pass"]) == 1
    ]
    strict.sort(key=lambda row: (
        -float(row["mean_iptm"]), -float(row["min_iptm"]), row["system"]
    ))
    if args.skip < 0:
        raise ValueError("--skip must be non-negative")
    strict = strict[args.skip: args.skip + args.top]

    args.out.mkdir(parents=True, exist_ok=True)
    if not strict:
        (args.out / "NO_SEED1_HITS").touch()
        with (args.out / "selected_candidates.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow([
                "candidate_rank", "source_system", "system", "base_shard",
                "source_seed1_rank",
                "n_mutations", "surface_mutation_fraction", "seed1_mean_iptm",
                "seed1_min_iptm", "sequence",
            ])
        print("strict seed-1 hits=0; no seeds2-5 inputs prepared")
        return

    source_rows = read_tsv(source_manifest_path)
    by_name = {row["name"]: row for row in source_rows}
    reference = fasta_sequence(args.reference_fasta, args.reference_name)
    surface = positions(args.surface_positions)
    output_rows: list[dict[str, object]] = []
    selected_rows: list[dict[str, object]] = []
    manifest_fields = list(source_rows[0])
    extra_fields = [
        "source_system", "source_experiment", "independent_seed",
        "n_mutations", "surface_mutation_fraction",
    ]
    for rank, score in enumerate(strict, 1):
        source_name = score["system"]
        if source_name not in by_name:
            raise ValueError(f"strict hit absent from source manifest: {source_name}")
        source = by_name[source_name]
        sequence = source["sequence"]
        if len(sequence) != len(reference):
            raise ValueError(
                f"{source_name}: sequence length {len(sequence)} != {len(reference)}"
            )
        mutation_positions = {
            index for index, (old, new) in enumerate(zip(reference, sequence), 1)
            if old != new
        }
        n_mutations = len(mutation_positions)
        surface_fraction = (
            len(mutation_positions & surface) / n_mutations if n_mutations else 1.0
        )
        if n_mutations >= 25 or surface_fraction < 0.80:
            raise ValueError(
                f"{source_name}: strict hit violates sequence gates: "
                f"nmut={n_mutations} surface={surface_fraction:.3f}"
            )
        payload = json.loads(pathlib.Path(source["input_json"]).read_text())
        output_name = f"{source_name}_ms25_r{rank:02d}"
        base_shard = (rank - 1) * 4
        selected_rows.append({
            "candidate_rank": rank,
            "source_system": source_name,
            "system": output_name,
            "base_shard": base_shard,
            "source_seed1_rank": args.skip + rank,
            "n_mutations": n_mutations,
            "surface_mutation_fraction": f"{surface_fraction:.6f}",
            "seed1_mean_iptm": score["mean_iptm"],
            "seed1_min_iptm": score["min_iptm"],
            "sequence": sequence,
        })
        for offset, seed in enumerate(range(2, 6)):
            shard = base_shard + offset
            input_dir = args.out / f"in_s{shard}"
            input_dir.mkdir(exist_ok=True)
            current = json.loads(json.dumps(payload))
            current["name"] = output_name
            current["modelSeeds"] = [seed]
            input_json = input_dir / f"{output_name}.json"
            input_json.write_text(json.dumps(current, indent=2) + "\n")
            row: dict[str, object] = dict(source)
            row.update({
                "name": output_name,
                "shard": shard,
                "seeds": 1,
                "expected_samples": 5,
                "description": (
                    f"{source['description']}; independent AF3 seed {seed}; "
                    f"source_system={source_name}"
                ),
                "input_json": str(input_json),
                "source_system": source_name,
                "source_experiment": str(args.source.resolve()),
                "independent_seed": seed,
                "n_mutations": n_mutations,
                "surface_mutation_fraction": f"{surface_fraction:.6f}",
            })
            output_rows.append(row)

    with (args.out / "manifest.tsv").open("w", newline="") as handle:
        fields = manifest_fields + [field for field in extra_fields if field not in manifest_fields]
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    with (args.out / "selected_candidates.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, list(selected_rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(selected_rows)
    (args.out / "CANDIDATES_FROZEN").touch()
    print(
        f"prepared strict hits={len(selected_rows)} skip={args.skip} "
        f"shards={len(output_rows)} "
        f"in {args.out}"
    )


if __name__ == "__main__":
    main()
