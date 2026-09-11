#!/usr/bin/env python3
"""Seed AF3 interface targets with diverse active-learning sequences.

Multiple ``--target-manifest`` inputs are consumed round-robin.  This makes it
possible to expose one tied sequence-design wave to independently predicted C4
backbones from more than one accepted sequence instead of overfitting every
trajectory to a single structural ensemble.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import shutil


PROTECTED = {12, 31, 74, 121, 368, 383} | set(range(258, 266))


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
    raise ValueError(f"FASTA record {name!r} not found")


def positions(path: pathlib.Path) -> set[int]:
    return {
        int(value)
        for value in path.read_text().replace(",", " ").split()
        if value.strip()
    }


def hamming(left: str, right: str) -> int:
    return sum(a != b for a, b in zip(left, right))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--surrogate", type=pathlib.Path, required=True)
    parser.add_argument(
        "--target-manifest", type=pathlib.Path, action="append", required=True,
        help="target manifest; repeat to interleave independent backbone ensembles",
    )
    parser.add_argument("--surface-positions", type=pathlib.Path, required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--exclude-glob", action="append", default=[])
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--top", type=int, default=4)
    parser.add_argument("--minimum-distance", type=int, default=4)
    parser.add_argument("--ranking-column", default="acquisition_average_rank")
    parser.add_argument("--reported-rank-column", default="rank")
    args = parser.parse_args()

    if (args.out / "target_manifest.tsv").exists():
        raise SystemExit(f"experiment already prepared: {args.out}")
    reference = fasta_sequence(args.reference_fasta, args.reference_name)
    surface = positions(args.surface_positions)
    excluded: set[str] = set()
    for pattern in args.exclude_glob:
        for path in pathlib.Path.cwd().glob(pattern):
            if path.is_file() and path.suffix == ".tsv":
                excluded.update(
                    row.get("sequence", "") for row in read_tsv(path)
                    if len(row.get("sequence", "")) == len(reference)
                )

    candidates = read_tsv(args.surrogate)
    candidates.sort(key=lambda row: (
        float(row[args.ranking_column]),
        -float(row.get("predicted_af3", 0.0)),
        int(row["n_mutations"]),
    ))
    selected: list[dict[str, str]] = []
    for minimum in range(args.minimum_distance, 0, -1):
        for row in candidates:
            sequence = row["sequence"]
            if sequence in excluded or any(sequence == old["sequence"] for old in selected):
                continue
            if all(hamming(sequence, old["sequence"]) >= minimum for old in selected):
                selected.append(row)
            if len(selected) == args.top:
                break
        if len(selected) == args.top:
            break
    if len(selected) != args.top:
        raise ValueError(f"only {len(selected)} novel diverse surrogate seeds")

    target_groups = [read_tsv(path) for path in args.target_manifest]
    targets = [
        group[index]
        for index in range(max(map(len, target_groups)))
        for group in target_groups
        if index < len(group)
    ]
    if len(targets) < args.top:
        raise ValueError(f"need at least {args.top} target rows, found {len(targets)}")
    args.out.mkdir(parents=True)
    design_dir = args.out / "design_masks"
    design_dir.mkdir()
    shutil.copyfile(args.surface_positions, args.out / "surface_positions.txt")

    manifest: list[dict[str, object]] = []
    selected_out: list[dict[str, object]] = []
    with (args.out / "anchors.fasta").open("w") as fasta:
        for index, (candidate, target) in enumerate(zip(selected, targets), 1):
            sequence = candidate["sequence"]
            changed = {
                position for position, (old, new) in enumerate(zip(reference, sequence), 1)
                if old != new
            }
            if len(changed) >= 25:
                raise ValueError(f"seed {index}: mutation budget exceeded")
            surface_fraction = len(changed & surface) / len(changed)
            if surface_fraction < 0.80:
                raise ValueError(f"seed {index}: surface fraction {surface_fraction:.3f}")
            anchor = f"XM2-SEED-{index:02d}"
            fasta.write(f">{anchor}\n{sequence}\n")
            design = (positions(pathlib.Path(target["joint_design"])) | changed) - PROTECTED
            design_path = design_dir / f"{anchor.lower()}_joint_design.txt"
            design_path.write_text(",".join(map(str, sorted(design))) + "\n")
            record: dict[str, object] = dict(target)
            record.update({
                "anchor": anchor,
                "anchor_n_mutations": len(changed),
                "joint_design": str(design_path),
                "joint_design_positions": len(design),
                "surrogate_name": candidate["name"],
                "surrogate_rank": candidate[args.reported_rank_column],
                "surrogate_ranking_metric": args.ranking_column,
                "surrogate_ranking_value": candidate[args.ranking_column],
                "surrogate_predicted_af3": candidate.get("predicted_af3", ""),
                "surrogate_predicted_protenix": candidate.get("predicted_protenix", ""),
            })
            manifest.append(record)
            selected_out.append({
                "anchor": anchor,
                "target_anchor": target["anchor"],
                "surrogate_name": candidate["name"],
                "surrogate_rank": candidate[args.reported_rank_column],
                "surrogate_ranking_metric": args.ranking_column,
                "surrogate_ranking_value": candidate[args.ranking_column],
                "n_mutations": len(changed),
                "surface_mutation_fraction": f"{surface_fraction:.6f}",
                "nearest_selected_distance": min(
                    [hamming(sequence, old["sequence"]) for old in selected[:index - 1]]
                    or [len(sequence)]
                ),
                "predicted_af3": candidate.get("predicted_af3", ""),
                "predicted_protenix": candidate.get("predicted_protenix", ""),
                "predicted_opendde": candidate.get("predicted_opendde", ""),
                "sequence": sequence,
            })

    fields = list(manifest[0])
    with (args.out / "target_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest)
    with (args.out / "selected_surrogate_seeds.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(selected_out[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(selected_out)
    print(
        f"prepared {len(manifest)} active-learning seeds across four AF3 targets; "
        f"nmut={min(int(row['n_mutations']) for row in selected_out)}-"
        f"{max(int(row['n_mutations']) for row in selected_out)}"
    )


if __name__ == "__main__":
    main()
