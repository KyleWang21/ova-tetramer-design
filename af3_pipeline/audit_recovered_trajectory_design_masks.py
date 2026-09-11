#!/usr/bin/env python3
"""Audit E257 recovered candidates against their source-specific AC/AD masks."""

from __future__ import annotations

import argparse
import csv
import pathlib
import re


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
    raise ValueError(f"missing FASTA record {name}")


def positions(path: pathlib.Path) -> set[int]:
    return {int(value) for value in re.findall(r"[0-9]+", path.read_text())}


def source_experiment(row: dict[str, str]) -> pathlib.Path:
    match = re.search(r"/experiments/([^/]+)/", row.get("recovered_library", ""))
    if not match:
        raise ValueError(f"cannot parse recovered source: {row.get('recovered_library')}")
    return pathlib.Path("experiments") / match.group(1)


def task_index(source: str) -> int:
    match = re.search(r"/run_[^/]*_([0-9]{2})/", "/" + source)
    if not match:
        raise ValueError(f"cannot parse task index from {source}")
    return int(match.group(1))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shortlist", type=pathlib.Path, required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    reference = fasta_sequence(args.reference_fasta, args.reference_name)
    manifests: dict[pathlib.Path, list[dict[str, str]]] = {}
    output = []
    for row in read_tsv(args.shortlist):
        experiment = source_experiment(row)
        manifests.setdefault(experiment, read_tsv(experiment / "target_manifest.tsv"))
        index = task_index(row["source"])
        anchor_index = index // 2
        anchor = manifests[experiment][anchor_index]
        mask1 = positions(pathlib.Path(anchor["design1"]))
        mask2 = positions(pathlib.Path(anchor["design2"]))
        changed = {
            position for position, (old, new) in enumerate(zip(reference, row["sequence"]), 1)
            if old != new
        }
        in1, in2 = changed & mask1, changed & mask2
        output.append({
            "candidate": row["name"],
            "source_experiment": str(experiment),
            "source_trajectory": row["source"],
            "source_anchor": anchor["anchor"],
            "n_mutations": len(changed),
            "AC_designed_mutation_count": len(in1),
            "AC_designed_mutation_positions": ",".join(map(str, sorted(in1))),
            "AD_designed_mutation_count": len(in2),
            "AD_designed_mutation_positions": ",".join(map(str, sorted(in2))),
            "both_design_masks_ge3": int(len(in1) >= 3 and len(in2) >= 3),
            "all_mutations_in_joint_mask": int(changed <= (mask1 | mask2)),
            "sequence": row["sequence"],
        })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(output[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(output)
    print(
        f"candidates={len(output)} both_masks_ge3="
        f"{sum(int(row['both_design_masks_ge3']) for row in output)} "
        f"AC_min={min(int(row['AC_designed_mutation_count']) for row in output)} "
        f"AD_min={min(int(row['AD_designed_mutation_count']) for row in output)} "
        f"all_in_joint={sum(int(row['all_mutations_in_joint_mask']) for row in output)}"
    )


if __name__ == "__main__":
    main()
