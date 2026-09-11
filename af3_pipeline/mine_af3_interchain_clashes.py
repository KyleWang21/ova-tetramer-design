#!/usr/bin/env python3
"""Enumerate exact <2.4-A interchain heavy-atom contacts in AF3 C4 models."""

from __future__ import annotations

import argparse
import csv
import pathlib
import statistics
from collections import defaultdict

import numpy as np
from scipy.spatial import cKDTree

from screen_final20_af3_structures import structure_arrays


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: pathlib.Path, name_prefix: str) -> str:
    records: list[tuple[str, str]] = []
    name = ""
    sequence: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    records.append((name, "".join(sequence)))
                name = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)
    if name:
        records.append((name, "".join(sequence)))
    matches = [seq for record_name, seq in records if record_name.startswith(name_prefix)]
    if len(matches) != 1:
        raise ValueError(
            f"expected one FASTA record starting with {name_prefix!r}; found {len(matches)}"
        )
    return matches[0]


def mutation_map_from_manifest(
    manifest_path: pathlib.Path,
    reference_path: pathlib.Path,
    reference_name: str,
) -> dict[str, dict[int, str]]:
    reference = read_fasta(reference_path, reference_name)
    result: dict[str, dict[int, str]] = {}
    for row in read_tsv(manifest_path):
        sequence = row["sequence"]
        if len(sequence) != len(reference):
            raise ValueError(
                f"{row['name']}: sequence length {len(sequence)} != reference {len(reference)}"
            )
        result[row["name"]] = {
            position: f"{old}{position}{new}"
            for position, (old, new) in enumerate(zip(reference, sequence), start=1)
            if old != new
        }
    return result


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"empty table: {path}")
    fields = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--per-model", type=pathlib.Path, required=True)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--final", type=pathlib.Path)
    source.add_argument("--manifest", type=pathlib.Path)
    parser.add_argument("--reference-fasta", type=pathlib.Path)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--cutoff", type=float, default=2.4)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    if args.manifest:
        if not args.reference_fasta:
            parser.error("--reference-fasta is required with --manifest")
        mutation = mutation_map_from_manifest(
            args.manifest, args.reference_fasta, args.reference_name
        )
    else:
        final = {row["candidate"]: row for row in read_tsv(args.final)}
        mutation = {
            name: {
                int(item[1:-1]): item
                for item in row["mutations"].split(",")
                if item
            }
            for name, row in final.items()
        }
    events = []
    for model in read_tsv(args.per_model):
        if int(model["atomic_clash_count_2p4"]) == 0:
            continue
        candidate = model.get("candidate") or model.get("system")
        if not candidate:
            raise ValueError("per-model table requires a candidate or system column")
        if candidate not in mutation:
            raise KeyError(f"missing sequence/mutation metadata for {candidate}")
        _, chains = structure_arrays(pathlib.Path(model["cif_path"]))
        chain_ids = sorted(chains)
        for left_index, chain_a in enumerate(chain_ids):
            a = chains[chain_a]
            for chain_b in chain_ids[left_index + 1:]:
                b = chains[chain_b]
                hits = cKDTree(a["coords"]).query_ball_tree(cKDTree(b["coords"]), args.cutoff)
                for atom_a, partners in enumerate(hits):
                    for atom_b in partners:
                        distance = float(np.linalg.norm(a["coords"][atom_a] - b["coords"][atom_b]))
                        pos_a, pos_b = int(a["resids"][atom_a]), int(b["resids"][atom_b])
                        events.append({
                            "candidate": candidate,
                            "seed": model.get("seed", ""),
                            "sample": model.get("sample", ""),
                            "chain_pair": f"{chain_a}-{chain_b}",
                            "chain_a": chain_a,
                            "position_a": pos_a,
                            "resname_a": str(a["resnames"][atom_a]),
                            "atom_a": str(a["names"][atom_a]),
                            "mutation_a": mutation[candidate].get(pos_a, ""),
                            "chain_b": chain_b,
                            "position_b": pos_b,
                            "resname_b": str(b["resnames"][atom_b]),
                            "atom_b": str(b["names"][atom_b]),
                            "mutation_b": mutation[candidate].get(pos_b, ""),
                            "distance_a": distance,
                            "cif_path": model["cif_path"],
                        })

    write_tsv(args.out / "clash_events.tsv", events)
    grouped: dict[tuple, list[dict[str, object]]] = defaultdict(list)
    for row in events:
        left = (int(row["position_a"]), str(row["resname_a"]), str(row["atom_a"]))
        right = (int(row["position_b"]), str(row["resname_b"]), str(row["atom_b"]))
        key = tuple(sorted((left, right)))
        grouped[key].append(row)
    summary = []
    for (left, right), rows in grouped.items():
        distances = [float(row["distance_a"]) for row in rows]
        summary.append({
            "position_1": left[0], "resname_1": left[1], "atom_1": left[2],
            "position_2": right[0], "resname_2": right[1], "atom_2": right[2],
            "candidate_count": len({str(row["candidate"]) for row in rows}),
            "model_count": len({str(row["cif_path"]) for row in rows}),
            "event_count": len(rows),
            "min_distance_a": min(distances),
            "median_distance_a": statistics.median(distances),
            "candidates": ",".join(sorted({str(row["candidate"]) for row in rows})),
        })
    summary.sort(key=lambda row: (-int(row["model_count"]), -int(row["candidate_count"]), row["position_1"]))
    write_tsv(args.out / "clash_position_atom_summary.tsv", summary)
    print(
        f"enumerated {len(events)} short contacts across "
        f"{len({row['cif_path'] for row in events})} models; patterns={len(summary)}"
    )


if __name__ == "__main__":
    main()
