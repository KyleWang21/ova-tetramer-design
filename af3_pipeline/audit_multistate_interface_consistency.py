#!/usr/bin/env python3
"""Fail closed unless each interface class is contact-map consistent across backbones."""

from __future__ import annotations

import argparse
import csv
import itertools
import pathlib

import numpy as np


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def atoms(path: pathlib.Path) -> dict[str, list[tuple[int, np.ndarray]]]:
    chains: dict[str, list[tuple[int, np.ndarray]]] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        element = line[76:78].strip() or line[12:16].strip()[0]
        if element == "H":
            continue
        chain = line[21]
        residue = int(line[22:26])
        xyz = np.asarray(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])],
            dtype=np.float32,
        )
        chains.setdefault(chain, []).append((residue, xyz))
    return chains


def contact_map(path: pathlib.Path, pair: str, cutoff: float) -> set[tuple[int, int]]:
    chains = atoms(path)
    left, right = chains[pair[0]], chains[pair[1]]
    left_res = np.asarray([item[0] for item in left])
    right_res = np.asarray([item[0] for item in right])
    left_xyz = np.stack([item[1] for item in left])
    right_xyz = np.stack([item[1] for item in right])
    distance2 = np.square(left_xyz[:, None, :] - right_xyz[None, :, :]).sum(axis=2)
    li, ri = np.where(distance2 <= cutoff * cutoff)
    return {(int(left_res[i]), int(right_res[j])) for i, j in zip(li, ri, strict=True)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--state-manifest", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--cutoff", type=float, default=5.3)
    parser.add_argument("--minimum-jaccard", type=float, default=0.70)
    args = parser.parse_args()
    project = pathlib.Path(__file__).resolve().parents[1]
    states = read_tsv(args.state_manifest)
    groups: dict[int, list[tuple[dict[str, str], set[tuple[int, int]]]]] = {1: [], 2: []}
    for row in states:
        state_index = int(row["state_index"])
        class_index = 1 if state_index % 2 else 2
        path = pathlib.Path(row["target_pdb"])
        if not path.is_absolute():
            path = project / path
        groups[class_index].append((row, contact_map(path, row["interface"], args.cutoff)))
    audit: list[dict[str, object]] = []
    for class_index, current in groups.items():
        if len(current) != 4:
            raise RuntimeError(f"interface class {class_index} has {len(current)} states, expected 4")
        for (left, left_map), (right, right_map) in itertools.combinations(current, 2):
            union = left_map | right_map
            score = len(left_map & right_map) / len(union) if union else 0.0
            audit.append({
                "interface_class": class_index,
                "left_context": left["context_index"],
                "right_context": right["context_index"],
                "left_pair": left["interface"],
                "right_pair": right["interface"],
                "left_contacts": len(left_map),
                "right_contacts": len(right_map),
                "jaccard": score,
                "pass": int(score >= args.minimum_jaccard),
            })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(audit[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(audit)
    minimum = min(float(row["jaccard"]) for row in audit)
    print(f"comparisons={len(audit)} minimum_jaccard={minimum:.3f}")
    if any(not int(row["pass"]) for row in audit):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
