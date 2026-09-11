#!/usr/bin/env python3
"""Summarize optional OpenDDE ProLIF chemistry and AF3 interaction consensus."""

from __future__ import annotations

import argparse
import csv
import pathlib
from collections import defaultdict


NONVDW = {"Hydrophobic", "HydrogenBond", "SaltBridge", "PiStacking", "CationPi"}


def read_tsv(path: pathlib.Path | None) -> list[dict[str, str]]:
    if path is None or not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def normalized_key(row: dict[str, str]) -> tuple[int, int, str]:
    positions = sorted((int(row["position_a"]), int(row["position_b"])))
    return positions[0], positions[1], row["category"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--final", type=pathlib.Path, required=True)
    parser.add_argument("--opendde-summary", type=pathlib.Path, required=True)
    parser.add_argument("--opendde-per-seed", type=pathlib.Path, required=True)
    parser.add_argument("--opendde-unique", type=pathlib.Path)
    parser.add_argument("--crossmodel-unique", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    final_rows = read_tsv(args.final)
    final = {row["candidate"]: row for row in final_rows}
    opendde = {row["candidate"]: row for row in read_tsv(args.opendde_summary)}
    selected_seed = {
        name: int(row["selected_seed"])
        for name, row in opendde.items()
        if name in final and row.get("opendde_evidence_available") == "1" and row.get("selected_seed")
    }
    selected_geometry = {}
    for row in read_tsv(args.opendde_per_seed):
        name = row["candidate"]
        if name in selected_seed and int(row["seed"]) == selected_seed[name]:
            selected_geometry[name] = row

    open_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(args.opendde_unique):
        if row.get("candidate") in final and row.get("source", "").endswith("__OpenDDE"):
            open_rows[row["candidate"]].append(row)
    af3_keys: dict[str, set[tuple[int, int, str]]] = defaultdict(set)
    for row in read_tsv(args.crossmodel_unique):
        if (row.get("candidate") in final and row.get("source", "").endswith("__AF3")
                and row.get("category") in NONVDW):
            af3_keys[row["candidate"]].add(normalized_key(row))

    output = []
    for name, base in final.items():
        mutations = base["mutations"].split(",")
        rows = open_rows[name]
        any_mutations: set[str] = set()
        nonvdw_mutations: set[str] = set()
        chain_support: dict[str, set[str]] = defaultdict(set)
        strength: dict[str, int] = defaultdict(int)
        for row in rows:
            for suffix in ("a", "b"):
                mutation = row.get(f"mutation_{suffix}", "")
                if not mutation:
                    continue
                any_mutations.add(mutation)
                if row["category"] in NONVDW:
                    nonvdw_mutations.add(mutation)
                    chain_support[mutation].add(row[f"chain_{suffix}"])
                    strength[mutation] += 1

        interface_groups = []
        geometry = selected_geometry.get(name, {})
        for group in geometry.get("interface_type_edges", "").split("|"):
            edges = {edge for edge in group.split(",") if edge}
            if edges:
                interface_groups.append(edges)
        type_mutations = []
        for edges in interface_groups[:2]:
            found = set()
            for row in rows:
                if row["chain_pair"] not in edges or row["category"] not in NONVDW:
                    continue
                found.update(filter(None, (row.get("mutation_a", ""), row.get("mutation_b", ""))))
            type_mutations.append(found)

        major = sorted(nonvdw_mutations, key=lambda item: (-strength[item], item))[:4]
        open_keys = {normalized_key(row) for row in rows if row["category"] in NONVDW}
        shared = open_keys & af3_keys[name]
        union = open_keys | af3_keys[name]
        evidence = name in selected_seed
        gates = {
            "gate_opendde_prolif_any_fraction_ge0p50": evidence and len(any_mutations) / len(mutations) >= 0.50,
            "gate_opendde_prolif_nonvdw_fraction_ge0p25": evidence and len(nonvdw_mutations) / len(mutations) >= 0.25,
            "gate_opendde_prolif_nonvdw_count_ge6": evidence and len(nonvdw_mutations) >= 6,
            "gate_opendde_prolif_two_interfaces_nonvdw_ge2": (
                evidence and len(type_mutations) >= 2 and min(map(len, type_mutations)) >= 2
            ),
            "gate_opendde_prolif_major_hotspots_symmetry_ge3of4": (
                evidence and bool(major) and all(len(chain_support[item]) >= 3 for item in major)
            ),
            "gate_opendde_prolif_shared_af3_nonvdw_ge1": evidence and len(shared) >= 1,
        }
        output.append({
            "candidate": name,
            "opendde_prolif_evidence_available": int(evidence),
            "opendde_prolif_design_any_count": len(any_mutations),
            "opendde_prolif_design_any_fraction": len(any_mutations) / len(mutations),
            "opendde_prolif_design_nonvdw_count": len(nonvdw_mutations),
            "opendde_prolif_design_nonvdw_fraction": len(nonvdw_mutations) / len(mutations),
            "opendde_prolif_interface_type_nonvdw_counts": ",".join(map(str, map(len, type_mutations))),
            "opendde_prolif_major_hotspots": ",".join(major),
            "opendde_prolif_major_hotspot_chain_support": ",".join(
                f"{item}:{len(chain_support[item])}" for item in major
            ),
            "opendde_nonvdw_position_type_count": len(open_keys),
            "af3_opendde_shared_nonvdw_position_type_count": len(shared),
            "af3_opendde_nonvdw_position_type_jaccard": len(shared) / len(union) if union else 0.0,
            "af3_opendde_shared_nonvdw_position_types": ";".join(
                f"{left}-{right}:{category}" for left, right, category in sorted(shared)
            ),
            **{key: int(value) for key, value in gates.items()},
            "opendde_prolif_all_optional_gates_pass": int(all(gates.values())),
        })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = list(output[0])
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(output)
    print(
        f"OpenDDE ProLIF evidence={sum(row['opendde_prolif_evidence_available'] for row in output)}/20; "
        f"all_optional_gates={sum(row['opendde_prolif_all_optional_gates_pass'] for row in output)}/20"
    )


if __name__ == "__main__":
    main()
