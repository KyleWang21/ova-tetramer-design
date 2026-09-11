#!/usr/bin/env python3
"""Recompute fail-closed AF3 geometry and topology metrics for the final OVA C4 set.

Run with the AlphaFold3 virtualenv.  All cutoffs match screening-plan v1.1:
2.4 A inter-chain heavy-atom clash, 2.5 A inter-chain CA overlap, and
5.3 A heavy-atom residue contact. Homomer topology comparison explicitly
optimizes over all 24 chain-label permutations.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import pathlib
import statistics
from collections import defaultdict

import numpy as np
from alphafold3.structure import from_mmcif
from scipy.spatial import cKDTree


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def jaccard(left: set, right: set) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 1.0


def components(nodes: list[str], edges: list[tuple[str, str]]) -> list[set[str]]:
    adjacent = {node: set() for node in nodes}
    for left, right in edges:
        adjacent[left].add(right)
        adjacent[right].add(left)
    output = []
    unseen = set(nodes)
    while unseen:
        start = min(unseen)
        group = {start}
        stack = [start]
        unseen.remove(start)
        while stack:
            node = stack.pop()
            for other in adjacent[node]:
                if other in unseen:
                    unseen.remove(other)
                    group.add(other)
                    stack.append(other)
        output.append(group)
    return output


def graph_clusters(similarity: np.ndarray, cutoff: float) -> list[list[int]]:
    n = similarity.shape[0]
    adjacent = {idx: set() for idx in range(n)}
    for left in range(n):
        for right in range(left + 1, n):
            if similarity[left, right] >= cutoff:
                adjacent[left].add(right)
                adjacent[right].add(left)
    output = []
    unseen = set(range(n))
    while unseen:
        start = min(unseen)
        unseen.remove(start)
        group = [start]
        stack = [start]
        while stack:
            node = stack.pop()
            for other in sorted(adjacent[node]):
                if other in unseen:
                    unseen.remove(other)
                    group.append(other)
                    stack.append(other)
        output.append(sorted(group))
    return sorted(output, key=lambda group: (-len(group), group))


def kabsch_rmsd(reference: np.ndarray, query: np.ndarray) -> float:
    reference = reference - reference.mean(axis=0)
    query = query - query.mean(axis=0)
    covariance = query.T @ reference
    left, _, right = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(left @ right))
    rotation = left @ correction @ right
    aligned = query @ rotation
    return float(np.sqrt(np.mean(np.sum((aligned - reference) ** 2, axis=1))))


def structure_arrays(cif: pathlib.Path) -> tuple[object, dict[str, dict[str, object]]]:
    structure = from_mmcif(cif.read_text())
    chains = {}
    coords = np.asarray(structure.coords, dtype=np.float32)
    for chain in structure.chains:
        mask = np.asarray(structure.chain_id == chain)
        elements = np.asarray(structure.atom_element[mask]).astype(str)
        heavy = np.asarray([value.upper() not in {"H", "D"} for value in elements])
        atom_names = np.asarray(structure.atom_name[mask]).astype(str)
        chains[str(chain)] = {
            "coords": coords[mask][heavy],
            "resids": np.asarray(structure.res_id[mask], dtype=int)[heavy],
            "names": atom_names[heavy],
            "resnames": np.asarray(structure.res_name[mask]).astype(str)[heavy],
            "plddt": np.asarray(structure.atom_b_factor[mask], dtype=float)[heavy],
        }
    return structure, chains


def residue_pair_contacts(
    chain_a: dict[str, object], chain_b: dict[str, object], cutoff: float
) -> tuple[set[tuple[int, int]], int]:
    tree_a = cKDTree(chain_a["coords"])
    tree_b = cKDTree(chain_b["coords"])
    neighbors = tree_a.query_ball_tree(tree_b, cutoff)
    pairs = set()
    atom_pairs = 0
    residues_a = chain_a["resids"]
    residues_b = chain_b["resids"]
    for atom_a, hits in enumerate(neighbors):
        atom_pairs += len(hits)
        pairs.update((int(residues_a[atom_a]), int(residues_b[atom_b])) for atom_b in hits)
    return pairs, atom_pairs


def close_atom_pairs(chain_a: dict[str, object], chain_b: dict[str, object], cutoff: float) -> int:
    return sum(len(hits) for hits in cKDTree(chain_a["coords"]).query_ball_tree(
        cKDTree(chain_b["coords"]), cutoff
    ))


def close_atom_pair_residues(
    chain_a: dict[str, object], chain_b: dict[str, object], cutoff: float
) -> list[tuple[int, int]]:
    neighbors = cKDTree(chain_a["coords"]).query_ball_tree(
        cKDTree(chain_b["coords"]), cutoff
    )
    residues_a = np.asarray(chain_a["resids"])
    residues_b = np.asarray(chain_b["resids"])
    return [
        (int(residues_a[index_a]), int(residues_b[index_b]))
        for index_a, hits in enumerate(neighbors)
        for index_b in hits
    ]


def ca_coords(chain: dict[str, object]) -> tuple[np.ndarray, np.ndarray]:
    mask = np.asarray(chain["names"] == "CA")
    return np.asarray(chain["resids"])[mask], np.asarray(chain["coords"])[mask]


def sg_coords(chain: dict[str, object]) -> dict[int, np.ndarray]:
    mask = (np.asarray(chain["names"]) == "SG") & (np.asarray(chain["resnames"]) == "CYS")
    return {
        int(residue): coord
        for residue, coord in zip(np.asarray(chain["resids"])[mask], np.asarray(chain["coords"])[mask])
    }


def orientless_pair_similarity(left: set[tuple[int, int]], right: set[tuple[int, int]]) -> float:
    swapped = {(b, a) for a, b in left}
    return max(jaccard(left, right), jaccard(swapped, right))


def within_model_interface_types(
    pair_maps: dict[tuple[str, str], set[tuple[int, int]]], effective: list[tuple[str, str]]
) -> list[list[tuple[str, str]]]:
    if not effective:
        return []
    n = len(effective)
    similarity = np.eye(n)
    for left in range(n):
        for right in range(left + 1, n):
            value = orientless_pair_similarity(pair_maps[effective[left]], pair_maps[effective[right]])
            similarity[left, right] = similarity[right, left] = value
    return [[effective[index] for index in cluster] for cluster in graph_clusters(similarity, 0.45)]


def relabel_map(
    contact_map: set[tuple[str, str, int, int]], mapping: dict[str, str]
) -> set[tuple[str, str, int, int]]:
    output = set()
    for chain_a, chain_b, residue_a, residue_b in contact_map:
        new_a, new_b = mapping[chain_a], mapping[chain_b]
        if new_a < new_b:
            output.add((new_a, new_b, residue_a, residue_b))
        else:
            output.add((new_b, new_a, residue_b, residue_a))
    return output


def permutation_jaccard(
    left: set[tuple[str, str, int, int]], right: set[tuple[str, str, int, int]], chains: list[str]
) -> float:
    best = 0.0
    for permutation in itertools.permutations(chains):
        mapping = dict(zip(chains, permutation, strict=True))
        best = max(best, jaccard(relabel_map(left, mapping), right))
    return best


def reference_ca(path: pathlib.Path, chain_id: str) -> dict[int, np.ndarray]:
    _, chains = structure_arrays(path)
    residues, coords = ca_coords(chains[chain_id])
    return {int(residue): coord for residue, coord in zip(residues, coords)}


def score_one(
    cif: pathlib.Path,
    mutation_positions: set[int],
    ref_ca: dict[int, np.ndarray],
) -> tuple[dict[str, object], set[tuple[str, str, int, int]]]:
    structure, chains_data = structure_arrays(cif)
    chains = sorted(chains_data)
    if len(chains) != 4:
        raise ValueError(f"{cif}: expected four chains, found {chains}")
    pair_maps: dict[tuple[str, str], set[tuple[int, int]]] = {}
    pair_atom_contacts = {}
    atomic_clashes = 0
    terminal_atomic_clashes = 0
    design_involved_atomic_clashes = 0
    ca_overlaps = 0
    for chain_a, chain_b in itertools.combinations(chains, 2):
        pair = (chain_a, chain_b)
        pair_maps[pair], pair_atom_contacts[pair] = residue_pair_contacts(
            chains_data[chain_a], chains_data[chain_b], 5.3
        )
        clash_residues = close_atom_pair_residues(
            chains_data[chain_a], chains_data[chain_b], 2.4
        )
        atomic_clashes += len(clash_residues)
        terminal_atomic_clashes += sum(
            (residue_a <= 3 or residue_a >= 384)
            and (residue_b <= 3 or residue_b >= 384)
            for residue_a, residue_b in clash_residues
        )
        design_involved_atomic_clashes += sum(
            residue_a in mutation_positions or residue_b in mutation_positions
            for residue_a, residue_b in clash_residues
        )
        _, ca_a = ca_coords(chains_data[chain_a])
        _, ca_b = ca_coords(chains_data[chain_b])
        ca_overlaps += sum(len(hits) for hits in cKDTree(ca_a).query_ball_tree(cKDTree(ca_b), 2.5))

    interface_residues_by_pair = {
        pair: {a for a, _ in contacts} | {b for _, b in contacts}
        for pair, contacts in pair_maps.items()
    }
    # The plan's "20 residue contacts" is implemented as 20 distinct
    # residue-residue contact pairs (not the union of positions on both sides).
    effective = [pair for pair in pair_maps if len(pair_maps[pair]) >= 20]
    component_sizes = sorted((len(group) for group in components(chains, effective)), reverse=True)
    degree = {chain: 0 for chain in chains}
    for chain_a, chain_b in effective:
        degree[chain_a] += 1
        degree[chain_b] += 1
    types = within_model_interface_types(pair_maps, effective)
    type_design_counts = []
    for group in types:
        positions = set()
        for pair in group:
            positions |= interface_residues_by_pair[pair]
        type_design_counts.append(len(positions & mutation_positions))

    all_interface_positions = set().union(*(interface_residues_by_pair.values()))
    designed_interface = mutation_positions & all_interface_positions
    native_ss = 0
    interchain_ss = 0
    sg_by_chain = {chain: sg_coords(chains_data[chain]) for chain in chains}
    for chain in chains:
        sg = sg_by_chain[chain]
        if 74 in sg and 121 in sg and np.linalg.norm(sg[74] - sg[121]) <= 2.5:
            native_ss += 1
    for chain_a, chain_b in itertools.combinations(chains, 2):
        for coord_a in sg_by_chain[chain_a].values():
            for coord_b in sg_by_chain[chain_b].values():
                if np.linalg.norm(coord_a - coord_b) <= 2.5:
                    interchain_ss += 1

    rmsds = []
    plddts = []
    for chain in chains:
        residues, coords = ca_coords(chains_data[chain])
        query = {int(residue): coord for residue, coord in zip(residues, coords)}
        common = sorted(set(query) & set(ref_ca))
        rmsds.append(kabsch_rmsd(
            np.asarray([ref_ca[position] for position in common]),
            np.asarray([query[position] for position in common]),
        ))
        ca_mask = np.asarray(chains_data[chain]["names"] == "CA")
        plddts.append(float(np.mean(np.asarray(chains_data[chain]["plddt"])[ca_mask])))

    full_map = {
        (chain_a, chain_b, residue_a, residue_b)
        for (chain_a, chain_b), contacts in pair_maps.items()
        for residue_a, residue_b in contacts
    }
    row = {
        "atomic_clash_count_2p4": atomic_clashes,
        "terminal_atomic_clash_count_2p4": terminal_atomic_clashes,
        "nonterminal_atomic_clash_count_2p4": atomic_clashes - terminal_atomic_clashes,
        "design_involved_atomic_clash_count_2p4": design_involved_atomic_clashes,
        "severe_ca_overlap_count_2p5": ca_overlaps,
        "effective_interface_edges": len(effective),
        "effective_edge_list": ",".join("-".join(pair) for pair in effective),
        "largest_component": component_sizes[0],
        "min_chain_degree": min(degree.values()),
        "valid_c4_network": int(component_sizes[0] == 4 and min(degree.values()) >= 2),
        "interface_type_count": len(types),
        "interface_type_multiplicities": ",".join(str(len(group)) for group in types),
        "interface_type_edges": "|".join(
            ",".join("-".join(pair) for pair in group) for group in types
        ),
        "interface_type_design_counts": ",".join(str(value) for value in type_design_counts),
        "two_interfaces_each_design_ge3": int(
            sum(value >= 3 for value in type_design_counts) >= 2
        ),
        "contact_residue_pairs": len(full_map),
        "interface_positions": len(all_interface_positions),
        "design_interface_count": len(designed_interface),
        # WT/P3 calibration controls intentionally have no design positions.
        # Their geometry is still scored, but this design-specific ratio is 0
        # rather than an undefined division by zero.
        "design_interface_fraction": (
            len(designed_interface) / len(mutation_positions)
            if mutation_positions else 0.0
        ),
        "design_interface_positions": ",".join(map(str, sorted(designed_interface))),
        "native_c74_c121_bonds": native_ss,
        "interchain_ss_bonds": interchain_ss,
        "mean_chain_plddt": statistics.mean(plddts),
        "min_chain_plddt": min(plddts),
        "mean_chain_ca_rmsd_to_1ova": statistics.mean(rmsds),
        "max_chain_ca_rmsd_to_1ova": max(rmsds),
        "cif_path": str(cif),
    }
    return row, full_map


def mean(values) -> float:
    values = list(values)
    return statistics.mean(values) if values else math.nan


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--final", type=pathlib.Path, required=True)
    parser.add_argument("--per-model", type=pathlib.Path, required=True)
    parser.add_argument("--reference", type=pathlib.Path, required=True)
    parser.add_argument("--reference-chain", default="A")
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    final = {row["candidate"]: row for row in read_tsv(args.final)}
    model_input = read_tsv(args.per_model)
    ref_ca = reference_ca(args.reference, args.reference_chain)
    scored = []
    contact_maps: dict[str, list[set[tuple[str, str, int, int]]]] = defaultdict(list)
    for index, source_row in enumerate(model_input, 1):
        candidate = source_row["candidate"]
        mutation_positions = {int(value[1:-1]) for value in final[candidate]["mutations"].split(",")}
        metrics, contact_map = score_one(pathlib.Path(source_row["cif_path"]), mutation_positions, ref_ca)
        row = {
            "candidate": candidate,
            "seed": source_row["seed"],
            "model": source_row["model"],
            "iptm": source_row["iptm"],
            "min_incident_iptm": source_row["min_incident_iptm"],
            "interface_pae": source_row["contact_pair_pae_median"],
            **metrics,
        }
        scored.append(row)
        contact_maps[candidate].append(contact_map)
        if index % 25 == 0:
            print(f"scored {index}/{len(model_input)} models ({candidate})", flush=True)
    grouped = defaultdict(list)
    for row in scored:
        grouped[row["candidate"]].append(row)
    summaries = []
    topology_rows = []
    for candidate, base in final.items():
        rows = grouped[candidate]
        maps = contact_maps[candidate]
        if len(rows) != 25 or len(maps) != 25:
            raise ValueError(f"{candidate}: expected 25 models, got {len(rows)}")
        similarity = np.eye(25)
        for left in range(25):
            for right in range(left + 1, 25):
                value = permutation_jaccard(maps[left], maps[right], list("ABCD"))
                similarity[left, right] = similarity[right, left] = value
                topology_rows.append({
                    "candidate": candidate,
                    "model_a": rows[left]["model"],
                    "model_b": rows[right]["model"],
                    "best_chain_permutation_contact_jaccard": value,
                })
        clusters = graph_clusters(similarity, 0.60)
        main_cluster = set(clusters[0])
        medoid = max(
            main_cluster,
            key=lambda index: statistics.mean(similarity[index, other] for other in main_cluster),
        )
        for cluster_id, cluster in enumerate(clusters, 1):
            for index in cluster:
                rows[index]["topology_cluster"] = cluster_id
                rows[index]["in_main_topology_cluster"] = int(cluster_id == 1)
                rows[index]["main_topology_medoid"] = int(index == medoid)
        design_recurrence = []
        mutation_positions = {int(value[1:-1]) for value in base["mutations"].split(",")}
        for position in mutation_positions:
            hits = sum(position in {
                residue
                for _, _, residue_a, residue_b in contact_map
                for residue in (residue_a, residue_b)
            } for contact_map in maps)
            design_recurrence.append(hits / 25)
        iptms = [float(row["iptm"]) for row in rows]
        weakest = [float(row["min_incident_iptm"]) for row in rows]
        pae = [float(row["interface_pae"]) for row in rows]
        summary = {
            "candidate": candidate,
            "n_models": len(rows),
            "mean_iptm": mean(iptms),
            "median_iptm": statistics.median(iptms),
            "min_iptm": min(iptms),
            "min_seed_mean_iptm": min(mean(float(row["iptm"]) for row in rows if int(row["seed"]) == seed) for seed in range(1, 6)),
            "mean_weakest_interface_iptm": mean(weakest),
            "min_weakest_interface_iptm": min(weakest),
            "mean_interface_pae": mean(pae),
            "max_interface_pae": max(pae),
            "zero_atomic_clash_models": sum(int(row["atomic_clash_count_2p4"]) == 0 for row in rows),
            "terminal_only_clash_models": sum(
                int(row["atomic_clash_count_2p4"]) > 0
                and int(row["terminal_atomic_clash_count_2p4"])
                == int(row["atomic_clash_count_2p4"])
                for row in rows
            ),
            "nonterminal_clash_models": sum(
                int(row["nonterminal_atomic_clash_count_2p4"]) > 0 for row in rows
            ),
            "design_involved_clash_models": sum(
                int(row["design_involved_atomic_clash_count_2p4"]) > 0 for row in rows
            ),
            "zero_ca_overlap_models": sum(int(row["severe_ca_overlap_count_2p5"]) == 0 for row in rows),
            "valid_c4_network_models": sum(int(row["valid_c4_network"]) for row in rows),
            "zero_interchain_ss_models": sum(int(row["interchain_ss_bonds"]) == 0 for row in rows),
            "native_ss_complete_models": sum(int(row["native_c74_c121_bonds"]) == 4 for row in rows),
            "two_interfaces_design_ge3_models": sum(int(row["two_interfaces_each_design_ge3"]) for row in rows),
            "main_topology_cluster_size": len(main_cluster),
            "topology_cluster_sizes": ",".join(map(str, map(len, clusters))),
            "mean_pairwise_permutation_jaccard": mean(similarity[np.triu_indices(25, 1)]),
            "mean_design_interface_recurrence": mean(design_recurrence),
            "design_mutations_recurrence_ge80_fraction": sum(value >= 0.80 for value in design_recurrence) / len(design_recurrence),
            "mean_chain_plddt": mean(float(row["mean_chain_plddt"]) for row in rows),
            "min_model_chain_plddt": min(float(row["min_chain_plddt"]) for row in rows),
            "mean_chain_ca_rmsd_to_1ova": mean(float(row["mean_chain_ca_rmsd_to_1ova"]) for row in rows),
            "max_chain_ca_rmsd_to_1ova": max(float(row["max_chain_ca_rmsd_to_1ova"]) for row in rows),
        }
        gates = {
            "gate_af3_model_count_25": len(rows) == 25,
            "gate_af3_mean_iptm_gt_0p80": summary["mean_iptm"] > 0.80,
            "gate_af3_median_iptm_ge_0p80": summary["median_iptm"] >= 0.80,
            "gate_af3_min_iptm_ge_0p75": summary["min_iptm"] >= 0.75,
            "gate_af3_each_seed_mean_ge_0p78": summary["min_seed_mean_iptm"] >= 0.78,
            "gate_af3_mean_weakest_ge_0p75": summary["mean_weakest_interface_iptm"] >= 0.75,
            "gate_af3_min_weakest_ge_0p70": summary["min_weakest_interface_iptm"] >= 0.70,
            "gate_af3_mean_pae_le_4p0": summary["mean_interface_pae"] <= 4.0,
            "gate_af3_max_pae_le_5p0": summary["max_interface_pae"] <= 5.0,
            "gate_af3_25_zero_atomic_clash": summary["zero_atomic_clash_models"] == 25,
            "gate_af3_25_zero_ca_overlap": summary["zero_ca_overlap_models"] == 25,
            "gate_af3_25_valid_c4_network": summary["valid_c4_network_models"] == 25,
            "gate_af3_25_zero_interchain_ss": summary["zero_interchain_ss_models"] == 25,
            "gate_af3_25_native_ss_complete": summary["native_ss_complete_models"] == 25,
            "gate_af3_main_topology_cluster_ge20": summary["main_topology_cluster_size"] >= 20,
            "gate_af3_mean_design_recurrence_ge_0p65": summary["mean_design_interface_recurrence"] >= 0.65,
            "gate_af3_design_ge80_fraction_ge_0p60": summary["design_mutations_recurrence_ge80_fraction"] >= 0.60,
            "gate_af3_two_interfaces_design_ge3_in_ge20_models": summary["two_interfaces_design_ge3_models"] >= 20,
            "gate_af3_mean_plddt_ge_85": summary["mean_chain_plddt"] >= 85.0,
            "gate_af3_max_chain_rmsd_le_2p0": summary["max_chain_ca_rmsd_to_1ova"] <= 2.0,
        }
        summary.update({key: int(value) for key, value in gates.items()})
        failed = [key for key, value in gates.items() if not value]
        summary["accepted_af3_recomputed"] = int(not failed)
        summary["failed_af3_gates"] = ";".join(failed)
        summaries.append(summary)
    # Cluster labels are computed candidate-wise after all contact maps are scored,
    # so the per-model table is written only after those labels have been attached.
    write_tsv(args.out / "af3_per_model_recomputed.tsv", scored)
    write_tsv(args.out / "af3_topology_pairwise.tsv", topology_rows)
    write_tsv(args.out / "af3_candidate_recomputed_summary.tsv", summaries)
    print(f"wrote {len(scored)} model rows and {len(summaries)} candidate summaries to {args.out}")


if __name__ == "__main__":
    main()
