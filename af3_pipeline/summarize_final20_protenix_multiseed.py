#!/usr/bin/env python3
"""Score five Protenix-v2 seeds and compare their C4 contact maps with AF3."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import re
import statistics
from collections import defaultdict

import numpy as np

from screen_final20_af3_structures import (
    graph_clusters, permutation_jaccard, reference_ca, score_one,
)


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def candidate_token(candidate: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", candidate.lower()).strip("_")


def numeric(row: dict[str, str], field: str, default: float) -> float:
    try:
        return float(row.get(field, ""))
    except (TypeError, ValueError):
        return default


def af3_representative_rank(row: dict[str, str]) -> tuple[float, ...]:
    """Rank both legacy screen tables and generic strict 25-model tables."""
    zero_clash = int(numeric(row, "atomic_clash_count_2p4", 0.0) == 0.0)
    valid_c4 = int(numeric(row, "valid_c4_network", 1.0) == 1.0)
    two_interfaces = int(numeric(row, "two_interfaces_each_design_ge3", 1.0) == 1.0)
    weakest = numeric(
        row, "min_incident_iptm",
        numeric(row, "min_chain_plddt", numeric(row, "mean_chain_plddt", 0.0)) / 100.0,
    )
    geometry_penalty = numeric(
        row, "interface_pae", numeric(row, "max_chain_ca_rmsd_to_1ova", 999.0)
    )
    return (
        zero_clash, valid_c4, two_interfaces,
        numeric(row, "iptm", 0.0), weakest, -geometry_penalty,
    )


def find_outputs(roots: list[pathlib.Path], candidate: str) -> list[tuple[int, pathlib.Path, pathlib.Path]]:
    token = candidate_token(candidate)
    found = []
    for root in roots:
        for summary in root.glob(f"**/*{token}*summary_confidence_sample_0.json"):
            seed_dir = next((parent for parent in summary.parents if parent.name.startswith("seed_")), None)
            if seed_dir is None:
                continue
            seed = int(seed_dir.name.split("_", 1)[1])
            cif_name = summary.name.replace("_summary_confidence_sample_0.json", "_sample_0.cif")
            cif = summary.with_name(cif_name)
            if cif.exists():
                found.append((seed, summary, cif))
    unique = {seed: (seed, summary, cif) for seed, summary, cif in found}
    return [unique[seed] for seed in sorted(unique)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--final", type=pathlib.Path, required=True)
    parser.add_argument("--af3-models", type=pathlib.Path, required=True)
    parser.add_argument("--protenix-root", type=pathlib.Path, action="append", required=True)
    parser.add_argument("--reference", type=pathlib.Path, required=True)
    parser.add_argument("--shared-nonvdw", type=pathlib.Path,
                        help="optional candidate table with shared_nonvdw_position_type_count")
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    final = {row["candidate"]: row for row in read_tsv(args.final)}
    af3_models = defaultdict(list)
    for row in read_tsv(args.af3_models):
        af3_models[row["candidate"]].append(row)
    shared = {}
    if args.shared_nonvdw and args.shared_nonvdw.exists():
        shared = {row["candidate"]: int(row["shared_nonvdw_position_type_count"])
                  for row in read_tsv(args.shared_nonvdw)}
    ref = reference_ca(args.reference, "A")
    per_model = []
    summaries = []
    for candidate, base in final.items():
        mutation_positions = {int(value[1:-1]) for value in base["mutations"].split(",")}
        outputs = find_outputs(args.protenix_root, candidate)
        if not af3_models[candidate]:
            raise RuntimeError(f"no AF3 model rows found for {candidate}")
        af3_best = max(af3_models[candidate], key=af3_representative_rank)
        _, af3_map = score_one(pathlib.Path(af3_best["cif_path"]), mutation_positions, ref)
        maps = []
        rows = []
        for seed, summary_path, cif in outputs:
            confidence = json.loads(summary_path.read_text())
            geometry, contact_map = score_one(cif, mutation_positions, ref)
            maps.append(contact_map)
            chains = list("ABCD")
            edge_indices = []
            for edge in str(geometry["effective_edge_list"]).split(","):
                if not edge:
                    continue
                left, right = edge.split("-")
                edge_indices.append((chains.index(left), chains.index(right)))
            pair_iptm = confidence.get("chain_pair_iptm") or []
            pair_gpde = confidence.get("chain_pair_gpde") or []
            incident = [max(row[:index] + row[index + 1:]) for index, row in enumerate(pair_iptm)]
            interface_gpde = [float(pair_gpde[left][right]) for left, right in edge_indices] if pair_gpde else []
            row = {
                "candidate": candidate, "seed": seed,
                "iptm": confidence.get("iptm"), "ptm": confidence.get("ptm"),
                "plddt_normalized": float(confidence.get("plddt", 0.0)) / 100.0,
                "min_incident_iptm": min(incident) if incident else "",
                "mean_effective_edge_gpde": statistics.mean(interface_gpde) if interface_gpde else "",
                **geometry,
                "summary_path": str(summary_path), "cif_path": str(cif),
            }
            rows.append(row); per_model.append(row)
        if maps:
            similarity = np.eye(len(maps))
            for left in range(len(maps)):
                for right in range(left + 1, len(maps)):
                    value = permutation_jaccard(maps[left], maps[right], list("ABCD"))
                    similarity[left, right] = similarity[right, left] = value
            clusters = graph_clusters(similarity, 0.60)
            main_indices = clusters[0]
            medoid = max(main_indices, key=lambda index: statistics.mean(similarity[index, other] for other in main_indices))
            crossmodel = permutation_jaccard(maps[medoid], af3_map, list("ABCD"))
            for cluster_id, cluster in enumerate(clusters, 1):
                for index in cluster:
                    rows[index]["topology_cluster"] = cluster_id
                    rows[index]["in_main_topology_cluster"] = int(cluster_id == 1)
                    rows[index]["main_topology_medoid"] = int(index == medoid)
        else:
            clusters, crossmodel, medoid = [], None, None
        n = len(rows)
        get = lambda field: [float(row[field]) for row in rows]
        summary = {
            "candidate": candidate, "n_seeds": n,
            "seeds": ",".join(str(row["seed"]) for row in rows),
            "mean_iptm": statistics.mean(get("iptm")) if rows else "",
            "median_iptm": statistics.median(get("iptm")) if rows else "",
            "min_iptm": min(get("iptm")) if rows else "",
            "mean_min_incident_iptm": statistics.mean(get("min_incident_iptm")) if rows else "",
            "mean_ptm": statistics.mean(get("ptm")) if rows else "",
            "mean_plddt_normalized": statistics.mean(get("plddt_normalized")) if rows else "",
            "mean_interface_gpde": statistics.mean(get("mean_effective_edge_gpde")) if rows else "",
            "full_connected_zero_clash_clean_models": sum(
                int(row["valid_c4_network"]) and int(row["atomic_clash_count_2p4"]) == 0
                and int(row["interchain_ss_bonds"]) == 0 and int(row["native_c74_c121_bonds"]) == 4
                for row in rows
            ),
            "main_topology_cluster_size": len(clusters[0]) if clusters else 0,
            "topology_cluster_sizes": ",".join(str(len(cluster)) for cluster in clusters),
            "af3_protenix_contact_jaccard": crossmodel if crossmodel is not None else "",
            "shared_nonvdw_position_type_count": shared.get(candidate, ""),
            "representative_seed": rows[medoid]["seed"] if medoid is not None else "",
            "representative_cif_path": rows[medoid]["cif_path"] if medoid is not None else "",
        }
        gates = {
            "gate_protenix_5_seeds": n == 5,
            "gate_protenix_5of5_geometry_clean": summary["full_connected_zero_clash_clean_models"] == 5,
            "gate_protenix_mean_iptm_ge0p70": n == 5 and summary["mean_iptm"] >= 0.70,
            "gate_protenix_median_iptm_ge0p70": n == 5 and summary["median_iptm"] >= 0.70,
            "gate_protenix_min_iptm_ge0p60": n == 5 and summary["min_iptm"] >= 0.60,
            "gate_protenix_mean_weakest_ge0p65": n == 5 and summary["mean_min_incident_iptm"] >= 0.65,
            "gate_protenix_mean_ptm_ge0p82": n == 5 and summary["mean_ptm"] >= 0.82,
            "gate_protenix_mean_plddt_ge0p85": n == 5 and summary["mean_plddt_normalized"] >= 0.85,
            "gate_protenix_mean_interface_pae_le7p5": n == 5 and summary["mean_interface_gpde"] <= 7.5,
            "gate_protenix_main_cluster_ge4": summary["main_topology_cluster_size"] >= 4,
            "gate_protenix_af3_jaccard_ge0p40": crossmodel is not None and crossmodel >= 0.40,
            "gate_protenix_shared_nonvdw_ge8": shared.get(candidate, 0) >= 8,
        }
        summary.update({key: int(value) for key, value in gates.items()})
        summary["crossmodel_validated"] = int(all(gates.values()))
        summary["failed_protenix_gates"] = ";".join(key for key, value in gates.items() if not value)
        summaries.append(summary)
    if per_model:
        write_tsv(args.out / "protenix_per_seed_recomputed.tsv", per_model)
    write_tsv(args.out / "protenix_candidate_summary.tsv", summaries)
    print(f"scored {len(per_model)} Protenix models; complete candidates={sum(int(row['n_seeds']) == 5 for row in summaries)}")


if __name__ == "__main__":
    main()
