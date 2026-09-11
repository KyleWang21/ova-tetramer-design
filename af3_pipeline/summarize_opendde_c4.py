#!/usr/bin/env python3
"""Score official OpenDDE-v1 C4 seeds and select evidence representatives.

DockQ is deliberately reported as DockQ_to_design because the native argument is
an AF3 design target, not an experimental tetramer structure.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
import re
import statistics
import subprocess
import tempfile
from collections import defaultdict

import numpy as np

from screen_final20_af3_structures import (
    graph_clusters, permutation_jaccard, reference_ca, score_one,
)


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("empty table")
    fields = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def ranks(rows: list[dict[str, object]], field: str) -> dict[int, float]:
    """Return descending average ranks (midranks) with deterministic ties."""
    ordered = sorted((float(row[field]) for row in rows), reverse=True)
    positions: dict[float, list[int]] = defaultdict(list)
    for index, value in enumerate(ordered, 1):
        positions[value].append(index)
    rank = {value: statistics.mean(indices) for value, indices in positions.items()}
    return {int(row["seed"]): rank[float(row[field])] for row in rows}


def dockq_to_design(binary: pathlib.Path, model: pathlib.Path, design: pathlib.Path) -> tuple[float, str]:
    with tempfile.TemporaryDirectory(prefix="dockq_to_design_") as directory:
        output = pathlib.Path(directory) / "dockq.json"
        command = [str(binary), str(model), str(design), "--json", str(output), "--n_cpu", "1", "--short"]
        completed = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        if completed.returncode or not output.exists():
            return float("nan"), completed.stderr[-1000:]
        data = json.loads(output.read_text())
        return float(data["GlobalDockQ"]), str(data.get("best_mapping_str", ""))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--out-root", type=pathlib.Path, required=True)
    parser.add_argument("--final", type=pathlib.Path, required=True)
    parser.add_argument("--af3-models", type=pathlib.Path, required=True)
    parser.add_argument("--reference", type=pathlib.Path, required=True)
    parser.add_argument("--dockq-bin", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    final = {row["candidate"]: row for row in read_tsv(args.final)}
    mutation_positions = {
        candidate: {int(value[1:-1]) for value in row["mutations"].split(",")}
        for candidate, row in final.items()
    }
    af3 = defaultdict(list)
    for row in read_tsv(args.af3_models):
        af3[row["candidate"]].append(row)
    design_target = {}
    for candidate, rows in af3.items():
        # Generic v1 staging preserves the strict AF3 geometry columns but
        # older full25_per_model tables do not carry the topology labels that
        # are computed by screen_final20_af3_structures.  In that case there
        # is no evidence for a separate cluster, so select the best AF3 model
        # from the complete set; this does not alter any OpenDDE gate.
        topology_values = {row.get("topology_cluster", "") for row in rows} - {""}
        if topology_values:
            cluster = max(
                topology_values,
                key=lambda value: sum(row.get("topology_cluster", "") == value for row in rows),
            )
            eligible = [row for row in rows if row.get("topology_cluster", "") == cluster]
        else:
            eligible = rows
        def optional_float(row: dict[str, str], field: str, default: float) -> float:
            try:
                return float(row.get(field, default))
            except (TypeError, ValueError):
                return default

        best = max(eligible, key=lambda row: (
            int(row["atomic_clash_count_2p4"]) == 0, int(row["valid_c4_network"]),
            optional_float(row, "iptm", float("-inf")),
            optional_float(row, "min_incident_iptm", float("-inf")),
            -optional_float(row, "interface_pae", float("inf")),
        ))
        design_target[candidate] = pathlib.Path(best["cif_path"])

    systems = []
    for item in read_tsv(args.manifest):
        input_path = args.manifest.parent / item["input_json"]
        job_name = json.loads(input_path.read_text())[0]["name"]
        systems.append((item, job_name))
    ref = reference_ca(args.reference, "A")
    all_rows = []
    summaries = []
    for item, job_name in systems:
        name = item["name"]
        positions = mutation_positions.get(name, set())
        outputs = []
        for summary_path in args.out_root.glob(f"**/{job_name}/seed_*/predictions/*_summary_confidence_sample_0.json"):
            match = re.search(r"seed_(\d+)", str(summary_path))
            if not match:
                continue
            seed = int(match.group(1))
            cif = summary_path.with_name(summary_path.name.replace(
                "_summary_confidence_sample_0.json", "_sample_0.cif"
            ))
            if cif.exists():
                outputs.append((seed, summary_path, cif))
        outputs = sorted({seed: (seed, summary, cif) for seed, summary, cif in outputs}.values())
        rows = []
        maps = []
        for seed, summary_path, cif in outputs:
            confidence = json.loads(summary_path.read_text())
            geometry, contact_map = score_one(cif, positions, ref)
            pair = confidence.get("chain_pair_iptm") or []
            incident = [max(values[:index] + values[index + 1:]) for index, values in enumerate(pair)] if pair else []
            row = {
                "candidate": name, "kind": item["kind"], "seed": seed,
                "iptm": confidence.get("iptm", ""), "ptm": confidence.get("ptm", ""),
                "plddt_normalized": float(confidence.get("plddt", 0.0)) / 100.0,
                "min_incident_iptm": min(incident) if incident else "",
                **geometry, "summary_path": str(summary_path), "cif_path": str(cif),
            }
            rows.append(row); maps.append(contact_map); all_rows.append(row)
        if maps:
            similarity = np.eye(len(maps))
            for left in range(len(maps)):
                for right in range(left + 1, len(maps)):
                    value = permutation_jaccard(maps[left], maps[right], list("ABCD"))
                    similarity[left, right] = similarity[right, left] = value
            clusters = graph_clusters(similarity, 0.60)
            for cluster_id, cluster_rows in enumerate(clusters, 1):
                for index in cluster_rows:
                    rows[index]["topology_cluster"] = cluster_id
                    rows[index]["in_main_topology_cluster"] = int(cluster_id == 1)
        else:
            clusters = []

        geometry_selectable = [row for row in rows if (
            float(row["iptm"]) > 0.80
            and int(row["in_main_topology_cluster"]) == 1
            and int(row["valid_c4_network"]) == 1
            and int(row["atomic_clash_count_2p4"]) == 0
            and int(row["severe_ca_overlap_count_2p5"]) == 0
            and int(row["interchain_ss_bonds"]) == 0
        )]
        selectable = list(geometry_selectable)
        if name in design_target:
            for row in selectable:
                dockq, mapping = dockq_to_design(
                    args.dockq_bin, pathlib.Path(str(row["cif_path"])), design_target[name]
                )
                row["dockq_to_design"] = dockq
                row["dockq_best_mapping"] = mapping
        selectable = [row for row in selectable if (
            row.get("dockq_to_design") is not None
            and math.isfinite(float(row["dockq_to_design"]))
        )]
        selected = None
        # A partial run must never be promoted to evidence. Selection is made only
        # after all ten frozen seeds are present, even if an early seed looks good.
        complete_10_seeds = len(rows) == 10 and len({int(row["seed"]) for row in rows}) == 10
        if complete_10_seeds and selectable and name in design_target:
            iptm_rank = ranks(selectable, "iptm")
            dockq_rank = ranks(selectable, "dockq_to_design")
            for row in selectable:
                row["opendde_borda"] = 0.5 * iptm_rank[int(row["seed"])] + 0.5 * dockq_rank[int(row["seed"])]
            selected = min(selectable, key=lambda row: (
                float(row["opendde_borda"]),
                -min(float(row["iptm"]), float(row["dockq_to_design"])),
                -float(row["dockq_to_design"]), -float(row["iptm"]), int(row["seed"]),
            ))
            selected["selected_opendde_representative"] = 1
        mean = lambda field: statistics.mean(float(row[field]) for row in rows) if rows else ""
        summaries.append({
            "candidate": name, "kind": item["kind"], "n_seeds": len(rows),
            "mean_iptm": mean("iptm"), "median_iptm": statistics.median(float(row["iptm"]) for row in rows) if rows else "",
            "max_iptm": max((float(row["iptm"]) for row in rows), default=""),
            "mean_ptm": mean("ptm"), "mean_plddt_normalized": mean("plddt_normalized"),
            "main_topology_cluster_size": len(clusters[0]) if clusters else 0,
            "gate_opendde_10_seeds": int(complete_10_seeds),
            "geometry_selectable_seed_count": len(geometry_selectable),
            "selectable_seed_count": len(selectable),
            "selected_seed": selected["seed"] if selected else "",
            "selected_iptm": selected["iptm"] if selected else "",
            "selected_dockq_to_design": selected.get("dockq_to_design", "") if selected else "",
            "selected_cif_path": selected["cif_path"] if selected else "",
            "opendde_evidence_available": int(complete_10_seeds and selected is not None),
        })

    write_tsv(args.out / "opendde_per_seed_recomputed.tsv", all_rows)
    write_tsv(args.out / "opendde_candidate_summary.tsv", summaries)
    controls = [row for row in summaries if row["kind"] == "negative_background"]
    positives = [row for row in summaries if row["kind"] == "process_positive"]
    calibration = {
        "negative_control_count": len(controls), "process_positive_count": len(positives),
        "negative_max_mean_iptm": max((float(row["mean_iptm"]) for row in controls), default=""),
        "positive_min_mean_iptm": min((float(row["mean_iptm"]) for row in positives), default=""),
        "negative_geometry_selectable_seed_count": sum(int(row["geometry_selectable_seed_count"]) for row in controls),
        "positive_geometry_selectable_seed_count": sum(int(row["geometry_selectable_seed_count"]) for row in positives),
        "calibration_all_systems_10_seeds": int(
            len(controls) == 2 and len(positives) == 2
            and all(int(row["gate_opendde_10_seeds"]) for row in controls + positives)
        ),
        "positive_min_minus_negative_max_mean_iptm": (
            min(float(row["mean_iptm"]) for row in positives)
            - max(float(row["mean_iptm"]) for row in controls)
        ) if controls and positives else "",
        "calibration_review_required": 1,
        "promoted_to_required_filter": 0,
    }
    write_tsv(args.out / "opendde_calibration_summary.tsv", [calibration])
    print(f"scored {len(all_rows)} OpenDDE models; representatives={sum(int(row['opendde_evidence_available']) for row in summaries)}")


if __name__ == "__main__":
    main()
