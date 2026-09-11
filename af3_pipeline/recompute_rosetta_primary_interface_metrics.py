#!/usr/bin/env python3
"""Recompute Rosetta chemistry gates on the fixed primary design interfaces.

This is a post-hoc migration for an already completed 200-relax ensemble.  It
does not rerun Rosetta: all per-pair InterfaceAnalyzer values are read from
``pair_metrics_json``.  Geometry/topology fields are retained from the
original run, while the four chemistry summaries are recomputed over the
primary edges recorded in the v1.1 manifest.  Missing primary-pair data is a
fail-closed error rather than an implicit pass.
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import statistics


def edges(value: str) -> list[str]:
    out = []
    seen = set()
    for token in value.replace("|", ",").split(","):
        token = token.strip()
        if not token:
            continue
        left, right = token.split("-", 1)
        name = "-".join(sorted((left, right)))
        if name not in seen:
            seen.add(name)
            out.append(name)
    return out


def f(value: object) -> float:
    return float(value)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicates", type=pathlib.Path, required=True)
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    with args.manifest.open() as handle:
        manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(manifest_rows) != 1:
        raise RuntimeError(f"expected one manifest row, found {len(manifest_rows)}")
    manifest = manifest_rows[0]
    candidate = manifest["candidate"]
    source_model = manifest.get("source_model", "AF3")
    effective_edges = edges(manifest.get("effective_edges", ""))
    design_edges = edges(manifest.get("design_interface_edges", "") or manifest.get("effective_edges", ""))
    if not design_edges:
        raise RuntimeError("manifest has no primary design interface edges")
    if not set(design_edges).issubset(set(effective_edges)):
        raise RuntimeError("primary design edges are not a subset of effective edges")

    with args.replicates.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise RuntimeError("empty replicate table")
    for row in rows:
        if row.get("candidate") != candidate or row.get("source_model") != source_model:
            raise RuntimeError("replicate candidate/source does not match manifest")
        pair = json.loads(row["pair_metrics_json"])
        missing = [edge for edge in design_edges if edge not in pair]
        if missing:
            raise RuntimeError(
                f"replicate {row.get('replicate')} is missing primary interface metrics: {','.join(missing)}"
            )
        all_values = list(pair.values())
        selected = [pair[edge] for edge in design_edges]
        row["design_interface_edges"] = ",".join(design_edges)
        row["background_edges"] = ",".join(edge for edge in effective_edges if edge not in design_edges)
        row["weakest_sc"] = str(min(f(item["sc_value"]) for item in selected))
        row["weakest_hbonds"] = str(min(f(item["hbonds_int"]) for item in selected))
        row["worst_dG_dSASA_x100"] = str(max(f(item["dG_dSASA_ratio"]) for item in selected))
        row["worst_unsat_per_1000A2"] = str(max(f(item["delta_unsatHbonds"]) / f(item["dSASA_int"]) * 1000.0 for item in selected))
        row["all_weakest_sc"] = str(min(f(item["sc_value"]) for item in all_values))
        row["all_weakest_hbonds"] = str(min(f(item["hbonds_int"]) for item in all_values))
        row["all_worst_dG_dSASA_x100"] = str(max(f(item["dG_dSASA_ratio"]) for item in all_values))
        row["all_worst_unsat_per_1000A2"] = str(max(f(item["delta_unsatHbonds"]) / f(item["dSASA_int"]) * 1000.0 for item in all_values))

    args.out.mkdir(parents=True, exist_ok=True)
    stem = f"{candidate}_{source_model}"
    replicate_out = args.out / f"{stem}_rosetta_replicates.tsv"
    fields = list(rows[0])
    with replicate_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)

    def mean(field: str) -> float:
        return statistics.mean(int(float(row[field])) for row in rows)

    summary: dict[str, object] = {
        "candidate": candidate,
        "source_model": source_model,
        "nrelax": len(rows),
        "production_ensemble_complete": int(len(rows) >= 100),
        "design_interface_edges": ",".join(design_edges),
        "background_edges": ",".join(edge for edge in effective_edges if edge not in design_edges),
        "topology_retained_fraction": mean("valid_c4_network"),
        "zero_clash_fraction": statistics.mean(int(float(row["atomic_clash_count_2p4"])) == 0 for row in rows),
        "zero_overlap_fraction": statistics.mean(int(float(row["severe_ca_overlap_count_2p5"])) == 0 for row in rows),
        "median_weakest_sc": statistics.median(f(row["weakest_sc"]) for row in rows),
        "median_weakest_hbonds": statistics.median(f(row["weakest_hbonds"]) for row in rows),
        "median_worst_dG_dSASA_x100": statistics.median(f(row["worst_dG_dSASA_x100"]) for row in rows),
        "median_worst_unsat_per_1000A2": statistics.median(f(row["worst_unsat_per_1000A2"]) for row in rows),
        "max_chain_ca_rmsd_to_input": max(f(row["max_chain_ca_rmsd_to_input"]) for row in rows),
        "median_all_weakest_sc": statistics.median(f(row["all_weakest_sc"]) for row in rows),
        "median_all_weakest_hbonds": statistics.median(f(row["all_weakest_hbonds"]) for row in rows),
        "median_all_worst_dG_dSASA_x100": statistics.median(f(row["all_worst_dG_dSASA_x100"]) for row in rows),
        "median_all_worst_unsat_per_1000A2": statistics.median(f(row["all_worst_unsat_per_1000A2"]) for row in rows),
    }
    gates = {
        "gate_rosetta_ensemble_100": summary["production_ensemble_complete"] == 1,
        "gate_rosetta_topology_ge0p90": summary["topology_retained_fraction"] >= 0.90,
        "gate_rosetta_zero_clash_fraction_ge0p90": summary["zero_clash_fraction"] >= 0.90,
        "gate_rosetta_zero_overlap_fraction_ge0p90": summary["zero_overlap_fraction"] >= 0.90,
        "gate_rosetta_weakest_sc_ge0p60": summary["median_weakest_sc"] >= 0.60,
        "gate_rosetta_weakest_hbonds_ge3": summary["median_weakest_hbonds"] >= 3.0,
        "gate_rosetta_dG_dSASA_le_minus1": summary["median_worst_dG_dSASA_x100"] <= -1.0,
        "gate_rosetta_unsat_per_1000A2_le2": summary["median_worst_unsat_per_1000A2"] <= 2.0,
        "gate_rosetta_chain_rmsd_le2": summary["max_chain_ca_rmsd_to_input"] <= 2.0,
    }
    summary.update({key: int(value) for key, value in gates.items()})
    summary["rosetta_pass"] = int(all(gates.values()))
    summary["failed_rosetta_gates"] = ";".join(key for key, value in gates.items() if not value)
    for suffix, payload in (("json", json.dumps(summary, indent=2) + "\n"), ("tsv", None)):
        path = args.out / f"{stem}_rosetta_summary.{suffix}"
        if suffix == "json":
            path.write_text(payload)
        else:
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(summary), delimiter="\t", lineterminator="\n")
                writer.writeheader(); writer.writerow(summary)
    print(f"wrote v1.1 primary-interface metrics for {candidate}: rosetta_pass={summary['rosetta_pass']}")


if __name__ == "__main__":
    main()
