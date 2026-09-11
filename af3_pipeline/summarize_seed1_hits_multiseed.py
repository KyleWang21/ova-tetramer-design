#!/usr/bin/env python3
"""Aggregate strict seed1 plus independent seeds2--5 into 25-model gates."""

from __future__ import annotations

import argparse
import csv
import pathlib
import re
import statistics
from collections import Counter


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing empty output: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def seed_number(model: str) -> int:
    match = re.match(r"seed-(\d+)_sample-\d+$", model)
    if match is None:
        raise ValueError(f"unexpected AF3 model name: {model}")
    return int(match.group(1))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=pathlib.Path, required=True)
    parser.add_argument("--multiseed", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    selected = read_tsv(args.multiseed / "selected_candidates.tsv")
    if not selected:
        args.out.mkdir(parents=True, exist_ok=True)
        with (args.out / "full25_summary.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow([
                "system", "source_system", "n_models", "mean_iptm",
                "accepted_af3", "failed_gates", "sequence",
            ])
        (args.out / "NO_SEED1_HITS").touch()
        print("no selected seed-1 hits")
        return

    source_rows = read_tsv(
        args.source / "seed1_strict_final" / "seed1_strict_per_model.tsv"
    )
    source_by_system: dict[str, list[dict[str, str]]] = {}
    for row in source_rows:
        source_by_system.setdefault(row["system"], []).append(row)
    reseed_by_system: dict[str, list[dict[str, str]]] = {}
    for path in sorted(args.multiseed.glob("strict_s*/seed1_strict_per_model.tsv")):
        for row in read_tsv(path):
            reseed_by_system.setdefault(row["system"], []).append(row)

    summaries: list[dict[str, object]] = []
    all_models: list[dict[str, object]] = []
    for candidate in selected:
        system = candidate["system"]
        source_system = candidate["source_system"]
        models = source_by_system.get(source_system, []) + reseed_by_system.get(system, [])
        models.sort(key=lambda row: (seed_number(row["model"]), row["model"]))
        seed_counts: dict[int, int] = {}
        for row in models:
            seed = seed_number(row["model"])
            seed_counts[seed] = seed_counts.get(seed, 0) + 1
            all_models.append({
                **row,
                "system": system,
                "source_system": source_system,
                "seed": seed,
            })
        iptm = [float(row["iptm"]) for row in models]
        design_fractions = [float(row["design_interface_fraction"]) for row in models]
        design_counts = [int(row["design_interface_count"]) for row in models]
        design_position_counts: Counter[int] = Counter()
        for row in models:
            design_position_counts.update(
                int(value) for value in row.get("design_interface_positions", "").split(",")
                if value
            )
        union_design_positions = sorted(design_position_counts)
        consensus_design_positions = sorted(
            position for position, count in design_position_counts.items()
            if count == len(models)
        )
        gates = {
            "gate_25_models_five_seeds": seed_counts == {1: 5, 2: 5, 3: 5, 4: 5, 5: 5},
            "gate_nmut_lt25": int(candidate["n_mutations"]) < 25,
            "gate_surface_fraction_ge0p80": float(candidate["surface_mutation_fraction"]) >= 0.80,
            "gate_mean_iptm_ge0p80": len(iptm) == 25 and statistics.fmean(iptm) >= 0.80,
            "gate_min_iptm_ge0p75": len(iptm) == 25 and min(iptm) >= 0.75,
            "gate_25_zero_atomic_clash": len(models) == 25 and all(
                int(row["atomic_clash_count_2p4"]) == 0 for row in models
            ),
            "gate_25_zero_ca_overlap": len(models) == 25 and all(
                int(row["severe_ca_overlap_count_2p5"]) == 0 for row in models
            ),
            "gate_25_valid_c4_network": len(models) == 25 and all(
                int(row["valid_c4_network"]) == 1 for row in models
            ),
            "gate_25_native_ss_complete": len(models) == 25 and all(
                int(row["native_c74_c121_bonds"]) == 4 for row in models
            ),
            "gate_25_zero_interchain_ss": len(models) == 25 and all(
                int(row["interchain_ss_bonds"]) == 0 for row in models
            ),
            "gate_25_two_interfaces_design_ge3": len(models) == 25 and all(
                int(row["two_interfaces_each_design_ge3"]) == 1 for row in models
            ),
            "gate_25_max_chain_rmsd_le2p0": len(models) == 25 and all(
                float(row["max_chain_ca_rmsd_to_1ova"]) <= 2.0 for row in models
            ),
        }
        failed = [name for name, passed in gates.items() if not passed]
        summaries.append({
            "system": system,
            "source_system": source_system,
            "candidate_rank": candidate["candidate_rank"],
            "n_mutations": candidate["n_mutations"],
            "surface_mutation_fraction": candidate["surface_mutation_fraction"],
            "n_models": len(models),
            "seed_counts": ";".join(f"{seed}:{seed_counts.get(seed, 0)}" for seed in range(1, 6)),
            "mean_iptm": f"{statistics.fmean(iptm):.6f}" if iptm else "",
            "median_iptm": f"{statistics.median(iptm):.6f}" if iptm else "",
            "min_iptm": f"{min(iptm):.6f}" if iptm else "",
            "max_iptm": f"{max(iptm):.6f}" if iptm else "",
            "target_mean_iptm_ge0p90": int(len(iptm) == 25 and statistics.fmean(iptm) >= 0.90),
            "models_two_interfaces_design_ge3": sum(
                int(row["two_interfaces_each_design_ge3"]) for row in models
            ),
            "mean_design_interface_recurrence": (
                f"{statistics.fmean(design_fractions):.6f}" if design_fractions else ""
            ),
            "min_design_interface_recurrence": (
                f"{min(design_fractions):.6f}" if design_fractions else ""
            ),
            "mean_design_interface_count": (
                f"{statistics.fmean(design_counts):.6f}" if design_counts else ""
            ),
            "union_design_interface_position_count": len(union_design_positions),
            "union_design_interface_positions": ",".join(map(str, union_design_positions)),
            "all_model_design_interface_position_count": len(consensus_design_positions),
            "all_model_design_interface_positions": ",".join(map(str, consensus_design_positions)),
            **{name: int(passed) for name, passed in gates.items()},
            "accepted_af3": int(not failed),
            "failed_gates": ";".join(failed),
            "sequence": candidate["sequence"],
        })

    summaries.sort(key=lambda row: (
        -int(row["accepted_af3"]),
        -float(row["mean_iptm"] or -1.0),
        int(row["candidate_rank"]),
    ))
    args.out.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out / "full25_summary.tsv", summaries)
    write_tsv(args.out / "full25_per_model.tsv", all_models)
    accepted = sum(int(row["accepted_af3"]) for row in summaries)
    print(f"candidates={len(summaries)} accepted_af3={accepted}")
    for row in summaries:
        print(
            row["system"], f"models={row['n_models']}",
            f"mean={row['mean_iptm'] or 'NA'}", f"pass={row['accepted_af3']}"
        )


if __name__ == "__main__":
    main()
