#!/usr/bin/env python3
"""Build one sequence-deduplicated registry of historical and new Accepted_AF3 hits."""

from __future__ import annotations

import argparse
import csv
import itertools
import pathlib
import statistics
from collections import Counter


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(row: dict[str, str], *fields: str) -> str:
    for field in fields:
        value = row.get(field, "")
        if value not in {"", "NA", "nan"}:
            return value
    return ""


def representative_cif(per_model: list[dict[str, str]], system: str, target: float) -> str:
    rows = [
        row for row in per_model
        if (row.get("system") or row.get("candidate")) == system and row.get("cif_path")
    ]
    if not rows:
        return ""
    medoids = [row for row in rows if row.get("main_topology_medoid") == "1"]
    if medoids:
        rows = medoids
    return min(rows, key=lambda row: (
        abs(float(row["iptm"]) - target), -float(row["iptm"]), row.get("model", "")
    ))["cif_path"]


def design_interface_coverage(per_model: list[dict[str, str]], system: str) -> dict[str, object]:
    rows = [
        row for row in per_model
        if (row.get("system") or row.get("candidate")) == system
    ]
    if not rows:
        return {}
    counts: Counter[int] = Counter()
    fractions = []
    for row in rows:
        value = row.get("design_interface_fraction", "")
        if value not in {"", "NA", "nan"}:
            fractions.append(float(value))
        counts.update(
            int(position) for position in row.get("design_interface_positions", "").split(",")
            if position
        )
    union = sorted(counts)
    consensus = sorted(position for position, count in counts.items() if count == len(rows))
    return {
        "models_two_interfaces_design_ge3": sum(
            int(row.get("two_interfaces_each_design_ge3", "0")) for row in rows
        ),
        "mean_design_interface_recurrence": (
            statistics.fmean(fractions) if fractions else ""
        ),
        "min_design_interface_recurrence": min(fractions) if fractions else "",
        "union_design_interface_position_count": len(union),
        "union_design_interface_positions": ",".join(map(str, union)),
        "all_model_design_interface_position_count": len(consensus),
        "all_model_design_interface_positions": ",".join(map(str, consensus)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=pathlib.Path, default=pathlib.Path.cwd())
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    candidates: list[dict[str, object]] = []
    screened_sequences: dict[str, dict[str, str]] = {}
    screen_paths = itertools.chain(
        args.project.glob("experiments/e*/unified*/final20_failclosed_screen.tsv"),
        args.project.glob(
            "experiments/current_pending_v1_screen/*/unified*/final20_failclosed_screen.tsv"
        ),
    )
    def screen_priority(path: pathlib.Path) -> tuple[int, float, str]:
        name = str(path)
        # Prefer the explicitly regenerated v1.1 no-OpenDDE table over legacy
        # v1.0 tables, then prefer a finalized table over an initial partial one.
        if "unified_v11_no_opendde" in name:
            rank = 0
        elif "unified_final" in name:
            rank = 1
        elif "unified_current_after_opendde" in name:
            rank = 2
        elif "unified_current" in name:
            rank = 3
        elif "unified_initial" in name:
            rank = 4
        else:
            rank = 5
        return (rank, -path.stat().st_mtime, name)

    for path in screen_paths:
        for row in read_tsv(path):
            sequence = row.get("sequence", "")
            if len(sequence) != 386:
                continue
            current = screened_sequences.get(sequence)
            current_path = pathlib.Path(current["screen_table"]) if current else None
            if current is None or screen_priority(path) < screen_priority(current_path):
                screened_sequences[sequence] = {**row, "screen_table": str(path)}

    historical_per_model: list[dict[str, str]] = []
    for path in args.project.glob("experiments/e*/final25_screen/af3_per_model_recomputed.tsv"):
        historical_per_model.extend(read_tsv(path))
    for path in args.project.glob("experiments/e*/full25_screen/af3_per_model_recomputed.tsv"):
        historical_per_model.extend(read_tsv(path))

    # New generic promotion outputs.
    for summary_path in args.project.glob("experiments/e*/full25/full25_summary.tsv"):
        experiment = summary_path.parents[1]
        per_model = read_tsv(summary_path.with_name("full25_per_model.tsv"))
        for row in read_tsv(summary_path):
            if int(row.get("accepted_af3", "0")) != 1:
                continue
            mean = float(row["mean_iptm"])
            candidates.append({
                "candidate": row["system"],
                "source_experiment": str(experiment),
                "source_kind": "generic_25model_promotion",
                "n_mutations": row["n_mutations"],
                "surface_mutation_fraction": row["surface_mutation_fraction"],
                "n_models": row["n_models"],
                "mean_iptm": row["mean_iptm"],
                "median_iptm": row["median_iptm"],
                "min_iptm": row["min_iptm"],
                "models_two_interfaces_design_ge3": row.get(
                    "models_two_interfaces_design_ge3", ""
                ),
                "mean_design_interface_recurrence": row.get(
                    "mean_design_interface_recurrence", ""
                ),
                "min_design_interface_recurrence": row.get(
                    "min_design_interface_recurrence", ""
                ),
                "union_design_interface_position_count": row.get(
                    "union_design_interface_position_count", ""
                ),
                "union_design_interface_positions": row.get(
                    "union_design_interface_positions", ""
                ),
                "all_model_design_interface_position_count": row.get(
                    "all_model_design_interface_position_count", ""
                ),
                "all_model_design_interface_positions": row.get(
                    "all_model_design_interface_positions", ""
                ),
                "representative_af3_cif": representative_cif(per_model, row["system"], mean),
                "sequence": row["sequence"],
            })

    # Historical Accepted_AF3 entries already present in a v1.0 unified table.
    for sequence, row in screened_sequences.items():
        if int(row.get("af3_accepted_af3_recomputed", "0")) != 1:
            continue
        coverage = design_interface_coverage(historical_per_model, row["candidate"])
        candidates.append({
            "candidate": row["candidate"],
            "source_experiment": row["screen_table"].split("/unified", 1)[0],
            "source_kind": "historical_v1_screen",
            "n_mutations": row.get("n_mutations", ""),
            "surface_mutation_fraction": row.get("surface_mutation_fraction_rsasa_ge0p20", ""),
            "n_models": row.get("af3_n_models", ""),
            "mean_iptm": row.get("af3_mean_iptm", ""),
            "median_iptm": row.get("af3_median_iptm", ""),
            "min_iptm": row.get("af3_min_iptm", ""),
            "models_two_interfaces_design_ge3": coverage.get(
                "models_two_interfaces_design_ge3",
                row.get("af3_two_interfaces_design_ge3_models", ""),
            ),
            "mean_design_interface_recurrence": coverage.get(
                "mean_design_interface_recurrence",
                row.get("af3_mean_design_interface_recurrence", ""),
            ),
            "min_design_interface_recurrence": coverage.get(
                "min_design_interface_recurrence",
                row.get("af3_min_design_interface_recurrence", ""),
            ),
            "union_design_interface_position_count": coverage.get(
                "union_design_interface_position_count", ""
            ),
            "union_design_interface_positions": coverage.get(
                "union_design_interface_positions", ""
            ),
            "all_model_design_interface_position_count": coverage.get(
                "all_model_design_interface_position_count", ""
            ),
            "all_model_design_interface_positions": coverage.get(
                "all_model_design_interface_positions", ""
            ),
            "representative_af3_cif": representative_cif(
                historical_per_model, row["candidate"], float(row["af3_mean_iptm"])
            ),
            "sequence": sequence,
        })

    by_sequence: dict[str, dict[str, object]] = {}
    for row in candidates:
        sequence = str(row["sequence"])
        previous = by_sequence.get(sequence)
        if previous is None or float(row["mean_iptm"]) > float(previous["mean_iptm"]):
            by_sequence[sequence] = row

    output: list[dict[str, object]] = []
    for row in by_sequence.values():
        screen = screened_sequences.get(str(row["sequence"]))
        row["v1_screen_available"] = int(screen is not None)
        # OpenDDE is retained as an optional cross-model diagnostic and is not
        # required for v1 completion.  ProLIF is a required chemistry stage,
        # but its legacy table uses ``prolif_pass`` rather than a status field.
        stage_available = {
            "af3": bool(screen and screen.get("af3_status", "MISSING") != "MISSING"),
            "prolif": bool(screen and screen.get("prolif_pass", "") not in {"", "MISSING"}),
            "protenix": bool(screen and screen.get("protenix_status", "MISSING") != "MISSING"),
            "stoichiometry": bool(screen and screen.get("stoichiometry_status", "MISSING") != "MISSING"),
            "rosetta": bool(screen and screen.get("rosetta_status", "MISSING") != "MISSING"),
        }
        screen_complete = bool(screen) and all(stage_available.values())
        row["v1_screen_complete"] = int(screen_complete)
        row["accepted_v1"] = int(screen is not None and screen.get("all_required_stages_pass") == "1")
        row["failed_required_stages"] = screen.get("failed_required_stages", "") if screen else ""
        row["prolif_design_nonvdw_count"] = (
            screen.get("prolif_design_nonvdw_count", "") if screen else ""
        )
        row["prolif_design_nonvdw_fraction"] = (
            screen.get("prolif_design_nonvdw_fraction", "") if screen else ""
        )
        row["v1_screen_table"] = screen.get("screen_table", "") if screen else ""
        row["needs_v1_screen"] = int(not screen_complete)
        output.append(row)
    output.sort(key=lambda row: (
        -int(row["accepted_v1"]), -float(row["mean_iptm"]), int(row["n_mutations"])
    ))
    for rank, row in enumerate(output, 1):
        row["rank"] = rank
    fields = ["rank"] + [field for field in output[0] if field != "rank"]
    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "accepted_af3_registry.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    with (args.out / "accepted_af3_sequences.fasta").open("w") as handle:
        for row in output:
            handle.write(
                f">RANK{int(row['rank']):02d}_{row['candidate']} "
                f"nmut={row['n_mutations']} mean_iptm={float(row['mean_iptm']):.4f} "
                f"accepted_v1={row['accepted_v1']}\n{row['sequence']}\n"
            )
    print(
        f"Accepted_AF3={len(output)} Accepted_v1={sum(int(row['accepted_v1']) for row in output)} "
        f"needs_v1={sum(int(row['needs_v1_screen']) for row in output)}"
    )


if __name__ == "__main__":
    main()
