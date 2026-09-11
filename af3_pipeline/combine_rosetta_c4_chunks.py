#!/usr/bin/env python3
"""Combine independently seeded Rosetta C4 FastRelax chunks fail-closed.

Each chunk is produced by ``rosetta_c4_relax_screen.py``.  This combiner
requires the exact expected seed set, recomputes the production gates over the
union, and copies the best relaxed structure selected with the same ordering as
the worker.  A partial or duplicated ensemble never receives a production
summary.
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import shutil
import statistics


def as_int(row: dict[str, str], field: str) -> int:
    return int(float(row[field]))


def as_float(row: dict[str, str], field: str) -> float:
    return float(row[field])


def selection_key(row: dict[str, str]) -> tuple[int, int, float, float]:
    return (
        as_int(row, "valid_c4_network"),
        -as_int(row, "atomic_clash_count_2p4"),
        as_float(row, "weakest_sc"),
        -as_float(row, "worst_dG_dSASA_x100"),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--chunks-root", type=pathlib.Path, required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--source-model", choices=("AF3", "Protenix"), required=True)
    parser.add_argument("--expected", type=int, default=200)
    parser.add_argument("--expected-chunks", type=int, default=8)
    parser.add_argument("--seed-start", type=int, default=7001)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    stem = f"{args.candidate}_{args.source_model}"
    replicate_name = f"{stem}_rosetta_replicates.tsv"
    chunk_tables = sorted(args.chunks_root.glob(f"chunk*/{replicate_name}"))
    if len(chunk_tables) != args.expected_chunks:
        raise RuntimeError(
            f"expected {args.expected_chunks} completed chunk tables, found {len(chunk_tables)}"
        )

    rows: list[dict[str, str]] = []
    rows_by_chunk: dict[pathlib.Path, list[dict[str, str]]] = {}
    fields: list[str] | None = None
    for table in chunk_tables:
        with table.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise RuntimeError(f"missing header in {table}")
            if fields is None:
                fields = list(reader.fieldnames)
            elif list(reader.fieldnames) != fields:
                raise RuntimeError(f"inconsistent columns in {table}")
            chunk_rows = list(reader)
        if not chunk_rows:
            raise RuntimeError(f"empty chunk table {table}")
        for row in chunk_rows:
            if row["candidate"] != args.candidate or row["source_model"] != args.source_model:
                raise RuntimeError(f"candidate/source mismatch in {table}")
        rows_by_chunk[table.parent] = chunk_rows
        rows.extend(chunk_rows)

    expected_seeds = list(range(args.seed_start, args.seed_start + args.expected))
    observed_seeds = sorted(as_int(row, "seed") for row in rows)
    if len(rows) != args.expected:
        raise RuntimeError(f"expected {args.expected} replicates, found {len(rows)}")
    if observed_seeds != expected_seeds:
        missing = sorted(set(expected_seeds) - set(observed_seeds))
        duplicate_count = len(observed_seeds) - len(set(observed_seeds))
        unexpected = sorted(set(observed_seeds) - set(expected_seeds))
        raise RuntimeError(
            f"seed-set mismatch: missing={missing} unexpected={unexpected} duplicates={duplicate_count}"
        )

    rows.sort(key=lambda row: as_int(row, "seed"))
    for replicate, row in enumerate(rows, 1):
        row["replicate"] = str(replicate)

    args.out.mkdir(parents=True, exist_ok=True)
    combined_table = args.out / replicate_name
    assert fields is not None
    with combined_table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    best_chunk = max(
        rows_by_chunk,
        key=lambda chunk: max(selection_key(row) for row in rows_by_chunk[chunk]),
    )
    best_pdb = best_chunk / f"{stem}_relaxed_best.pdb"
    clean_pdb = best_chunk / f"{stem}_clean.pdb"
    if not best_pdb.is_file() or not clean_pdb.is_file():
        raise RuntimeError(f"best chunk is missing a required PDB: {best_chunk}")
    shutil.copy2(best_pdb, args.out / best_pdb.name)
    shutil.copy2(clean_pdb, args.out / clean_pdb.name)

    summary: dict[str, object] = {
        "candidate": args.candidate,
        "source_model": args.source_model,
        "nrelax": len(rows),
        "production_ensemble_complete": 1,
        "parallel_chunk_count": len(chunk_tables),
        "seed_start": args.seed_start,
        "seed_end": args.seed_start + args.expected - 1,
        "design_interface_edges": rows[0].get("design_interface_edges", ""),
        "background_edges": rows[0].get("background_edges", ""),
        "topology_retained_fraction": statistics.mean(
            as_int(row, "valid_c4_network") for row in rows
        ),
        "zero_clash_fraction": statistics.mean(
            as_int(row, "atomic_clash_count_2p4") == 0 for row in rows
        ),
        "zero_overlap_fraction": statistics.mean(
            as_int(row, "severe_ca_overlap_count_2p5") == 0 for row in rows
        ),
        "median_weakest_sc": statistics.median(as_float(row, "weakest_sc") for row in rows),
        "median_weakest_hbonds": statistics.median(
            as_float(row, "weakest_hbonds") for row in rows
        ),
        "median_worst_dG_dSASA_x100": statistics.median(
            as_float(row, "worst_dG_dSASA_x100") for row in rows
        ),
        "median_worst_unsat_per_1000A2": statistics.median(
            as_float(row, "worst_unsat_per_1000A2") for row in rows
        ),
        "max_chain_ca_rmsd_to_input": max(
            as_float(row, "max_chain_ca_rmsd_to_input") for row in rows
        ),
    }
    if "all_weakest_sc" in rows[0]:
        summary.update({
            "median_all_weakest_sc": statistics.median(as_float(row, "all_weakest_sc") for row in rows),
            "median_all_weakest_hbonds": statistics.median(as_float(row, "all_weakest_hbonds") for row in rows),
            "median_all_worst_dG_dSASA_x100": statistics.median(as_float(row, "all_worst_dG_dSASA_x100") for row in rows),
            "median_all_worst_unsat_per_1000A2": statistics.median(as_float(row, "all_worst_unsat_per_1000A2") for row in rows),
        })
    gates = {
        "gate_rosetta_ensemble_100": len(rows) >= 100,
        "gate_rosetta_topology_ge0p90": summary["topology_retained_fraction"] >= 0.90,
        "gate_rosetta_zero_clash_fraction_ge0p90": summary["zero_clash_fraction"] >= 0.90,
        "gate_rosetta_zero_overlap_fraction_ge0p90": summary["zero_overlap_fraction"] >= 0.90,
        "gate_rosetta_weakest_sc_ge0p60": summary["median_weakest_sc"] >= 0.60,
        "gate_rosetta_weakest_hbonds_ge3": summary["median_weakest_hbonds"] >= 3.0,
        "gate_rosetta_dG_dSASA_le_minus1": summary["median_worst_dG_dSASA_x100"] <= -1.0,
        "gate_rosetta_unsat_per_1000A2_le2": summary["median_worst_unsat_per_1000A2"] <= 2.0,
        "gate_rosetta_chain_rmsd_le2": summary["max_chain_ca_rmsd_to_input"] <= 2.0,
    }
    summary.update({name: int(passed) for name, passed in gates.items()})
    summary["rosetta_pass"] = int(all(gates.values()))
    summary["failed_rosetta_gates"] = ";".join(
        name for name, passed in gates.items() if not passed
    )

    summary_json = args.out / f"{stem}_rosetta_summary.json"
    summary_tsv = args.out / f"{stem}_rosetta_summary.tsv"
    summary_json.write_text(json.dumps(summary, indent=2) + "\n")
    with summary_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(summary), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(summary)
    print(f"combined {len(rows)} Rosetta replicates from {len(chunk_tables)} chunks")
    print(f"rosetta_pass={summary['rosetta_pass']} failed={summary['failed_rosetta_gates']}")


if __name__ == "__main__":
    main()
