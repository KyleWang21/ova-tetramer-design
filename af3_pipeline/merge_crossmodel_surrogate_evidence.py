#!/usr/bin/env python3
"""Merge v1.0 screen evidence into compact active-learning input tables."""

from __future__ import annotations

import argparse
import csv
import pathlib


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing empty output: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--unified", type=pathlib.Path, action="append", required=True)
    parser.add_argument("--prolif", type=pathlib.Path, action="append", required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    unified: dict[str, dict[str, str]] = {}
    for path in args.unified:
        for row in read_tsv(path):
            name = row["candidate"]
            if name in unified and unified[name].get("sequence") != row.get("sequence"):
                raise ValueError(f"conflicting sequences for {name}")
            unified[name] = row

    prolif: dict[str, dict[str, str]] = {}
    for path in args.prolif:
        for row in read_tsv(path):
            prolif[row["candidate"]] = row

    final_rows: list[dict[str, object]] = []
    geometry_rows: list[dict[str, object]] = []
    protenix_rows: list[dict[str, object]] = []
    opendde_rows: list[dict[str, object]] = []
    prolif_rows: list[dict[str, object]] = []
    for name, row in unified.items():
        sequence = row.get("sequence", "")
        if len(sequence) != 386:
            raise ValueError(f"{name}: expected 386-aa sequence, found {len(sequence)}")
        final_rows.append({"candidate": name, "sequence": sequence})
        geometry_rows.append({
            "candidate": name,
            "n_models": row.get("af3_n_models", "0"),
            "zero_atomic_clash_models": row.get("af3_zero_atomic_clash_models", "0"),
        })
        protenix_rows.append({
            "candidate": name,
            "mean_iptm": row.get("protenix_mean_iptm", ""),
        })
        opendde_rows.append({
            "candidate": name,
            "gate_opendde_10_seeds": row.get("opendde_gate_opendde_10_seeds", "0"),
            "mean_iptm": row.get("opendde_mean_iptm", ""),
        })
        interaction = prolif.get(name)
        if interaction is not None:
            prolif_rows.append({
                "candidate": name,
                "shared_nonvdw_position_type_count": interaction[
                    "shared_nonvdw_position_type_count"
                ],
            })

    args.out.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out / "final_candidates.tsv", final_rows)
    write_tsv(args.out / "geometry_summary.tsv", geometry_rows)
    write_tsv(args.out / "protenix_summary.tsv", protenix_rows)
    write_tsv(args.out / "opendde_summary.tsv", opendde_rows)
    write_tsv(args.out / "crossmodel_prolif.tsv", prolif_rows)
    print(
        f"merged candidates={len(final_rows)} crossmodel_prolif={len(prolif_rows)} "
        f"from {len(args.unified)} screen tables"
    )


if __name__ == "__main__":
    main()
