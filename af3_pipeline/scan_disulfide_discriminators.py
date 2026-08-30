#!/usr/bin/env python3
"""Scan surface residue pairs that geometrically distinguish C4 rings from failed states."""

from __future__ import annotations

import argparse
import csv
import itertools
import pathlib
import statistics

import numpy as np
from Bio.PDB import MMCIFParser


def cycle_geometry(cb, chains: list[str], left: int, right: int) -> tuple[float, float]:  # noqa: ANN001
    anchor = chains[0]
    choices = []
    for order in itertools.permutations(chains[1:]):
        cycle = (anchor, *order)
        distances = [
            float(np.linalg.norm(cb[(cycle[index], left)] - cb[(cycle[(index + 1) % 4], right)]))
            for index in range(4)
        ]
        choices.append((max(distances), statistics.mean(distances)))
    return min(choices)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ring-table", type=pathlib.Path, action="append", required=True)
    ap.add_argument("--interface-table", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--recurrence-min", type=float, default=.7)
    ap.add_argument("--sasa-min", type=float, default=20.0)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    positions = []
    for row in csv.DictReader(args.interface_table.open(), delimiter="\t"):
        if (float(row["contact_chain_recurrence"]) >= args.recurrence_min
                and float(row["isolated_chain_sasa_a2"]) >= args.sasa_min
                and float(row["mean_plddt"]) >= 70
                and not int(row["excluded_functional_region"])):
            positions.append((int(row["position"]), row["wt_aa"]))

    model_rows = []
    for table in args.ring_table:
        model_rows.extend(csv.DictReader(table.open(), delimiter="\t"))
    selected = [row for row in model_rows if row["system"] in {"ova_body_c4_04", "ova_body_double_ring_reseed"}]
    parser = MMCIFParser(QUIET=True)
    models = []
    for row in selected:
        model = next(parser.get_structure(row["model"], row["cif_path"]).get_models())
        chains = [str(chain.id) for chain in model]
        cb = {}
        for chain in model:
            for residue in chain:
                if residue.id[0] == " " and "CA" in residue:
                    atom = residue["CB"] if "CB" in residue else residue["CA"]
                    cb[(str(chain.id), int(residue.id[1]))] = np.asarray(atom.coord)
        models.append((bool(int(row["full_ss_ring"])), chains, cb))

    output = []
    for (left, left_aa), (right, right_aa) in itertools.permutations(positions, 2):
        positive_max, negative_max, intra = [], [], []
        for is_ring, chains, cb in models:
            maximum, _ = cycle_geometry(cb, chains, left, right)
            (positive_max if is_ring else negative_max).append(maximum)
            intra.extend(float(np.linalg.norm(cb[(chain, left)] - cb[(chain, right)])) for chain in chains)
        success_median = statistics.median(positive_max)
        failure_median = statistics.median(negative_max)
        intra_median = statistics.median(intra)
        output.append({
            "left": left, "left_aa": left_aa, "right": right, "right_aa": right_aa,
            "success_median_max_cb_a": success_median,
            "success_all4_within_6_fraction": sum(value <= 6 for value in positive_max) / len(positive_max),
            "failure_median_max_cb_a": failure_median,
            "failure_all4_within_6_fraction": sum(value <= 6 for value in negative_max) / len(negative_max),
            "failure_minus_success_a": failure_median - success_median,
            "intrachain_median_cb_a": intra_median,
            "discriminator_score": failure_median - success_median - max(0, 8 - intra_median),
        })
    output.sort(key=lambda row: (
        -float(row["success_all4_within_6_fraction"]),
        float(row["failure_all4_within_6_fraction"]),
        -float(row["discriminator_score"]),
    ))
    fields = list(output[0])
    with (args.out / "disulfide_discriminator_scan.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t"); writer.writeheader(); writer.writerows(output)
    for row in output[:20]:
        print(f"{row['left_aa']}{row['left']}->{row['right_aa']}{row['right']}",
              f"success={float(row['success_median_max_cb_a']):.2f}",
              f"failure={float(row['failure_median_max_cb_a']):.2f}",
              f"intra={float(row['intrachain_median_cb_a']):.2f}")


if __name__ == "__main__":
    main()
