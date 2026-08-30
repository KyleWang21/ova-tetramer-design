#!/usr/bin/env python3
"""Find two self-symmetric Cys sites that encode an even-chain D2 tetramer.

Each site must support two disjoint Cys--Cys bonds in a C4 model.  A pair of
sites is retained only when their best perfect matchings differ, so the union
of the four bonds connects all four chains.  Three identical chains cannot
saturate either self-symmetric Cys site, providing explicit odd-stoichiometry
negative design.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import statistics
from collections import Counter
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser


MATCHINGS = (
    ((0, 1), (2, 3)),
    ((0, 2), (1, 3)),
    ((0, 3), (1, 2)),
)


def read_rows(path: Path, iptm_min: float, weak_min: float) -> list[dict[str, str]]:
    rows = []
    for row in csv.DictReader(path.open(), delimiter="\t"):
        weak = float(row.get("min_incident_iptm", row.get("weakest_chain", 0)))
        if float(row["iptm"]) >= iptm_min and weak >= weak_min and int(float(row.get("has_clash", 0))) == 0:
            rows.append(row)
    if not rows:
        raise SystemExit("no passing C4 rows")
    return rows


def best_matching(points: list[np.ndarray]) -> tuple[int, float, float]:
    choices = []
    for index, matching in enumerate(MATCHINGS):
        distances = [float(np.linalg.norm(points[left] - points[right])) for left, right in matching]
        choices.append((max(distances), statistics.mean(distances), index))
    maximum, mean, index = min(choices)
    return index, maximum, mean


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c4-scores", type=Path, required=True)
    parser.add_argument("--hotspots", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--iptm-min", type=float, default=0.50)
    parser.add_argument("--weak-min", type=float, default=0.45)
    parser.add_argument("--recurrence-min", type=float, default=0.50)
    parser.add_argument("--sasa-min", type=float, default=20.0)
    parser.add_argument("--distance-cutoff", type=float, default=7.0)
    parser.add_argument("--site-success-min", type=float, default=0.70)
    parser.add_argument("--pair-success-min", type=float, default=0.60)
    parser.add_argument("--exclude", default="14,69,74,77,79,94,96,121,155,190-193,196,213,257-264,289-295")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    excluded = set()
    for part in args.exclude.split(","):
        if "-" in part:
            lo, hi = (int(value) for value in part.split("-", 1))
            excluded.update(range(lo, hi + 1))
        else:
            excluded.add(int(part))
    positions = []
    amino_acids = {}
    for row in csv.DictReader(args.hotspots.open(), delimiter="\t"):
        position = int(row["position"])
        amino_acid = row.get("aa", row.get("wt_aa", "X"))
        amino_acids[position] = amino_acid
        if (
            position not in excluded
            and amino_acid not in {"C", "G", "P", "X"}
            and float(row["c4_contact_recurrence"]) >= args.recurrence_min
            and float(row["isolated_sasa_a2"]) >= args.sasa_min
            and float(row["c4_mean_plddt"]) >= 70
        ):
            positions.append(position)

    score_rows = read_rows(args.c4_scores, args.iptm_min, args.weak_min)
    parser_cif = MMCIFParser(QUIET=True)
    models = []
    for row in score_rows:
        model = next(parser_cif.get_structure(row.get("model", "model"), row["cif_path"]).get_models())
        chains = list(model.get_chains())
        if len(chains) != 4:
            continue
        coordinates = {}
        for chain in chains:
            for residue in chain:
                if residue.id[0] == " " and int(residue.id[1]) in positions and "CA" in residue:
                    atom = residue["CB"] if "CB" in residue else residue["CA"]
                    coordinates[(str(chain.id), int(residue.id[1]))] = np.asarray(atom.coord)
        models.append(([str(chain.id) for chain in chains], coordinates))

    per_site = {}
    site_rows = []
    for position in positions:
        values = []
        for chains, coordinates in models:
            points = [coordinates[(chain, position)] for chain in chains]
            values.append(best_matching(points))
        success = sum(maximum <= args.distance_cutoff for _, maximum, _ in values) / len(values)
        per_site[position] = values
        site_rows.append({
            "position": position,
            "aa": amino_acids[position],
            "success_fraction": success,
            "median_best_max_cb_a": statistics.median(maximum for _, maximum, _ in values),
            "median_best_mean_cb_a": statistics.median(mean for _, _, mean in values),
            "dominant_matching_fraction": Counter(index for index, _, _ in values).most_common(1)[0][1] / len(values),
        })
    site_rows.sort(key=lambda row: (-row["success_fraction"], row["median_best_max_cb_a"]))
    with (args.out / "self_cys_sites.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(site_rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(site_rows)

    pair_rows = []
    eligible = [int(row["position"]) for row in site_rows if row["success_fraction"] >= args.site_success_min]
    for first, second in itertools.combinations(eligible, 2):
        connected = []
        first_ok, second_ok = [], []
        for first_value, second_value in zip(per_site[first], per_site[second]):
            first_match, first_max, _ = first_value
            second_match, second_max, _ = second_value
            ok1 = first_max <= args.distance_cutoff
            ok2 = second_max <= args.distance_cutoff
            first_ok.append(ok1); second_ok.append(ok2)
            connected.append(ok1 and ok2 and first_match != second_match)
        intrachain = []
        for chains, coordinates in models:
            intrachain.extend(float(np.linalg.norm(coordinates[(chain, first)] - coordinates[(chain, second)])) for chain in chains)
        pair_rows.append({
            "site_1": first,
            "aa_1": amino_acids[first],
            "site_2": second,
            "aa_2": amino_acids[second],
            "connected_d2_fraction": sum(connected) / len(connected),
            "site_1_success_fraction": sum(first_ok) / len(first_ok),
            "site_2_success_fraction": sum(second_ok) / len(second_ok),
            "median_intrachain_cb_a": statistics.median(intrachain),
            "topology_score": sum(connected) / len(connected) - max(0.0, 8.0 - statistics.median(intrachain)) / 8.0,
        })
    pair_rows.sort(key=lambda row: (-row["topology_score"], -row["connected_d2_fraction"]))
    with (args.out / "d2_self_cys_pairs.tsv").open("w", newline="") as handle:
        if pair_rows:
            writer = csv.DictWriter(handle, list(pair_rows[0]), delimiter="\t")
            writer.writeheader(); writer.writerows(pair_rows)
    retained = [row for row in pair_rows if row["connected_d2_fraction"] >= args.pair_success_min]
    print(f"C4 models={len(models)} surface positions={len(positions)} eligible sites={len(eligible)} D2 pairs={len(retained)}")
    for row in retained[:30]:
        print(f"{row['aa_1']}{row['site_1']}C/{row['aa_2']}{row['site_2']}C",
              f"D2={row['connected_d2_fraction']:.2f}",
              f"site={row['site_1_success_fraction']:.2f}/{row['site_2_success_fraction']:.2f}",
              f"intra={row['median_intrachain_cb_a']:.1f}")


if __name__ == "__main__":
    main()
