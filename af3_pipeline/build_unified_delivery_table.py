#!/usr/bin/env python3
"""Build one candidate-per-row table joining AF3, Protenix-v2, and ProLIF results."""

from __future__ import annotations

import argparse
import csv
import pathlib
from collections import defaultdict


PAIRS = ["A-B", "A-C", "A-D", "B-C", "B-D", "C-D"]
PROLIF_COUNTS = [
    "unique_residue_pairs", "unique_interaction_types", "atom_level_events",
    "non_vdw_interactions", "Hydrophobic", "HydrogenBond", "SaltBridge",
    "PiStacking", "CationPi", "VdWContact",
]


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def model_sort(row: dict[str, str]) -> tuple[int, str]:
    return int(row["seed"]), row["model"]


def compact(value: float) -> str:
    return f"{value:.4f}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--final", type=pathlib.Path, required=True)
    parser.add_argument("--per-seed", type=pathlib.Path, required=True)
    parser.add_argument("--per-model", type=pathlib.Path, required=True)
    parser.add_argument("--mutation-recurrence", type=pathlib.Path, required=True)
    parser.add_argument("--protenix", type=pathlib.Path, required=True)
    parser.add_argument("--prolif-pairs", type=pathlib.Path, required=True)
    parser.add_argument("--prolif-hotspots", type=pathlib.Path, required=True)
    parser.add_argument("--crossmodel-consensus", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    final = read_tsv(args.final)
    seed_rows: dict[str, dict[int, dict[str, str]]] = defaultdict(dict)
    for row in read_tsv(args.per_seed):
        seed_rows[row["candidate"]][int(row["seed"])] = row
    model_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(args.per_model):
        model_rows[row["candidate"]].append(row)
    recurrence_rows: dict[str, dict[str, str]] = defaultdict(dict)
    for row in read_tsv(args.mutation_recurrence):
        recurrence_rows[row["candidate"]][row["mutation"]] = row["interface_recurrence"]
    protenix = {row["candidate"]: row for row in read_tsv(args.protenix)}
    pair_rows: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in read_tsv(args.prolif_pairs):
        pair_rows[row["source"]][row["chain_pair"]] = row
    hotspot_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(args.prolif_hotspots):
        hotspot_rows[row["candidate"]].append(row)
    crossmodel = read_tsv(args.crossmodel_consensus)

    output: list[dict[str, object]] = []
    for base in final:
        candidate = base["candidate"]
        row: dict[str, object] = {
            "af3_rank": base["rank"],
            "candidate": candidate,
            "source_system": base["source_system"],
            "n_mutations": base["n_mutations"],
            "mutations": base["mutations"],
            "af3_n_models": base["n_models"],
            "af3_mean_iptm": base["mean_iptm"],
            "af3_median_iptm": base["median_iptm"],
            "af3_min_iptm": base["min_iptm"],
            "af3_mean_weakest_chain": base["mean_weakest_chain"],
            "af3_best_model_iptm": base["best_model_iptm"],
            "af3_full_connected_models": base["full_connected_models"],
            "af3_unclashed_models": base["unclashed_models"],
            "af3_no_interchain_ss_models": base["no_interchain_ss_models"],
            "af3_chemically_clean_models": base["chemically_clean_models"],
            "af3_all_native_c74_c121_models": base["all_native_c74_c121_models"],
            "siinfecl_intact": base["siinfecl_intact"],
            "af3_cif": f"cif/{int(base['rank']):02d}_{candidate}.cif",
        }
        for seed in range(1, 6):
            seed_row = seed_rows[candidate][seed]
            row[f"af3_seed{seed}_mean_iptm"] = seed_row["mean_iptm"]
            row[f"af3_seed{seed}_min_iptm"] = seed_row["min_iptm"]
            row[f"af3_seed{seed}_max_iptm"] = seed_row["max_iptm"]

        models = sorted(model_rows[candidate], key=model_sort)
        if len(models) != 25:
            raise ValueError(f"{candidate}: expected 25 AF3 models, found {len(models)}")
        row["af3_iptm_25_models"] = ",".join(model["iptm"] for model in models)
        row["af3_weakest_interface_25_models"] = ",".join(
            model["min_incident_iptm"] for model in models
        )
        row["af3_interface_pae_25_models"] = ",".join(
            model["contact_pair_pae_median"] for model in models
        )
        mutations = base["mutations"].split(",")
        mutation_recurrences = {
            mutation: float(recurrence_rows[candidate][mutation])
            for mutation in mutations
        }
        row["mutation_interface_recurrence"] = ";".join(
            f"{mutation}:{mutation_recurrences[mutation]:.2f}"
            for mutation in mutations
        )
        robust_interface_mutations = [
            mutation for mutation in mutations
            if mutation_recurrences[mutation] >= 0.80
        ]
        row["design_interface_match_score_25models"] = compact(
            sum(mutation_recurrences.values()) / len(mutations)
        )
        row["design_mutations_interface_ge80_count_25models"] = len(
            robust_interface_mutations
        )
        row["design_mutations_interface_ge80_fraction_25models"] = compact(
            len(robust_interface_mutations) / len(mutations)
        )
        row["design_mutations_interface_ge80_list_25models"] = ",".join(
            robust_interface_mutations
        )

        ptx = protenix[candidate]
        row.update({
            "protenix_rank": ptx["protenix_rank"],
            "protenix_seed": ptx["protenix_seed"],
            "protenix_iptm": ptx["protenix_iptm"],
            "protenix_ptm": ptx["protenix_ptm"],
            "protenix_plddt": ptx["protenix_plddt"],
            "protenix_min_incident_iptm": ptx["protenix_min_incident_iptm"],
            "protenix_largest_component": ptx["largest_component"],
            "protenix_full_connected": ptx["full_connected"],
            "protenix_has_clash": ptx["has_clash"],
            "protenix_native_c74_c121_bonds": ptx["native_c74_c121_bonds"],
            "protenix_interchain_ss_bonds": ptx["interchain_ss_bonds"],
            "protenix_contact_edges": ptx["contact_edges"],
        })

        pairs = pair_rows[candidate]
        for count_name in PROLIF_COUNTS:
            row[f"prolif_total_{count_name}"] = sum(
                int(pairs[pair][count_name]) for pair in PAIRS
            )
        for pair in PAIRS:
            compact_pair = pair.replace("-", "")
            row[f"prolif_{compact_pair}_unique_residue_pairs"] = pairs[pair]["unique_residue_pairs"]
            row[f"prolif_{compact_pair}_non_vdw"] = pairs[pair]["non_vdw_interactions"]
        row["prolif_dominant_chain_pair"] = max(
            PAIRS, key=lambda pair: int(pairs[pair]["non_vdw_interactions"])
        )

        residue_summary: dict[int, dict[str, object]] = {}
        mutation_summary: dict[str, dict[str, int]] = {}
        for hotspot in hotspot_rows[candidate]:
            position = int(hotspot["position"])
            residue = residue_summary.setdefault(position, {"types": 0, "nonvdw": 0})
            residue["types"] = int(residue["types"]) + int(hotspot["interaction_types"])
            residue["nonvdw"] = int(residue["nonvdw"]) + sum(
                int(hotspot[name]) for name in [
                    "Hydrophobic", "HydrogenBond", "SaltBridge", "PiStacking", "CationPi"
                ]
            )
            mutation = hotspot["mutation"]
            if mutation:
                item = mutation_summary.setdefault(mutation, {"types": 0, "nonvdw": 0})
                item["types"] += int(hotspot["interaction_types"])
                item["nonvdw"] += sum(
                    int(hotspot[name]) for name in [
                        "Hydrophobic", "HydrogenBond", "SaltBridge", "PiStacking", "CationPi"
                    ]
                )
        top_positions = sorted(
            residue_summary.items(), key=lambda item: (int(item[1]["types"]), int(item[1]["nonvdw"])),
            reverse=True,
        )[:10]
        row["prolif_top_residue_hotspots"] = ";".join(
            f"{position}:types={values['types']},nonvdw={values['nonvdw']}"
            for position, values in top_positions
        )
        top_mutations = sorted(
            mutation_summary.items(), key=lambda item: (item[1]["types"], item[1]["nonvdw"]),
            reverse=True,
        )
        row["prolif_mutation_hotspots"] = ";".join(
            f"{mutation}:types={values['types']},nonvdw={values['nonvdw']}"
            for mutation, values in top_mutations
        )
        designed_any_contact = [
            mutation for mutation in mutations if mutation in mutation_summary
        ]
        designed_nonvdw_contact = [
            mutation for mutation in mutations
            if mutation in mutation_summary and mutation_summary[mutation]["nonvdw"] > 0
        ]
        row["best_af3_prolif_design_any_contact_count"] = len(designed_any_contact)
        row["best_af3_prolif_design_any_contact_fraction"] = compact(
            len(designed_any_contact) / len(mutations)
        )
        row["best_af3_prolif_design_any_contact_list"] = ",".join(designed_any_contact)
        row["best_af3_prolif_design_nonvdw_count"] = len(designed_nonvdw_contact)
        row["best_af3_prolif_design_nonvdw_fraction"] = compact(
            len(designed_nonvdw_contact) / len(mutations)
        )
        row["best_af3_prolif_design_nonvdw_list"] = ",".join(designed_nonvdw_contact)

        designed_positions = {int(mutation[1:-1]) for mutation in mutations}
        interface_positions = set(residue_summary)
        designed_interface_positions = designed_positions & interface_positions
        row["best_af3_prolif_interface_position_count"] = len(interface_positions)
        row["best_af3_prolif_interface_designed_position_count"] = len(
            designed_interface_positions
        )
        row["best_af3_prolif_interface_designed_position_fraction"] = compact(
            len(designed_interface_positions) / len(interface_positions)
        )

        if candidate == "OVA-C4-IPTM80-20A":
            row["prolif_af3_protenix_consensus_count"] = len(crossmodel)
            row["prolif_af3_protenix_consensus"] = ";".join(
                f"{item['residue_1']}-{item['residue_2']}:{item['category']}"
                for item in crossmodel
            )
        else:
            row["prolif_af3_protenix_consensus_count"] = ""
            row["prolif_af3_protenix_consensus"] = ""
        row["sequence"] = base["sequence"]
        output.append(row)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = list(output[0])
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader(); writer.writerows(output)
    print(f"wrote {len(output)} candidates and {len(fields)} columns to {args.out}")


if __name__ == "__main__":
    main()
