#!/usr/bin/env python3
"""Apply C1/C2/C3/C4/C5/C6 specificity gates to the final OVA C4 set."""

from __future__ import annotations

import argparse
import csv
import pathlib
import re
import statistics
from collections import defaultdict

import numpy as np

from screen_final20_af3_structures import (
    ca_coords, kabsch_rmsd, orientless_pair_similarity, reference_ca,
    residue_pair_contacts, structure_arrays, within_model_interface_types,
)


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def pair_maps(path: pathlib.Path):
    _, chains = structure_arrays(path)
    ids = sorted(chains)
    maps = {}
    for left_index, left in enumerate(ids):
        for right in ids[left_index + 1:]:
            contacts, _ = residue_pair_contacts(chains[left], chains[right], 5.3)
            maps[(left, right)] = contacts
    effective = [pair for pair, contacts in maps.items() if len(contacts) >= 20]
    return maps, effective


def monomer_rmsd(path: pathlib.Path, reference: dict[int, np.ndarray]) -> float:
    _, chains = structure_arrays(path)
    if len(chains) != 1:
        raise ValueError(f"{path}: expected one chain, found {sorted(chains)}")
    residues, coords = ca_coords(next(iter(chains.values())))
    mobile = {int(residue): coord for residue, coord in zip(residues, coords)}
    common = sorted(set(reference) & set(mobile))
    if len(common) < 350:
        raise ValueError(f"{path}: only {len(common)} common CA positions")
    return kabsch_rmsd(
        np.stack([reference[position] for position in common]),
        np.stack([mobile[position] for position in common]),
    )


def canonical_candidate_key(name: str) -> str:
    """Normalize legacy system labels only for an otherwise-unmatched join."""
    name = re.sub(r"_ms$", "", name.lower())
    return re.sub(r"\d+", lambda match: str(int(match.group())), name)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--stoich-models", type=pathlib.Path, required=True)
    parser.add_argument("--c4-models", type=pathlib.Path, required=True)
    parser.add_argument("--reference", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    manifest = {row["name"]: row for row in read_tsv(args.manifest)}
    grouped = defaultdict(lambda: defaultdict(list))
    for row in read_tsv(args.stoich_models):
        meta = manifest[row["system"]]
        candidate = meta.get("candidate", "")
        if not candidate:
            # Single-candidate stoichiometry manifests prepared by the generic
            # AF3 helper encode the parent as `<candidate>_n{1,2,3,5,6}`.
            # Keep the explicit column authoritative when present, and only
            # use this exact suffix convention as a compatibility fallback.
            candidate = re.sub(r"_n(?:1|2|3|5|6)$", "", meta["name"])
        grouped[candidate][int(meta["n_chains"])].append(row)
    c4 = defaultdict(list)
    for row in read_tsv(args.c4_models):
        c4[row["candidate"]].append(row)
    ref_ca = reference_ca(args.reference, "A")
    summaries = []
    for candidate, c4_rows in c4.items():
        states = grouped.get(candidate)
        if states is None:
            matches = [
                name for name in grouped
                if canonical_candidate_key(name) == canonical_candidate_key(candidate)
            ]
            if len(matches) > 1:
                raise ValueError(f"ambiguous stoichiometry join for {candidate}: {matches}")
            states = grouped[matches[0]] if matches else {}
        state_metrics = {}
        for n in (1, 2, 3, 5, 6):
            rows = states.get(n, [])
            state_metrics[n] = {
                "n": len(rows),
                "mean_iptm": statistics.mean(float(row["iptm"]) for row in rows) if rows else None,
                "mean_ptm": statistics.mean(float(row["ptm"]) for row in rows) if rows else None,
                "mean_plddt": statistics.mean(float(row["min_ova_plddt"]) for row in rows) if rows else None,
                "high_connected_fraction": (
                    sum(float(row["iptm"]) >= 0.75 and int(row["full_connected"]) for row in rows) / len(rows)
                    if rows else None
                ),
                "high_any_fraction": (
                    sum(float(row["iptm"]) >= 0.75 for row in rows) / len(rows) if rows else None
                ),
            }
        c4_mean = statistics.mean(float(row["iptm"]) for row in c4_rows)
        c1_rmsds = [monomer_rmsd(pathlib.Path(row["cif_path"]), ref_ca)
                    for row in states.get(1, [])]
        complete = all(state_metrics[n]["n"] == 15 for n in (1, 2, 3, 5, 6))
        negative_means = [state_metrics[n]["mean_iptm"] for n in (3, 5, 6)]
        specificity = c4_mean - max(negative_means) if complete else None

        # C2-vs-C4 contact-map negative design check. A true second interface is
        # a repeated C4 interface type that is dissimilar to every C2 contact map.
        c2_contact_sets = []
        for row in states.get(2, []):
            maps, effective = pair_maps(pathlib.Path(row["cif_path"]))
            c2_contact_sets.extend(maps[pair] for pair in effective)
        second_unique_models = 0
        for row in c4_rows:
            maps, effective = pair_maps(pathlib.Path(row["cif_path"]))
            types = within_model_interface_types(maps, effective)
            type_novelty = []
            for group in types:
                similarities = []
                for pair in group:
                    similarities.extend(
                        orientless_pair_similarity(maps[pair], c2_map) for c2_map in c2_contact_sets
                    )
                # If C2 has no effective contact interface, every bona fide C4
                # interface type is absent from C2 and is therefore maximally novel.
                type_novelty.append(max(similarities) if similarities else 0.0)
            if len(types) >= 2 and min(type_novelty) < 0.50:
                second_unique_models += 1

        row = {
            "candidate": candidate,
            "c4_n_models": len(c4_rows), "c4_mean_iptm": c4_mean,
            "c1_n_models": state_metrics[1]["n"], "c1_mean_ptm": state_metrics[1]["mean_ptm"],
            "c1_mean_plddt": state_metrics[1]["mean_plddt"],
            "c1_mean_ca_rmsd_to_1ova": statistics.mean(c1_rmsds) if c1_rmsds else None,
            "c1_max_ca_rmsd_to_1ova": max(c1_rmsds) if c1_rmsds else None,
            "c2_n_models": state_metrics[2]["n"], "c2_mean_iptm": state_metrics[2]["mean_iptm"],
            "c3_n_models": state_metrics[3]["n"], "c3_mean_iptm": state_metrics[3]["mean_iptm"],
            "c5_n_models": state_metrics[5]["n"], "c5_mean_iptm": state_metrics[5]["mean_iptm"],
            "c6_n_models": state_metrics[6]["n"], "c6_mean_iptm": state_metrics[6]["mean_iptm"],
            "c3_high_connected_fraction": state_metrics[3]["high_connected_fraction"],
            "c5_high_connected_fraction": state_metrics[5]["high_connected_fraction"],
            "c6_high_connected_fraction": state_metrics[6]["high_connected_fraction"],
            "c3_high_any_fraction": state_metrics[3]["high_any_fraction"],
            "c5_high_any_fraction": state_metrics[5]["high_any_fraction"],
            "c6_high_any_fraction": state_metrics[6]["high_any_fraction"],
            "c4_specificity_margin": specificity,
            "c4_second_interface_unique_vs_c2_models": second_unique_models,
        }
        gates = {
            "gate_stoich_complete_15_models_each": complete and len(c4_rows) == 25,
            "gate_stoich_specificity_ge0p15": specificity is not None and specificity >= 0.15,
            "gate_stoich_c3_c5_c6_mean_below0p65": complete and max(negative_means) < 0.65,
            "gate_stoich_c3_c5_c6_high_connected_le0p20": complete and max(
                state_metrics[n]["high_connected_fraction"] for n in (3, 5, 6)
            ) <= 0.20,
            "gate_stoich_no_disconnected_high_conf_competitor": complete and max(
                state_metrics[n]["high_any_fraction"] for n in (3, 5, 6)
            ) <= 0.20,
            "gate_stoich_c1_ptm_ge0p88": complete and state_metrics[1]["mean_ptm"] >= 0.88,
            "gate_stoich_c1_plddt_ge85": complete and state_metrics[1]["mean_plddt"] >= 85.0,
            "gate_stoich_c1_max_rmsd_le2": complete and bool(c1_rmsds) and max(c1_rmsds) <= 2.0,
            "gate_stoich_second_interface_unique_ge20of25": second_unique_models >= 20,
        }
        row.update({key: int(value) for key, value in gates.items()})
        row["stoichiometry_pass"] = int(all(gates.values()))
        row["failed_stoichiometry_gates"] = ";".join(key for key, value in gates.items() if not value)
        summaries.append(row)
    with (args.out / "stoichiometry_candidate_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(summaries[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(summaries)
    print(f"wrote {len(summaries)} candidates; passed={sum(row['stoichiometry_pass'] for row in summaries)}")


if __name__ == "__main__":
    main()
