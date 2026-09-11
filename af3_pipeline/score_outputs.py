#!/usr/bin/env python3
"""Score AF3 oligomer outputs using confidence, contacts, and graph connectivity.

Run this script with the AlphaFold3 virtualenv so its mmCIF parser is available.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import pathlib
import statistics
from collections import defaultdict
from typing import Iterable

import numpy as np
from alphafold3.structure import from_mmcif


def median(values: Iterable[float]) -> float:
    vals = [float(x) for x in values]
    return statistics.median(vals) if vals else float("nan")


def representative_coords(structure) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Return one CB (or CA for Gly/missing CB) coordinate per residue."""
    ca: dict[tuple[str, int], np.ndarray] = {}
    cb: dict[tuple[str, int], np.ndarray] = {}
    for chain, residue, atom, x, y, z in zip(
        structure.chain_id,
        structure.res_id,
        structure.atom_name,
        structure.atom_x,
        structure.atom_y,
        structure.atom_z,
        strict=True,
    ):
        key = (str(chain), int(residue))
        coord = np.array([x, y, z], dtype=np.float32)
        if atom == "CA":
            ca[key] = coord
        elif atom == "CB":
            cb[key] = coord
    by_chain: dict[str, list[tuple[int, np.ndarray]]] = {}
    for (chain, residue), coord in ca.items():
        by_chain.setdefault(chain, []).append((residue, cb.get((chain, residue), coord)))
    output = {}
    for chain, items in by_chain.items():
        ordered = sorted(items)
        output[chain] = (
            np.array([residue for residue, _ in ordered], dtype=np.int32),
            np.stack([coord for _, coord in ordered]),
        )
    return output


def contact_details(
    a: np.ndarray, b: np.ndarray, cutoff: float = 8.0
) -> tuple[int, np.ndarray, np.ndarray, np.ndarray]:
    delta = a[:, None, :] - b[None, :, :]
    mask = np.einsum("ijk,ijk->ij", delta, delta) < cutoff * cutoff
    return int(np.count_nonzero(mask)), mask, np.any(mask, axis=1), np.any(mask, axis=0)


def largest_component(nodes: list[str], edges: list[tuple[str, str]]) -> int:
    adjacent = {node: set() for node in nodes}
    for a, b in edges:
        adjacent[a].add(b)
        adjacent[b].add(a)
    best = 0
    seen: set[str] = set()
    for start in nodes:
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        size = 0
        while stack:
            node = stack.pop()
            size += 1
            for other in adjacent[node]:
                if other not in seen:
                    seen.add(other)
                    stack.append(other)
        best = max(best, size)
    return best


def parse_ranges(spec: str, expected_chains: int) -> list[tuple[int, int]]:
    ranges = []
    for item in spec.split(";"):
        lo, hi = item.split("-", 1)
        ranges.append((int(lo), int(hi)))
    if len(ranges) != expected_chains:
        raise ValueError(f"expected {expected_chains} ranges, found {spec}")
    return ranges


def range_plddt(structure, chains: list[str], ranges: list[tuple[int, int]]) -> float:
    """Return the weakest chain's mean CA pLDDT within the given ranges."""
    chain_values: list[float] = []
    for chain, (lo, hi) in zip(chains, ranges, strict=True):
        if lo == 0 and hi == 0:
            continue
        values = [
            float(b)
            for cid, rid, atom, b in zip(
                structure.chain_id,
                structure.res_id,
                structure.atom_name,
                structure.atom_b_factor,
                strict=True,
            )
            if str(cid) == chain and atom == "CA" and lo <= int(rid) <= hi
        ]
        if not values:
            raise ValueError(f"no CA atoms for chain {chain}, range {lo}-{hi}")
        chain_values.append(statistics.mean(values))
    return min(chain_values) if chain_values else float("nan")


def score_model(
    summary_path: pathlib.Path,
    system: str,
    expected_chains: int,
    ova_ranges: list[tuple[int, int]],
    module_ranges: list[tuple[int, int]],
    binder_ranges: list[tuple[int, int]],
) -> dict[str, object]:
    summary = json.loads(summary_path.read_text())
    cif_path = pathlib.Path(str(summary_path).replace("_summary_confidences.json", "_model.cif"))
    structure = from_mmcif(cif_path.read_text())
    chains = list(structure.chains)
    if len(chains) != expected_chains:
        raise ValueError(f"{cif_path}: expected {expected_chains} chains, found {chains}")
    coords = representative_coords(structure)
    pair_iptm = summary.get("chain_pair_iptm") or []
    pair_pae = summary.get("chain_pair_pae_min") or []
    contacts: dict[str, int] = {}
    contact_edges: list[tuple[str, str]] = []
    contact_iptm: list[float] = []
    contact_pae: list[float] = []
    interface_residues: set[int] = set()
    module_edges: list[tuple[str, str]] = []
    module_contacts = 0
    ova_contacts = 0
    max_partner_iptm = {chain: 0.0 for chain in chains}
    for i, chain_a in enumerate(chains):
        for j in range(i + 1, len(chains)):
            chain_b = chains[j]
            residues_a, xyz_a = coords[chain_a]
            residues_b, xyz_b = coords[chain_b]
            count, pair_mask, mask_a, mask_b = contact_details(xyz_a, xyz_b)
            ova_a = (residues_a >= ova_ranges[i][0]) & (residues_a <= ova_ranges[i][1])
            ova_b = (residues_b >= ova_ranges[j][0]) & (residues_b <= ova_ranges[j][1])
            if ova_ranges[i] != (0, 0) and ova_ranges[j] != (0, 0):
                ova_contacts += int(np.count_nonzero(pair_mask[np.ix_(ova_a, ova_b)]))
            mod_a = (residues_a >= module_ranges[i][0]) & (residues_a <= module_ranges[i][1])
            mod_b = (residues_b >= module_ranges[j][0]) & (residues_b <= module_ranges[j][1])
            if module_ranges[i] != (0, 0) and module_ranges[j] != (0, 0):
                mod_count = int(np.count_nonzero(pair_mask[np.ix_(mod_a, mod_b)]))
                module_contacts += mod_count
                if mod_count >= 5:
                    module_edges.append((chain_a, chain_b))
            contacts[f"{chain_a}-{chain_b}"] = count
            iptm = float(pair_iptm[i][j]) if pair_iptm else 0.0
            pae = float(min(pair_pae[i][j], pair_pae[j][i])) if pair_pae else float("nan")
            max_partner_iptm[chain_a] = max(max_partner_iptm[chain_a], iptm)
            max_partner_iptm[chain_b] = max(max_partner_iptm[chain_b], iptm)
            if count >= 5:
                contact_edges.append((chain_a, chain_b))
                contact_iptm.append(iptm)
                contact_pae.append(pae)
                interface_residues.update(int(x) for x in residues_a[mask_a])
                interface_residues.update(int(x) for x in residues_b[mask_b])
    connected = largest_component(chains, contact_edges)
    return {
        "system": system,
        "model": summary_path.parent.name,
        "n_chains": len(chains),
        "iptm": float(summary.get("iptm") or 0.0),
        "ptm": float(summary.get("ptm") or 0.0),
        "min_ova_plddt": range_plddt(structure, chains, ova_ranges),
        "min_module_plddt": range_plddt(structure, chains, module_ranges),
        "min_binder_plddt": range_plddt(structure, chains, binder_ranges),
        "ranking_score": float(summary.get("ranking_score") or 0.0),
        "has_clash": float(summary.get("has_clash") or 0.0),
        "contact_edges": len(contact_edges),
        "largest_component": connected,
        "full_connected": int(len(chains) > 1 and connected == len(chains)),
        "module_full_connected": int(
            len(chains) > 1 and largest_component(chains, module_edges) == len(chains)
        ),
        "total_contacts": sum(contacts.values()),
        "ova_contacts": ova_contacts,
        "module_contacts": module_contacts,
        "min_incident_iptm": min(max_partner_iptm.values()) if len(chains) > 1 else 0.0,
        "contact_pair_iptm_median": median(contact_iptm),
        "contact_pair_pae_median": median(contact_pae),
        "contacts_json": json.dumps(contacts, sort_keys=True),
        "interface_residues_json": json.dumps(sorted(interface_residues)),
        "summary_path": str(summary_path),
        "cif_path": str(cif_path),
    }


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def aggregate(rows: list[dict[str, object]], manifest: dict[str, dict[str, str]]) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault(str(row["system"]), []).append(row)
    output: list[dict[str, object]] = []
    for system in manifest:
        models = grouped.get(system, [])
        if not models:
            continue
        expected = int(manifest[system]["expected_samples"])
        best = max(models, key=lambda x: (
            int(x["full_connected"]), -float(x["has_clash"]), float(x["iptm"]),
            float(x["min_incident_iptm"]), int(x["total_contacts"]),
        ))
        residue_sets = [set(json.loads(str(x["interface_residues_json"]))) for x in models]
        jaccards = []
        for i, first in enumerate(residue_sets):
            for second in residue_sets[i + 1:]:
                union = first | second
                jaccards.append(len(first & second) / len(union) if union else 0.0)
        output.append({
            "system": system,
            "kind": manifest[system]["kind"],
            "n_chains": int(manifest[system]["n_chains"]),
            "n_models": len(models),
            "expected_models": expected,
            "complete": int(len(models) == expected),
            "median_iptm": round(median(float(x["iptm"]) for x in models), 4),
            "max_iptm": round(max(float(x["iptm"]) for x in models), 4),
            "median_ptm": round(median(float(x["ptm"]) for x in models), 4),
            "median_min_ova_plddt": round(median(float(x["min_ova_plddt"]) for x in models), 2),
            "median_min_module_plddt": round(median(float(x["min_module_plddt"]) for x in models), 2),
            "median_min_binder_plddt": round(median(float(x["min_binder_plddt"]) for x in models), 2),
            "median_min_incident_iptm": round(median(float(x["min_incident_iptm"]) for x in models), 4),
            "median_contact_pair_iptm": round(median(float(x["contact_pair_iptm_median"]) for x in models), 4),
            "median_contact_pair_pae": round(median(float(x["contact_pair_pae_median"]) for x in models), 4),
            "median_total_contacts": round(median(int(x["total_contacts"]) for x in models), 1),
            "median_ova_contacts": round(median(int(x["ova_contacts"]) for x in models), 1),
            "median_module_contacts": round(median(int(x["module_contacts"]) for x in models), 1),
            "full_connected_fraction": round(sum(int(x["full_connected"]) for x in models) / len(models), 4),
            "module_connected_fraction": round(
                sum(int(x["module_full_connected"]) for x in models) / len(models), 4
            ),
            "unclashed_fraction": round(sum(float(x["has_clash"]) == 0.0 for x in models) / len(models), 4),
            "median_interface_jaccard": round(median(jaccards), 4),
            "best_model": best["model"],
            "best_iptm": best["iptm"],
            "best_min_incident_iptm": best["min_incident_iptm"],
            "best_cif": best["cif_path"],
        })
    return output


def svg_plot(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    width = 1180
    row_h = 92
    height = 105 + row_h * len(rows)
    left = 245
    scale = 700
    colors = {
        "median_iptm": "#2563eb",
        "median_min_incident_iptm": "#f59e0b",
        "full_connected_fraction": "#16a34a",
        "median_min_ova_plddt": "#7c3aed",
        "median_interface_jaccard": "#0891b2",
    }
    labels = {
        "median_iptm": "median ipTM",
        "median_min_incident_iptm": "weak-chain ipTM",
        "full_connected_fraction": "all-chain fraction",
        "median_min_ova_plddt": "OVA pLDDT / 100",
        "median_interface_jaccard": "interface consistency",
    }
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#172033}.label{font-size:16px}.small{font-size:13px}.title{font-size:22px;font-weight:700}</style>',
        '<text x="28" y="34" class="title">AF3 oligomer screen</text>',
    ]
    legend_x = left
    for key in colors:
        parts.append(f'<rect x="{legend_x}" y="18" width="16" height="16" fill="{colors[key]}" rx="3"/>')
        parts.append(f'<text x="{legend_x + 22}" y="31" class="small">{labels[key]}</text>')
        legend_x += 185
    for tick in range(6):
        x = left + scale * tick / 5
        parts.append(f'<line x1="{x}" y1="58" x2="{x}" y2="{height - 25}" stroke="#e5e7eb"/>')
        parts.append(f'<text x="{x - 7}" y="76" class="small">{tick / 5:.1f}</text>')
    for row_idx, row in enumerate(rows):
        y = 92 + row_idx * row_h
        parts.append(f'<text x="28" y="{y + 28}" class="label">{html.escape(str(row["system"]))}</text>')
        for bar_idx, key in enumerate(colors):
            value = float(row[key])
            if key == "median_min_ova_plddt":
                value /= 100.0
            if not math.isfinite(value):
                value = 0.0
            value = max(0.0, min(1.0, value))
            bar_y = y + 3 + bar_idx * 17
            parts.append(f'<rect x="{left}" y="{bar_y}" width="{scale * value:.1f}" height="12" fill="{colors[key]}" rx="3"/>')
            parts.append(f'<text x="{left + scale * value + 6:.1f}" y="{bar_y + 11}" class="small">{value:.2f}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", type=pathlib.Path, required=True)
    args = parser.parse_args()
    manifest_rows = list(csv.DictReader((args.experiment / "manifest.tsv").open(), delimiter="\t"))
    rows_by_system: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in manifest_rows:
        rows_by_system[row["name"]].append(row)
    manifest: dict[str, dict[str, str]] = {}
    for system, rows in rows_by_system.items():
        if len({row.get("sequence", "") for row in rows}) != 1:
            raise ValueError(f"{system}: duplicate manifest rows have different sequences")
        meta = dict(rows[0])
        expected = [int(row["expected_samples"]) for row in rows]
        # prepare_reseed_shards records the combined total on every row, while
        # other splitters may record a per-shard total. Support both conventions.
        meta["expected_samples"] = str(
            max(expected) if max(expected) >= 5 * len(rows) else sum(expected)
        )
        manifest[system] = meta
    scored: list[dict[str, object]] = []
    seen_paths: set[pathlib.Path] = set()
    for system, metas in rows_by_system.items():
        for meta in metas:
            # Use only manifest-assigned shards. This excludes outputs from a
            # stopped duplicate task while allowing intentional reseed shards
            # that repeat the same system name in multiple manifest rows.
            shard_pattern = f"out_s{meta['shard']}" if meta.get("shard", "") else "out_s*"
            # AF3 writes either directly below <shard>/<system>/seed-* or,
            # when the JSON name is retained as an output subdirectory, below
            # <shard>/<system>/<system>/seed-*.  Accept both layouts.
            candidates = []
            for pattern in (
                f"{shard_pattern}/{system}/seed-*/*_summary_confidences.json",
                f"{shard_pattern}/{system}/{system}/seed-*/*_summary_confidences.json",
            ):
                candidates.extend(args.experiment.glob(pattern))
            paths = [
                path for path in sorted(set(candidates))
                if path.parent.name.startswith("seed-") and path not in seen_paths
            ]
            seen_paths.update(paths)
            for path in paths:
                n_chains = int(meta["n_chains"])
                scored.append(score_model(
                    path,
                    system,
                    n_chains,
                    parse_ranges(meta["ova_ranges"], n_chains),
                    parse_ranges(meta["module_ranges"], n_chains),
                    parse_ranges(meta.get("binder_ranges", meta["module_ranges"]), n_chains),
                ))
    if not scored:
        raise SystemExit("no AF3 summary outputs found")
    write_tsv(args.experiment / "model_scores.tsv", scored)
    aggregated = aggregate(scored, manifest)
    write_tsv(args.experiment / "system_scores.tsv", aggregated)
    svg_plot(args.experiment / "oligomer_screen.svg", aggregated)
    print("system\tmodels\tmedian_ipTM\tweakest_chain\tconnected\tcontacts\tbest")
    for row in aggregated:
        print(
            f"{row['system']}\t{row['n_models']}/{row['expected_models']}\t{row['median_iptm']}\t"
            f"{row['median_min_incident_iptm']}\t{row['full_connected_fraction']}\t"
            f"{row['median_total_contacts']}\t{row['best_model']}"
        )


if __name__ == "__main__":
    main()
