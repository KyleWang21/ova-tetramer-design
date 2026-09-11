#!/usr/bin/env python3
"""Summarize Protenix oligomer confidence, residue contacts, and disulfides."""

from __future__ import annotations

import argparse
import itertools
import json
import pathlib

import numpy as np
from Bio.PDB import MMCIFParser


def component_size(chains: list[str], edges: list[tuple[str, str]]) -> int:
    graph = {chain: set() for chain in chains}
    for left, right in edges:
        graph[left].add(right)
        graph[right].add(left)
    best = 0
    for start in chains:
        seen, stack = {start}, [start]
        while stack:
            for other in graph[stack.pop()]:
                if other not in seen:
                    seen.add(other)
                    stack.append(other)
        best = max(best, len(seen))
    return best


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=pathlib.Path, required=True)
    parser.add_argument("--cif", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--contact-cutoff", type=float, default=8.0)
    parser.add_argument("--edge-min-contacts", type=int, default=10)
    args = parser.parse_args()

    confidence = json.loads(args.summary.read_text())
    model = next(MMCIFParser(QUIET=True).get_structure(args.cif.stem, args.cif).get_models())
    chains = [str(chain.id) for chain in model]
    coords: dict[str, list[tuple[int, np.ndarray]]] = {}
    sulfur: dict[tuple[str, int], np.ndarray] = {}
    for chain in model:
        chain_id = str(chain.id)
        coords[chain_id] = []
        for residue in chain:
            if residue.id[0] != " " or "CA" not in residue:
                continue
            atom = residue["CB"] if "CB" in residue else residue["CA"]
            coords[chain_id].append((int(residue.id[1]), np.asarray(atom.coord)))
            if "SG" in residue:
                sulfur[(chain_id, int(residue.id[1]))] = np.asarray(residue["SG"].coord)

    contacts = {}
    edges = []
    for left, right in itertools.combinations(chains, 2):
        xyz_left = np.stack([xyz for _, xyz in coords[left]])
        xyz_right = np.stack([xyz for _, xyz in coords[right]])
        count = int((np.linalg.norm(xyz_left[:, None] - xyz_right[None, :], axis=-1) <= args.contact_cutoff).sum())
        contacts[f"{left}-{right}"] = count
        if count >= args.edge_min_contacts:
            edges.append((left, right))

    close_ss = []
    for (key_a, xyz_a), (key_b, xyz_b) in itertools.combinations(sorted(sulfur.items()), 2):
        distance = float(np.linalg.norm(xyz_a - xyz_b))
        if distance <= 2.6:
            close_ss.append({"left": key_a, "right": key_b, "distance_a": distance})
    interchain_ss = [bond for bond in close_ss if bond["left"][0] != bond["right"][0]]
    native_ss = [
        bond for bond in close_ss
        if bond["left"][0] == bond["right"][0]
        and {bond["left"][1], bond["right"][1]} == {74, 121}
    ]
    pair = confidence.get("chain_pair_iptm", [])
    incident = [max(row[:idx] + row[idx + 1 :]) for idx, row in enumerate(pair)] if pair else []
    output = {
        "plddt": confidence.get("plddt"),
        "ptm": confidence.get("ptm"),
        "iptm": confidence.get("iptm"),
        "has_clash": confidence.get("has_clash"),
        "chain_pair_iptm": pair,
        "min_incident_iptm": min(incident) if incident else None,
        "contacts_8a": contacts,
        "contact_edges": [f"{a}-{b}" for a, b in edges],
        "largest_component": component_size(chains, edges),
        "full_connected": component_size(chains, edges) == len(chains),
        "native_c74_c121_bonds": len(native_ss),
        "interchain_ss_bonds": len(interchain_ss),
        "close_ss": close_ss,
        "summary": str(args.summary),
        "cif": str(args.cif),
    }
    args.out.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
