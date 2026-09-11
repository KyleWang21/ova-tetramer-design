#!/usr/bin/env python3
"""Run ProLIF protein-protein fingerprints on homotetramer chain interfaces."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import pathlib
from collections import Counter, defaultdict

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import prolif as plf
import seaborn as sns
from Bio.PDB import MMCIFParser, PDBIO
from Bio.Data.PDBData import protein_letters_3to1
from rdkit import Chem


INTERACTIONS = [
    "Hydrophobic", "HBDonor", "HBAcceptor", "Cationic", "Anionic",
    "PiStacking", "CationPi", "PiCation",
    "VdWContact",
]


def category(name: str) -> str:
    if name in {"HBDonor", "HBAcceptor"}:
        return "HydrogenBond"
    if name in {"Cationic", "Anionic"}:
        return "SaltBridge"
    if name in {"PiStacking", "FaceToFace", "EdgeToFace"}:
        return "PiStacking"
    if name in {"CationPi", "PiCation"}:
        return "CationPi"
    return name


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def to_pdb(cif: pathlib.Path, pdb: pathlib.Path) -> None:
    structure = MMCIFParser(QUIET=True).get_structure(cif.stem, str(cif))
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(pdb))


def protonated_chain_molecules(pdb: pathlib.Path) -> tuple[dict[str, plf.Molecule], dict[str, Chem.Mol]]:
    # Parse each chain separately. Parsing the complete assembly with RDKit's
    # proximity bonding can incorrectly create covalent bonds across a tight
    # noncovalent interface.
    pdb_lines = pdb.read_text().splitlines()
    chain_ids = sorted({line[21].strip() for line in pdb_lines
                        if line.startswith(("ATOM  ", "HETATM")) and line[21].strip()})
    rdkit_chains: dict[str, Chem.Mol] = {}
    for chain in chain_ids:
        chain_pdb = pdb.with_name(f"{pdb.stem}_chain_{chain}.pdb")
        selected = [line for line in pdb_lines
                    if line.startswith(("ATOM  ", "HETATM")) and line[21].strip() == chain]
        chain_pdb.write_text("\n".join(selected + ["TER", "END"]) + "\n")
        raw = Chem.MolFromPDBFile(str(chain_pdb), sanitize=False, removeHs=False, proximityBonding=False)
        if raw is None:
            raise ValueError(f"RDKit could not read coordinates from {chain_pdb}")
        residues: list[tuple[int, str]] = []
        seen_residues = set()
        coordinate_atoms: dict[tuple[int, str], int] = {}
        for atom in raw.GetAtoms():
            info = atom.GetPDBResidueInfo()
            if not info:
                continue
            key = (info.GetResidueNumber(), info.GetResidueName().strip())
            if key not in seen_residues:
                seen_residues.add(key)
                residues.append(key)
            coordinate_atoms[(info.GetResidueNumber(), info.GetName().strip())] = atom.GetIdx()
        sequence = "".join(protein_letters_3to1[name.upper()] for _, name in residues)
        mol = Chem.MolFromSequence(sequence)
        if mol is None:
            raise ValueError(f"RDKit could not build sequence topology for {chain_pdb}")
        conformer = Chem.Conformer(mol.GetNumAtoms())
        source_conformer = raw.GetConformer()
        for atom in mol.GetAtoms():
            info = atom.GetPDBResidueInfo()
            info.SetChainId(chain)
            key = (info.GetResidueNumber(), info.GetName().strip())
            if key not in coordinate_atoms:
                raise ValueError(f"{chain_pdb}: missing coordinate for residue atom {key}")
            conformer.SetAtomPosition(atom.GetIdx(), source_conformer.GetAtomPosition(coordinate_atoms[key]))
        mol.AddConformer(conformer, assignId=True)
        residue_numbers = {number for number, _ in residues}
        terminus = (min(residue_numbers), max(residue_numbers))
        rw = Chem.RWMol(mol)
        disulfide_atoms: dict[int, int] = {}
        for atom in rw.GetAtoms():
            info = atom.GetPDBResidueInfo()
            if not info:
                continue
            resname = info.GetResidueName().strip()
            atomname = info.GetName().strip()
            resid = info.GetResidueNumber()
            if resname == "CYS" and resid in {74, 121} and atomname == "SG":
                disulfide_atoms[resid] = atom.GetIdx()
                atom.SetNoImplicit(True)
                atom.SetNumExplicitHs(0)
            if (resname, atomname) in {("LYS", "NZ"), ("ARG", "NH2")}:
                atom.SetFormalCharge(1)
            elif (resname, atomname) in {("ASP", "OD2"), ("GLU", "OE2")}:
                atom.SetFormalCharge(-1)
                atom.SetNoImplicit(True)
                atom.SetNumExplicitHs(0)
            elif resid == terminus[0] and atomname == "N":
                atom.SetFormalCharge(1)
            elif resid == terminus[1] and atomname == "OXT":
                atom.SetFormalCharge(-1)
                atom.SetNoImplicit(True)
                atom.SetNumExplicitHs(0)
        if set(disulfide_atoms) == {74, 121}:
            rw.AddBond(disulfide_atoms[74], disulfide_atoms[121], Chem.BondType.SINGLE)
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
        rdkit_chains[chain] = Chem.AddHs(mol, addCoords=True, addResidueInfo=True)
    return {chain: plf.Molecule(chain_mol) for chain, chain_mol in rdkit_chains.items()}, rdkit_chains


def atom_names(mol: Chem.Mol, indices: tuple[int, ...] | list[int]) -> str:
    names = []
    for index in indices:
        atom = mol.GetAtomWithIdx(int(index))
        info = atom.GetPDBResidueInfo()
        names.append(info.GetName().strip() if info else atom.GetSymbol())
    return ",".join(names)


def jsonable(value):
    if hasattr(value, "item"):
        return value.item()
    if isinstance(value, tuple):
        return list(value)
    raise TypeError(type(value).__name__)


def parse_structures(values: list[str]) -> dict[str, pathlib.Path]:
    parsed = {}
    for value in values:
        if "=" not in value:
            raise ValueError("--structure must be LABEL=PATH")
        label, path = value.split("=", 1)
        parsed[label] = pathlib.Path(path)
    return parsed


def read_structure_manifest(path: pathlib.Path) -> dict[str, pathlib.Path]:
    rows = read_tsv(path)
    required = {"source", "cif_path"}
    if rows and not required.issubset(rows[0]):
        raise ValueError(f"{path}: expected columns {sorted(required)}")
    return {row["source"]: pathlib.Path(row["cif_path"]) for row in rows}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--structure", action="append", default=[], help="LABEL=CIF_PATH")
    parser.add_argument("--structure-manifest", type=pathlib.Path,
                        help="optional TSV with source and cif_path columns")
    parser.add_argument("--final-candidates", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    pdb_dir = args.out / "pdb"
    pdb_dir.mkdir(exist_ok=True)

    candidates = {row["candidate"]: row for row in read_tsv(args.final_candidates)}
    mutation_maps = {
        name: {int(mutation[1:-1]): mutation for mutation in row["mutations"].split(",")}
        for name, row in candidates.items()
    }
    structures = parse_structures(args.structure)
    if args.structure_manifest:
        structures.update(read_structure_manifest(args.structure_manifest))
    if not structures:
        parser.error("at least one --structure or --structure-manifest is required")
    fingerprint = plf.Fingerprint(INTERACTIONS, count=True, vicinity_cutoff=6.0)
    event_rows: list[dict[str, object]] = []

    for source, cif in structures.items():
        candidate = next((name for name in candidates if source.startswith(name)), "")
        pdb = pdb_dir / f"{source}.pdb"
        to_pdb(cif, pdb)
        chains, rdkit_chains = protonated_chain_molecules(pdb)
        if len(chains) != 4:
            raise ValueError(f"{source}: expected four chains, found {sorted(chains)}")
        for chain_a, chain_b in itertools.combinations(sorted(chains), 2):
            result = fingerprint.generate(chains[chain_a], chains[chain_b], metadata=True)
            for item in result.interactions():
                metadata = item.metadata
                parent = metadata.get("parent_indices", {})
                left_indices = tuple(parent.get("ligand", ()))
                right_indices = tuple(parent.get("protein", ()))
                event_rows.append({
                    "source": source,
                    "candidate": candidate,
                    "chain_pair": f"{chain_a}-{chain_b}",
                    "chain_a": chain_a,
                    "resname_a": item.ligand.name,
                    "position_a": item.ligand.number,
                    "residue_a": str(item.ligand),
                    "mutation_a": mutation_maps.get(candidate, {}).get(item.ligand.number, ""),
                    "atoms_a": atom_names(rdkit_chains[chain_a], left_indices),
                    "chain_b": chain_b,
                    "resname_b": item.protein.name,
                    "position_b": item.protein.number,
                    "residue_b": str(item.protein),
                    "mutation_b": mutation_maps.get(candidate, {}).get(item.protein.number, ""),
                    "atoms_b": atom_names(rdkit_chains[chain_b], right_indices),
                    "interaction": item.interaction,
                    "category": category(item.interaction),
                    "distance_a": metadata.get("distance", ""),
                    "angle_deg": metadata.get("DHA_angle", metadata.get("angle", "")),
                    "metadata_json": json.dumps(metadata, sort_keys=True, default=jsonable),
                    "input_cif": str(cif),
                })
        print(f"{source}: {sum(row['source'] == source for row in event_rows)} ProLIF events")

    event_fields = [
        "source", "candidate", "chain_pair", "chain_a", "resname_a", "position_a", "residue_a",
        "mutation_a", "atoms_a", "chain_b", "resname_b", "position_b", "residue_b", "mutation_b",
        "atoms_b", "interaction", "category", "distance_a", "angle_deg", "metadata_json", "input_cif",
    ]
    write_tsv(args.out / "interaction_events.tsv", event_rows, event_fields)

    unique = {}
    for row in event_rows:
        key = (row["source"], row["chain_pair"], row["residue_a"], row["residue_b"], row["interaction"])
        unique[key] = row
    unique_rows = list(unique.values())
    write_tsv(args.out / "unique_residue_interactions.tsv", unique_rows, event_fields)

    category_order = ["Hydrophobic", "HydrogenBond", "SaltBridge", "PiStacking", "CationPi", "VdWContact"]
    pair_rows = []
    type_rows = []
    for source in structures:
        for chain_a, chain_b in itertools.combinations("ABCD", 2):
            pair = f"{chain_a}-{chain_b}"
            rows = [row for row in unique_rows if row["source"] == source and row["chain_pair"] == pair]
            event_subset = [row for row in event_rows if row["source"] == source and row["chain_pair"] == pair]
            counts = Counter(row["category"] for row in rows)
            residue_pairs = {(row["residue_a"], row["residue_b"]) for row in rows}
            output = {
                "source": source,
                "chain_pair": pair,
                "unique_residue_pairs": len(residue_pairs),
                "unique_interaction_types": len(rows),
                "atom_level_events": len(event_subset),
                "non_vdw_interactions": sum(count for name, count in counts.items() if name != "VdWContact"),
            }
            output.update({name: counts[name] for name in category_order})
            pair_rows.append(output)
            for name in category_order:
                type_rows.append({
                    "source": source,
                    "chain_pair": pair,
                    "category": name,
                    "unique_residue_pair_types": counts[name],
                    "atom_level_events": sum(1 for row in event_subset if row["category"] == name),
                })
    pair_fields = ["source", "chain_pair", "unique_residue_pairs", "unique_interaction_types",
                   "atom_level_events", "non_vdw_interactions"] + category_order
    write_tsv(args.out / "chain_pair_summary.tsv", pair_rows, pair_fields)
    write_tsv(args.out / "interaction_type_summary.tsv", type_rows,
              ["source", "chain_pair", "category", "unique_residue_pair_types", "atom_level_events"])

    hotspot_counts: dict[tuple[str, str, str, int, str], Counter] = defaultdict(Counter)
    for row in unique_rows:
        for suffix in ("a", "b"):
            key = (str(row["source"]), str(row["candidate"]), str(row[f"chain_{suffix}"]),
                   int(row[f"position_{suffix}"]), str(row[f"resname_{suffix}"]))
            hotspot_counts[key][str(row["category"])] += 1
    hotspot_rows = []
    for (source, candidate, chain, position, resname), counts in hotspot_counts.items():
        hotspot_rows.append({
            "source": source, "candidate": candidate, "chain": chain, "resname": resname,
            "position": position, "mutation": mutation_maps.get(candidate, {}).get(position, ""),
            "interaction_types": sum(counts.values()),
            **{name: counts[name] for name in category_order},
        })
    hotspot_rows.sort(key=lambda row: (row["source"], -int(row["interaction_types"]), row["chain"], row["position"]))
    write_tsv(args.out / "interface_residue_hotspots.tsv", hotspot_rows,
              ["source", "candidate", "chain", "resname", "position", "mutation", "interaction_types"] + category_order)

    # Cross-model consensus for the top candidate, ignoring arbitrary homomer chain labels.
    top = next(iter(candidates))
    top_sources = [source for source in structures if source.startswith(top)]
    consensus_rows = []
    if len(top_sources) >= 2:
        support: dict[str, Counter] = {}
        exemplar: dict[tuple, dict[str, object]] = {}
        for source in top_sources:
            support[source] = Counter()
            for row in unique_rows:
                if row["source"] != source:
                    continue
                left = (row["resname_a"], int(row["position_a"]))
                right = (row["resname_b"], int(row["position_b"]))
                key = tuple(sorted((left, right))) + (row["category"],)
                support[source][key] += 1
                exemplar[key] = row
        common = set.intersection(*(set(values) for values in support.values()))
        for key in sorted(common, key=str):
            row = exemplar[key]
            consensus_rows.append({
                "residue_1": f"{key[0][0]}{key[0][1]}", "residue_2": f"{key[1][0]}{key[1][1]}",
                "category": key[2],
                **{f"chain_pair_support_{source}": support[source][key] for source in top_sources},
                "mutation_1": mutation_maps[top].get(key[0][1], ""),
                "mutation_2": mutation_maps[top].get(key[1][1], ""),
            })
        consensus_fields = ["residue_1", "residue_2", "category"] + [
            f"chain_pair_support_{source}" for source in top_sources
        ] + ["mutation_1", "mutation_2"]
        write_tsv(args.out / "top_candidate_cross_model_consensus.tsv", consensus_rows, consensus_fields)

    sns.set_theme(style="whitegrid")
    heat = pd.DataFrame(type_rows).pivot_table(
        index=["source", "chain_pair"], columns="category", values="unique_residue_pair_types", fill_value=0
    ).reindex(columns=category_order, fill_value=0)
    fig, ax = plt.subplots(figsize=(9.0, max(6.0, 0.27 * len(heat))))
    sns.heatmap(heat, cmap="mako", annot=True, fmt="g", linewidths=0.3, ax=ax)
    ax.set(xlabel="ProLIF interaction category", ylabel="Structure and chain pair",
           title="Homotetramer interface fingerprints")
    fig.tight_layout()
    fig.savefig(args.out / "interaction_type_heatmap.png", dpi=220)
    fig.savefig(args.out / "interaction_type_heatmap.svg")
    plt.close(fig)

    totals = pd.DataFrame(type_rows).groupby(["source", "category"])["unique_residue_pair_types"].sum().unstack(fill_value=0)
    totals = totals.reindex(columns=category_order, fill_value=0)
    ax = totals.plot(kind="bar", stacked=True, figsize=(10, 5.5),
                     color=["#0072B2", "#56B4E9", "#D55E00", "#CC79A7", "#E69F00", "#999999"])
    ax.set(xlabel="Structure", ylabel="Unique residue-pair interaction types",
           title="ProLIF tetramer interface chemistry")
    ax.legend(frameon=False, ncol=3, fontsize=8)
    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    plt.savefig(args.out / "interaction_category_totals.png", dpi=220)
    plt.savefig(args.out / "interaction_category_totals.svg")
    plt.close()

    for source in structures:
        graph = nx.Graph()
        graph.add_nodes_from("ABCD")
        for row in pair_rows:
            if row["source"] != source:
                continue
            left, right = str(row["chain_pair"]).split("-")
            graph.add_edge(left, right, total=int(row["unique_interaction_types"]),
                           nonvdw=int(row["non_vdw_interactions"]))
        pos = {"A": (-1.0, 1.0), "B": (1.0, 1.0), "C": (-1.0, -1.0), "D": (1.0, -1.0)}
        fig, ax = plt.subplots(figsize=(5.4, 5.2))
        widths = [0.8 + 0.22 * graph[u][v]["nonvdw"] for u, v in graph.edges]
        nx.draw_networkx_nodes(graph, pos, node_color="#56B4E9", node_size=1800, ax=ax)
        nx.draw_networkx_labels(graph, pos, font_size=15, font_weight="bold", ax=ax)
        nx.draw_networkx_edges(graph, pos, width=widths, edge_color="#444444", alpha=0.8, ax=ax)
        offsets = {"A-B": (0.0, 0.12), "A-C": (-0.15, 0.0), "A-D": (0.0, 0.18),
                   "B-C": (0.0, -0.18), "B-D": (0.15, 0.0), "C-D": (0.0, -0.12)}
        for u, v in graph.edges:
            pair = "-".join(sorted((u, v)))
            midpoint = ((pos[u][0] + pos[v][0]) / 2, (pos[u][1] + pos[v][1]) / 2)
            dx, dy = offsets[pair]
            ax.text(midpoint[0] + dx, midpoint[1] + dy,
                    f"{graph[u][v]['total']} ({graph[u][v]['nonvdw']})",
                    fontsize=9, ha="center", va="center", backgroundcolor="white")
        ax.set_title(f"{source}\nedge label: all unique types (non-vdW)")
        ax.axis("off")
        fig.tight_layout()
        fig.savefig(args.out / f"chain_network_{source}.png", dpi=220)
        fig.savefig(args.out / f"chain_network_{source}.svg")
        plt.close(fig)

    mutation_rows = [row for row in hotspot_rows if row["mutation"]]
    if mutation_rows:
        mutation_frame = pd.DataFrame(mutation_rows)
        mutation_frame["mutation_position"] = mutation_frame["mutation"]
        mutation_heat = mutation_frame.groupby(["mutation_position", "source"])["interaction_types"].sum().unstack(fill_value=0)
        fig, ax = plt.subplots(figsize=(max(8.0, 1.5 * len(mutation_heat.columns)), max(5.5, 0.33 * len(mutation_heat))))
        sns.heatmap(mutation_heat, cmap="rocket_r", annot=True, fmt="g", linewidths=0.3, ax=ax)
        ax.set(xlabel="Structure", ylabel="Designed mutation",
               title="Designed residues participating in ProLIF interface fingerprints")
        fig.tight_layout()
        fig.savefig(args.out / "mutation_interaction_hotspots.png", dpi=220)
        fig.savefig(args.out / "mutation_interaction_hotspots.svg")
        plt.close(fig)

    summary = {
        "prolif_version": plf.__version__,
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "interactions": INTERACTIONS,
        "structures": {source: str(path) for source, path in structures.items()},
        "atom_level_events": len(event_rows),
        "unique_residue_interaction_types": len(unique_rows),
        "top_candidate_consensus_interactions": len(consensus_rows),
        "charge_model": "Lys/Arg +1; Asp/Glu -1; charged termini; neutral His; RDKit explicit H",
    }
    (args.out / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
