#!/usr/bin/env python3
"""Constrained full-C4 FastRelax and weakest-edge InterfaceAnalyzer screening.

This is the executable Rosetta stage for screening-plan v1.1. It relaxes the complete
tetramer (never isolated subunits), keeps the docking pose with coordinate constraints,
then evaluates the fixed primary design-interface copies for chemical gates and all
effective chain pairs for topology/clash diagnostics on the relaxed coordinates. The
script supports the production 100-replicate AF3 or Protenix ensemble; `--nrelax 1`
is only a timing/diagnostic run and is labeled incomplete in the output.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import pathlib
import statistics

import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Select
from scipy.spatial import cKDTree

import pyrosetta


METRICS = [
    "dG_separated", "dSASA_int", "sc_value", "hbonds_int",
    "delta_unsatHbonds", "nres_int", "packstat", "dG_dSASA_ratio",
]


class ProteinOnly(Select):
    def accept_residue(self, residue):
        return int(residue.id[0] == " ")

    def accept_atom(self, atom):
        return int(atom.element != "H" and atom.altloc in {" ", "", "A"})


def convert_to_pdb(source: pathlib.Path, target: pathlib.Path) -> None:
    if source.suffix.lower() == ".pdb":
        target.write_text(source.read_text())
        return
    structure = MMCIFParser(QUIET=True).get_structure(source.stem, str(source))
    io = PDBIO(); io.set_structure(structure); io.save(str(target), ProteinOnly())


def chain_letters(pose) -> list[str]:
    info = pose.pdb_info()
    return [info.chain(pose.chain_begin(index)) for index in range(1, pose.num_chains() + 1)]


def interface_residues(pose, cutoff: float = 10.0) -> set[int]:
    output = set()
    for left in range(1, pose.total_residue() + 1):
        residue_left = pose.residue(left)
        if not residue_left.is_protein():
            continue
        xyz_left = residue_left.nbr_atom_xyz()
        for right in range(left + 1, pose.total_residue() + 1):
            if pose.chain(left) == pose.chain(right):
                continue
            residue_right = pose.residue(right)
            if residue_right.is_protein() and xyz_left.distance(residue_right.nbr_atom_xyz()) <= cutoff:
                output.update((left, right))
    return output


def constrained_relax(clean_pose, scorefxn, seed: int):
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.core.pack.task import TaskFactory, operation
    from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector, NotResidueSelector
    from pyrosetta.rosetta.protocols.relax import FastRelax

    pose = clean_pose.clone()
    pyrosetta.rosetta.numeric.random.rg().set_seed(seed)
    selected = interface_residues(pose, 10.0)
    selector = ResidueIndexSelector(",".join(map(str, sorted(selected))))
    tf = TaskFactory()
    tf.push_back(operation.InitializeFromCommandline())
    tf.push_back(operation.RestrictToRepacking())
    tf.push_back(operation.OperateOnResidueSubset(
        operation.PreventRepackingRLT(), NotResidueSelector(selector)
    ))
    movemap = MoveMap(); movemap.set_bb(False); movemap.set_chi(False)
    for index in selected:
        movemap.set_bb(index, True); movemap.set_chi(index, True)
    relax = FastRelax(scorefxn, 1)
    relax.set_task_factory(tf)
    relax.set_movemap(movemap)
    relax.constrain_relax_to_start_coords(True)
    relax.ramp_down_constraints(False)
    relax.apply(pose)
    return pose


def pair_pose(pose, left_chain: int, right_chain: int):
    split = pose.split_by_chain()
    output = split[left_chain].clone()
    pyrosetta.rosetta.core.pose.append_pose_to_pose(output, split[right_chain], True)
    info = output.pdb_info()
    for index in range(output.chain_begin(1), output.chain_end(1) + 1):
        info.chain(index, "A")
    for index in range(output.chain_begin(2), output.chain_end(2) + 1):
        info.chain(index, "B")
    # split_by_chain preserves residue connection metadata that still points to
    # the original full-pose indices. Re-detect the four native intrachain Cys
    # pairs after appending, otherwise scoring can report a stale partner index.
    output.conformation().detect_disulfides()
    return output


def analyze_pair(pose, left_chain: int, right_chain: int, scorefxn) -> dict[str, float | None]:
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
    from pyrosetta.rosetta.std import set_int_t
    pair = pair_pose(pose, left_chain, right_chain)
    fixed = set_int_t(); fixed.add(1)
    mover = InterfaceAnalyzerMover(fixed, True, scorefxn, True, True, True)
    mover.set_pack_separated(True); mover.set_compute_packstat(True); mover.apply(pair)
    scores = pair.scores
    get = lambda name: float(scores[name]) if name in scores else None
    return {
        "dG_separated": get("dG_separated"),
        "dSASA_int": get("dSASA_int"),
        "sc_value": get("sc_value"),
        "hbonds_int": get("hbonds_int"),
        "delta_unsatHbonds": get("delta_unsatHbonds"),
        "nres_int": get("nres_int"),
        "packstat": get("packstat"),
        "dG_dSASA_ratio": get("dG_separated/dSASAx100"),
    }


def pose_atoms(pose, chain_index: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    coords, residue_numbers, atom_names = [], [], []
    info = pose.pdb_info()
    for index in range(pose.chain_begin(chain_index), pose.chain_end(chain_index) + 1):
        residue = pose.residue(index)
        if not residue.is_protein():
            continue
        number = int(info.number(index))
        for atom_index in range(1, residue.nheavyatoms() + 1):
            xyz = residue.xyz(atom_index)
            coords.append((xyz.x, xyz.y, xyz.z))
            residue_numbers.append(number)
            atom_names.append(residue.atom_name(atom_index).strip())
    return np.asarray(coords), np.asarray(residue_numbers), np.asarray(atom_names)


def geometry(pose) -> dict[str, object]:
    chains = chain_letters(pose)
    arrays = {chain: pose_atoms(pose, index + 1) for index, chain in enumerate(chains)}
    clashes = overlaps = 0
    effective = []
    contact_counts = {}
    for left, right in itertools.combinations(chains, 2):
        left_coords, left_res, left_names = arrays[left]
        right_coords, right_res, right_names = arrays[right]
        neighbors = cKDTree(left_coords).query_ball_tree(cKDTree(right_coords), 5.3)
        residue_pairs = {(int(left_res[i]), int(right_res[j]))
                         for i, hits in enumerate(neighbors) for j in hits}
        contact_counts[f"{left}-{right}"] = len(residue_pairs)
        if len(residue_pairs) >= 20:
            effective.append(f"{left}-{right}")
        clashes += sum(len(hits) for hits in cKDTree(left_coords).query_ball_tree(cKDTree(right_coords), 2.4))
        ca_left = left_coords[left_names == "CA"]
        ca_right = right_coords[right_names == "CA"]
        overlaps += sum(len(hits) for hits in cKDTree(ca_left).query_ball_tree(cKDTree(ca_right), 2.5))
    degree = {chain: 0 for chain in chains}
    for edge in effective:
        left, right = edge.split("-"); degree[left] += 1; degree[right] += 1
    return {
        "atomic_clash_count_2p4": clashes,
        "severe_ca_overlap_count_2p5": overlaps,
        "effective_edges": ",".join(effective),
        "effective_edge_count": len(effective),
        "min_chain_degree": min(degree.values()),
        "valid_c4_network": int(len(effective) >= 4 and min(degree.values()) >= 2),
        "contact_counts_json": json.dumps(contact_counts, sort_keys=True),
    }


def ca_rmsd(reference, mobile) -> float:
    reference = reference - reference.mean(axis=0)
    mobile = mobile - mobile.mean(axis=0)
    left, _, right = np.linalg.svd(mobile.T @ reference)
    correction = np.eye(3); correction[-1, -1] = np.sign(np.linalg.det(left @ right))
    aligned = mobile @ (left @ correction @ right)
    return float(np.sqrt(np.mean(np.sum((aligned - reference) ** 2, axis=1))))


def chain_rmsds(reference_pose, mobile_pose) -> list[float]:
    output = []
    for chain_index in range(1, reference_pose.num_chains() + 1):
        reference, _, reference_names = pose_atoms(reference_pose, chain_index)
        mobile, _, mobile_names = pose_atoms(mobile_pose, chain_index)
        output.append(ca_rmsd(reference[reference_names == "CA"], mobile[mobile_names == "CA"]))
    return output


def summarize(values: list[float]) -> dict[str, object]:
    return {
        "mean": statistics.mean(values), "median": statistics.median(values),
        "min": min(values), "max": max(values), "values": values,
    }


def parse_edges(value: str) -> list[tuple[str, str]]:
    """Parse a comma-separated chain-pair list and reject malformed input."""
    output = []
    seen = set()
    for token in value.split(","):
        token = token.strip()
        if not token:
            continue
        parts = token.split("-")
        if len(parts) != 2 or any(len(part) != 1 for part in parts):
            raise ValueError(f"invalid chain edge: {token!r}")
        edge = (parts[0], parts[1])
        if edge[0] == edge[1]:
            raise ValueError(f"self edge is not allowed: {token!r}")
        canonical = "-".join(sorted(edge))
        if canonical in seen:
            continue
        seen.add(canonical)
        output.append(tuple(canonical.split("-")))
    if not output:
        raise ValueError("at least one chain edge is required")
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--in", dest="source", type=pathlib.Path, required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--effective-edges", required=True, help="comma list such as A-B,A-C,B-D,C-D")
    parser.add_argument(
        "--design-interface-edges",
        default="",
        help=(
            "comma list of the primary designed interface edges. Chemistry gates "
            "(SC, H-bonds, dG/dSASA, unsatisfied polar atoms) use only these edges; "
            "all effective edges remain in the geometry/topology diagnostics. "
            "If omitted, all effective edges are used for backwards compatibility."
        ),
    )
    parser.add_argument("--source-model", choices=("AF3", "Protenix"), required=True)
    parser.add_argument("--nrelax", type=int, default=100)
    parser.add_argument("--seed-start", type=int, default=7001)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    pdb = args.out / f"{args.candidate}_{args.source_model}_clean.pdb"
    convert_to_pdb(args.source, pdb)

    pyrosetta.init(" ".join([
        "-mute all", "-use_input_sc", "-ignore_unrecognized_res", "-ignore_zero_occupancy false",
        "-load_PDB_components false", "-no_fconfig", "-constant_seed",
    ]))
    clean = pyrosetta.pose_from_pdb(str(pdb))
    if clean.num_chains() != 4:
        raise ValueError(f"expected 4 chains, got {clean.num_chains()}")
    letters = chain_letters(clean)
    chain_index = {letter: index + 1 for index, letter in enumerate(letters)}
    edges = parse_edges(args.effective_edges)
    design_edges = parse_edges(args.design_interface_edges or args.effective_edges)
    effective_edge_names = {"-".join(edge) for edge in edges}
    design_edge_names = {"-".join(edge) for edge in design_edges}
    if not design_edge_names.issubset(effective_edge_names):
        missing = sorted(design_edge_names - effective_edge_names)
        raise ValueError(
            "design-interface-edges must be a subset of effective-edges; "
            f"missing={','.join(missing)}"
        )
    scorefxn = pyrosetta.create_score_function("ref2015")
    rows = []
    best = None
    for replicate in range(args.nrelax):
        seed = args.seed_start + replicate
        relaxed = constrained_relax(clean, scorefxn, seed)
        geo = geometry(relaxed)
        rmsds = chain_rmsds(clean, relaxed)
        pair_metrics = {}
        for left, right in edges:
            pair_metrics[f"{left}-{right}"] = analyze_pair(
                relaxed, chain_index[left], chain_index[right], scorefxn
            )
        design_pair_metrics = {
            edge: pair_metrics[edge] for edge in sorted(design_edge_names)
        }
        all_weakest_sc = min(value["sc_value"] for value in pair_metrics.values())
        all_weakest_hbonds = min(value["hbonds_int"] for value in pair_metrics.values())
        all_worst_ratio = max(value["dG_dSASA_ratio"] for value in pair_metrics.values())
        all_worst_unsat_density = max(
            value["delta_unsatHbonds"] / value["dSASA_int"] * 1000.0
            for value in pair_metrics.values()
        )
        # The production chemistry gates deliberately use only the primary
        # design interface class (normally two symmetry-related copies).  A
        # diagonal/background contact can be geometrically valid while having
        # no hydrogen bond, so including it in a minimum would reject an
        # otherwise well-supported design.  Geometry/topology below still uses
        # the complete effective-edge set.
        weakest_sc = min(value["sc_value"] for value in design_pair_metrics.values())
        weakest_hbonds = min(value["hbonds_int"] for value in design_pair_metrics.values())
        worst_ratio = max(value["dG_dSASA_ratio"] for value in design_pair_metrics.values())
        worst_unsat_density = max(
            value["delta_unsatHbonds"] / value["dSASA_int"] * 1000.0
            for value in design_pair_metrics.values()
        )
        row = {
            "candidate": args.candidate, "source_model": args.source_model,
            "replicate": replicate + 1, "seed": seed,
            **geo,
            "design_interface_edges": ",".join(sorted(design_edge_names)),
            "background_edges": ",".join(sorted(effective_edge_names - design_edge_names)),
            "max_chain_ca_rmsd_to_input": max(rmsds),
            "weakest_sc": weakest_sc,
            "weakest_hbonds": weakest_hbonds,
            "worst_dG_dSASA_x100": worst_ratio,
            "worst_unsat_per_1000A2": worst_unsat_density,
            "all_weakest_sc": all_weakest_sc,
            "all_weakest_hbonds": all_weakest_hbonds,
            "all_worst_dG_dSASA_x100": all_worst_ratio,
            "all_worst_unsat_per_1000A2": all_worst_unsat_density,
            "pair_metrics_json": json.dumps(pair_metrics, sort_keys=True),
        }
        rows.append(row)
        key = (int(row["valid_c4_network"]), -int(row["atomic_clash_count_2p4"]),
               float(row["weakest_sc"]), -float(row["worst_dG_dSASA_x100"]))
        if best is None or key > best[0]:
            best = (key, relaxed.clone())
        print(args.candidate, args.source_model, f"relax {replicate + 1}/{args.nrelax}",
              f"clash={row['atomic_clash_count_2p4']}", f"SCmin={weakest_sc:.3f}", flush=True)
    fields = list(rows[0])
    with (args.out / f"{args.candidate}_{args.source_model}_rosetta_replicates.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    best[1].dump_pdb(str(args.out / f"{args.candidate}_{args.source_model}_relaxed_best.pdb"))
    complete = args.nrelax >= 100
    summary = {
        "candidate": args.candidate, "source_model": args.source_model,
        "nrelax": args.nrelax, "production_ensemble_complete": int(complete),
        "design_interface_edges": ",".join(sorted(design_edge_names)),
        "background_edges": ",".join(sorted(effective_edge_names - design_edge_names)),
        "topology_retained_fraction": statistics.mean(int(row["valid_c4_network"]) for row in rows),
        "zero_clash_fraction": statistics.mean(int(row["atomic_clash_count_2p4"]) == 0 for row in rows),
        "zero_overlap_fraction": statistics.mean(int(row["severe_ca_overlap_count_2p5"]) == 0 for row in rows),
        "median_weakest_sc": statistics.median(float(row["weakest_sc"]) for row in rows),
        "median_weakest_hbonds": statistics.median(float(row["weakest_hbonds"]) for row in rows),
        "median_worst_dG_dSASA_x100": statistics.median(float(row["worst_dG_dSASA_x100"]) for row in rows),
        "median_worst_unsat_per_1000A2": statistics.median(float(row["worst_unsat_per_1000A2"]) for row in rows),
        "max_chain_ca_rmsd_to_input": max(float(row["max_chain_ca_rmsd_to_input"]) for row in rows),
        "median_all_weakest_sc": statistics.median(float(row["all_weakest_sc"]) for row in rows),
        "median_all_weakest_hbonds": statistics.median(float(row["all_weakest_hbonds"]) for row in rows),
        "median_all_worst_dG_dSASA_x100": statistics.median(float(row["all_worst_dG_dSASA_x100"]) for row in rows),
        "median_all_worst_unsat_per_1000A2": statistics.median(float(row["all_worst_unsat_per_1000A2"]) for row in rows),
    }
    gates = {
        "gate_rosetta_ensemble_100": complete,
        "gate_rosetta_topology_ge0p90": summary["topology_retained_fraction"] >= 0.90,
        "gate_rosetta_zero_clash_fraction_ge0p90": summary["zero_clash_fraction"] >= 0.90,
        "gate_rosetta_zero_overlap_fraction_ge0p90": summary["zero_overlap_fraction"] >= 0.90,
        "gate_rosetta_weakest_sc_ge0p60": summary["median_weakest_sc"] >= 0.60,
        "gate_rosetta_weakest_hbonds_ge3": summary["median_weakest_hbonds"] >= 3.0,
        "gate_rosetta_dG_dSASA_le_minus1": summary["median_worst_dG_dSASA_x100"] <= -1.0,
        "gate_rosetta_unsat_per_1000A2_le2": summary["median_worst_unsat_per_1000A2"] <= 2.0,
        "gate_rosetta_chain_rmsd_le2": summary["max_chain_ca_rmsd_to_input"] <= 2.0,
    }
    summary.update({key: int(value) for key, value in gates.items()})
    summary["rosetta_pass"] = int(all(gates.values()))
    summary["failed_rosetta_gates"] = ";".join(key for key, value in gates.items() if not value)
    (args.out / f"{args.candidate}_{args.source_model}_rosetta_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    with (args.out / f"{args.candidate}_{args.source_model}_rosetta_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(summary), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerow(summary)


if __name__ == "__main__":
    main()
