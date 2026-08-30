#!/usr/bin/env python3
"""Generate one de novo C4-symmetric OVA docking trajectory with Rosetta SymDock."""

from __future__ import annotations

import argparse
import csv
import pathlib
import time

import numpy as np
import pyrosetta


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--monomer", type=pathlib.Path, required=True)
    parser.add_argument("--symmdef", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--max-repeats", type=int, default=1)
    args = parser.parse_args()

    pyrosetta.init(
        f"-mute all -ex1 -ex2aro -constant_seed -jran {args.seed} "
        "-docking:dock_mcm_trans_magnitude 3.0 -docking:dock_mcm_rot_magnitude 8.0"
    )
    from pyrosetta import rosetta

    args.out.mkdir(parents=True, exist_ok=True)
    pose = pyrosetta.pose_from_pdb(str(args.monomer))
    rosetta.protocols.symmetry.SetupForSymmetryMover(str(args.symmdef)).apply(pose)
    mover = rosetta.protocols.symmetric_docking.SymDockProtocol()
    mover.set_fullatom(True)
    mover.set_local_refine(False)
    mover.set_max_repeats(args.max_repeats)
    start = time.time()
    mover.apply(pose)
    elapsed = time.time() - start

    name = f"c4_symdock_seed{args.seed}"
    pdb_path = args.out / f"{name}.pdb"
    pose.dump_pdb(str(pdb_path))

    # Reload as an ordinary four-chain pose so total-minus-separated energies are
    # directly comparable and do not include virtual symmetry residues.
    ordinary = pyrosetta.pose_from_pdb(str(pdb_path))
    scorefxn = pyrosetta.get_fa_scorefxn()
    total = float(scorefxn(ordinary))
    split = ordinary.split_by_chain()
    separated = sum(float(scorefxn(split[chain])) for chain in range(1, len(split) + 1))

    pair_stats: list[tuple[int, int, float, int, int, int]] = []
    for chain_a, chain_b in [(1, 2), (2, 3), (3, 4), (4, 1), (1, 3), (2, 4)]:
        a, b = split[chain_a], split[chain_b]
        xyz_a = np.asarray([
            list(a.residue(i).xyz("CA")) for i in range(1, a.total_residue() + 1)
        ])
        xyz_b = np.asarray([
            list(b.residue(i).xyz("CA")) for i in range(1, b.total_residue() + 1)
        ])
        distances = np.sqrt(np.square(xyz_a[:, None] - xyz_b[None]).sum(-1))
        pair_stats.append((
            chain_a,
            chain_b,
            float(distances.min()),
            int((distances < 8.0).sum()),
            int((distances.min(1) < 10.0).sum()),
            int((distances.min(0) < 10.0).sum()),
        ))
    adjacent = pair_stats[:4]
    opposite = pair_stats[4:]
    row = {
        "name": name,
        "seed": args.seed,
        "elapsed_seconds": elapsed,
        "mover_status": str(mover.get_last_move_status()),
        "total_score": total,
        "separated_score": separated,
        "interface_energy": total - separated,
        "mean_adjacent_ca_pairs_lt8": float(np.mean([item[3] for item in adjacent])),
        "mean_adjacent_interface_residues": float(np.mean([item[4] + item[5] for item in adjacent]) / 2),
        "min_adjacent_ca_distance": float(min(item[2] for item in adjacent)),
        "min_opposite_ca_distance": float(min(item[2] for item in opposite)),
        "pdb": str(pdb_path),
    }
    with (args.out / f"{name}.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    print(
        name,
        f"seconds={elapsed:.1f}",
        f"dG={row['interface_energy']:.2f}",
        f"CApairs={row['mean_adjacent_ca_pairs_lt8']:.1f}",
        f"iface_res={row['mean_adjacent_interface_residues']:.1f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
