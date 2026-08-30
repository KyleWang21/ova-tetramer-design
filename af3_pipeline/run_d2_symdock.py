#!/usr/bin/env python3
"""Generate one de novo D2-symmetric OVA docking trajectory with Rosetta SymDock."""

from __future__ import annotations

import argparse
import csv
import pathlib
import time

import numpy as np
import pyrosetta


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--monomer", type=pathlib.Path, required=True)
    ap.add_argument("--symmdef", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--seed", type=int, required=True)
    args = ap.parse_args()
    pyrosetta.init(
        f"-mute all -ex1 -ex2aro -constant_seed -jran {args.seed} "
        "-docking:dock_mcm_trans_magnitude 3.0 -docking:dock_mcm_rot_magnitude 8.0"
    )
    from pyrosetta import rosetta

    args.out.mkdir(parents=True, exist_ok=True)
    pose = pyrosetta.pose_from_pdb(str(args.monomer))
    rosetta.protocols.symmetry.SetupForSymmetryMover(str(args.symmdef)).apply(pose)
    if pose.num_chains() < 4:
        raise RuntimeError(f"symmetry setup produced {pose.num_chains()} chains")
    mover = rosetta.protocols.symmetric_docking.SymDockProtocol()
    mover.set_fullatom(True); mover.set_local_refine(False); mover.set_max_repeats(1)
    start = time.time(); mover.apply(pose); elapsed = time.time() - start
    name = f"d2_symdock_seed{args.seed}"
    pdb_path = args.out / f"{name}.pdb"; pose.dump_pdb(str(pdb_path))

    ordinary = pyrosetta.pose_from_pdb(str(pdb_path))
    scorefxn = pyrosetta.get_fa_scorefxn(); total = float(scorefxn(ordinary))
    split = ordinary.split_by_chain()
    separated = sum(float(scorefxn(split[index])) for index in range(1, len(split) + 1))
    pair_rows = []
    for chain_a in range(1, 5):
        for chain_b in range(chain_a + 1, 5):
            a, b = split[chain_a], split[chain_b]
            xyz_a = np.asarray([list(a.residue(i).xyz("CA")) for i in range(1, a.total_residue() + 1)])
            xyz_b = np.asarray([list(b.residue(i).xyz("CA")) for i in range(1, b.total_residue() + 1)])
            distance = np.sqrt(np.square(xyz_a[:, None] - xyz_b[None]).sum(-1))
            pair_rows.append((chain_a, chain_b, float(distance.min()), int((distance < 8).sum()),
                              int((distance.min(1) < 10).sum()), int((distance.min(0) < 10).sum())))
    counts = sorted((row[3] for row in pair_rows), reverse=True)
    interface_sizes = sorted(((row[4] + row[5]) / 2 for row in pair_rows), reverse=True)
    row = {
        "name": name, "seed": args.seed, "elapsed_seconds": elapsed,
        "mover_status": str(mover.get_last_move_status()), "total_score": total,
        "separated_score": separated, "interface_energy": total - separated,
        "contacting_chain_pairs": sum(count > 0 for count in counts),
        "ca_pairs_lt8_sorted": ",".join(map(str, counts)),
        "interface_residues_sorted": ",".join(f"{value:.1f}" for value in interface_sizes),
        "minimum_ca_distance": min(item[2] for item in pair_rows), "pdb": str(pdb_path),
    }
    with (args.out / f"{name}.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(row), delimiter="\t"); writer.writeheader(); writer.writerow(row)
    print(name, f"seconds={elapsed:.1f}", f"dG={row['interface_energy']:.2f}",
          f"pairs={row['contacting_chain_pairs']}", f"contacts={row['ca_pairs_lt8_sorted']}")


if __name__ == "__main__":
    main()
