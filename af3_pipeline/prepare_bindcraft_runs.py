#!/usr/bin/env python3
"""Create reproducible stock-binder BindCraft jobs for ranked OVA hotspots."""

from __future__ import annotations

import argparse
import json
import pathlib


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--patches", type=pathlib.Path, required=True)
    ap.add_argument("--bindcraft", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--patch-ranks", default="1,2")
    ap.add_argument("--max-trajectories", type=int, default=4)
    ap.add_argument("--length-min", type=int, default=65)
    ap.add_argument("--length-max", type=int, default=85)
    ap.add_argument("--target-pdb", type=pathlib.Path,
                    help="override the full target with a hotspot-centred structural crop")
    args = ap.parse_args()

    patch_data = json.loads(args.patches.read_text())
    wanted = {int(x) for x in args.patch_ranks.split(",")}
    defaults = json.loads(
        (args.bindcraft / "settings_advanced/default_4stage_multimer.json").read_text()
    )
    filters = json.loads(
        (args.bindcraft / "settings_filters/no_filters.json").read_text()
    )
    target_pdb = (args.target_pdb or pathlib.Path(patch_data["target_pdb"])).resolve()
    params = args.bindcraft.parent / "params"

    manifest = []
    for patch in patch_data["patches"]:
        if int(patch["rank"]) not in wanted:
            continue
        name = str(patch["name"])
        root = (args.out / name).resolve()
        settings = root / "settings"
        settings.mkdir(parents=True, exist_ok=True)
        target = {
            "design_path": str(root / "design"),
            "binder_name": f"ova_p3_{name}",
            "starting_pdb": str(target_pdb),
            "chains": "A",
            "target_hotspot_residues": patch["bindcraft_hotspots"],
            "lengths": [args.length_min, args.length_max],
            "number_of_final_designs": 2,
        }
        advanced = dict(defaults)
        advanced.update({
            # A bounded pilot: real AF2 backpropagation + ProteinMPNN, but enough
            # iterations to be meaningful and small enough to learn from quickly.
            "design_algorithm": "4stage",
            "soft_iterations": 75,
            "temporary_iterations": 45,
            "hard_iterations": 5,
            "greedy_iterations": 15,
            "greedy_percentage": 1,
            "save_design_animations": False,
            "save_design_trajectory_plots": True,
            "zip_animations": False,
            "zip_plots": False,
            "save_trajectory_pickle": False,
            "max_trajectories": args.max_trajectories,
            "enable_rejection_check": False,
            "enable_relax_scoring": False,
            "enable_mpnn": True,
            # The vendored mask-cap extension calls its optional initializer in
            # binder mode; make the stock random initialization explicit.
            "init": "random",
            "num_seqs": 8,
            "max_mpnn_sequences": 2,
            "sampling_temp": 0.15,
            "save_mpnn_fasta": True,
            "af_params_dir": str(params.resolve()),
            "dssp_path": str((args.bindcraft / "functions/dssp").resolve()),
            "dalphaball_path": str((args.bindcraft / "functions/DAlphaBall.gcc").resolve()),
        })
        for filename, value in (
            ("target.json", target), ("advanced.json", advanced), ("filters.json", filters)
        ):
            (settings / filename).write_text(json.dumps(value, indent=2) + "\n")
        manifest.append({
            "rank": patch["rank"], "name": name,
            "hotspots": patch["bindcraft_hotspots"],
            "settings": str(settings), "design_path": target["design_path"],
        })
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
