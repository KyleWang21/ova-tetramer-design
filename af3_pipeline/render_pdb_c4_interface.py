#!/usr/bin/env python3
"""Render two CA-trace views of a four-chain PDB and highlight design positions."""

from __future__ import annotations

import argparse
import pathlib

import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import PDBParser


COLORS = ["#2563eb", "#dc2626", "#16a34a", "#9333ea"]


def equal_axes(ax, xyz: np.ndarray, pad: float = 4.0) -> None:  # noqa: ANN001
    lo, hi = xyz.min(0), xyz.max(0)
    center = (lo + hi) / 2
    radius = max(hi - lo) / 2 + pad
    ax.set_xlim(center[0] - radius, center[0] + radius)
    ax.set_ylim(center[1] - radius, center[1] + radius)
    ax.set_zlim(center[2] - radius, center[2] + radius)
    ax.set_box_aspect((1, 1, 1)); ax.set_axis_off()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", type=pathlib.Path, required=True)
    ap.add_argument("--positions", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--title", default="Noncovalent C4 design backbone")
    args = ap.parse_args()
    positions = {int(value) for value in args.positions.read_text().strip().split(",")}
    model = next(PDBParser(QUIET=True).get_structure("c4", str(args.pdb)).get_models())
    chains = []
    for chain in model:
        points = [(int(res.id[1]), np.asarray(res["CA"].coord))
                  for res in chain if res.id[0] == " " and "CA" in res and int(res.id[1]) <= 386]
        if len(points) == 386:
            chains.append((str(chain.id), points))
    all_xyz = np.concatenate([np.stack([xyz for _, xyz in points]) for _, points in chains])

    fig = plt.figure(figsize=(12, 6))
    for view, (elev, azim) in enumerate(((88, 0), (16, 35)), 1):
        ax = fig.add_subplot(1, 2, view, projection="3d")
        for index, (chain_id, points) in enumerate(chains):
            color = COLORS[index % len(COLORS)]
            xyz = np.stack([point for _, point in points])
            interface = np.stack([point for residue, point in points if residue in positions])
            ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], color=color, lw=1.6, alpha=0.62,
                    label=f"chain {chain_id}")
            ax.scatter(interface[:, 0], interface[:, 1], interface[:, 2], color=color,
                       s=20, edgecolor="black", linewidth=0.25, depthshade=False)
        equal_axes(ax, all_xyz)
        ax.view_init(elev=elev, azim=azim)
    fig.suptitle(args.title + "\nblack-edged points: ProteinMPNN design surface; no engineered Cys")
    fig.legend(loc="lower center", ncol=4, frameon=False)
    fig.tight_layout(rect=(0, 0.06, 1, 0.92))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=220, bbox_inches="tight")
    fig.savefig(args.out.with_suffix(".svg"), bbox_inches="tight")
    print(args.out)


if __name__ == "__main__":
    main()
