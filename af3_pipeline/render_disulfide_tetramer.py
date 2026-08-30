#!/usr/bin/env python3
"""Render an AF3 OVA tetramer and its predicted inter-chain disulfides."""

from __future__ import annotations

import argparse
import pathlib

import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import MMCIFParser


COLORS = ["#2563eb", "#dc2626", "#16a34a", "#9333ea", "#ea580c", "#0891b2"]


def equal_axes(ax, xyz: np.ndarray, pad: float = 4.0) -> None:
    lo, hi = xyz.min(0), xyz.max(0)
    center = (lo + hi) / 2
    radius = max(hi - lo) / 2 + pad
    ax.set_xlim(center[0] - radius, center[0] + radius)
    ax.set_ylim(center[1] - radius, center[1] + radius)
    ax.set_zlim(center[2] - radius, center[2] + radius)
    ax.set_box_aspect((1, 1, 1))
    ax.set_axis_off()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--pairs", default="96-190,155-196")
    ap.add_argument("--cutoff", type=float, default=2.6)
    ap.add_argument("--title", default="AF3 OVA-BODY-SS-04 tetramer")
    ap.add_argument("--annotate-bonds", action="store_true")
    args = ap.parse_args()
    intended = [tuple(map(int, item.split("-", 1))) for item in args.pairs.split(",")]

    model = next(MMCIFParser(QUIET=True).get_structure("complex", str(args.cif)).get_models())
    chains = []
    sg: dict[tuple[str, int], np.ndarray] = {}
    for chain in model:
        ca = []
        for residue in chain:
            if residue.id[0] != " ":
                continue
            position = int(residue.id[1])
            if "CA" in residue:
                ca.append(np.asarray(residue["CA"].coord))
            if "SG" in residue:
                sg[(str(chain.id), position)] = np.asarray(residue["SG"].coord)
        if ca:
            chains.append((str(chain.id), np.stack(ca)))
    all_xyz = np.concatenate([xyz for _, xyz in chains])

    # Greedily match each intended residue pair exactly as in the quantitative scorer.
    bonds = []
    chain_ids = [cid for cid, _ in chains]
    for left, right in intended:
        candidates = []
        for a in chain_ids:
            for b in chain_ids:
                if a == b or (a, left) not in sg or (b, right) not in sg:
                    continue
                distance = float(np.linalg.norm(sg[(a, left)] - sg[(b, right)]))
                if distance <= args.cutoff:
                    candidates.append((distance, (a, left), (b, right)))
        used = set()
        for distance, first, second in sorted(candidates):
            if first in used or second in used:
                continue
            used.update((first, second))
            bonds.append((distance, first, second))

    fig = plt.figure(figsize=(12, 6.2))
    for view, azim in enumerate((28, 118), 1):
        ax = fig.add_subplot(1, 2, view, projection="3d")
        for idx, (cid, xyz) in enumerate(chains):
            ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], color=COLORS[idx % len(COLORS)],
                    lw=1.8, alpha=.82, label=f"chain {cid}")
        for distance, first, second in bonds:
            xyz = np.stack([sg[first], sg[second]])
            ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], color="#f4b400", lw=5.0,
                    solid_capstyle="round")
            if args.annotate_bonds:
                mid = xyz.mean(0)
                ax.text(*mid, f"{first[0]}{first[1]}–{second[0]}{second[1]}\n{distance:.2f} Å",
                        fontsize=6, color="#6b4d00")
        equal_axes(ax, all_xyz)
        ax.view_init(elev=18, azim=azim)
    distance_text = ""
    if bonds:
        distances = [bond[0] for bond in bonds]
        distance_text = f"; S–S {np.median(distances):.2f} Å median ({min(distances):.2f}–{max(distances):.2f})"
    fig.suptitle(
        f"{args.title}\n386-aa OVA body only; gold = predicted inter-chain S–S bonds ({len(bonds)}){distance_text}"
    )
    handles, labels = fig.axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(chains), frameon=False)
    fig.tight_layout(rect=(0, .06, 1, .93))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=240, bbox_inches="tight")
    fig.savefig(args.out.with_suffix(".svg"), bbox_inches="tight")
    print(f"{args.out}\tbonds={len(bonds)}")


if __name__ == "__main__":
    main()
