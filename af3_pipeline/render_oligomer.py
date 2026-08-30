#!/usr/bin/env python3
"""Render a compact two-view CA trace of an AF3 oligomer and its module."""

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
    ax.set_box_aspect((1, 1, 1)); ax.set_axis_off()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--ova-end", type=int, default=386)
    ap.add_argument("--module-start", type=int, required=True)
    ap.add_argument("--module-end", type=int, required=True)
    ap.add_argument("--title", default="AF3 oligomer")
    ap.add_argument("--highlight-label", default="oligomerization module")
    args = ap.parse_args()
    model = next(MMCIFParser(QUIET=True).get_structure("complex", str(args.cif)).get_models())
    chains = []
    for chain in model:
        points = [(int(r.id[1]), np.asarray(r["CA"].coord)) for r in chain if r.id[0] == " " and "CA" in r]
        if points:
            chains.append((chain.id, points))
    all_xyz = np.concatenate([np.stack([p for _, p in pts]) for _, pts in chains])

    fig = plt.figure(figsize=(12, 6))
    for view, azim in enumerate((35, 215), 1):
        ax = fig.add_subplot(1, 2, view, projection="3d")
        for idx, (cid, pts) in enumerate(chains):
            color = COLORS[idx % len(COLORS)]
            for lo, hi, lw, alpha in (
                (1, args.ova_end, 1.6, .72),
                (args.ova_end + 1, args.module_start - 1, 1.2, .35),
                (args.module_start, args.module_end, 5.0, 1.0),
            ):
                xyz = np.stack([p for r, p in pts if lo <= r <= hi]) if any(lo <= r <= hi for r, _ in pts) else None
                if xyz is not None and len(xyz) > 1:
                    ax.plot(xyz[:,0], xyz[:,1], xyz[:,2], color=color, lw=lw, alpha=alpha,
                            label=f"chain {cid}" if view == 1 and lo == args.module_start else None)
        equal_axes(ax, all_xyz)
        ax.view_init(elev=18, azim=azim)
    fig.suptitle(args.title + f"\nthick segments: {args.highlight_label}")
    fig.legend(loc="lower center", ncol=len(chains), frameon=False)
    fig.tight_layout(rect=(0, .06, 1, .93))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=220, bbox_inches="tight")
    fig.savefig(args.out.with_suffix(".svg"), bbox_inches="tight")
    print(args.out)


if __name__ == "__main__":
    main()
