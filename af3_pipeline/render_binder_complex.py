#!/usr/bin/env python3
"""Render two views of a full-length OVA/BindCraft AF3 complex."""

from __future__ import annotations

import argparse
import pathlib

import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import MMCIFParser


def equal_axes(ax, xyz: np.ndarray) -> None:
    lo, hi = xyz.min(0), xyz.max(0); center = (lo + hi) / 2
    radius = max(hi - lo) / 2 + 4
    ax.set_xlim(center[0]-radius, center[0]+radius)
    ax.set_ylim(center[1]-radius, center[1]+radius)
    ax.set_zlim(center[2]-radius, center[2]+radius)
    ax.set_box_aspect((1, 1, 1)); ax.set_axis_off()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--hotspots", default="99,155,157,182,184,334,336,338")
    args = ap.parse_args()
    hotspots = {int(x) for x in args.hotspots.split(",")}
    model = next(MMCIFParser(QUIET=True).get_structure("complex", str(args.cif)).get_models())
    traces = {}
    for chain in model:
        traces[chain.id] = [(int(r.id[1]), np.asarray(r["CA"].coord))
                            for r in chain if r.id[0] == " " and "CA" in r]
    all_xyz = np.concatenate([np.stack([p for _, p in pts]) for pts in traces.values()])
    fig = plt.figure(figsize=(12, 6))
    for panel, azim in enumerate((35, 215), 1):
        ax = fig.add_subplot(1, 2, panel, projection="3d")
        for cid, color, width, label in (("A", "#94a3b8", 1.8, "P3-13R"),
                                          ("B", "#dc2626", 5.0, "AI binder")):
            xyz = np.stack([p for _, p in traces[cid]])
            ax.plot(xyz[:,0], xyz[:,1], xyz[:,2], color=color, lw=width, alpha=.92, label=label)
        hxyz = np.stack([p for r, p in traces["A"] if r in hotspots])
        ax.scatter(hxyz[:,0], hxyz[:,1], hxyz[:,2], color="#facc15", edgecolor="#713f12",
                   s=55, depthshade=False, label="8-residue hotspot")
        equal_axes(ax, all_xyz); ax.view_init(elev=18, azim=azim)
    fig.suptitle("Independent AF3 validation of BindCraft patch-1 binder\n"
                 "median ipTM 0.75 · interface PAE 2.57 Å · 15/15 complexes")
    handles, labels = fig.axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False)
    fig.tight_layout(rect=(0, .06, 1, .92))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=220, bbox_inches="tight")
    fig.savefig(args.out.with_suffix(".svg"), bbox_inches="tight")
    print(args.out)


if __name__ == "__main__":
    main()
