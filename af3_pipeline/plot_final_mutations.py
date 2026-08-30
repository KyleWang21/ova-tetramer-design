#!/usr/bin/env python3
"""Plot mutation classes along the current pure-OVA tetramer candidate."""

from __future__ import annotations

import argparse
import pathlib

import matplotlib.pyplot as plt


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--label", default="OVA-BODY-C4-SS3-AI03")
    ap.add_argument("--p3", default="14,69,77,79,191,192,193")
    ap.add_argument("--disulfide", default="94,96,155,190,196,213")
    ap.add_argument("--ai", default="93,97,99,101,152,195,237,239,340")
    args = ap.parse_args()
    p3 = [int(value) for value in args.p3.split(",") if value]
    disulfide = [int(value) for value in args.disulfide.split(",") if value]
    ai_lock = [int(value) for value in args.ai.split(",") if value]
    fig, ax = plt.subplots(figsize=(11, 3.8))
    ax.hlines(0, 1, 386, color="#777777", lw=3)
    ax.scatter(p3, [.22] * len(p3), marker="v", s=75, color="#6b7280",
               label=f"original P3-13R ({len(p3)})")
    ax.scatter(disulfide, [0] * len(disulfide), marker="o", s=85, color="#e6a000",
               label=f"three directed Cys pairs ({len(disulfide)})", zorder=3)
    ax.scatter(ai_lock, [-.22] * len(ai_lock), marker="^", s=75, color="#2563eb",
               label=f"ProteinMPNN lock-ring surface ({len(ai_lock)})")
    for index, position in enumerate(disulfide):
        ax.text(position, .08 + .08 * (index % 2), str(position), ha="center", va="bottom",
                fontsize=8, rotation=45)
    ax.axvspan(258, 265, color="#2ca02c", alpha=.18, label="SIINFEKL preserved")
    ax.axvspan(289, 295, color="#9333ea", alpha=.12, label="N292 neighbourhood preserved")
    ax.set_xlim(1, 386); ax.set_ylim(-.45, .45)
    ax.set_yticks([]); ax.set_xlabel("OVA sequence position (source Excel numbering)")
    ax.set_title(f"{args.label}: all changes remain within the 386-aa OVA body")
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(.5, -.22))
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=220, bbox_inches="tight")
    print(args.out)


if __name__ == "__main__":
    main()
