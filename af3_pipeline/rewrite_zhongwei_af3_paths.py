#!/usr/bin/env python3
"""Point copied AF3 model tables at their candidate-local CIF copies."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True,
                        help="current_pending_v1_screen root")
    parser.add_argument("--model-root", type=Path, required=True,
                        help="root containing <candidate>/<row>.cif copies")
    args = parser.parse_args()
    changed = 0
    for table in sorted(args.root.glob("ova_body_c4_*/af3_per_model_recomputed.tsv")):
        candidate = table.parent.name
        with table.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if not rows:
            continue
        model_dir = args.model_root / candidate
        for index, row in enumerate(rows):
            cif = model_dir / f"{index}.cif"
            if not cif.is_file():
                raise FileNotFoundError(cif)
            row["cif_path"] = str(cif)
        temp = table.with_suffix(".tsv.tmp")
        with temp.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t",
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temp, table)
        changed += 1
    print(f"rewrote {changed} AF3 model tables")


if __name__ == "__main__":
    main()
