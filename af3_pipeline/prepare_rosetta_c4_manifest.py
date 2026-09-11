#!/usr/bin/env python3
"""Select one AF3 representative per final candidate for Rosetta C4 screening."""

import argparse
import csv
import pathlib
from collections import defaultdict


def number(row: dict[str, str], *keys: str, default: float = 0.0) -> float:
    for key in keys:
        value = row.get(key, "")
        if value not in {"", None}:
            return float(value)
    return default


def edge_list(value: str) -> list[str]:
    """Return canonical, de-duplicated chain-pair names from a class/list."""
    output = []
    seen = set()
    for token in value.replace("|", ",").split(","):
        token = token.strip()
        if not token or "-" not in token:
            continue
        left, right = (part.strip() for part in token.split("-", 1))
        if len(left) != 1 or len(right) != 1 or left == right:
            continue
        edge = "-".join(sorted((left, right)))
        if edge not in seen:
            seen.add(edge)
            output.append(edge)
    return output


def select_primary_interfaces(row: dict[str, str]) -> tuple[str, str, str, str]:
    """Select the highest-design-count interface class and its two copies.

    AF3 reports interface classes as symmetry-related edge groups, e.g.
    ``A-B,C-D|A-C,B-D|A-D,B-C``.  The class with the fewest designed-residue
    The selected class normally contains two symmetry-related chain pairs,
    such as ``A-B,C-D``.  Every other class remains in the effective-edge
    geometry set, but is excluded from chemistry minima.
    """
    effective = edge_list(row.get("effective_edge_list", ""))
    classes = [edge_list(group) for group in row.get("interface_type_edges", "").split("|")]
    classes = [group for group in classes if group]
    raw_counts = [part.strip() for part in row.get("interface_type_design_counts", "").split(",")]
    try:
        counts = [int(float(part)) for part in raw_counts if part != ""]
    except ValueError:
        counts = []
    if not classes or len(counts) != len(classes):
        # Older recomputed tables did not carry interface classes.  Preserve a
        # safe compatibility path while recording that no primary split was
        # available; new manifests always use the explicit class split.
        primary = effective
        status = "FALLBACK_ALL_EFFECTIVE"
        class_text = ""
        count_text = ""
    else:
        order = sorted(range(len(classes)), key=lambda i: (-counts[i], i))[:1]
        primary = edge_list(",".join(edge for i in order for edge in classes[i]))
        # Only pass edges that are present in the selected representative's
        # effective geometry set.  This is important when one class is absent
        # in a particular AF3 sample.
        primary = [edge for edge in primary if edge in set(effective)]
        status = "TOP1_PRIMARY_DESIGN_CLASS"
        class_text = "|".join(",".join(classes[i]) for i in order)
        count_text = ",".join(str(counts[i]) for i in order)
    if not primary:
        primary = effective
        status = "FALLBACK_ALL_EFFECTIVE_EMPTY_PRIMARY"
    background = [edge for edge in effective if edge not in set(primary)]
    return ",".join(effective), ",".join(primary), ",".join(background), status + ":" + class_text + ":" + count_text


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--models", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()
    grouped = defaultdict(list)
    with args.models.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            grouped[row["candidate"]].append(row)
    rows = []
    for candidate, models in sorted(grouped.items()):
        best = max(models, key=lambda row: (
            int(number(row, "atomic_clash_count_2p4", default=1.0) == 0),
            int(number(row, "valid_c4_network") == 1),
            int(number(row, "two_interfaces_each_design_ge3") == 1),
            number(row, "iptm"),
            number(row, "min_incident_iptm", "min_chain_plddt"),
            -number(row, "interface_pae", "max_chain_ca_rmsd_to_1ova", default=999.0),
        ))
        effective_edges, design_edges, background_edges, selection = select_primary_interfaces(best)
        rows.append({
            "candidate": candidate,
            "source_model": "AF3",
            "iptm": best["iptm"],
            "effective_edges": effective_edges,
            "design_interface_edges": design_edges,
            "background_edges": background_edges,
            "design_interface_selection": selection,
            "cif_path": str(pathlib.Path(best["cif_path"]).resolve()),
        })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    print(f"wrote {len(rows)} Rosetta representatives to {args.out}")


if __name__ == "__main__":
    main()
