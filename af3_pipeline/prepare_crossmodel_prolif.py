#!/usr/bin/env python3
"""Select AF3 and Protenix-v2 representative CIFs for normalized ProLIF comparison."""

from __future__ import annotations

import argparse
import csv
import pathlib


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def numeric(row: dict[str, str], field: str, default: float) -> float:
    try:
        return float(row.get(field, ""))
    except (TypeError, ValueError):
        return default


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--final", type=pathlib.Path, required=True)
    parser.add_argument("--af3-models", type=pathlib.Path, required=True)
    parser.add_argument("--protenix-summary", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    candidates = [row["candidate"] for row in read_tsv(args.final)]
    af3_rows: dict[str, list[dict[str, str]]] = {candidate: [] for candidate in candidates}
    for row in read_tsv(args.af3_models):
        if row["candidate"] in af3_rows:
            af3_rows[row["candidate"]].append(row)
    protenix = {row["candidate"]: row for row in read_tsv(args.protenix_summary)}

    output = []
    for candidate in candidates:
        # Representative must be in the dominant AF3 topology. Within that set,
        # prefer an exactly clean model, then confidence and weakest-chain support.
        models = af3_rows[candidate]
        if not models:
            raise ValueError(f"{candidate}: missing AF3 model rows")
        clusters = {row.get("topology_cluster", "") for row in models} - {""}
        if clusters:
            main_cluster = max(
                clusters,
                key=lambda cluster: sum(row.get("topology_cluster") == cluster for row in models),
            )
            eligible = [row for row in models if row.get("topology_cluster") == main_cluster]
        else:
            # Generic strict 25-model tables have already required a valid C4
            # network and two designed interfaces for every model, but do not
            # carry the legacy topology-cluster annotation.
            eligible = models
        best = max(eligible, key=lambda row: (
            numeric(row, "atomic_clash_count_2p4", 1.0) == 0.0,
            numeric(row, "valid_c4_network", 0.0) == 1.0,
            numeric(row, "two_interfaces_each_design_ge3", 0.0) == 1.0,
            numeric(row, "iptm", 0.0),
            numeric(
                row, "min_incident_iptm",
                numeric(row, "min_chain_plddt", numeric(row, "mean_chain_plddt", 0.0)) / 100.0,
            ),
            -numeric(
                row, "interface_pae", numeric(row, "max_chain_ca_rmsd_to_1ova", 999.0)
            ),
        ))
        output.append({
            "source": f"{candidate}__AF3", "candidate": candidate,
            "model_system": "AF3", "cif_path": best["cif_path"],
        })
        pt = protenix.get(candidate)
        if not pt or not pt.get("representative_cif_path"):
            raise ValueError(f"{candidate}: missing Protenix representative")
        output.append({
            "source": f"{candidate}__Protenix", "candidate": candidate,
            "model_system": "Protenix-v2", "cif_path": pt["representative_cif_path"],
        })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(output)
    print(f"wrote {len(output)} structures to {args.out}")


if __name__ == "__main__":
    main()
