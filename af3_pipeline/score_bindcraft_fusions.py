#!/usr/bin/env python3
"""Score whether fused BindCraft domains bind OVA in trans across four chains."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import statistics

from alphafold3.structure import from_mmcif

from score_outputs import contact_details, largest_component, parse_ranges, representative_coords


HOTSPOTS = {99, 155, 157, 182, 184, 334, 336, 338}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--experiment", type=pathlib.Path, required=True)
    args = ap.parse_args()
    manifest = list(csv.DictReader((args.experiment / "manifest.tsv").open(), delimiter="\t"))
    rows = []
    for meta in manifest:
        n = int(meta["n_chains"]); ova = parse_ranges(meta["ova_ranges"], n)
        binder = parse_ranges(meta["binder_ranges"], n)
        for cif in sorted(args.experiment.glob(f"out_s*/{meta['name']}/seed-*/*_model.cif")):
            s = from_mmcif(cif.read_text()); chains = list(s.chains); coords = representative_coords(s)
            trans_edges = []; binder_engaged = set(); trans_contacts = 0; cis_contacts = 0
            contacted_hotspots = set(); directed = {}
            for i, ci in enumerate(chains):
                ri, xi = coords[ci]; bi = (ri >= binder[i][0]) & (ri <= binder[i][1])
                for j, cj in enumerate(chains):
                    rj, xj = coords[cj]; oj = (rj >= ova[j][0]) & (rj <= ova[j][1])
                    count, mask, _, ova_mask = contact_details(xi[bi], xj[oj])
                    directed[f"{ci}->{cj}"] = count
                    if i == j:
                        cis_contacts += count
                    else:
                        trans_contacts += count
                        if count >= 5:
                            trans_edges.append((ci, cj)); binder_engaged.add(ci)
                        canonical = rj[oj][ova_mask] - ova[j][0] + 1
                        contacted_hotspots.update(HOTSPOTS & set(map(int, canonical)))
            undirected = [(a, b) for a, b in trans_edges]
            rows.append({
                "system": meta["name"], "model": cif.parent.name,
                "trans_binder_contacts": trans_contacts, "cis_binder_contacts": cis_contacts,
                "all_binders_trans_engaged": int(len(binder_engaged) == n),
                "trans_largest_component": largest_component(chains, undirected),
                "hotspot_coverage": len(contacted_hotspots),
                "hotspots_json": json.dumps(sorted(contacted_hotspots)),
                "directed_contacts_json": json.dumps(directed, sort_keys=True),
                "cif_path": str(cif),
            })
    if not rows:
        raise SystemExit("no fusion models found")
    with (args.experiment / "bindcraft_fusion_model_scores.tsv").open("w", newline="") as h:
        w = csv.DictWriter(h, list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    print("system\tmodels\ttrans-engaged\ttrans-contacts\tcis-contacts\thotspots")
    for meta in manifest:
        z = [x for x in rows if x["system"] == meta["name"]]
        if not z: continue
        print(meta["name"], len(z),
              round(sum(x["all_binders_trans_engaged"] for x in z) / len(z), 3),
              statistics.median(x["trans_binder_contacts"] for x in z),
              statistics.median(x["cis_binder_contacts"] for x in z),
              statistics.median(x["hotspot_coverage"] for x in z), sep="\t")


if __name__ == "__main__":
    main()
