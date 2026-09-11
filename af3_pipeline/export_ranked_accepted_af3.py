#!/usr/bin/env python3
"""Export rank-prefixed AF3 CIFs, one FASTA and one compact table from the registry."""

from __future__ import annotations

import argparse
import csv
import pathlib
import re
import shutil


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    with args.registry.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("Accepted_AF3 registry is empty")
    args.out.mkdir(parents=True, exist_ok=True)
    cif_dir = args.out / "cif"
    cif_dir.mkdir(exist_ok=True)
    compact: list[dict[str, object]] = []
    for row in rows:
        source = pathlib.Path(row["representative_af3_cif"])
        if not source.is_absolute():
            source = pathlib.Path.cwd() / source
        if not source.is_file():
            raise FileNotFoundError(f"missing representative AF3 CIF for {row['candidate']}: {source}")
        rank = int(row["rank"])
        filename = (
            f"RANK{rank:02d}_{safe_name(row['candidate'])}_"
            f"AF3_meaniptm{float(row['mean_iptm']):.4f}.cif"
        )
        destination = cif_dir / filename
        shutil.copyfile(source, destination)
        compact.append({
            "rank": rank,
            "candidate": row["candidate"],
            "n_mutations": row["n_mutations"],
            "surface_mutation_fraction": row["surface_mutation_fraction"],
            "af3_n_models": row["n_models"],
            "af3_mean_iptm": row["mean_iptm"],
            "af3_median_iptm": row["median_iptm"],
            "af3_min_iptm": row["min_iptm"],
            "af3_models_two_interfaces_each_contact_ge3_designed_mutations": row[
                "models_two_interfaces_design_ge3"
            ],
            "af3_mean_fraction_mutations_at_predicted_interface": row[
                "mean_design_interface_recurrence"
            ],
            "af3_min_fraction_mutations_at_predicted_interface": row[
                "min_design_interface_recurrence"
            ],
            "af3_union_interface_mutation_count": row[
                "union_design_interface_position_count"
            ],
            "af3_union_interface_mutation_positions": row[
                "union_design_interface_positions"
            ],
            "af3_all25_interface_mutation_count": row[
                "all_model_design_interface_position_count"
            ],
            "af3_all25_interface_mutation_positions": row[
                "all_model_design_interface_positions"
            ],
            "af3_predicted_interfaces_match_designed_mutations": int(
                int(row["models_two_interfaces_design_ge3"]) == int(row["n_models"])
            ),
            "prolif_design_nonvdw_count": row["prolif_design_nonvdw_count"],
            "prolif_design_nonvdw_fraction": row["prolif_design_nonvdw_fraction"],
            "accepted_af3": 1,
            "v1_screen_available": row["v1_screen_available"],
            "accepted_v1": row["accepted_v1"],
            "failed_required_stages": row["failed_required_stages"],
            # Keep the delivery self-contained.  ``args.out`` is an atomic
            # staging directory during refresh and must never leak into the
            # published table.
            "ranked_cif": str(pathlib.Path("cif") / filename),
            "source_af3_cif": row["representative_af3_cif"],
            "sequence": row["sequence"],
        })
    with (args.out / "all_sequences.fasta").open("w") as handle:
        for row in compact:
            handle.write(
                f">RANK{int(row['rank']):02d}_{row['candidate']} "
                f"nmut={row['n_mutations']} AF3_mean_ipTM={float(row['af3_mean_iptm']):.4f}\n"
                f"{row['sequence']}\n"
            )
    with (args.out / "unified_candidates.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(compact[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(compact)
    accepted_v1 = sum(int(row["accepted_v1"]) for row in compact)
    summary_lines = [
        "# 当前 OVA 非共价同源四聚体交付摘要",
        "",
        f"- Accepted_AF3：{len(compact)}",
        f"- Accepted_v1.1：{accepted_v1}",
        "- 排名依据：先按 Accepted_v1.1 状态，再按无模板 AF3 五 seed×五 sample 算术平均 ipTM。",
        "- 每个候选仅交付一个代表性 AF3 CIF；文件名以 `RANKNN_` 开头。",
        "",
        "| 排名 | 候选 | 突变数 | AF3平均ipTM | AF3最低ipTM | v1.1 | 失败阶段 |",
        "|---:|---|---:|---:|---:|---:|---|",
    ]
    for row in compact:
        failed = str(row["failed_required_stages"]) or "—"
        summary_lines.append(
            f"| {int(row['rank'])} | {row['candidate']} | {row['n_mutations']} | "
            f"{float(row['af3_mean_iptm']):.4f} | {float(row['af3_min_iptm']):.4f} | "
            f"{row['accepted_v1']} | {failed} |"
        )
    summary_lines.extend([
        "",
        "`unified_candidates.tsv` 包含突变—AF3预测界面复现、ProLIF非纯VdW作用、"
        "v1.1状态及AF3坐标来源等字段。",
    ])
    (args.out / "SUMMARY.md").write_text("\n".join(summary_lines) + "\n")
    print(f"exported candidates={len(compact)} CIFs={len(compact)} to {args.out}")


if __name__ == "__main__":
    main()
