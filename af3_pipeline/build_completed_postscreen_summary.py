#!/usr/bin/env python3
"""Build a compact report for every completed full post-screen row."""

from __future__ import annotations

import argparse
import csv
import pathlib
from collections import Counter


def read_rows(path: pathlib.Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_one(path: pathlib.Path) -> dict[str, str]:
    rows = read_rows(path)
    return rows[0] if rows else {}


def first(row: dict[str, str], *keys: str) -> str:
    for key in keys:
        if row.get(key, "") != "":
            return row[key]
    return ""


def passed(value: str) -> bool:
    return str(value).strip().upper() in {"1", "PASS"}


def compact(label: str, row: dict[str, str], rosetta_v11: dict[str, str] | None = None) -> dict[str, str]:
    rosetta_v11 = rosetta_v11 or {}
    zero_clash = row.get("af3_zero_atomic_clash_models", "")
    if not zero_clash and row.get("af3_gate_25_zero_atomic_clash") == "1":
        zero_clash = "25"
    valid_c4 = row.get("af3_valid_c4_network_models", "")
    if not valid_c4 and row.get("af3_gate_25_valid_c4_network") == "1":
        valid_c4 = "25"
    core_statuses = [
        ("sequence", row.get("sequence_pass", "")),
        ("af3", row.get("af3_status", "")),
        ("prolif", row.get("prolif_pass", "")),
        ("protenix", row.get("protenix_status", "")),
        ("stoichiometry", row.get("stoichiometry_status", "")),
        ("rosetta", "PASS" if rosetta_v11.get("rosetta_pass") == "1" else row.get("rosetta_status", "")),
    ]
    failed_core = [stage for stage, value in core_statuses if not passed(value)]
    return {
        "batch": label,
        "sample_id": row.get("candidate", ""),
        "n_mutations": row.get("n_mutations", ""),
        "surface_fraction": row.get("surface_mutation_fraction_rsasa_ge0p20", ""),
        "sequence_pass": row.get("sequence_pass", ""),
        "af3_mean_iptm": row.get("af3_mean_iptm", ""),
        "af3_min_iptm": row.get("af3_min_iptm", ""),
        "af3_zero_clash_models": zero_clash,
        "af3_valid_c4_models": valid_c4,
        "af3_status": row.get("af3_status", ""),
        "prolif_pass": row.get("prolif_pass", ""),
        "prolif_nonvdw_count": row.get("prolif_design_nonvdw_count", ""),
        "prolif_nonvdw_fraction": row.get("prolif_design_nonvdw_fraction", ""),
        "prolif_interface_nonvdw_counts": row.get("prolif_interface_type_nonvdw_counts", ""),
        "prolif_major_hotspots": row.get("prolif_major_hotspots", ""),
        "protenix_status": row.get("protenix_status", ""),
        "protenix_mean_iptm": row.get("protenix_mean_iptm", ""),
        "protenix_min_iptm": row.get("protenix_min_iptm", ""),
        "protenix_zero_clash_clean_models": row.get("protenix_full_connected_zero_clash_clean_models", ""),
        "protenix_contact_jaccard": row.get("protenix_af3_protenix_contact_jaccard", ""),
        "protenix_failed_gates": row.get("protenix_failed_protenix_gates", ""),
        "stoichiometry_status": row.get("stoichiometry_status", ""),
        "stoich_specificity_margin": first(row, "stoich_c4_specificity_margin"),
        "stoich_second_interface_unique": first(row, "stoich_c4_second_interface_unique_vs_c2_models"),
        "rosetta_status": ("PASS" if rosetta_v11.get("rosetta_pass") == "1" else "FAIL") if rosetta_v11 else row.get("rosetta_status", ""),
        "rosetta_screen_version": "v1.1_primary_design_interface" if rosetta_v11 else "v1.0_all_effective_edges",
        "rosetta_nrelax": rosetta_v11.get("nrelax", row.get("rosetta_nrelax", "")),
        "rosetta_zero_clash_fraction": rosetta_v11.get("zero_clash_fraction", row.get("rosetta_zero_clash_fraction", "")),
        "rosetta_weakest_sc": rosetta_v11.get("median_weakest_sc", row.get("rosetta_median_weakest_sc", "")),
        "rosetta_weakest_hbonds": rosetta_v11.get("median_weakest_hbonds", row.get("rosetta_median_weakest_hbonds", "")),
        "rosetta_unsat_per_1000A2": rosetta_v11.get("median_worst_unsat_per_1000A2", row.get("rosetta_median_worst_unsat_per_1000A2", "")),
        "rosetta_design_interface_edges": rosetta_v11.get("design_interface_edges", ""),
        "rosetta_background_edges": rosetta_v11.get("background_edges", ""),
        "rosetta_all_weakest_hbonds": rosetta_v11.get("median_all_weakest_hbonds", ""),
        "rosetta_all_unsat_per_1000A2": rosetta_v11.get("median_all_worst_unsat_per_1000A2", ""),
        "rosetta_failed_gates": rosetta_v11.get("failed_rosetta_gates", row.get("rosetta_failed_rosetta_gates", "")),
        "opendde_status": row.get("opendde_status", ""),
        "opendde_n_seeds": row.get("opendde_n_seeds", ""),
        "opendde_mean_iptm": row.get("opendde_mean_iptm", ""),
        "opendde_max_iptm": row.get("opendde_max_iptm", ""),
        "opendde_selectable_seeds": first(row, "opendde_selectable_seed_count", "opendde_geometry_selectable_seed_count"),
        # OpenDDE is intentionally excluded from the current required set.  Its
        # status and metrics remain in the report as optional diagnostics.
        "all_required_stages_pass": "1" if not failed_core else "0",
        "failed_required_stages": ";".join(failed_core),
    }


def fmt(value: str, digits: int = 3) -> str:
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return value or "—"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=pathlib.Path, default=pathlib.Path.cwd())
    parser.add_argument("--out-tsv", type=pathlib.Path, required=True)
    parser.add_argument("--out-md", type=pathlib.Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    inputs = [
        ("E151 Final20", root / "experiments/e151_final20_screen_v1/unified/final20_failclosed_screen.tsv"),
        ("AI-MULTIHIT-11", root / "experiments/e190_c4_multihit11_screen_v1/unified_current_after_opendde/final20_failclosed_screen.tsv"),
        ("OVA-C4-073", root / "experiments/e213_c4_073_screen_v1/unified_final/final20_failclosed_screen.tsv"),
        ("c06", root / "experiments/current_pending_v1_screen/ova_body_c4_06_ms25_r04/unified_final/final20_failclosed_screen.tsv"),
    ]
    v11_paths = {
        "AI-MULTIHIT-11": root / "experiments/e190_c4_multihit11_screen_v1/rosetta_primary_v11/ova_body_c4_multihit11_ms_AF3_rosetta_summary.tsv",
        "OVA-C4-073": root / "experiments/e213_c4_073_screen_v1/rosetta_primary_v11/ova_body_c4_073_ms_AF3_rosetta_summary.tsv",
        "c06": root / "experiments/current_pending_v1_screen/ova_body_c4_06_ms25_r04/rosetta_primary_v11/ova_body_c4_06_ms25_r04_AF3_rosetta_summary.tsv",
    }
    rows: list[dict[str, str]] = []
    source_files: list[str] = []
    for label, path in inputs:
        loaded = read_rows(path)
        if loaded:
            rosetta_v11 = read_one(v11_paths.get(label, pathlib.Path("/nonexistent")))
            rows.extend(compact(label, row, rosetta_v11) for row in loaded)
            source_files.append(str(path.relative_to(root)))
    fields = list(rows[0]) if rows else []
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    by_batch = Counter(row["batch"] for row in rows)
    md: list[str] = [
        "# 已完成完整后筛链的 OVA 四聚体样本",
        "",
        "本表只收录已经完成核心 fail-closed 后筛的结果；未完成的候选不按缺失值算通过。OpenDDE从当前required集合移除，仅作为可选诊断列保留。三条单候选若存在v1.1结果，Rosetta栏采用“主设计界面”口径；E151保留v1.0全有效接触口径。",
        "",
        f"共 {len(rows)} 条记录：E151 Final20 的 {by_batch.get('E151 Final20', 0)} 条，以及 3 条单独完成样本。",
        "",
        "## 批次级结果",
        "",
        "| 批次 | 数量 | sequence | AF3 | ProLIF | Protenix-v2 | 化学计量 | Rosetta | OpenDDE（可选） | 全部 required（不含OpenDDE） |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for label in ("E151 Final20", "AI-MULTIHIT-11", "OVA-C4-073", "c06"):
        subset = [r for r in rows if r["batch"] == label]
        if not subset:
            continue
        def npass(field: str) -> str:
            return f"{sum(r.get(field) in ('1', 'PASS') for r in subset)}/{len(subset)}"
        complete_opendde = sum(bool(r.get("opendde_n_seeds")) for r in subset)
        md.append(
            f"| {label} | {len(subset)} | {npass('sequence_pass') if 'sequence_pass' in subset[0] else ('—')} | "
            f"{npass('af3_status')} | {npass('prolif_pass')} | {npass('protenix_status')} | "
            f"{npass('stoichiometry_status')} | {npass('rosetta_status')} | "
            f"{complete_opendde}/{len(subset)} complete | "
            f"{npass('all_required_stages_pass')} |"
        )
    md += [
        "",
        "说明：OpenDDE列仅表示可选诊断是否完成/有证据，不参与全部required判定。E151的OpenDDE已完成10 seeds，但0条有可选代表结构。",
        "",
        "## E151 Final20 逐样本摘要",
        "",
        "| 样本 | 突变/表面比例 | AF3均值/最低 | 零clash/C4 | ProLIF | Protenix均值/最低 | Rosetta SC/H-bond/unsat | OpenDDE均值/最高（可选） | 失败 required（不含OpenDDE） |",
        "|---|---:|---:|---:|---:|---:|---|---:|---|",
    ]
    for row in rows:
        if row["batch"] != "E151 Final20":
            continue
        md.append(
            f"| {row['sample_id']} | {row['n_mutations']}/{fmt(row['surface_fraction'])} | "
            f"{fmt(row['af3_mean_iptm'])}/{fmt(row['af3_min_iptm'])} | "
            f"{row['af3_zero_clash_models'] or '—'}/25; {row['af3_valid_c4_models'] or '—'}/25 | "
            f"{row['prolif_pass']} | {fmt(row['protenix_mean_iptm'])}/{fmt(row['protenix_min_iptm'])} | "
            f"{fmt(row['rosetta_weakest_sc'])}/{row['rosetta_weakest_hbonds'] or '—'}/{fmt(row['rosetta_unsat_per_1000A2'], 2)} | "
            f"{fmt(row['opendde_mean_iptm'])}/{fmt(row['opendde_max_iptm'])} | {row['failed_required_stages']} |"
        )
    md += [
        "",
        "## 三条单独完成样本",
        "",
        "| 样本 | 突变 | AF3均值/最低 | AF3零clash | ProLIF | Protenix均值/最低 | 计量特异性 | Rosetta SC/H-bond/unsat | OpenDDE均值/最高（可选） | 失败 required（不含OpenDDE） |",
        "|---|---:|---:|---:|---|---:|---:|---|---:|---|",
    ]
    for row in rows:
        if row["batch"] == "E151 Final20":
            continue
        md.append(
            f"| {row['sample_id']} | {row['n_mutations']} | {fmt(row['af3_mean_iptm'])}/{fmt(row['af3_min_iptm'])} | "
            f"{row['af3_zero_clash_models']}/25 | {row['prolif_pass']} "
            f"({row['prolif_nonvdw_count']} non-VdW) | {fmt(row['protenix_mean_iptm'])}/{fmt(row['protenix_min_iptm'])} | "
            f"{fmt(row['stoich_specificity_margin'])}; unique {row['stoich_second_interface_unique'] or '—'} | "
            f"{fmt(row['rosetta_weakest_sc'])}/{row['rosetta_weakest_hbonds'] or '—'}/{fmt(row['rosetta_unsat_per_1000A2'], 2)} | "
            f"{fmt(row['opendde_mean_iptm'])}/{fmt(row['opendde_max_iptm'])} | {row['failed_required_stages']} |"
        )
    md += [
        "",
        "## 解释",
        "",
        "- `COMPLETE_NO_ELIGIBLE_REP` 表示 OpenDDE的10个seed都跑完了，但没有同时满足ipTM>0.80和几何条件的代表结构；该状态只作诊断，不再导致required失败。",
        "- Protenix-v2 的 `0/5 geometry-clean`、低均值/最低 ipTM 是这几条候选最一致的独立模型失败模式。",
        "- Rosetta v1.1 的化学门槛只在一类主设计界面的两个对称副本上计算；其他界面类型及对角背景链对仍计入完整四聚体的clash、穿插和拓扑检查，但不因没有氢键单独淘汰。",
        "",
        "逐样本机器可读结果见 `COMPLETED_POSTSCREEN_RESULTS.tsv`；原始宽表仍保留在各实验目录。",
        "",
        "### 原始来源",
        "",
    ]
    md.extend(f"- `{path}`" for path in source_files)
    args.out_md.write_text("\n".join(md) + "\n")
    print(f"wrote {len(rows)} rows to {args.out_tsv} and {args.out_md}")


if __name__ == "__main__":
    main()
