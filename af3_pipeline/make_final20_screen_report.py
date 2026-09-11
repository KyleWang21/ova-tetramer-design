#!/usr/bin/env python3
"""Create one compact table, stage matrix, and Markdown summary from the unified screen."""

from __future__ import annotations

import argparse
import csv
import html
import pathlib


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def field(row: dict[str, str], name: str) -> str:
    return row.get(name, "")


def formatted_float(value: str, digits: int = 3) -> str:
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return "NA"


def compact(row: dict[str, str]) -> dict[str, str]:
    return {
        "evidence_rank_not_acceptance": field(row, "provisional_rank_not_final"),
        "candidate": field(row, "candidate"),
        "n_mutations": field(row, "n_mutations"),
        "sequence_pass": field(row, "sequence_pass"),
        "surface_mutation_fraction": field(row, "surface_mutation_fraction_rsasa_ge0p20"),
        "af3_mean_iptm": field(row, "af3_mean_iptm"),
        "af3_min_iptm": field(row, "af3_min_iptm"),
        "af3_zero_clash_models_of25": field(row, "af3_zero_atomic_clash_models"),
        "af3_main_cluster_of25": field(row, "af3_main_topology_cluster_size"),
        "af3_status": field(row, "af3_status"),
        "prolif_nonvdw_design_count": field(row, "prolif_design_nonvdw_count"),
        "prolif_status": "PASS" if field(row, "prolif_pass") == "1" else "FAIL",
        "protenix_mean_iptm": field(row, "protenix_mean_iptm"),
        "protenix_min_iptm": field(row, "protenix_min_iptm"),
        "protenix_main_cluster_of5": field(row, "protenix_main_topology_cluster_size"),
        "af3_protenix_contact_jaccard": field(row, "protenix_af3_protenix_contact_jaccard"),
        "af3_protenix_shared_nonvdw": field(row, "protenix_shared_nonvdw_position_type_count"),
        "protenix_status": field(row, "protenix_status"),
        "c4_specificity_margin": field(row, "stoich_c4_specificity_margin"),
        "stoichiometry_status": field(row, "stoichiometry_status"),
        "rosetta_diagnostic_zero_clash_fraction": field(row, "rosetta_zero_clash_fraction"),
        "rosetta_diagnostic_weakest_sc": field(row, "rosetta_median_weakest_sc"),
        "rosetta_diagnostic_weakest_hbonds": field(row, "rosetta_median_weakest_hbonds"),
        "rosetta_diagnostic_unsat_per_1000A2": field(row, "rosetta_median_worst_unsat_per_1000A2"),
        "rosetta_status": field(row, "rosetta_status"),
        "opendde_mean_iptm_optional": field(row, "opendde_mean_iptm"),
        "opendde_selected_dockq_to_design": field(row, "opendde_selected_dockq_to_design"),
        "opendde_geometry_selectable_seeds_of10": field(row, "opendde_geometry_selectable_seed_count"),
        "opendde_prolif_design_nonvdw_count": field(row, "opendde_prolif_design_nonvdw_count"),
        "af3_opendde_shared_nonvdw": field(row, "af3_opendde_shared_nonvdw_position_type_count"),
        "opendde_prolif_optional_pass": field(row, "opendde_prolif_all_optional_gates_pass"),
        "opendde_status": field(row, "opendde_status"),
        "all_required_stages_pass": field(row, "all_required_stages_pass"),
        "failed_required_stages": field(row, "failed_required_stages"),
        "mutations": field(row, "mutations"),
        "af3_candidate_cif": field(row, "af3_candidate_cif"),
        "sequence": field(row, "sequence"),
    }


def stage_state(row: dict[str, str], stage: str) -> str:
    if stage == "Sequence":
        return "PASS" if row.get("sequence_pass") == "1" else "FAIL"
    if stage == "ProLIF":
        return "PASS" if row.get("prolif_pass") == "1" else "FAIL"
    if stage.startswith("OpenDDE"):
        if row.get("opendde_status") == "EVIDENCE":
            return "EVIDENCE"
        if row.get("opendde_status") == "COMPLETE_NO_ELIGIBLE_REP":
            return "NO_HIT"
        return "MISSING"
    return row.get({
        "AF3": "af3_status", "Protenix": "protenix_status",
        "Stoichiometry": "stoichiometry_status", "Rosetta": "rosetta_status",
    }[stage], "MISSING")


def draw_svg(path: pathlib.Path, rows: list[dict[str, str]]) -> None:
    stages = ["Sequence", "AF3", "ProLIF", "Protenix", "Stoichiometry", "Rosetta", "OpenDDE (optional)"]
    cell_w, cell_h, left, top = 128, 30, 235, 68
    width, height = left + len(stages) * cell_w + 25, top + len(rows) * cell_h + 45
    colors = {
        "PASS": "#22c55e", "FAIL": "#ef4444", "MISSING": "#94a3b8",
        "EVIDENCE": "#2563eb", "NO_HIT": "#f59e0b",
    }
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#172033}.h{font-size:13px;font-weight:700}.r{font-size:12px}.s{font-size:11px;fill:white;font-weight:700}</style>',
        '<text x="18" y="28" font-size="20" font-weight="700">OVA C4 fail-closed screening matrix</text>',
    ]
    for index, stage in enumerate(stages):
        parts.append(f'<text class="h" x="{left + index * cell_w + 8}" y="{top - 12}">{stage}</text>')
    for row_index, row in enumerate(rows):
        y = top + row_index * cell_h
        parts.append(f'<text class="r" x="18" y="{y + 20}">{html.escape(row["candidate"])}</text>')
        for stage_index, stage in enumerate(stages):
            state = stage_state(row, stage)
            state = state if state in colors else "MISSING"
            x = left + stage_index * cell_w
            parts.append(f'<rect x="{x}" y="{y + 3}" width="{cell_w - 6}" height="{cell_h - 6}" rx="4" fill="{colors[state]}"/>')
            text_x = x + (28 if state in {"EVIDENCE", "NO_HIT"} else 43)
            parts.append(f'<text class="s" x="{text_x}" y="{y + 20}">{state}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--unified", type=pathlib.Path, required=True)
    parser.add_argument("--opendde-calibration", type=pathlib.Path)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rows = read_tsv(args.unified)
    rows.sort(key=lambda row: int(row.get("provisional_rank_not_final", 10**9)))
    compact_rows = [compact(row) for row in rows]
    write_tsv(args.out / "final20_screen_compact.tsv", compact_rows)
    gate_rows = []
    for name in sorted({key for row in rows for key in row if "gate_" in key}):
        present = [row[name] for row in rows if row.get(name, "") != ""]
        gate_rows.append({
            "gate": name, "pass_count": str(sum(value == "1" for value in present)),
            "evaluated_count": str(len(present)), "candidate_count": str(len(rows)),
        })
    if gate_rows:
        write_tsv(args.out / "gate_pass_counts.tsv", gate_rows)
    draw_svg(args.out / "final20_stage_matrix.svg", rows)
    stages = {
        "sequence": sum(row.get("sequence_pass") == "1" for row in rows),
        "AF3": sum(row.get("af3_status") == "PASS" for row in rows),
        "ProLIF": sum(row.get("prolif_pass") == "1" for row in rows),
        "Protenix-v2": sum(row.get("protenix_status") == "PASS" for row in rows),
        "stoichiometry": sum(row.get("stoichiometry_status") == "PASS" for row in rows),
        "Rosetta": sum(row.get("rosetta_status") == "PASS" for row in rows),
    }
    def numbers(name: str) -> list[float]:
        output = []
        for row in rows:
            try:
                output.append(float(row[name]))
            except (KeyError, TypeError, ValueError):
                pass
        return output

    ptx_top = sorted(
        (row for row in rows if row.get("protenix_mean_iptm", "") != ""),
        key=lambda row: float(row["protenix_mean_iptm"]), reverse=True,
    )[:5]
    opendde_top = sorted(
        (row for row in rows if row.get("opendde_mean_iptm", "") != ""),
        key=lambda row: float(row["opendde_mean_iptm"]), reverse=True,
    )[:5]
    calibration = None
    if args.opendde_calibration and args.opendde_calibration.exists():
        calibration = read_tsv(args.opendde_calibration)[0]
    zero_clash = numbers("af3_zero_atomic_clash_models")
    specificity = numbers("stoich_c4_specificity_margin")
    relaxed_clean = sum(row.get("rosetta_zero_clash_fraction") == "1" for row in rows)
    lines = [
        "# OVA非共价四聚体候选筛选总结", "",
        "本文件由当前v1.1 fail-closed筛选脚本生成；OpenDDE仅为可选诊断，证据排名不等于Accepted。", "",
        "## 硬门槛通过数", "",
    ]
    candidate_count = len(rows)
    lines.extend(f"- {stage}: {count}/{candidate_count}" for stage, count in stages.items())
    lines.extend([
        "", f"- 全部required阶段通过：{sum(row.get('all_required_stages_pass') == '1' for row in rows)}/{candidate_count}",
        "", "## 关键结果", "",
        f"- AF3按2.4 Å重原子规则的零clash模型数范围：{min(zero_clash):.0f}–{max(zero_clash):.0f}/25；因此Accepted_AF3为{stages['AF3']}/{candidate_count}。" if zero_clash else "- AF3严格几何：未评估。",
        f"- C4化学计量特异性裕量范围：{min(specificity):.3f}–{max(specificity):.3f}；通过{stages['stoichiometry']}/{candidate_count}。" if specificity else "- 化学计量：未评估。",
        f"- Rosetta中{relaxed_clean}/{candidate_count}候选的松弛集合达到零clash；正式通过数为{stages['Rosetta']}/{candidate_count}。",
        "- Protenix-v2五seed最高：" + ", ".join(
            f"{row['candidate']}={float(row['protenix_mean_iptm']):.3f}" for row in ptx_top
        ) + ".",
        ("- OpenDDE第三模型证据可用："
         f"{sum(row.get('opendde_status') == 'EVIDENCE' for row in rows)}/{candidate_count}；平均ipTM最高：" +
         ", ".join(f"{row['candidate']}={float(row['opendde_mean_iptm']):.3f}" for row in opendde_top) + ".")
        if opendde_top else
        "- OpenDDE在对照校准完成前只作第三模型证据；DockQ无实验四聚体参考时只报告为DockQ_to_design。",
        "- OpenDDE不参与required总通过判定；DockQ无实验四聚体参考时只解释为DockQ_to_design。",
        ("- OpenDDE代表结构完成全部可选ProLIF门槛："
         f"{sum(row.get('opendde_prolif_all_optional_gates_pass') == '1' for row in rows)}/{candidate_count}。"),
        "", "## 产物", "",
        "- `final20_screen_compact.tsv`：单一候选总表。",
        "- `final20_stage_matrix.svg`：各候选逐阶段PASS/FAIL/MISSING矩阵。",
        "- `gate_pass_counts.tsv`：每一冻结硬门槛的通过数和已评估数。",
        "- 原始逐模型、逐相互作用和逐Relax结果保留在上级实验目录。", "",
    ])
    if calibration:
        insert_at = lines.index("## 产物") - 1
        lines.insert(insert_at, (
            "- OpenDDE校准：阴性对照最大平均ipTM="
            f"{formatted_float(calibration.get('negative_max_mean_iptm', ''))}，过程阳性最小平均ipTM="
            f"{formatted_float(calibration.get('positive_min_mean_iptm', ''))}；promoted_to_required_filter="
            f"{calibration['promoted_to_required_filter']}。"
        ))
    (args.out / "筛选总结.md").write_text("\n".join(lines))
    print("; ".join(f"{key}={value}/{candidate_count}" for key, value in stages.items()))


if __name__ == "__main__":
    main()
