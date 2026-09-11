#!/usr/bin/env bash
# Join completed core stages into fail-closed v1.1 tables.
# OpenDDE is optional: if an old result exists it is copied into the table as a
# diagnostic, but its absence never blocks finalization.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
ROOT="$PROJECT/experiments/current_pending_v1_screen"
PY_AF3=/root/400083/alphafold3/repo/.venv/bin/python
PY_FILTER=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"

for screen in "$ROOT"/*; do
  [[ -f "$screen/V1_CORE_COMPLETE" ]] || continue
  [[ ! -f "$screen/SCREEN_V1_COMPLETE" ]] || continue
  candidate=$(awk -F '\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="candidate")c=i;next} NR==2{print $c}' "$screen/final_candidates.tsv")
  rosetta_summary="$screen/rosetta/${candidate}_AF3_rosetta_summary.tsv"
  open_args=()
  if [[ -f "$screen/opendde_summary/opendde_candidate_summary.tsv" ]]; then
    open_args+=(--opendde-summary "$screen/opendde_summary/opendde_candidate_summary.tsv")
  elif [[ -f "$screen/opendde/opendde_candidate_summary.tsv" ]]; then
    open_args+=(--opendde-summary "$screen/opendde/opendde_candidate_summary.tsv")
  fi
  if [[ -f "$screen/opendde_summary/opendde_prolif_summary.tsv" ]]; then
    open_args+=(--opendde-prolif-summary "$screen/opendde_summary/opendde_prolif_summary.tsv")
  elif [[ -f "$screen/opendde/opendde_prolif_summary.tsv" ]]; then
    open_args+=(--opendde-prolif-summary "$screen/opendde/opendde_prolif_summary.tsv")
  fi
  "$PY_FILTER" af3_pipeline/run_final20_failclosed_filter.py \
    --final "$screen/final_candidates.tsv" \
    --af3-summary "$screen/af3_candidate_recomputed_summary.tsv" \
    --af3-models "$screen/af3_per_model_recomputed.tsv" \
    --prolif-unique "$screen/prolif_af3/unique_residue_interactions.tsv" \
    --reference references/1OVA.cif \
    --protenix-summary "$screen/protenix_summary/protenix_candidate_summary.tsv" \
    --stoichiometry-summary "$screen/stoichiometry_summary/stoichiometry_candidate_summary.tsv" \
    --rosetta-summary "$rosetta_summary" "${open_args[@]}" \
    --out "$screen/unified_v11_no_opendde"
  report_args=()
  if [[ -f "$screen/opendde_summary/opendde_calibration_summary.tsv" ]]; then
    report_args+=(--opendde-calibration "$screen/opendde_summary/opendde_calibration_summary.tsv")
  elif [[ -f "$screen/opendde/opendde_calibration_summary.tsv" ]]; then
    report_args+=(--opendde-calibration "$screen/opendde/opendde_calibration_summary.tsv")
  fi
  "$PY_AF3" af3_pipeline/make_final20_screen_report.py \
    --unified "$screen/unified_final/final20_failclosed_screen.tsv" \
    "${report_args[@]}" \
    --out "$screen/report_v11"
  touch "$screen/SCREEN_V1_COMPLETE"
  echo "$candidate generic v1.1 core screen complete"
done

"$PY_AF3" af3_pipeline/update_accepted_af3_registry.py \
  --project "$PROJECT" --out "$PROJECT/experiments/current_accepted_af3_registry"
bash af3_pipeline/refresh_current_accepted_delivery.sh
