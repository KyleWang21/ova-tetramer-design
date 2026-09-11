#!/usr/bin/env bash
# Complete the remaining OVA-C4-073 stoichiometry/Rosetta stages and write one
# fail-closed v1.1 table. Existing OpenDDE output is retained only as optional
# diagnostic evidence.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E165="$PROJECT/experiments/e165_c4_iptm90_073_multiseed"
E171="$PROJECT/experiments/e171_c4_073_stoichiometry_v1"
E213="$PROJECT/experiments/e213_c4_073_screen_v1"
E214="$PROJECT/experiments/e214_c4_073_opendde_v1"
PY_AF3=/root/400083/alphafold3/repo/.venv/bin/python
PY_FILTER=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"

while true; do
  complete=$(find "$E171" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  [[ "$complete" -eq 75 ]] && break
  sleep 30
done
"$PY_AF3" af3_pipeline/score_outputs.py --experiment "$E171"
"$PY_AF3" af3_pipeline/summarize_final20_stoichiometry.py \
  --manifest "$E171/manifest.tsv" --stoich-models "$E171/model_scores.tsv" \
  --c4-models "$E165/final25_screen/af3_per_model_recomputed.tsv" \
  --reference references/1OVA.cif --out "$E213/stoichiometry"
touch "$E213/STOICHIOMETRY_COMPLETE"

while true; do
  rosetta_summary=$(find "$E213/rosetta" -maxdepth 1 -name '*_rosetta_summary.tsv' 2>/dev/null | head -1 || true)
  [[ -n "$rosetta_summary" ]] && break
  sleep 30
done

"$PY_FILTER" af3_pipeline/run_final20_failclosed_filter.py \
  --final "$E165/final25_screen_inputs/final_candidates.tsv" \
  --af3-summary "$E165/final25_screen/af3_candidate_recomputed_summary.tsv" \
  --af3-models "$E165/final25_screen/af3_per_model_recomputed.tsv" \
  --prolif-unique "$E165/prolif_af3_medoid/unique_residue_interactions.tsv" \
  --reference references/1OVA.cif \
  --protenix-summary "$E213/protenix_summary/protenix_candidate_summary.tsv" \
  --stoichiometry-summary "$E213/stoichiometry/stoichiometry_candidate_summary.tsv" \
  --rosetta-summary "$rosetta_summary" \
  --opendde-summary "$E214/summary/opendde_candidate_summary.tsv" \
  --out "$E213/unified_final"
"$PY_AF3" af3_pipeline/make_final20_screen_report.py \
  --unified "$E213/unified_final/final20_failclosed_screen.tsv" \
  --opendde-calibration "$E214/summary/opendde_calibration_summary.tsv" \
  --out "$E213/report"
touch "$E213/SCREEN_V1_COMPLETE"
echo "OVA-C4-073 screening v1.1 core stages complete (OpenDDE optional)"
