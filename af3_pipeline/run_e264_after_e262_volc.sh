#!/usr/bin/env bash
# Fill the one missing v1.0 stage for historical Accepted_AF3 AI-MULTIHIT-11.
# This remains in the single Fire chain and only submits when c20250508 has
# more than ten free GPUs.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
PREVIOUS="$PROJECT/experiments/e262_c4_e243_crossbackbone_continuation_af3_seed1"
EXPERIMENT="$PROJECT/experiments/e264_c4_multihit11_stoichiometry_v1"
E186="$PROJECT/experiments/e186_c4_multihit11_multiseed"
E190="$PROJECT/experiments/e190_c4_multihit11_screen_v1"
E203="$PROJECT/experiments/e203_c4_multihit11_opendde_msa1024_v2"
PY_AF3=/root/400083/alphafold3/repo/.venv/bin/python
PY_FILTER=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"

while [[ ! -f "$PREVIOUS/SEED1_STRICT_SCREEN_COMPLETE" || ! -f "$EXPERIMENT/INPUTS_FROZEN" ]]; do
  sleep 30
done

while [[ ! -f "$EXPERIMENT/VOLC_AF3_SHARDS_SUBMITTED" ]]; do
  "$PY_AF3" af3_pipeline/submit_af3_volc_if_capacity.py \
    --experiment "$EXPERIMENT" --job-number-start 1184
  [[ -f "$EXPERIMENT/VOLC_AF3_SHARDS_SUBMITTED" ]] && break
  sleep 30
done

expected_models=$(awk -F'\t' '
  NR==1 {for (i=1; i<=NF; i++) if ($i=="expected_samples") c=i; next}
  {sum += $c}
  END {print sum+0}
' "$EXPERIMENT/manifest.tsv")
while true; do
  complete=$(find "$EXPERIMENT" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  [[ "$complete" -eq "$expected_models" ]] && break
  sleep 30
done

"$PY_AF3" af3_pipeline/score_outputs.py --experiment "$EXPERIMENT"
"$PY_AF3" af3_pipeline/summarize_final20_stoichiometry.py \
  --manifest "$EXPERIMENT/manifest.tsv" \
  --stoich-models "$EXPERIMENT/model_scores.tsv" \
  --c4-models "$E186/full25_screen/af3_per_model_recomputed.tsv" \
  --reference references/1OVA.cif --out "$EXPERIMENT/summary"

"$PY_FILTER" af3_pipeline/run_final20_failclosed_filter.py \
  --final "$E186/full25_screen_inputs/final_candidates.tsv" \
  --af3-summary "$E186/full25_screen/af3_candidate_recomputed_summary.tsv" \
  --af3-models "$E186/full25_screen/af3_per_model_recomputed.tsv" \
  --prolif-unique "$E190/prolif_crossmodel/unique_residue_interactions.tsv" \
  --reference references/1OVA.cif \
  --protenix-summary "$E190/protenix_summary/protenix_candidate_summary.tsv" \
  --stoichiometry-summary "$EXPERIMENT/summary/stoichiometry_candidate_summary.tsv" \
  --rosetta-summary "$E190/rosetta/ova_body_c4_multihit11_ms_AF3_rosetta_summary.tsv" \
  --opendde-summary "$E203/summary/opendde_candidate_summary.tsv" \
  --out "$E190/unified_final"
"$PY_AF3" af3_pipeline/make_final20_screen_report.py \
  --unified "$E190/unified_final/final20_failclosed_screen.tsv" \
  --opendde-calibration "$E203/summary/opendde_calibration_summary.tsv" \
  --out "$E190/report_final"

touch "$EXPERIMENT/STOICHIOMETRY_COMPLETE" "$E190/SCREEN_V1_COMPLETE"
"$PY_AF3" af3_pipeline/update_accepted_af3_registry.py \
  --project "$PROJECT" --out "$PROJECT/experiments/current_accepted_af3_registry"
echo "AI-MULTIHIT-11 historical v1.0 stoichiometry backfill complete"
