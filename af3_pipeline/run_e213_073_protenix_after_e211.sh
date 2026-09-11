#!/usr/bin/env bash
# Give the local A100 to the missing OVA-C4-073 Protenix-v2 five-seed screen
# after all priority E168 AF3 follow-ups have released it.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
PREREQUISITE="$PROJECT/experiments/e211_c4_e168_rev001_multiseed/RESEED_DECISION_COMPLETE"
SCREEN="$PROJECT/experiments/e213_c4_073_screen_v1"
PROTENIX="$SCREEN/protenix_local"
FINAL="$PROJECT/experiments/e165_c4_iptm90_073_multiseed/final25_screen_inputs/final_candidates.tsv"
AF3_MODELS="$PROJECT/experiments/e165_c4_iptm90_073_multiseed/final25_screen/af3_per_model_recomputed.tsv"
SCORE_PY=/root/400083/alphafold3/repo/.venv/bin/python
PROLIF_PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"

while [[ ! -f "$PREREQUISITE" ]]; do
  sleep 30
done

for shard in 0 1 2 3 4; do
  seed=$((11 + shard * 10))
  found=$(find "$PROTENIX/out" -path "*/seed_${seed}/*_summary_confidence_sample_0.json" 2>/dev/null | wc -l || true)
  [[ "$found" -eq 1 ]] || bash af3_pipeline/run_protenix_multiseed_volc_shard.sh "$PROTENIX" "$shard"
done
complete=$(find "$PROTENIX/out" -path '*/seed_*/*_summary_confidence_sample_0.json' 2>/dev/null | wc -l || true)
[[ "$complete" -eq 5 ]] || { echo "expected five Protenix outputs, found $complete" >&2; exit 2; }

export PYTHONPATH=/root/400083/alphafold3/repo/src${PYTHONPATH:+:$PYTHONPATH}
"$SCORE_PY" af3_pipeline/summarize_final20_protenix_multiseed.py \
  --final "$FINAL" --af3-models "$AF3_MODELS" --protenix-root "$PROTENIX/out" \
  --reference references/1OVA.cif --out "$SCREEN/protenix_summary_preprolif"
"$SCORE_PY" af3_pipeline/prepare_crossmodel_prolif.py \
  --final "$FINAL" --af3-models "$AF3_MODELS" \
  --protenix-summary "$SCREEN/protenix_summary_preprolif/protenix_candidate_summary.tsv" \
  --out "$SCREEN/prolif_crossmodel/structure_manifest.tsv"
"$PROLIF_PY" af3_pipeline/analyze_prolif_interfaces.py \
  --structure-manifest "$SCREEN/prolif_crossmodel/structure_manifest.tsv" \
  --final-candidates "$FINAL" --out "$SCREEN/prolif_crossmodel"
"$SCORE_PY" af3_pipeline/summarize_crossmodel_prolif.py \
  --unique "$SCREEN/prolif_crossmodel/unique_residue_interactions.tsv" \
  --out "$SCREEN/prolif_crossmodel/crossmodel_nonvdw.tsv"
"$SCORE_PY" af3_pipeline/summarize_final20_protenix_multiseed.py \
  --final "$FINAL" --af3-models "$AF3_MODELS" --protenix-root "$PROTENIX/out" \
  --reference references/1OVA.cif \
  --shared-nonvdw "$SCREEN/prolif_crossmodel/crossmodel_nonvdw.tsv" \
  --out "$SCREEN/protenix_summary"
touch "$SCREEN/PROTENIX_FIVE_SEEDS_COMPLETE"
echo "OVA-C4-073 Protenix-v2 five-seed screen complete"
