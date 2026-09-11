#!/usr/bin/env bash
# After the final E174 reseed decision, validate AI11 with Protenix-v2 first,
# then hand the local A100 to E177 unless Volc already owns those shards.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E189="$PROJECT/experiments/e189_c4_multihit15_multiseed"
E177="$PROJECT/experiments/e177_c4_tied_proteinmpnn_af3_seed1"
PROTENIX="$PROJECT/experiments/e190_c4_multihit11_screen_v1/protenix_local"
MARKER="$PROTENIX/PROTENIX_FIVE_SEEDS_COMPLETE"
FINAL="$PROJECT/experiments/e186_c4_multihit11_multiseed/full25_screen_inputs/final_candidates.tsv"
AF3_MODELS="$PROJECT/experiments/e186_c4_multihit11_multiseed/full25_screen/af3_per_model_recomputed.tsv"
PRE="$PROJECT/experiments/e190_c4_multihit11_screen_v1/protenix_summary_preprolif"
CROSS="$PROJECT/experiments/e190_c4_multihit11_screen_v1/prolif_crossmodel"
SUMMARY="$PROJECT/experiments/e190_c4_multihit11_screen_v1/protenix_summary"
SCORE_PY=/root/400083/alphafold3/repo/.venv/bin/python
PROLIF_PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"

deadline=$((SECONDS + 259200))
while [[ ! -f "$E189/RESEED_DECISION_COMPLETE" ]]; do
  if (( SECONDS >= deadline )); then
    echo "timed out waiting for E189 reseed decision" >&2
    exit 2
  fi
  sleep 30
done

# When the eight E177 shards were claimed by Volc, Zhongwei+Volc already use
# the user-authorized 16 GPUs.  Do not turn the local card into a 17th worker.
if [[ -f "$E177/VOLC_AF3_SHARDS_SUBMITTED" ]]; then
  while true; do
    complete=$(find "$E177" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
    [[ "$complete" -eq 160 ]] && break
    if (( SECONDS >= deadline )); then
      echo "timed out waiting for Volc E177: found $complete/160" >&2
      exit 2
    fi
    sleep 30
  done
fi

if [[ ! -f "$MARKER" ]]; then
  for shard in 0 1 2 3 4; do
    seed=$((11 + shard * 10))
    found=$(find "$PROTENIX/out" -path "*/seed_${seed}/*_summary_confidence_sample_0.json" 2>/dev/null | wc -l || true)
    if [[ "$found" -eq 1 ]]; then
      echo "AI11 Protenix seed $seed already complete"
      continue
    fi
    bash af3_pipeline/run_protenix_multiseed_volc_shard.sh "$PROTENIX" "$shard"
  done
  complete=$(find "$PROTENIX/out" -path '*/seed_*/*_summary_confidence_sample_0.json' 2>/dev/null | wc -l || true)
  if [[ "$complete" -ne 5 ]]; then
    echo "expected five Protenix summaries, found $complete" >&2
    exit 2
  fi
  "$SCORE_PY" af3_pipeline/summarize_final20_protenix_multiseed.py \
    --final "$FINAL" --af3-models "$AF3_MODELS" \
    --protenix-root "$PROTENIX/out" --reference references/1OVA.cif --out "$PRE"
  "$SCORE_PY" af3_pipeline/prepare_crossmodel_prolif.py \
    --final "$FINAL" --af3-models "$AF3_MODELS" \
    --protenix-summary "$PRE/protenix_candidate_summary.tsv" \
    --out "$CROSS/structure_manifest.tsv"
  "$PROLIF_PY" af3_pipeline/analyze_prolif_interfaces.py \
    --structure-manifest "$CROSS/structure_manifest.tsv" \
    --final-candidates "$FINAL" --out "$CROSS"
  "$SCORE_PY" af3_pipeline/summarize_crossmodel_prolif.py \
    --unique "$CROSS/unique_residue_interactions.tsv" \
    --out "$CROSS/crossmodel_nonvdw.tsv"
  "$SCORE_PY" af3_pipeline/summarize_final20_protenix_multiseed.py \
    --final "$FINAL" --af3-models "$AF3_MODELS" \
    --protenix-root "$PROTENIX/out" --reference references/1OVA.cif \
    --shared-nonvdw "$CROSS/crossmodel_nonvdw.tsv" --out "$SUMMARY"
  touch "$MARKER"
fi

bash af3_pipeline/run_e177_after_e174_local.sh
python af3_pipeline/submit_af3_volc_if_capacity.py \
  --experiment experiments/e178_c4_073_af3ensemble_ai_af3_seed1 \
  --job-number-start 618
exec bash af3_pipeline/run_e178_after_e177_local.sh
