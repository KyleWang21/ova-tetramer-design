#!/usr/bin/env bash
# Refresh E311 anchors with the newest 25-model AF3 evidence before submission.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E308="$PROJECT/experiments/e308_e307_top4_af3_multiseed"
TARGET="$PROJECT/experiments/e311_c4_current_best_target_ensemble"
E311="$PROJECT/experiments/e311_c4_e294_detbest_continuation_tied_design"
ARCHIVE="$PROJECT/experiments/archive"
PY=/root/400083/alphafold3/repo/.venv/bin/python
cd "$PROJECT"

while [[ ! -f "$E308/MULTISEED_VALIDATION_COMPLETE" ]]; do sleep 30; done
if [[ -f "$E311/PARACLOUD_JOBS_SUBMITTED" ]]; then
  echo "refusing to refresh E311 after remote submission" >&2
  exit 2
fi
if [[ -f "$E311/E308_REFRESH_COMPLETE" ]]; then
  echo "E311 already refreshed after E308"
  exit 0
fi

stamp=$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p "$ARCHIVE"
best=$(awk -F'\t' 'NR==2{print $2}' experiments/current_accepted_af3_registry/accepted_af3_registry.tsv)
if [[ "$best" == ova_body_c4_25_ms25_r01 ]]; then
  mv "$TARGET" "$ARCHIVE/e311_target_pre_e308_refresh_$stamp"
  "$PY" af3_pipeline/prepare_strict_hit_ensemble_tied_design.py \
    --summary "$E308/full25/full25_summary.tsv" \
    --per-model "$E308/full25/full25_per_model.tsv" \
    --system "$best" \
    --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
    --reference-name D0-P3-13R \
    --surface-positions experiments/e283_c4_chainmapped_eightstate_tied_design/surface_positions.txt \
    --pairs AC,AD --backbones 4 --anchor-prefix E311-TARGET \
    --out "$TARGET"
fi

mv "$E311" "$ARCHIVE/e311_pre_e308_refresh_$stamp"
"$PY" af3_pipeline/prepare_e311_current_best_multistate.py \
  --target-ensemble "$TARGET" \
  --seed1-summary experiments/e307_c4_e303_label_reflow_af3_seed1/seed1_strict_final/seed1_strict_summary.tsv \
  --multiseed-summary experiments/e304_e303_top4_af3_multiseed/full25/full25_summary.tsv \
  --multiseed-summary experiments/e305_e303_rank05_12_af3_multiseed/full25/full25_summary.tsv \
  --multiseed-summary experiments/e306_e303_rank13_20_af3_multiseed/full25/full25_summary.tsv \
  --multiseed-summary experiments/e308_e307_top4_af3_multiseed/full25/full25_summary.tsv \
  --multiseed-summary experiments/e300_e299_top4_af3_multiseed/full25/full25_summary.tsv \
  --multiseed-summary experiments/e302_e301_top4_af3_multiseed/full25/full25_summary.tsv \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --reference-name D0-P3-13R \
  --surface-positions experiments/e283_c4_chainmapped_eightstate_tied_design/surface_positions.txt \
  --out "$E311"
"$PY" af3_pipeline/audit_multistate_interface_consistency.py \
  --state-manifest "$E311/state_manifest.tsv" \
  --out "$E311/contact_map_consistency_audit.tsv" \
  --cutoff 5.3 --minimum-jaccard 0.70
touch "$E311/E308_REFRESH_COMPLETE"
echo "E311 refreshed from E308 evidence; best=$best"
