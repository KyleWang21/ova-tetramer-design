#!/usr/bin/env bash
# Freeze candidates selected by measured dual-state tied-AF2 evidence.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
SOURCE="$PROJECT/experiments/e235_c4_crossmodel_active_tied_design/candidate_pool_novel/hamming_ai_ranking_current_weighted.tsv"
SELECTED="$PROJECT/experiments/e259_c4_direct_tied_af2_rescue_library"
E259="$PROJECT/experiments/e259_c4_direct_tied_af2_rescue_af3_seed1"
cd "$PROJECT"
if [[ -f "$E259/CANDIDATES_FROZEN" ]]; then
  echo "E259 already frozen"
  exit 0
fi

"$PY" af3_pipeline/select_direct_tied_af2_rescues.py \
  --ranking "$SOURCE" --out "$SELECTED" --top 8 \
  --max-mutations 24 --min-surface-fraction 0.8
"$PY" af3_pipeline/prepare_ova_body_af3.py \
  --fasta "$SELECTED/direct_tied_af2_rescues.fasta" \
  --ova-msa experiments/e202_c4_ai11v2_ai_af3_seed1/msas/ova_body_c4_01.a3m \
  --out "$E259" --top 8 --seeds 1 --seed-start 1 --seed-offset 0 \
  --shards 8 --n-chains 4
touch "$E259/CANDIDATES_FROZEN"
echo "E259 direct tied-AF2 rescue candidates frozen"
