#!/usr/bin/env bash
# Validate E303 strict seed-1 ranks 5--20 in two <=8-candidate batches.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
SOURCE="$PROJECT/experiments/e303_c4_cleanbridge_allele_ai_af3_seed1"
TOP4="$PROJECT/experiments/e304_e303_top4_af3_multiseed"
NEXT8="$PROJECT/experiments/e305_e303_rank05_12_af3_multiseed"
LAST8="$PROJECT/experiments/e306_e303_rank13_20_af3_multiseed"
cd "$PROJECT"

while [[ ! -f "$TOP4/MULTISEED_VALIDATION_COMPLETE" ]]; do sleep 30; done
bash af3_pipeline/run_remaining_seed1_hits_multiseed_batched_volc.sh \
  "$SOURCE" "$NEXT8" 1784 4 8
bash af3_pipeline/run_remaining_seed1_hits_multiseed_batched_volc.sh \
  "$SOURCE" "$LAST8" 1816 12 8
touch "$SOURCE/TOP20_MULTISEED_VALIDATION_COMPLETE"
echo "E303 strict seed-1 Top20 multiseed validation complete"
