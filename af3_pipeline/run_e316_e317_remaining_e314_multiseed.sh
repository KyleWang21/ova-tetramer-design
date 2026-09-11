#!/usr/bin/env bash
# Validate E314 strict seed-1 ranks 5--20 after the first Top4, in two
# fail-fast batches of at most eight concurrent Fire GPUs.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
SOURCE="$PROJECT/experiments/e314_c4_e307_epistasis_reflow_af3_seed1"
TOP4="$PROJECT/experiments/e315_e314_top4_af3_multiseed"
NEXT8="$PROJECT/experiments/e316_e314_rank05_12_af3_multiseed"
LAST8="$PROJECT/experiments/e317_e314_rank13_20_af3_multiseed"
cd "$PROJECT"

while [[ ! -f "$TOP4/MULTISEED_VALIDATION_COMPLETE" ]]; do sleep 30; done
bash af3_pipeline/run_remaining_seed1_hits_multiseed_batched_volc.sh \
  "$SOURCE" "$NEXT8" 3000 4 8
bash af3_pipeline/run_remaining_seed1_hits_multiseed_batched_volc.sh \
  "$SOURCE" "$LAST8" 3032 12 8
touch "$SOURCE/TOP20_MULTISEED_VALIDATION_COMPLETE"
echo "E314 strict seed-1 Top20 multiseed validation complete"
