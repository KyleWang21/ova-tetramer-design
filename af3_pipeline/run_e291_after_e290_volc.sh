#!/usr/bin/env bash
# Validate the five E289 strict seed1 hits not included in E290 Top4.
set -euo pipefail
PROJECT=/root/400083/ova_p3_tetramer_design
cd "$PROJECT"
while [[ ! -f experiments/e290_e289_top4_af3_multiseed/MULTISEED_VALIDATION_COMPLETE ]]; do
  sleep 30
done
bash af3_pipeline/run_remaining_seed1_hits_multiseed_batched_volc.sh \
  experiments/e289_c4_e252_clashrepair_ai_af3_seed1 \
  experiments/e291_e289_remaining5_af3_multiseed 1480 4 5
