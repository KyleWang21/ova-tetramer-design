#!/usr/bin/env bash
# Materialize generic v1.1 inputs whenever a new 25-model AF3 hit is registered.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
MARKER="$PROJECT/experiments/current_accepted_af3_registry/NEW_ACCEPTED_AF3_PENDING_V1"
LOCK="$PROJECT/experiments/VOLC_PRIORITY_BATCH.lock"
DESIGN_BARRIER="$PROJECT/experiments/e314_c4_e307_epistasis_reflow_af3_seed1/TOP20_MULTISEED_VALIDATION_COMPLETE"
cd "$PROJECT"
while true; do
  # Keep the Fire queue on the high-value E303 -> E307 -> E314 AF3 label loop,
  # including independent validation of all twenty strict E314 hits, first.
  # The v1.1 backlog is retained by MARKER and resumes after that chain.
  while [[ ! -f "$MARKER" || ! -f "$DESIGN_BARRIER" ]]; do sleep 30; done
  while [[ -e "$LOCK" ]]; do sleep 30; done
  printf '%s\t%s\n' "$$" generic_pending_v1 > "$LOCK"
  cleanup() {
    if [[ -f "$LOCK" ]] && grep -q "^$$" "$LOCK"; then rm -f "$LOCK"; fi
  }
  trap cleanup EXIT
  if bash af3_pipeline/run_pending_accepted_v1_staging.sh && \
     bash af3_pipeline/run_pending_v1_gpu_screens.sh; then
    rm -f "$MARKER"
  else
    echo "pending v1.1 run failed; retaining marker for retry" >&2
    cleanup
    trap - EXIT
    sleep 60
    continue
  fi
  cleanup
  trap - EXIT
done
