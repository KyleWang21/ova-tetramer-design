#!/usr/bin/env bash
# Prepare and launch E294 as soon as the complete E286 archive is recovered.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
cd "$PROJECT"
while [[ ! -f experiments/e286_c4_chainmapped_bottomk_tied_design/REMOTE_ARCHIVE_RECOVERED ]]; do sleep 30; done
python af3_pipeline/prepare_e294_best_ensemble_from_e286.py
bash af3_pipeline/submit_e294_zhongwei.sh
setsid -f bash af3_pipeline/recover_zhongwei_archive.sh \
  experiments/e294_c4_detbest_ensemble_bottomk_tied_design/e294_detbest_results.tgz \
  experiments/e294_c4_detbest_ensemble_bottomk_tied_design/REMOTE_ARCHIVE_RECOVERED \
  > logs/e294_recovery.log 2>&1
setsid -f bash af3_pipeline/finalize_e294_generation.sh > logs/e294_finalizer.log 2>&1
echo "E294 submitted and recovery/finalizer watchers started"
