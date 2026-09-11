#!/usr/bin/env bash
# Keep Zhongwei's eight A800s occupied by continuing immediately after E294.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
cd "$PROJECT"
while [[ ! -f experiments/e294_c4_detbest_ensemble_bottomk_tied_design/REMOTE_ARCHIVE_RECOVERED ]]; do
  sleep 30
done
while [[ ! -f experiments/e308_e307_top4_af3_multiseed/MULTISEED_VALIDATION_COMPLETE ]]; do
  sleep 30
done
while [[ ! -f experiments/e311_c4_e294_detbest_continuation_tied_design/E308_REFRESH_COMPLETE ]]; do
  sleep 30
done
[[ -f experiments/e311_c4_e294_detbest_continuation_tied_design/INPUT_FROZEN ]]
python af3_pipeline/audit_multistate_interface_consistency.py \
  --state-manifest experiments/e311_c4_e294_detbest_continuation_tied_design/state_manifest.tsv \
  --out experiments/e311_c4_e294_detbest_continuation_tied_design/contact_map_consistency_audit.tsv \
  --cutoff 5.3 --minimum-jaccard 0.70
bash af3_pipeline/submit_e311_zhongwei.sh
setsid -f bash af3_pipeline/recover_zhongwei_archive.sh \
  experiments/e311_c4_e294_detbest_continuation_tied_design/e311_detbest_results.tgz \
  experiments/e311_c4_e294_detbest_continuation_tied_design/REMOTE_ARCHIVE_RECOVERED \
  > logs/e311_recovery.log 2>&1
echo "E311 submitted and recovery watcher started"
