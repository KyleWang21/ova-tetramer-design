#!/usr/bin/env bash
# Keep Zhongwei useful after E265 by running robust exact-sequence MPNN scoring.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
PREVIOUS="$PROJECT/experiments/e265_c4_crossbackbone_consensus_tied_design"
EXPERIMENT="$PROJECT/experiments/e271_c4_multibackbone_exact_mpnn_score"
cd "$PROJECT"
# E265's eight jobs must enter Slurm first, but E271 can then queue directly
# behind them rather than waiting for their archive round trip.  This preserves
# priority while avoiding a GPU-idle gap if local recovery is delayed.
while [[ ! -f "$PREVIOUS/PARACLOUD_JOBS_SUBMITTED" ]]; do sleep 30; done
bash af3_pipeline/submit_e271_zhongwei.sh
exec bash af3_pipeline/recover_zhongwei_archive.sh \
  experiments/e271_c4_multibackbone_exact_mpnn_score/e271_exact_scores.tgz \
  experiments/e271_c4_multibackbone_exact_mpnn_score/REMOTE_ARCHIVE_RECOVERED
