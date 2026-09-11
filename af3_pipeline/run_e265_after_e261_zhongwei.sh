#!/usr/bin/env bash
# Keep Zhongwei full after E261 using a <=20-mutation cross-backbone consensus.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
SOURCE="$PROJECT/experiments/e261_c4_e243_crossbackbone_continuation_tied_design"
EXPERIMENT="$PROJECT/experiments/e265_c4_crossbackbone_consensus_tied_design"
PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"
while [[ ! -f "$SOURCE/REMOTE_ARCHIVE_RECOVERED" ]]; do sleep 30; done
"$PY" af3_pipeline/build_e265_crossbackbone_consensus_continuation.py \
  --source "$SOURCE" --out "$EXPERIMENT" \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --reference-name D0-P3-13R --max-source-mutations 24 \
  --max-anchor-mutations 20 --min-anchor-mutations 12 \
  --min-surface-fraction 0.8
bash af3_pipeline/submit_e265_zhongwei.sh
exec bash af3_pipeline/recover_zhongwei_archive.sh \
  experiments/e265_c4_crossbackbone_consensus_tied_design/e265_cons2_results.tgz \
  experiments/e265_c4_crossbackbone_consensus_tied_design/REMOTE_ARCHIVE_RECOVERED
