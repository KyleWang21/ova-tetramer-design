#!/usr/bin/env bash
# Continue the best fully discrete E268 hit across four independent AF3 backbones.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
SOURCE="$PROJECT/experiments/e268_c4_e251_softbasin_stabilization_tied_design"
EXPERIMENT="$PROJECT/experiments/e277_c4_e268_hardhit_crossbackbone_tied_design"
PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"

while [[ ! -f "$SOURCE/REMOTE_ARCHIVE_RECOVERED" ]]; do sleep 30; done
"$PY" af3_pipeline/build_e261_e243_crossbackbone_continuation.py \
  --source "$SOURCE" --out "$EXPERIMENT" \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --reference-name D0-P3-13R --max-mutations 20 --min-surface-fraction 0.8 \
  --run-glob 'run_softstab_*' --anchor-prefix E268X \
  --source-label E268-BEST-FULLY-DISCRETE-HARD
bash af3_pipeline/submit_e277_zhongwei.sh
exec bash af3_pipeline/recover_zhongwei_archive.sh \
  experiments/e277_c4_e268_hardhit_crossbackbone_tied_design/e277_e268x_results.tgz \
  experiments/e277_c4_e268_hardhit_crossbackbone_tied_design/REMOTE_ARCHIVE_RECOVERED
