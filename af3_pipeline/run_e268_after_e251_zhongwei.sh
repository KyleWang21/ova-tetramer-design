#!/usr/bin/env bash
# Queue soft-basin stabilization behind the already-submitted E261 jobs.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
SOURCE="$PROJECT/experiments/e251_c4_epistasis_seeded_tied_design"
EXPERIMENT="$PROJECT/experiments/e268_c4_e251_softbasin_stabilization_tied_design"
PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"
while [[ ! -f "$SOURCE/REMOTE_ARCHIVE_RECOVERED" ]]; do sleep 30; done
"$PY" af3_pipeline/build_e268_e251_softbasin_stabilization.py \
  --source "$SOURCE" --out "$EXPERIMENT" \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --reference-name D0-P3-13R --min-soft 0.65 \
  --max-mutations 20 --min-surface-fraction 0.8
bash af3_pipeline/submit_e268_zhongwei.sh
exec bash af3_pipeline/recover_zhongwei_archive.sh \
  experiments/e268_c4_e251_softbasin_stabilization_tied_design/e268_softstab_results.tgz \
  experiments/e268_c4_e251_softbasin_stabilization_tied_design/REMOTE_ARCHIVE_RECOVERED
