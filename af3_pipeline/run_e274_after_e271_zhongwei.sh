#!/usr/bin/env bash
# Queue the E261 soft-basin rescue after E271 has entered Slurm.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
SOURCE="$PROJECT/experiments/e261_c4_e243_crossbackbone_continuation_tied_design"
PREVIOUS="$PROJECT/experiments/e271_c4_multibackbone_exact_mpnn_score"
EXPERIMENT="$PROJECT/experiments/e274_c4_e261_softbasin_stabilization_tied_design"
PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"
while [[ ! -f "$SOURCE/REMOTE_ARCHIVE_RECOVERED" || \
         ! -f "$PREVIOUS/PARACLOUD_JOBS_SUBMITTED" ]]; do sleep 30; done
"$PY" af3_pipeline/build_e268_e251_softbasin_stabilization.py \
  --source "$SOURCE" --out "$EXPERIMENT" \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --reference-name D0-P3-13R --min-soft 0.65 \
  --max-mutations 20 --min-surface-fraction 0.8 \
  --run-glob run_cb1 --anchor-prefix CB1SOFT \
  --source-label E261-BEST-CROSSBACKBONE-NEARDISCRETE-BASIN
bash af3_pipeline/submit_e274_zhongwei.sh
exec bash af3_pipeline/recover_zhongwei_archive.sh \
  experiments/e274_c4_e261_softbasin_stabilization_tied_design/e274_cbsoft_results.tgz \
  experiments/e274_c4_e261_softbasin_stabilization_tied_design/REMOTE_ARCHIVE_RECOVERED
