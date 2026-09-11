#!/usr/bin/env bash
# Prepare, submit, recover, and freeze E286 once E283 is completely recovered.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
cd "$PROJECT"
while [[ ! -f experiments/e283_c4_chainmapped_eightstate_tied_design/REMOTE_ARCHIVE_RECOVERED ]]; do
  sleep 30
done
progress_logs=experiments/e283_c4_chainmapped_eightstate_tied_design/live_progress_logs
mkdir -p "$progress_logs"
copied=1
while IFS=$'\t' read -r task_index task_name job_id method status; do
  [[ "$task_index" == "task_index" ]] && continue
  if ! scp -q -o ControlPath=/tmp/ova-zhongwei-control.sock \
      -o 'User=scwb671@NMCC-N46A8' -P 2222 \
      "ssh.cn-zhongwei-1.paracloud.com:/data/run01/scwb671/wky/ova_p3_tetramer_design/slurm-${job_id}.out" \
      "$progress_logs/"; then
    copied=0
  fi
done < experiments/e283_c4_chainmapped_eightstate_tied_design/paracloud_jobs.tsv
if [[ "$copied" -eq 1 ]]; then
  /root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python \
    af3_pipeline/plot_e283_chainmapping_progress.py \
    --experiment experiments/e283_c4_chainmapped_eightstate_tied_design \
    --logs "$progress_logs" \
    --out experiments/e283_c4_chainmapped_eightstate_tied_design/figures || true
fi
python af3_pipeline/prepare_e286_bottomk_from_e283.py
bash af3_pipeline/submit_e286_zhongwei.sh
setsid -f bash af3_pipeline/recover_zhongwei_archive.sh \
  experiments/e286_c4_chainmapped_bottomk_tied_design/e286_bottomk_results.tgz \
  experiments/e286_c4_chainmapped_bottomk_tied_design/REMOTE_ARCHIVE_RECOVERED \
  > logs/e286_recovery.log 2>&1
setsid -f bash af3_pipeline/finalize_e286_generation.sh > logs/e286_finalizer.log 2>&1
echo "E286 submitted and recovery/finalizer watchers started"
