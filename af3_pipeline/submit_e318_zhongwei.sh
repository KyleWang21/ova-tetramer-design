#!/usr/bin/env bash
# Submit E318's eight independent tied-AF2 trajectories to Zhongwei A800.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
EXPERIMENT=experiments/e318_c4_e311_diverse_tied_design
REMOTE_PROJECT=/data/home/scwb671/run/wky/ova_p3_tetramer_design
HOST=ssh.cn-zhongwei-1.paracloud.com
USER_NAME='scwb671@NMCC-N46A8'
PORT=2222
CONTROL=/tmp/ova-zhongwei-control.sock
JOB_BASE=2300
cd "$PROJECT"
[[ -s "$EXPERIMENT/anchor_manifest.tsv" && -f "$EXPERIMENT/INPUT_FROZEN" ]] || exit 2
if [[ -s "$EXPERIMENT/paracloud_jobs.tsv" ]] && \
   [[ "$(tail -n +2 "$EXPERIMENT/paracloud_jobs.tsv" | wc -l)" -ge 8 ]]; then
  echo "E318 already has eight recorded submissions"
  exit 0
fi
rsync -az -e "ssh -S $CONTROL -o User=$USER_NAME -p $PORT" \
  "$EXPERIMENT/" "$HOST:$REMOTE_PROJECT/$EXPERIMENT/"
rsync -az -e "ssh -S $CONTROL -o User=$USER_NAME -p $PORT" \
  af3_pipeline/af2_multistate_tied_lowmut_design.py \
  af3_pipeline/run_e318_diverse_worker_paracloud.sh \
  af3_pipeline/archive_zhongwei_tied_batch.sh \
  "$HOST:$REMOTE_PROJECT/af3_pipeline/"
date_tag=$(TZ=Asia/Shanghai date +%Y%m%d)
record="$EXPERIMENT/paracloud_jobs.tsv"
printf 'task_index\ttask_name\tslurm_job_id\tmethod\tstatus\n' > "$record"
for index in $(seq 0 7); do
  number=$((JOB_BASE + index))
  task_name=$(printf 'ky-%s-%03d' "$date_tag" "$number")
  job_id=$(ssh -S "$CONTROL" -o "User=$USER_NAME" -p "$PORT" "$HOST" \
    "cd '$REMOTE_PROJECT' && sbatch --parsable --job-name='$task_name' \
      --export=ALL,TASK_INDEX='$index',OVA_E318_EXPERIMENT='$REMOTE_PROJECT/$EXPERIMENT' \
      af3_pipeline/run_e318_diverse_worker_paracloud.sh")
  printf '%s\t%s\t%s\t%s\t%s\n' "$index" "$task_name" "$job_id" \
    e318_diverse_eight_state_tied_af2 SUBMITTED >> "$record"
  echo "$task_name -> Slurm $job_id"
done
touch "$EXPERIMENT/PARACLOUD_JOBS_SUBMITTED"
job_ids=$(tail -n +2 "$record" | cut -f3 | paste -sd, -)
ssh -S "$CONTROL" -o "User=$USER_NAME" -p "$PORT" "$HOST" \
  "cd '$REMOTE_PROJECT' && setsid -f bash af3_pipeline/archive_zhongwei_tied_batch.sh \
    '$EXPERIMENT' 'run_e318_*' '$job_ids' e318_diverse_results.tgz \
    > '$EXPERIMENT/archive.log' 2>&1 < /dev/null"
echo "E318 archive watcher started for $job_ids"
