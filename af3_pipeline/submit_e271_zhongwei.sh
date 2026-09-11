#!/usr/bin/env bash
# Submit eight exact-sequence ProteinMPNN backbone scoring jobs on Zhongwei.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
EXPERIMENT=experiments/e271_c4_multibackbone_exact_mpnn_score
REMOTE_PROJECT=/data/home/scwb671/run/wky/ova_p3_tetramer_design
HOST=ssh.cn-zhongwei-1.paracloud.com
USER_NAME='scwb671@NMCC-N46A8'
PORT=2222
CONTROL=/tmp/ova-zhongwei-control.sock
JOB_BASE=1256
cd "$PROJECT"
[[ -s "$EXPERIMENT/candidates.tsv" && -s "$EXPERIMENT/backbone_manifest.tsv" && \
   -s "$EXPERIMENT/score_positions.txt" && -f "$EXPERIMENT/INPUT_FROZEN" ]] || {
  echo "missing frozen E271 inputs" >&2
  exit 2
}
if [[ -s "$EXPERIMENT/paracloud_jobs.tsv" ]] && \
   [[ "$(tail -n +2 "$EXPERIMENT/paracloud_jobs.tsv" | wc -l)" -ge 8 ]]; then
  echo "E271 already has eight recorded submissions"
  exit 0
fi
rsync -az -e "ssh -S $CONTROL -o User=$USER_NAME -p $PORT" \
  "$EXPERIMENT/" "$HOST:$REMOTE_PROJECT/$EXPERIMENT/"
rsync -az -e "ssh -S $CONTROL -o User=$USER_NAME -p $PORT" \
  af3_pipeline/score_multibackbone_exact_mpnn.py \
  af3_pipeline/run_e271_exact_mpnn_score_worker.sh \
  af3_pipeline/run_c4_tied_mpnn_paracloud.sh \
  af3_pipeline/archive_zhongwei_exact_mpnn_scores.sh \
  "$HOST:$REMOTE_PROJECT/af3_pipeline/"
date_tag=$(TZ=Asia/Shanghai date +%Y%m%d)
record="$EXPERIMENT/paracloud_jobs.tsv"
printf 'task_index\ttask_name\tslurm_job_id\tmethod\tstatus\n' > "$record"
for index in $(seq 0 7); do
  number=$((JOB_BASE + index))
  task_name=$(printf 'ky-%s-%03d' "$date_tag" "$number")
  job_id=$(ssh -S "$CONTROL" -o "User=$USER_NAME" -p "$PORT" "$HOST" \
    "cd '$REMOTE_PROJECT' && sbatch --parsable --job-name='$task_name' \
      --output='$REMOTE_PROJECT/$EXPERIMENT/logs/slurm-%j.out' \
      --error='$REMOTE_PROJECT/$EXPERIMENT/logs/slurm-%j.err' \
      --export=ALL,TASK_INDEX='$index',OVA_MPNN_EXPERIMENT='$REMOTE_PROJECT/$EXPERIMENT',OVA_MPNN_WORKER_SCRIPT='run_e271_exact_mpnn_score_worker.sh' \
      af3_pipeline/run_c4_tied_mpnn_paracloud.sh")
  printf '%s\t%s\t%s\t%s\t%s\n' "$index" "$task_name" "$job_id" \
    exact_tied_c4_ProteinMPNN_eight_backbone_NLL SUBMITTED >> "$record"
  echo "$task_name -> Slurm $job_id"
done
touch "$EXPERIMENT/PARACLOUD_JOBS_SUBMITTED"
job_ids=$(tail -n +2 "$record" | cut -f3 | paste -sd, -)
ssh -S "$CONTROL" -o "User=$USER_NAME" -p "$PORT" "$HOST" \
  "cd '$REMOTE_PROJECT' && setsid -f bash af3_pipeline/archive_zhongwei_exact_mpnn_scores.sh \
    '$EXPERIMENT' '$job_ids' e271_exact_scores.tgz \
    > '$EXPERIMENT/archive.log' 2>&1 < /dev/null"
echo "E271 archive watcher started for $job_ids"
