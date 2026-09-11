#!/usr/bin/env bash
# Submit cross-backbone continuation of the best discrete E243 sequence.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
EXPERIMENT=experiments/e261_c4_e243_crossbackbone_continuation_tied_design
REMOTE_PROJECT=/data/home/scwb671/run/wky/ova_p3_tetramer_design
HOST=ssh.cn-zhongwei-1.paracloud.com
USER_NAME='scwb671@NMCC-N46A8'
PORT=2222
CONTROL=/tmp/ova-zhongwei-control.sock
JOB_BASE=1152
cd "$PROJECT"

[[ -s "$EXPERIMENT/target_manifest.tsv" && -f "$EXPERIMENT/SOURCE_SELECTION_FROZEN" ]] || {
  echo "missing frozen E261 target manifest" >&2
  exit 2
}
if [[ -s "$EXPERIMENT/paracloud_jobs.tsv" ]] && \
   [[ "$(tail -n +2 "$EXPERIMENT/paracloud_jobs.tsv" | wc -l)" -ge 8 ]]; then
  echo "E261 already has eight recorded Slurm submissions"
  exit 0
fi

rsync -az -e "ssh -S $CONTROL -o User=$USER_NAME -p $PORT" \
  "$EXPERIMENT/" "$HOST:$REMOTE_PROJECT/$EXPERIMENT/"
rsync -az -e "ssh -S $CONTROL -o User=$USER_NAME -p $PORT" \
  af3_pipeline/af2_multistate_tied_lowmut_design.py \
  af3_pipeline/run_iptm90_tied_design_worker.sh \
  af3_pipeline/run_iptm90_tied_design_paracloud.sh \
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
      --export=ALL,TASK_INDEX='$index',OVA_IPTM90_EXPERIMENT='$REMOTE_PROJECT/$EXPERIMENT',OVA_IPTM90_RUN_PREFIX='run_cb1',OVA_MUTATION_BUDGET='22',OVA_MIN_SURFACE_MUTATION_FRACTION='0.80',OVA_STATE1_WEIGHT_R0='11',OVA_STATE2_WEIGHT_R0='11',OVA_STATE1_WEIGHT_R1='13',OVA_STATE2_WEIGHT_R1='9',OVA_LOGITS_ITERS='32',OVA_SOFT_ITERS='16',OVA_HARD_ITERS='8',OVA_RECYCLES='1',OVA_W_BUDGET='28',OVA_W_SPARSE='1.8',OVA_W_BURIED_MUTATION='14',OVA_WORST_STATE_GRADIENT_BOOST='1.75',OVA_SEED_OFFSET='110000' \
      af3_pipeline/run_iptm90_tied_design_paracloud.sh")
  printf '%s\t%s\t%s\t%s\t%s\n' \
    "$index" "$task_name" "$job_id" e243_best_crossbackbone_continuation_tied_af2 SUBMITTED >> "$record"
  echo "$task_name -> Slurm $job_id"
done
touch "$EXPERIMENT/PARACLOUD_JOBS_SUBMITTED"

job_ids=$(tail -n +2 "$record" | cut -f3 | paste -sd, -)
ssh -S "$CONTROL" -o "User=$USER_NAME" -p "$PORT" "$HOST" \
  "cd '$REMOTE_PROJECT' && setsid -f bash af3_pipeline/archive_zhongwei_tied_batch.sh \
    '$EXPERIMENT' 'run_cb1_*' '$job_ids' e261_cb1_results.tgz \
    > '$EXPERIMENT/archive.log' 2>&1 < /dev/null"
echo "E261 archive watcher started for $job_ids"
