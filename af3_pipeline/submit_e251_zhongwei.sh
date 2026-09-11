#!/usr/bin/env bash
# Submit epistasis-ranked, cross-ensemble tied-AF2 designs to Zhongwei.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
EXPERIMENT=experiments/e251_c4_epistasis_seeded_tied_design
REMOTE_PROJECT=/data/home/scwb671/run/wky/ova_p3_tetramer_design
HOST=ssh.cn-zhongwei-1.paracloud.com
USER_NAME='scwb671@NMCC-N46A8'
PORT=2222
CONTROL=/tmp/ova-zhongwei-control.sock
JOB_BASE=1024
cd "$PROJECT"

[[ -s "$EXPERIMENT/target_manifest.tsv" ]] || {
  echo "missing frozen E251 target manifest" >&2
  exit 2
}
if [[ -s "$EXPERIMENT/paracloud_jobs.tsv" ]] && \
   [[ "$(tail -n +2 "$EXPERIMENT/paracloud_jobs.tsv" | wc -l)" -ge 8 ]]; then
  echo "E251 already has eight recorded Slurm submissions"
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
      --export=ALL,TASK_INDEX='$index',OVA_IPTM90_EXPERIMENT='$REMOTE_PROJECT/$EXPERIMENT',OVA_IPTM90_RUN_PREFIX='run_epi1',OVA_MUTATION_BUDGET='23',OVA_MIN_SURFACE_MUTATION_FRACTION='0.80',OVA_STATE1_WEIGHT_R0='11',OVA_STATE2_WEIGHT_R0='11',OVA_STATE1_WEIGHT_R1='12',OVA_STATE2_WEIGHT_R1='10',OVA_LOGITS_ITERS='28',OVA_SOFT_ITERS='14',OVA_HARD_ITERS='8',OVA_RECYCLES='1',OVA_W_BUDGET='26',OVA_W_SPARSE='1.8',OVA_W_BURIED_MUTATION='14',OVA_WORST_STATE_GRADIENT_BOOST='1.25',OVA_SEED_OFFSET='100000' \
      af3_pipeline/run_iptm90_tied_design_paracloud.sh")
  printf '%s\t%s\t%s\t%s\t%s\n' \
    "$index" "$task_name" "$job_id" epistasis_seeded_crossensemble_tied_af2 SUBMITTED >> "$record"
  echo "$task_name -> Slurm $job_id"
done
touch "$EXPERIMENT/PARACLOUD_JOBS_SUBMITTED"

job_ids=$(tail -n +2 "$record" | cut -f3 | paste -sd, -)
ssh -S "$CONTROL" -o "User=$USER_NAME" -p "$PORT" "$HOST" \
  "cd '$REMOTE_PROJECT' && setsid -f bash af3_pipeline/archive_zhongwei_tied_batch.sh \
    '$EXPERIMENT' 'run_epi1_*' '$job_ids' e251_epi1_slim_results.tgz \
    > '$EXPERIMENT/archive.log' 2>&1 < /dev/null"
echo "E251 archive watcher started for $job_ids"
