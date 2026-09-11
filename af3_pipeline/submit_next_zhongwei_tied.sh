#!/usr/bin/env bash
# Start the next Zhongwei eight-state tied-AF2 batch after E318 freezes.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
SOURCE="$PROJECT/experiments/e318_c4_e311_diverse_tied_design"
EXPERIMENT_REL="experiments/e321_c4_e318_next_tied_design"
EXPERIMENT="$PROJECT/$EXPERIMENT_REL"
REMOTE_PROJECT=/data/home/scwb671/run/wky/ova_p3_tetramer_design
HOST=ssh.cn-zhongwei-1.paracloud.com
USER_NAME='scwb671@NMCC-N46A8'
PORT=2222
CONTROL=/tmp/ova-zhongwei-control.sock
JOB_BASE=2500
PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"

if [[ -f "$EXPERIMENT/PARACLOUD_JOBS_SUBMITTED" ]]; then
  echo "E321 already submitted"
  exit 0
fi

# Start from seven completed E318 shards when they are available so the seven
# released A800s do not sit idle behind the slower eighth shard.  This partial
# pool is used only to choose E321 generation anchors; E319 AF3 still waits for
# the complete eight-shard E318 pool and its own CANDIDATES_FROZEN marker.
PARTIAL_POOL="$SOURCE/candidate_pool_partial7"
while true; do
  if [[ -f "$SOURCE/CANDIDATES_FROZEN" ]]; then
    POOL_ARGS=()
    break
  fi
  partial_tables=$(find "$SOURCE" -path '*/run_e318_*/joint_*/trajectory.tsv' 2>/dev/null | wc -l || true)
  if [[ "$partial_tables" -ge 7 && -s "$PARTIAL_POOL/hamming_ai_ranking_current_weighted.tsv" ]]; then
    POOL_ARGS=(--pool-dir "$PARTIAL_POOL")
    break
  fi
  sleep 30
done
if [[ ! -f "$EXPERIMENT/INPUT_FROZEN" ]]; then
  "$PY" af3_pipeline/prepare_next_zhongwei_tied.py \
    --source "$SOURCE" --out "$EXPERIMENT" --anchor-prefix E321 \
    "${POOL_ARGS[@]}" \
    --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
    --reference-name D0-P3-13R \
    --registry experiments/current_accepted_af3_registry/accepted_af3_registry.tsv
fi

rsync -az -e "ssh -S $CONTROL -o User=$USER_NAME -p $PORT" \
  "$EXPERIMENT/" "$HOST:$REMOTE_PROJECT/$EXPERIMENT_REL/"
rsync -az -e "ssh -S $CONTROL -o User=$USER_NAME -p $PORT" \
  af3_pipeline/af2_multistate_tied_lowmut_design.py \
  af3_pipeline/run_next_zhongwei_tied_worker.sh \
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
      --export=ALL,TASK_INDEX='$index',BATCH_EXPERIMENT='$REMOTE_PROJECT/experiments/e321_c4_e318_next_tied_design' \
      af3_pipeline/run_next_zhongwei_tied_worker.sh")
  printf '%s\t%s\t%s\t%s\t%s\n' "$index" "$task_name" "$job_id" \
    next_eight_state_tied_af2 SUBMITTED >> "$record"
  echo "$task_name -> Slurm $job_id"
done
touch "$EXPERIMENT/PARACLOUD_JOBS_SUBMITTED"
job_ids=$(tail -n +2 "$record" | cut -f3 | paste -sd, -)
ssh -S "$CONTROL" -o "User=$USER_NAME" -p "$PORT" "$HOST" \
  "cd '$REMOTE_PROJECT' && setsid -f bash af3_pipeline/archive_zhongwei_tied_batch.sh \
    'experiments/e321_c4_e318_next_tied_design' 'run_next_*' '$job_ids' e321_next_results.tgz \
    > 'experiments/e321_c4_e318_next_tied_design/archive.log' 2>&1 < /dev/null"
echo "E321 archive watcher started for $job_ids"
