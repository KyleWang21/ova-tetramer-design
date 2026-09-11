#!/usr/bin/env bash
# Corrected full-condition ProteinMPNN resampling on Zhongwei.
# Comma-valued temperature/native-bias settings are intentionally not passed
# through sbatch --export; both workers use their audited 2x3 defaults.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
REMOTE_PROJECT=/data/home/scwb671/run/wky/ova_p3_tetramer_design
HOST=ssh.cn-zhongwei-1.paracloud.com
USER_NAME='scwb671@NMCC-N46A8'
PORT=2222
CONTROL=/tmp/ova-zhongwei-control.sock
cd "$PROJECT"

rsync -az -e "ssh -S $CONTROL -o User=$USER_NAME -p $PORT" \
  af3_pipeline/mpnn_tied_c4_lowmut_design.py \
  af3_pipeline/run_c4_tied_mpnn_crossmodel_worker.sh \
  af3_pipeline/run_c4_073_rosetta_polar_mpnn_worker.sh \
  af3_pipeline/run_c4_tied_mpnn_paracloud.sh \
  af3_pipeline/archive_zhongwei_mpnn_batch.sh \
  "$HOST:$REMOTE_PROJECT/af3_pipeline/"

submit_batch() {
  local experiment="$1" run_prefix="$2" worker="$3" job_base="$4" seed_offset="$5" method="$6" archive="$7"
  local record="$experiment/paracloud_jobs.tsv"
  mkdir -p "$experiment"
  if [[ -s "$record" ]] && [[ "$(tail -n +2 "$record" | wc -l)" -ge 8 ]]; then
    echo "$experiment already has eight recorded Slurm submissions"
    return 0
  fi
  ssh -S "$CONTROL" -o "User=$USER_NAME" -p "$PORT" "$HOST" \
    "mkdir -p '$REMOTE_PROJECT/$experiment'"
  printf 'task_index\ttask_name\tslurm_job_id\tmethod\tconditions\texpected_raw_samples\tstatus\n' > "$record"
  local date_tag
  date_tag=$(TZ=Asia/Shanghai date +%Y%m%d)
  for index in $(seq 0 7); do
    local number task_name job_id
    number=$((job_base + index))
    task_name=$(printf 'ky-%s-%03d' "$date_tag" "$number")
    job_id=$(ssh -S "$CONTROL" -o "User=$USER_NAME" -p "$PORT" "$HOST" \
      "cd '$REMOTE_PROJECT' && sbatch --parsable --job-name='$task_name' \
        --export=ALL,TASK_INDEX='$index',OVA_MPNN_EXPERIMENT='$REMOTE_PROJECT/$experiment',OVA_MPNN_RUN_PREFIX='$run_prefix',OVA_MPNN_WORKER_SCRIPT='$worker',OVA_MPNN_SAMPLES_PER_CONDITION='512',OVA_MPNN_SEED_OFFSET='$seed_offset' \
        af3_pipeline/run_c4_tied_mpnn_paracloud.sh")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$index" "$task_name" "$job_id" "$method" '2_temperatures_x_3_native_biases' 3072 SUBMITTED >> "$record"
    echo "$task_name -> Slurm $job_id"
  done
  local job_ids
  job_ids=$(tail -n +2 "$record" | cut -f3 | paste -sd, -)
  ssh -S "$CONTROL" -o "User=$USER_NAME" -p "$PORT" "$HOST" \
    "cd '$REMOTE_PROJECT' && setsid -f bash af3_pipeline/archive_zhongwei_mpnn_batch.sh \
      '$experiment' '${run_prefix}_*' '$job_ids' '$archive' \
      > '$experiment/archive.log' 2>&1 < /dev/null"
  touch "$experiment/PARACLOUD_JOBS_SUBMITTED"
  echo "$experiment archive watcher started for $job_ids"
}

submit_batch \
  experiments/e246_c4_ai11_crossmodel_tied_mpnn_full \
  run_ai11xm_mpnn2 run_c4_tied_mpnn_crossmodel_worker.sh \
  1008 96000 crossmodel_tied_ProteinMPNN_fullconditions e246_ai11xm_mpnn2_results.tgz

submit_batch \
  experiments/e247_c4_073_rosetta_polar_tied_mpnn_full \
  run_073polar2 run_c4_073_rosetta_polar_mpnn_worker.sh \
  1016 98000 rosetta_polar_tied_ProteinMPNN_fullconditions e247_073polar2_results.tgz
