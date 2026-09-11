#!/usr/bin/env bash
# Run one staged candidate's ten-seed OpenDDE diagnostic in the verified
# Zhongwei environment. This one-shot legacy runner is not part of v1.1
# acceptance and is no longer called by the watcher or staging pipeline.
set -euo pipefail

SCREEN="${1:?usage: run_one_pending_v1_opendde_zhongwei.sh SCREEN_DIR}"
PROJECT=/root/400083/ova_p3_tetramer_design
REMOTE_PROJECT=/data/home/scwb671/run/wky/ova_p3_tetramer_design
HOST=ssh.cn-zhongwei-1.paracloud.com
USER_NAME='scwb671@NMCC-N46A8'
PORT=2222
CONTROL=/tmp/ova-zhongwei-control.sock
PY_AF3=/root/400083/alphafold3/repo/.venv/bin/python
PY_FILTER=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
DOCKQ="$PROJECT/experiments/e153_dockq/env/bin/DockQ"
cd "$PROJECT"

[[ -f "$SCREEN/V1_CORE_COMPLETE" && -f "$SCREEN/v1_job_base.tsv" ]] || exit 2
relative=${SCREEN#"$PROJECT/"}
remote_screen="$REMOTE_PROJECT/$relative"
base=$(awk -F '\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="job_base")c=i;next} NR==2{print $c}' "$SCREEN/v1_job_base.tsv")
number=$((base + 16))
record="$SCREEN/opendde/paracloud_jobs.tsv"

if [[ ! -s "$record" ]]; then
  ssh -S "$CONTROL" -o "User=$USER_NAME" -p "$PORT" "$HOST" \
    "mkdir -p '$remote_screen/opendde' '$REMOTE_PROJECT/af3_pipeline'"
  rsync -az -e "ssh -S $CONTROL -o User=$USER_NAME -p $PORT" \
    "$SCREEN/opendde/" "$HOST:$remote_screen/opendde/"
  rsync -az -e "ssh -S $CONTROL -o User=$USER_NAME -p $PORT" \
    af3_pipeline/run_opendde_c4_shard.sh af3_pipeline/run_opendde_paracloud.sh \
    "$HOST:$REMOTE_PROJECT/af3_pipeline/"
  date_tag=$(TZ=Asia/Shanghai date +%Y%m%d)
  task_name=$(printf 'ky-%s-%03d' "$date_tag" "$number")
  job_id=$(ssh -S "$CONTROL" -o "User=$USER_NAME" -p "$PORT" "$HOST" \
    "cd '$REMOTE_PROJECT' && sbatch --parsable --job-name='$task_name' \
      --export=ALL,EXP='$remote_screen/opendde',SHARD='0' \
      af3_pipeline/run_opendde_paracloud.sh")
  printf 'task_name\tslurm_job_id\tmethod\tstatus\n%s\t%s\t%s\t%s\n' \
    "$task_name" "$job_id" opendde_10seed SUBMITTED > "$record"
  echo "$task_name -> Slurm $job_id"
fi
job_id=$(awk -F '\t' 'NR==2{print $2}' "$record")
while ssh -S "$CONTROL" -o "User=$USER_NAME" -p "$PORT" "$HOST" \
  "squeue -h -j '$job_id'" | grep -q .; do sleep 30; done
state=$(ssh -S "$CONTROL" -o "User=$USER_NAME" -p "$PORT" "$HOST" \
  "sacct -j '$job_id' --format=State -n -X | head -1" | xargs)
[[ "$state" == COMPLETED* ]] || { echo "OpenDDE Slurm $job_id failed: $state" >&2; exit 3; }

mkdir -p "$SCREEN/opendde/out" "$SCREEN/opendde/logs"
rsync -az -e "ssh -S $CONTROL -o User=$USER_NAME -p $PORT" \
  "$HOST:$remote_screen/opendde/out/" "$SCREEN/opendde/out/"
rsync -az -e "ssh -S $CONTROL -o User=$USER_NAME -p $PORT" \
  "$HOST:$remote_screen/opendde/logs/" "$SCREEN/opendde/logs/"
complete=$(find "$SCREEN/opendde/out" -name '*_summary_confidence_sample_0.json' | wc -l)
[[ "$complete" -eq 10 ]] || { echo "expected 10 OpenDDE models, found $complete" >&2; exit 4; }

"$PY_AF3" af3_pipeline/summarize_opendde_c4.py \
  --manifest "$SCREEN/opendde/manifest.tsv" --out-root "$SCREEN/opendde/out" \
  --final "$SCREEN/final_candidates.tsv" --af3-models "$SCREEN/af3_per_model_recomputed.tsv" \
  --reference references/1OVA.cif --dockq-bin "$DOCKQ" --out "$SCREEN/opendde_summary"
"$PY_AF3" af3_pipeline/prepare_opendde_prolif.py \
  --final "$SCREEN/final_candidates.tsv" \
  --opendde-summary "$SCREEN/opendde_summary/opendde_candidate_summary.tsv" \
  --out "$SCREEN/opendde_summary/prolif_structure_manifest.tsv"
representatives=$(($(wc -l < "$SCREEN/opendde_summary/prolif_structure_manifest.tsv") - 1))
open_args=()
if [[ "$representatives" -gt 0 ]]; then
  "$PY_FILTER" af3_pipeline/analyze_prolif_interfaces.py \
    --structure-manifest "$SCREEN/opendde_summary/prolif_structure_manifest.tsv" \
    --final-candidates "$SCREEN/final_candidates.tsv" --out "$SCREEN/opendde_summary/prolif_results"
  open_args=(--opendde-unique "$SCREEN/opendde_summary/prolif_results/unique_residue_interactions.tsv")
fi
"$PY_AF3" af3_pipeline/summarize_opendde_prolif.py \
  --final "$SCREEN/final_candidates.tsv" \
  --opendde-summary "$SCREEN/opendde_summary/opendde_candidate_summary.tsv" \
  --opendde-per-seed "$SCREEN/opendde_summary/opendde_per_seed_recomputed.tsv" \
  "${open_args[@]}" --crossmodel-unique "$SCREEN/prolif_crossmodel/unique_residue_interactions.tsv" \
  --out "$SCREEN/opendde_summary/opendde_prolif_summary.tsv"
touch "$SCREEN/OPENDDE_COMPLETE"
bash af3_pipeline/finalize_pending_v1_screens.sh
