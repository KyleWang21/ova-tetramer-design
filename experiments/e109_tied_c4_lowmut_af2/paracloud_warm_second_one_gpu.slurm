#!/bin/bash
#SBATCH --partition=vip_gpu_a800_scwb671
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00

set -euo pipefail
PROJECT=/data/home/scwb671/run/wky/ova_p3_tetramer_design
export OVA_TETRA_ROOT="$PROJECT"
export BINDCRAFT_MOD=/data/home/scwb671/run/wky/antibody-design/pipeline_stage6_maskcap_bindcraft
mkdir -p "$PROJECT/experiments/e109_tied_c4_lowmut_af2/logs"
exec bash "$PROJECT/af3_pipeline/run_c4_lowmut_af2_warm_second_worker.sh" "$TASK_INDEX"
