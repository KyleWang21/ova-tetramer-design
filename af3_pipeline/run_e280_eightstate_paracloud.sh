#!/usr/bin/env bash
#SBATCH --partition=vip_gpu_a800_scwb671
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00

set -euo pipefail
TASK_INDEX="${TASK_INDEX:?submit with TASK_INDEX=0..7}"
ROOT=/data/home/scwb671/run/wky
export OVA_TETRA_ROOT="$ROOT/ova_p3_tetramer_design"
export BINDCRAFT_MOD="$ROOT/antibody-design/pipeline_stage6_maskcap_bindcraft"
export OVA_E280_EXPERIMENT="${OVA_E280_EXPERIMENT:-$OVA_TETRA_ROOT/experiments/e280_c4_eightstate_multibackbone_tied_design}"
echo "$(date -u '+%F %T UTC') E280 eight-state task=$TASK_INDEX host=$(hostname)"
nvidia-smi -L
exec bash "$OVA_TETRA_ROOT/af3_pipeline/run_e280_eightstate_worker.sh" "$TASK_INDEX"
