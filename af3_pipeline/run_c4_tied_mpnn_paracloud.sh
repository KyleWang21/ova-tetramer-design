#!/usr/bin/env bash
#SBATCH --partition=vip_gpu_a800_scwb671
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=04:00:00

set -euo pipefail
TASK_INDEX="${TASK_INDEX:?submit with --export=ALL,TASK_INDEX=0..7}"
ROOT=/data/home/scwb671/run/wky
export OVA_TETRA_ROOT="$ROOT/ova_p3_tetramer_design"
export BINDCRAFT_MOD="$ROOT/antibody-design/pipeline_stage6_maskcap_bindcraft"
export OVA_MPNN_EXPERIMENT="${OVA_MPNN_EXPERIMENT:-$OVA_TETRA_ROOT/experiments/e166_c4_tied_proteinmpnn}"
mkdir -p "$OVA_MPNN_EXPERIMENT/logs"
echo "$(date -u '+%F %T UTC') tied ProteinMPNN task=$TASK_INDEX host=$(hostname)"
nvidia-smi -L
worker="${OVA_MPNN_WORKER_SCRIPT:-run_c4_tied_mpnn_worker.sh}"
exec bash "$OVA_TETRA_ROOT/af3_pipeline/$worker" "$TASK_INDEX"
