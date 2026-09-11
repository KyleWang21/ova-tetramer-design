#!/usr/bin/env bash
# Summarize Zhongwei Protenix-v2 and stoichiometry outputs after both arrays.
#SBATCH --partition=vip_gpu_a800_scwb671
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
set -euo pipefail
PROJECT="${ZH_PROJECT:-/data/home/scwb671/run/wky/ova_p3_tetramer_design}"
ROOT="$PROJECT/experiments/current_pending_v1_screen"
PY="/data/home/scwb671/run/wky/alphafold3/repo/.venv/bin/python"
REFERENCE="${ZH_REFERENCE:-1OVA.cif}"
cd "$PROJECT"
while IFS=$'\t' read -r candidate directory; do
  [[ "$candidate" == candidate ]] && continue
  screen="$ROOT/$directory"
  pro_out="$screen/zh_protenix_out"
  sto_out="$screen/zh_stoich_out"
  pro_summary="$screen/protenix_summary_zhongwei"
  sto_summary="$screen/stoichiometry_summary_zhongwei"
  mkdir -p "$pro_summary" "$sto_summary"
  if [[ ! -s "$pro_summary/protenix_candidate_summary.tsv" ]]; then
    "$PY" af3_pipeline/summarize_final20_protenix_multiseed.py \
      --final "$screen/final_candidates.tsv" \
      --af3-models "$screen/af3_per_model_recomputed.tsv" \
      --protenix-root "$pro_out" --reference "$REFERENCE" \
      --out "$pro_summary"
  fi
  # score_outputs expects a candidate-local manifest and out_s{shard}/system layout.
  if [[ ! -s "$sto_out/model_scores.tsv" ]]; then
    "$PY" af3_pipeline/score_outputs.py --experiment "$sto_out"
  fi
  if [[ ! -s "$sto_summary/stoichiometry_candidate_summary.tsv" ]]; then
    "$PY" af3_pipeline/summarize_final20_stoichiometry.py \
      --manifest "$sto_out/manifest.tsv" --stoich-models "$sto_out/model_scores.tsv" \
      --c4-models "$screen/af3_per_model_recomputed.tsv" \
      --reference "$REFERENCE" --out "$sto_summary"
  fi
  touch "$screen/ZHONGWEI_GPU_CORES_COMPLETE"
  echo "$candidate complete"
done < "$ROOT/zh_candidates.tsv"
