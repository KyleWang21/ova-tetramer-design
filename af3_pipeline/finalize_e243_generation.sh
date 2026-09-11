#!/usr/bin/env bash
# Merge E243 cross-ensemble tied-AF2 trajectories and freeze E244 for AF3.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E234="$PROJECT/experiments/e234_crossmodel_active_learning_v2"
E243="$PROJECT/experiments/e243_c4_crossensemble_active_tied_design"
E244="$PROJECT/experiments/e244_c4_crossensemble_active_ai_af3_seed1"
POOL="$E243/candidate_pool_novel"
PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"

if [[ -f "$E244/CANDIDATES_FROZEN" ]]; then
  echo "E244 already frozen; nothing to do"
  exit 0
fi
while [[ ! -f "$E243/REMOTE_ARCHIVE_RECOVERED" ]]; do sleep 30; done

"$PY" af3_pipeline/collect_iptm90_crossmodel_candidates.py \
  --experiment "$E243" --run-glob 'run_xm3_*' \
  --surrogate-tsv "$E234/surrogate_v2/crossmodel_surrogate_candidates.tsv" \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
  --exclude-glob 'experiments/e*/manifest.tsv' \
  --exclude-glob 'experiments/e*/*/manifest.tsv' \
  --exclude-glob 'experiments/e*/*/*/manifest.tsv' \
  --exclude-glob 'experiments/e193_c4_multihit11_af3ensemble_tied_design/candidate_pool_ai11v*_novel/*.fasta' \
  --exclude-glob 'experiments/e205_c4_276_strictensemble_tied_mpnn/candidate_pool_novel/*.fasta' \
  --exclude-glob 'experiments/e215_c4_073_balanced_tied_design/candidate_pool_novel/*.fasta' \
  --exclude-glob 'experiments/e217_c4_073_highgain_tied_design/candidate_pool_novel/*.fasta' \
  --exclude-glob 'experiments/e219_c4_073_recycle1_tied_design/candidate_pool_novel/*.fasta' \
  --exclude-glob 'experiments/e221_c4_ai11_af3_protenix_crossmodel_tied_design/candidate_pool_*/*.fasta' \
  --exclude-glob 'experiments/e225_c4_073_rosetta_polar_tied_mpnn/candidate_pool_novel/*.fasta' \
  --exclude-glob 'experiments/e235_c4_crossmodel_active_tied_design/candidate_pool_novel/*.fasta' \
  --exclude-glob 'deliverables/**/*.fasta' \
  --exclude-ignore-path "$E244/manifest.tsv" \
  --out "$POOL" --top-af2 256 --top-surrogate 128 --top-mpnn 0 \
  --minimum-distance 2 --min-mutations 1 --max-mutations 24

mapfile -t rankings < <(awk '
  /^rankings=\(/{inside=1; next}
  inside && /^\)/{exit}
  inside {gsub(/^[[:space:]]+|[[:space:]]+$/, ""); if (length) print}
' af3_pipeline/finalize_e217_generation.sh)
ranking_args=()
for file in "${rankings[@]}"; do ranking_args+=(--ranking "$file"); done
RANKED="$POOL/hamming_ai_ranking_current_weighted.tsv"
"$PY" af3_pipeline/score_c4_library_hamming_ai.py \
  "${ranking_args[@]}" --library "$POOL/generation_shortlist.tsv" \
  --out "$RANKED" --shortlist 32 --minimum-distance 3 \
  --allow-short-shortlist --duplicate-label-policy max_models_mean --weighted-iptm

"$PY" af3_pipeline/prepare_ova_body_af3.py \
  --fasta "$POOL/hamming_ai_ranking_current_weighted_shortlist.fasta" \
  --ova-msa experiments/e202_c4_ai11v2_ai_af3_seed1/msas/ova_body_c4_01.a3m \
  --out "$E244" --top 32 --seeds 1 --seed-start 1 --seed-offset 0 \
  --shards 8 --n-chains 4
touch "$E244/CANDIDATES_FROZEN"
echo "E243 cross-ensemble active-learning/tied-AF2 library frozen as E244"
