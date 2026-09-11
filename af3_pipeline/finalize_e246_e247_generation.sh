#!/usr/bin/env bash
# Merge corrected full-condition cross-model/polar MPNN libraries and freeze E248.
set -euo pipefail
shopt -s nullglob

PROJECT=/root/400083/ova_p3_tetramer_design
E221="$PROJECT/experiments/e221_c4_ai11_af3_protenix_crossmodel_tied_design"
E246="$PROJECT/experiments/e246_c4_ai11_crossmodel_tied_mpnn_full"
E247="$PROJECT/experiments/e247_c4_073_rosetta_polar_tied_mpnn_full"
E248="$PROJECT/experiments/e248_c4_fullcondition_mpnn_joint_ai_af3_seed1"
POOL="$E246/candidate_pool_joint_e247_novel"
PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"

if [[ -f "$E248/CANDIDATES_FROZEN" ]]; then
  echo "E248 already frozen; nothing to do"
  exit 0
fi
while [[ ! -f "$E221/REMOTE_ARCHIVE_RECOVERED" || \
         ! -f "$E246/REMOTE_ARCHIVE_RECOVERED" || \
         ! -f "$E247/REMOTE_ARCHIVE_RECOVERED" ]]; do
  sleep 30
done

mpnn_files=(
  "$E246"/run_ai11xm_mpnn2_*/mpnn_shortlist.tsv
  "$E247"/run_073polar2_*/mpnn_shortlist.tsv
)
[[ "${#mpnn_files[@]}" -eq 16 ]] || {
  echo "expected sixteen recovered full-condition MPNN shortlists, found ${#mpnn_files[@]}" >&2
  exit 2
}
mpnn_args=()
for file in "${mpnn_files[@]}"; do mpnn_args+=(--mpnn-tsv "$file"); done

"$PY" af3_pipeline/collect_iptm90_crossmodel_candidates.py \
  --experiment "$E221" --run-glob 'run_ai11xm1_*' \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
  "${mpnn_args[@]}" \
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
  --exclude-glob 'experiments/e243_c4_crossensemble_active_tied_design/candidate_pool_novel/*.fasta' \
  --exclude-glob 'deliverables/**/*.fasta' \
  --exclude-ignore-path "$E248/manifest.tsv" \
  --out "$POOL" --top-af2 128 --top-surrogate 0 --top-mpnn 512 \
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
  --out "$E248" --top 32 --seeds 1 --seed-start 1 --seed-offset 0 \
  --shards 8 --n-chains 4
touch "$E248/CANDIDATES_FROZEN"
echo "E246/E247 full-condition MPNN libraries frozen as E248"
