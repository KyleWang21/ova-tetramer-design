#!/usr/bin/env bash
# Freeze E284 from corrected chain-mapped eight-state tied-AF2 trajectories.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E283="$PROJECT/experiments/e283_c4_chainmapped_eightstate_tied_design"
E284="$PROJECT/experiments/e284_c4_chainmapped_eightstate_af3_seed1"
POOL="$E283/candidate_pool_novel"
PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"
if [[ -f "$E284/CANDIDATES_FROZEN" ]]; then echo "E284 already frozen"; exit 0; fi
while [[ ! -f "$E283/REMOTE_ARCHIVE_RECOVERED" ]]; do sleep 30; done
"$PY" af3_pipeline/collect_iptm90_crossmodel_candidates.py \
  --experiment "$E283" --run-glob 'run_chainmapped_*' \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
  --exclude-glob 'experiments/e*/manifest.tsv' \
  --exclude-glob 'experiments/e*/*/manifest.tsv' \
  --exclude-glob 'experiments/e*/*/*/manifest.tsv' \
  --exclude-glob 'experiments/e*/candidate_pool*/*.fasta' \
  --exclude-glob 'experiments/e*/*/candidate_pool*/*.fasta' \
  --exclude-glob 'deliverables/**/*.fasta' \
  --exclude-ignore-path "$E284/manifest.tsv" \
  --out "$POOL" --top-af2 256 --top-surrogate 0 --top-mpnn 0 \
  --minimum-distance 2 --min-mutations 1 --max-mutations 20
mapfile -t rankings < <(awk '
  /^rankings=\(/{inside=1; next}
  inside && /^\)/{exit}
  inside {gsub(/^[[:space:]]+|[[:space:]]+$/, ""); if (length) print}
' af3_pipeline/finalize_e217_generation.sh)
ranking_args=(); for file in "${rankings[@]}"; do ranking_args+=(--ranking "$file"); done
"$PY" af3_pipeline/score_c4_library_hamming_ai.py \
  "${ranking_args[@]}" --library "$POOL/generation_shortlist.tsv" \
  --out "$POOL/hamming_ai_ranking_current_weighted.tsv" \
  --shortlist 32 --minimum-distance 3 --allow-short-shortlist \
  --duplicate-label-policy max_models_mean --weighted-iptm \
  --generation-rank-weight 0.25
"$PY" af3_pipeline/prepare_ova_body_af3.py \
  --fasta "$POOL/hamming_ai_ranking_current_weighted_shortlist.fasta" \
  --ova-msa experiments/e202_c4_ai11v2_ai_af3_seed1/msas/ova_body_c4_01.a3m \
  --out "$E284" --top 32 --seeds 1 --seed-start 1 --seed-offset 0 \
  --shards 8 --n-chains 4
touch "$E284/CANDIDATES_FROZEN"
echo "E283 corrected trajectories frozen as E284"
