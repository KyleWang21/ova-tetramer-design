#!/usr/bin/env bash
# Build and freeze an AF3-driven clash-repair batch after E252 seed1 completes.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
SOURCE="$PROJECT/experiments/e252_c4_epistasis_ai_af3_seed1"
OUT="$PROJECT/experiments/e289_c4_e252_clashrepair_ai_af3_seed1"
WORK="$PROJECT/experiments/e289_c4_e252_clashrepair_generation"
AF3_PY=/root/400083/alphafold3/repo/.venv/bin/python
BIO_PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"
if [[ -f "$OUT/CANDIDATES_FROZEN" ]]; then echo "E289 already frozen"; exit 0; fi
while [[ ! -f "$SOURCE/SEED1_STRICT_RESULTS_COMPLETE" ]]; do sleep 30; done

mkdir -p "$WORK"
"$AF3_PY" af3_pipeline/audit_af3_interchain_clashes.py \
  --experiment "$SOURCE" --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --out "$WORK/clash_audit"
"$BIO_PY" af3_pipeline/plot_af3_seed1_iptm_clash_tradeoff.py \
  --summary "$SOURCE/seed1_strict_final/seed1_strict_summary.tsv" \
  --out "$SOURCE/figures/e252_seed1_iptm_clash_tradeoff.png" \
  --title "E252 AF3 seed1: ipTM vs 2.4 A zero-clash reproducibility"
"$BIO_PY" af3_pipeline/prepare_e289_clash_repair_library.py \
  --experiment "$SOURCE" \
  --strict-summary "$SOURCE/seed1_strict_final/seed1_strict_summary.tsv" \
  --clash-recurrence "$WORK/clash_audit/clash_residue_pair_recurrence.tsv" \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --accepted-fasta deliverables/current_accepted_af3_ranked_20260902/all_sequences.fasta \
  --surface-positions experiments/e283_c4_chainmapped_eightstate_tied_design/surface_positions.txt \
  --exclude-glob 'experiments/e*/manifest.tsv' \
  --exclude-glob 'experiments/e*/*/manifest.tsv' \
  --exclude-glob 'experiments/e*/*/*/manifest.tsv' \
  --out "$WORK/generation_library.tsv" --max-parents 20 --max-edits 3

mapfile -t rankings < <(awk '
  /^rankings=\(/{inside=1; next}
  inside && /^\)/{exit}
  inside {gsub(/^[[:space:]]+|[[:space:]]+$/, ""); if (length) print}
' af3_pipeline/finalize_e217_generation.sh)
ranking_args=(); for file in "${rankings[@]}"; do ranking_args+=(--ranking "$file"); done
ranking_args+=(--ranking "$SOURCE/seed1_strict_final/seed1_strict_summary.tsv")
"$BIO_PY" af3_pipeline/score_c4_library_hamming_ai.py \
  "${ranking_args[@]}" --library "$WORK/generation_library.tsv" \
  --out "$WORK/hamming_ai_ranking.tsv" --shortlist 32 --minimum-distance 2 \
  --allow-short-shortlist --duplicate-label-policy max_models_mean --weighted-iptm
"$BIO_PY" af3_pipeline/prepare_ova_body_af3.py \
  --fasta "$WORK/hamming_ai_ranking_shortlist.fasta" \
  --ova-msa experiments/e202_c4_ai11v2_ai_af3_seed1/msas/ova_body_c4_01.a3m \
  --out "$OUT" --top 32 --seeds 1 --seed-start 1 --seed-offset 0 \
  --shards 8 --n-chains 4
touch "$OUT/CANDIDATES_FROZEN"
echo "E289 AF3-clash-directed repair candidates frozen"
