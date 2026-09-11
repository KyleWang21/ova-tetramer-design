#!/usr/bin/env bash
# Recompute the strict seed-1 table whenever another five-sample system completes.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
EXPERIMENT="${1:?usage: watch_af3_seed1_partial_scoring.sh EXPERIMENT [OUT]}"
OUT="${2:-$EXPERIMENT/live_seed1_strict}"
PY=/root/400083/alphafold3/repo/.venv/bin/python
PLOT_PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"

expected_systems=$(( $(wc -l < "$EXPERIMENT/manifest.tsv") - 1 ))
last_complete=-1
while true; do
  if compgen -G "$EXPERIMENT/WORKER_FAILED_s*" >/dev/null; then
    echo "AF3 worker failure marker detected in $EXPERIMENT" >&2
    exit 1
  fi
  complete_systems=$(find "$EXPERIMENT" -path '*/seed-*/*_summary_confidences.json' \
    -printf '%h\n' 2>/dev/null | sed 's#/seed-[^/]*$##' | sort | uniq -c | \
    awk '$1 == 5 {n++} END {print n+0}')
  if [[ "$complete_systems" -gt "$last_complete" ]]; then
    if [[ "$complete_systems" -gt 0 ]]; then
      "$PY" af3_pipeline/summarize_iptm90_seed1_strict.py \
        --experiment "$EXPERIMENT" \
        --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
        --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
        --out "$OUT"
      if [[ -x "$PLOT_PY" ]]; then
        "$PLOT_PY" af3_pipeline/plot_af3_seed1_iptm_clash_tradeoff.py \
          --summary "$OUT/seed1_strict_summary.tsv" \
          --out "$OUT/af3_seed1_iptm_clash_tradeoff.png" \
          --title "$(basename "$EXPERIMENT"): live AF3 seed1 screen" || true
        if [[ -s "$PROJECT/experiments/e289_c4_e252_clashrepair_generation/hamming_ai_ranking_shortlist.tsv" && \
              "$(basename "$EXPERIMENT")" == e289_* ]]; then
          "$PLOT_PY" af3_pipeline/plot_clash_repair_effect.py \
            --ranking "$PROJECT/experiments/e289_c4_e252_clashrepair_generation/hamming_ai_ranking_shortlist.tsv" \
            --manifest "$EXPERIMENT/manifest.tsv" \
            --summary "$OUT/seed1_strict_summary.tsv" \
            --out "$OUT/clash_repair_effect.png" || true
        fi
      fi
    fi
    echo "complete_systems=$complete_systems/$expected_systems"
    last_complete=$complete_systems
  fi
  if [[ "$complete_systems" -eq "$expected_systems" ]]; then
    touch "$OUT/LIVE_PARTIAL_SCORING_COMPLETE"
    exit 0
  fi
  sleep 30
done
