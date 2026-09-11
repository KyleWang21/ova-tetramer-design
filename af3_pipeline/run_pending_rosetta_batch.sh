#!/usr/bin/env bash
# Run the current pending candidates through the Rosetta primary-interface screen.
# Rosetta is CPU-bound; this launcher uses eight candidate jobs at a time and
# eight independent relax chunks per candidate (200 relaxes per candidate).
set -uo pipefail

PROJECT="${PROJECT:-/root/400083/ova_p3_tetramer_design}"
ROOT="$PROJECT/experiments/current_pending_v1_screen"
LOG_ROOT="$ROOT/rosetta_batch_20260904"
mkdir -p "$LOG_ROOT"
TMP_ROOT="${TMPDIR:-/tmp}/ova_c4_rosetta_20260904"
mkdir -p "$TMP_ROOT"
export OVA_ROSETTA_CHUNKS="${OVA_ROSETTA_CHUNKS:-8}"
export OVA_ROSETTA_TOTAL="${OVA_ROSETTA_TOTAL:-200}"

run_one() {
  local manifest="$1"
  local screen candidate source_model iptm effective_edges design_edges background selection cif
  screen="${manifest%/rosetta_manifest.tsv}"
  candidate="$(basename "$screen")"
  IFS=$'\t' read -r _ _ _ effective_edges design_edges background selection cif < <(sed -n '2p' "$manifest")
  local out="$screen/rosetta"
  local tmp="$TMP_ROOT/$candidate"
  local tmp_out="$tmp/out"
  local chunks="$tmp/chunks8"
  local summary="$out/${candidate}_AF3_rosetta_summary.tsv"
  if [[ -s "$summary" ]]; then
    echo "$(date -u '+%F %T UTC') skip $candidate (summary exists)"
    return 0
  fi
  echo "$(date -u '+%F %T UTC') start $candidate" > "$LOG_ROOT/$candidate.log"
  mkdir -p "$tmp_out" "$chunks"
  # The source-model field is fixed to AF3 by prepare_rosetta_c4_manifest.py.
  bash "$PROJECT/af3_pipeline/run_rosetta_c4_parallel.sh" \
    "$cif" "$candidate" "$effective_edges" "$design_edges" AF3 "$chunks" "$tmp_out" \
    >> "$LOG_ROOT/$candidate.log" 2>&1
  local rc=$?
  if [[ "$rc" -eq 0 ]]; then
    mkdir -p "$out"
    cp "$tmp_out/${candidate}_AF3_rosetta_summary.tsv" "$out/"
    cp "$tmp_out/${candidate}_AF3_rosetta_summary.json" "$out/" 2>/dev/null || true
    cp "$tmp_out/${candidate}_AF3_rosetta_replicates.tsv" "$out/" 2>/dev/null || true
    cp "$tmp_out/${candidate}_AF3_relaxed_best.pdb" "$out/" 2>/dev/null || true
  fi
  echo "$(date -u '+%F %T UTC') finish $candidate rc=$rc" >> "$LOG_ROOT/$candidate.log"
  return "$rc"
}

active=0
failed=0
MAX_ACTIVE_CANDIDATES="${OVA_ROSETTA_CANDIDATES_PARALLEL:-2}"
if ! [[ "$MAX_ACTIVE_CANDIDATES" =~ ^[1-9][0-9]*$ ]]; then
  echo "invalid OVA_ROSETTA_CANDIDATES_PARALLEL=$MAX_ACTIVE_CANDIDATES" >&2
  exit 2
fi
for manifest in "$ROOT"/ova_body_c4_*/rosetta_manifest.tsv; do
  [[ -f "$manifest" ]] || continue
  run_one "$manifest" &
  active=$((active + 1))
  if (( active >= MAX_ACTIVE_CANDIDATES )); then
    wait -n || failed=1
    active=$((active - 1))
  fi
done
while (( active > 0 )); do
  wait -n || failed=1
  active=$((active - 1))
done
echo "$(date -u '+%F %T UTC') batch finished failed=$failed"
exit "$failed"
