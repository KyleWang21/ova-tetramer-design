#!/usr/bin/env bash
# Volcengine-only wrapper for the OVA-C4-073 AF3-medoid E169 experiment.
# Keep the experiment path fixed so a missing environment variable can never
# silently fall back to the older E155 19K/21B/21S/20A parent set.
set -euo pipefail

PROJECT_ROOT="${OVA_TETRA_ROOT:-/root/400083/ova_p3_tetramer_design}"
export OVA_IPTM90_EXPERIMENT="$PROJECT_ROOT/experiments/e169_c4_073_af3ensemble_tied_design"

manifest="$OVA_IPTM90_EXPERIMENT/target_manifest.tsv"
[[ -s "$manifest" ]] || { echo "missing E169 target manifest: $manifest" >&2; exit 2; }
mapfile -t anchors < <(awk -F '\t' 'NR>1 {print $1}' "$manifest")
expected=(OVA-C4-073-BB1 OVA-C4-073-BB2 OVA-C4-073-BB3 OVA-C4-073-BB4)
[[ "${anchors[*]}" == "${expected[*]}" ]] || {
  echo "refusing non-E169 anchors: ${anchors[*]}" >&2
  exit 2
}

exec /bin/bash "$PROJECT_ROOT/af3_pipeline/run_iptm90_tied_design_worker.sh" "$@"
