#!/usr/bin/env bash
# Run one manifest row; intended for local xargs -P or a CPU-bearing batch job.
set -uo pipefail

MANIFEST="${1:?usage: run_rosetta_c4_manifest_row.sh MANIFEST ROW NRELAX OUT}"
ROW="${2:?usage: run_rosetta_c4_manifest_row.sh MANIFEST ROW NRELAX OUT}"
NRELAX="${3:?usage: run_rosetta_c4_manifest_row.sh MANIFEST ROW NRELAX OUT}"
OUT="${4:?usage: run_rosetta_c4_manifest_row.sh MANIFEST ROW NRELAX OUT}"
PY=/root/400083/antibody-design/pipeline_stage1_interface_energetics/env/bin/python
manifest_values=$(awk -F '\t' -v row="$ROW" 'NR==row+1 {
  if (NF >= 8) printf "%s|%s|%s|%s|%s|%s|%s|%s", $1,$2,$3,$4,$5,$6,$7,$8;
  else printf "%s|%s|%s|%s|%s||||%s", $1,$2,$3,$4,$4,$5;
}' "$MANIFEST")
# Use a non-whitespace separator so an empty background_edges field is not
# collapsed by bash (which would shift the CIF path into the wrong variable).
IFS='|' read -r candidate source_model iptm effective_edges design_edges background_edges design_selection cif_path <<<"$manifest_values"
[ -n "$candidate" ] || exit 2
candidate_out="$OUT/$candidate"
summary="$candidate_out/${candidate}_${source_model}_rosetta_summary.tsv"
mkdir -p "$candidate_out" "$OUT/logs"
if [ -s "$summary" ]; then
  echo "$candidate already complete"
  exit 0
fi
exec > "$OUT/logs/${candidate}.log" 2>&1
echo "$(date -u '+%F %T UTC') start $candidate"
"$PY" af3_pipeline/rosetta_c4_relax_screen.py \
  --in "$cif_path" --candidate "$candidate" --effective-edges "$effective_edges" \
  --design-interface-edges "$design_edges" \
  --source-model "$source_model" --nrelax "$NRELAX" --out "$candidate_out"
rc=$?
echo "$(date -u '+%F %T UTC') finish $candidate rc=$rc"
exit "$rc"
