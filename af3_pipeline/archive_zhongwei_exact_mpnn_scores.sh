#!/usr/bin/env bash
# Run on Zhongwei: wait for eight exact-scoring shards and archive all score tables.
set -euo pipefail
shopt -s nullglob

EXPERIMENT_REL="${1:?project-relative experiment directory required}"
JOB_IDS="${2:?comma-separated Slurm job IDs required}"
ARCHIVE_NAME="${3:?archive filename required}"
PROJECT="${OVA_TETRA_ROOT:-/data/home/scwb671/run/wky/ova_p3_tetramer_design}"
EXPERIMENT="$PROJECT/$EXPERIMENT_REL"
while [[ -n "$(squeue -h -j "$JOB_IDS")" ]]; do sleep 30; done
cd "$PROJECT"
score_tables=("$EXPERIMENT"/task_*/candidate_scores.tsv)
order_tables=("$EXPERIMENT"/task_*/per_order.tsv)
metadata=("$EXPERIMENT"/task_*/metadata.json)
if [[ "${#score_tables[@]}" -ne 8 || "${#order_tables[@]}" -ne 8 || \
      "${#metadata[@]}" -ne 8 ]]; then
  printf 'incomplete E271 outputs: scores=%s orders=%s metadata=%s\n' \
    "${#score_tables[@]}" "${#order_tables[@]}" "${#metadata[@]}" >&2
  exit 2
fi
for table in "${score_tables[@]}"; do
  rows=$(tail -n +2 "$table" | wc -l)
  [[ "$rows" -eq 512 ]] || { echo "$table has $rows candidates" >&2; exit 2; }
done
for table in "${order_tables[@]}"; do
  rows=$(tail -n +2 "$table" | wc -l)
  [[ "$rows" -eq 2048 ]] || { echo "$table has $rows order scores" >&2; exit 2; }
done
files=()
for file in "${score_tables[@]}" "${order_tables[@]}" "${metadata[@]}"; do
  files+=("${file#"$PROJECT/"}")
done
logs=("$EXPERIMENT"/logs/*)
for file in "${logs[@]}"; do files+=("${file#"$PROJECT/"}"); done
archive_rel="$EXPERIMENT_REL/$ARCHIVE_NAME"
tar -czf "$archive_rel" "${files[@]}"
sha256sum "$archive_rel" > "$archive_rel.sha256"
echo "archived E271 eight-backbone exact scores to $archive_rel"
