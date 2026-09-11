#!/usr/bin/env bash
# Run on Zhongwei: wait for eight tied-ProteinMPNN shards and archive proposals.
set -euo pipefail
shopt -s nullglob

EXPERIMENT_REL="${1:?project-relative experiment directory required}"
RUN_GLOB="${2:?run prefix glob required}"
JOB_IDS="${3:?comma-separated Slurm job IDs required}"
ARCHIVE_NAME="${4:?archive filename required}"
PROJECT="${OVA_TETRA_ROOT:-/data/home/scwb671/run/wky/ova_p3_tetramer_design}"
EXPERIMENT="$PROJECT/$EXPERIMENT_REL"

while [[ -n "$(squeue -h -j "$JOB_IDS")" ]]; do
  sleep 30
done

cd "$PROJECT"
shortlists=("$EXPERIMENT"/$RUN_GLOB/mpnn_shortlist.tsv)
if [[ "${#shortlists[@]}" -ne 8 ]]; then
  printf 'expected 8 completed MPNN shortlists, found %s\n' "${#shortlists[@]}" >&2
  exit 2
fi
relative_files=()
for table in "${shortlists[@]}"; do
  directory=$(dirname "$table")
  for name in mpnn_shortlist.tsv mpnn_shortlist.fasta mpnn_all_valid.tsv rejection_counts.tsv run_metadata.json; do
    file="$directory/$name"
    [[ -s "$file" ]] && relative_files+=("${file#"$PROJECT/"}")
  done
done
log_files=("$EXPERIMENT"/logs/*)
for file in "${log_files[@]}"; do
  relative_files+=("${file#"$PROJECT/"}")
done

archive_rel="$EXPERIMENT_REL/$ARCHIVE_NAME"
tar -czf "$archive_rel" "${relative_files[@]}"
sha256sum "$archive_rel" > "$archive_rel.sha256"
printf 'archived %s MPNN shortlists to %s\n' "${#shortlists[@]}" "$archive_rel"
