#!/usr/bin/env bash
# Run on Zhongwei: wait for an eight-shard tied-AF2 batch and create a small,
# checksummed archive containing the sequence trajectories needed downstream.
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
candidate_files=("$EXPERIMENT"/$RUN_GLOB/joint_*/af2diff_candidates.tsv)
trajectory_files=("$EXPERIMENT"/$RUN_GLOB/joint_*/trajectory.tsv)
if [[ "${#candidate_files[@]}" -ne 8 ]]; then
  printf 'expected 8 completed candidate tables, found %s\n' "${#candidate_files[@]}" >&2
  exit 2
fi
if [[ "${#trajectory_files[@]}" -ne 8 ]]; then
  printf 'expected 8 completed trajectory tables, found %s\n' "${#trajectory_files[@]}" >&2
  exit 2
fi

relative_files=()
for file in "${candidate_files[@]}" "${trajectory_files[@]}"; do
  relative_files+=("${file#"$PROJECT/"}")
done
log_files=("$EXPERIMENT"/logs/*)
for file in "${log_files[@]}"; do
  relative_files+=("${file#"$PROJECT/"}")
done

archive_rel="$EXPERIMENT_REL/$ARCHIVE_NAME"
tar -czf "$archive_rel" "${relative_files[@]}"
sha256sum "$archive_rel" > "$archive_rel.sha256"
printf 'archived %s candidate and %s trajectory tables to %s\n' \
  "${#candidate_files[@]}" "${#trajectory_files[@]}" "$archive_rel"
