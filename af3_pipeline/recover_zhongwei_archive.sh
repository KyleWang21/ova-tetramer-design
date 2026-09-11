#!/usr/bin/env bash
# Wait for one Zhongwei result archive, copy it through the persistent SSH
# control socket, verify its recorded SHA-256, and extract it in project root.
set -euo pipefail

REMOTE_RELATIVE="${1:?remote project-relative archive path required}"
MARKER="${2:?local completion marker required}"
PROJECT=/root/400083/ova_p3_tetramer_design
REMOTE_PROJECT=/data/home/scwb671/run/wky/ova_p3_tetramer_design
HOST=ssh.cn-zhongwei-1.paracloud.com
USER_NAME='scwb671@NMCC-N46A8'
PORT=2222
CONTROL=/tmp/ova-zhongwei-control.sock
LOCAL_ARCHIVE="$PROJECT/$REMOTE_RELATIVE"
LOCAL_SHA="$LOCAL_ARCHIVE.sha256"

cd "$PROJECT"
while ! ssh -S "$CONTROL" -o "User=$USER_NAME" -p "$PORT" "$HOST" \
  "test -s '$REMOTE_PROJECT/$REMOTE_RELATIVE' && test -s '$REMOTE_PROJECT/$REMOTE_RELATIVE.sha256'"; do
  sleep 30
done

scp -o "ControlPath=$CONTROL" -o "User=$USER_NAME" -P "$PORT" \
  "$HOST:$REMOTE_PROJECT/$REMOTE_RELATIVE" \
  "$HOST:$REMOTE_PROJECT/$REMOTE_RELATIVE.sha256" \
  "$(dirname "$LOCAL_ARCHIVE")/"
sha256sum -c "$LOCAL_SHA"
archive_members=$(tar -tzf "$LOCAL_ARCHIVE")
tied_candidates=$(printf '%s\n' "$archive_members" | \
  grep -Ec '/joint_[^/]*/af2diff_candidates\.tsv$' || true)
tied_trajectories=$(printf '%s\n' "$archive_members" | \
  grep -Ec '/joint_[^/]*/trajectory\.tsv$' || true)
if [[ "$tied_candidates" -gt 0 && "$tied_trajectories" -ne "$tied_candidates" ]]; then
  printf 'incomplete tied-AF2 archive: candidates=%s trajectories=%s\n' \
    "$tied_candidates" "$tied_trajectories" >&2
  exit 2
fi
tar -xzf "$LOCAL_ARCHIVE"
touch "$MARKER"
echo "recovered and extracted $REMOTE_RELATIVE candidates=$tied_candidates trajectories=$tied_trajectories"
