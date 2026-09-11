#!/usr/bin/env python3
"""Submit eight AF3 shards only when c20250508 has more than ten free GPUs."""

from __future__ import annotations

import argparse
import csv
import json
import fcntl
import os
import pathlib
import subprocess
import sys


QUEUE = "q-20251224105744-qhrgl"
VE = "/home/wangkuiyu/.local/bin/ve"
PROJECT = pathlib.Path(__file__).resolve().parents[1]
PRIORITY_LOCK = PROJECT / "experiments" / "VOLC_PRIORITY_BATCH.lock"
SUBMIT_MUTEX = PROJECT / "experiments" / "VOLC_CAPACITY_SUBMIT.mutex"


def _run() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", type=pathlib.Path, required=True)
    parser.add_argument("--job-number-start", type=int, required=True)
    parser.add_argument("--record-file", default="volc_af3_jobs.tsv")
    parser.add_argument("--marker", default="VOLC_AF3_SHARDS_SUBMITTED")
    parser.add_argument("--shards", default="auto")
    parser.add_argument(
        "--expected", type=int, default=0,
        help="expected submitted jobs; 0 infers the number of manifest shards",
    )
    parser.add_argument(
        "--priority", action="store_true",
        help="allow this explicitly prioritized batch to run while the global priority lock exists",
    )
    args = parser.parse_args()
    if PRIORITY_LOCK.exists() and not args.priority:
        print(f"priority batch lock active: {PRIORITY_LOCK}; no ordinary submission")
        return
    record = args.experiment / args.record_file
    marker = args.experiment / args.marker
    manifest_rows = list(csv.DictReader(
        (args.experiment / "manifest.tsv").open(), delimiter="\t"
    ))
    manifest_shards = sorted({int(row["shard"]) for row in manifest_rows})
    if not manifest_shards:
        raise SystemExit("manifest contains no shards")
    shards = args.shards
    if shards == "auto":
        shards = ",".join(map(str, manifest_shards))
    expected = args.expected or len(manifest_shards)
    if record.exists():
        rows = list(csv.DictReader(record.open(), delimiter="\t"))
        successful = [
            row for row in rows
            if row.get("mode") == "submit" and row.get("status") == "ok"
        ]
        if len(successful) >= expected:
            marker.touch()
            print(f"already submitted {len(successful)} AF3 shards")
            return

    command = [VE, "mlplatform20240701", "GetResourceQueue", "--Id", QUEUE]
    if os.geteuid() == 0:
        command = ["runuser", "-u", "wangkuiyu", "--", *command]
    response = subprocess.run(
        command, check=True, text=True, stdout=subprocess.PIPE
    ).stdout
    result = json.loads(response)["Result"]
    capability = int(result["QuotaCapability"]["GpuCount"])
    allocated = int(result["QuotaAllocated"]["GpuCount"])
    free = capability - allocated
    print(f"c20250508 capability={capability} allocated={allocated} free={free}")
    if free <= 10:
        print("capacity gate not met: require free GPU >10; no submission")
        return

    submitter = pathlib.Path(__file__).with_name("submit_volc.py")
    subprocess.run([
        sys.executable, str(submitter),
        "--experiment", str(args.experiment),
        "--shards", shards,
        "--job-number-start", str(args.job_number_start),
        "--record-file", args.record_file,
        "--submit",
    ], check=True)
    rows = list(csv.DictReader(record.open(), delimiter="\t"))
    successful = [
        row for row in rows
        if row.get("mode") == "submit" and row.get("status") == "ok"
    ]
    if len(successful) != expected:
        raise SystemExit(
            f"expected {expected} successful submissions, found {len(successful)}"
        )
    marker.touch()
    print(f"submitted {expected} shards and wrote {marker}")


def main() -> None:
    # Several long-running watchers can become eligible at the same time.
    # Serialize only the queue-capacity check and CreateJob calls so two
    # watchers cannot both observe the same free-GPU count and over-submit.
    SUBMIT_MUTEX.parent.mkdir(parents=True, exist_ok=True)
    with SUBMIT_MUTEX.open("a+") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            _run()
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


if __name__ == "__main__":
    main()
