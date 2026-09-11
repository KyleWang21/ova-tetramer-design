#!/usr/bin/env python3
"""Submit eight final20 Protenix-v2 multiseed shards to the approved queue."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import fcntl
import os
import pathlib
import subprocess

from submit_volc import IMAGE, QUEUE, ROOT, VE, VEPFS, call, portable

SUBMIT_MUTEX = ROOT / "experiments" / "VOLC_CAPACITY_SUBMIT.mutex"


def _run() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", type=pathlib.Path, required=True)
    parser.add_argument("--job-number-start", type=int, required=True)
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args()
    rows = list(csv.DictReader((args.experiment / "volc_manifest.tsv").open(), delimiter="\t"))
    shards = sorted({int(row["shard"]) for row in rows})
    if not shards or len(shards) > 8:
        raise SystemExit(f"expected 1-8 populated shards, found {shards}")
    if args.submit:
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
    jobs_path = args.experiment / "volc_jobs.tsv"
    previous = []
    if jobs_path.exists() and jobs_path.stat().st_size:
        with jobs_path.open() as handle:
            previous = list(csv.DictReader(handle, delimiter="\t"))
    completed_shards = {
        int(row["shard"]) for row in previous
        if row.get("mode") == "submit" and row.get("status") == "ok" and row.get("job_id")
    }
    shards = [shard for shard in shards if shard not in completed_shards]
    if not shards:
        print(f"all Protenix shards already submitted: {sorted(completed_shards)}")
        return
    exp = portable(args.experiment)
    records = []
    failed = False
    for shard in shards:
        number = args.job_number_start + shard
        date = dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).strftime("%Y%m%d")
        name = f"ky-{date}-{number:03d}"
        body = {
            "Name": name,
            "ProjectName": "default",
            "DryRun": not args.submit,
            "ResourceConfig": {
                "ResourceQueueId": QUEUE,
                "MaxRuntimeSeconds": 14400,
                "Roles": [{"Name": "worker", "Replicas": 1,
                           "Resource": {"Type": "Native", "InstanceTypeId": "ml.pni2.3xlarge"}}],
            },
            "RuntimeConfig": {
                "Framework": "Custom",
                "Image": {"Type": "VolcEngine", "Url": IMAGE},
                "Command": (
                    "/bin/bash /root/400083/ova_p3_tetramer_design/"
                    f"af3_pipeline/run_protenix_multiseed_volc_shard.sh {exp} {shard}"
                ),
            },
            "StorageConfig": {"Storages": [{
                "Type": "Vepfs", "MountPath": "/root/400083", "ReadOnly": False,
                "Config": {"Vepfs": {"Id": VEPFS, "HostPath": f"/mnt/{VEPFS}",
                                      "SubPath": "c20250508/400083"}},
            }]},
        }
        rc, response = call(body)
        try:
            job_id = str((json.loads(response).get("Result") or {}).get("Id") or "")
        except json.JSONDecodeError:
            job_id = ""
        ok = (rc == 0) or (not args.submit and "DryRunOperation" in response)
        if args.submit and not job_id:
            ok = False
        failed |= not ok
        print(f"s{shard}\tname={name}\trc={rc}\tjob={job_id or '-'}\tok={ok}")
        if not ok:
            print(response)
        records.append({
            "submitted_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "shard": shard, "name": name, "job_id": job_id, "queue_id": QUEUE,
            "mode": "submit" if args.submit else "dry_run", "status": "ok" if ok else "error",
            "response": response.replace("\t", " ").replace("\n", " "),
        })
    if args.submit:
        exists = jobs_path.exists() and jobs_path.stat().st_size > 0
        with jobs_path.open("a", newline="") as handle:
            writer = csv.DictWriter(handle, list(records[0]), delimiter="\t")
            if not exists:
                writer.writeheader()
            writer.writerows(records)
    if failed:
        raise SystemExit(1)


def main() -> None:
    # Share the same short critical section as AF3 submission so concurrent
    # watchers cannot both pass the free-GPU gate on a stale queue snapshot.
    SUBMIT_MUTEX.parent.mkdir(parents=True, exist_ok=True)
    with SUBMIT_MUTEX.open("a+") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            _run()
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


if __name__ == "__main__":
    main()
