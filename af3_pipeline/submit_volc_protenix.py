#!/usr/bin/env python3
"""Submit eight final-candidate Protenix-v2 shards to the authorized queue."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import pathlib
import subprocess


QUEUE = "q-20251224105744-qhrgl"
IMAGE = "vemlp-cn-beijing.cr.volces.com/preset-images/cuda:12.8.1-py3.12-ubuntu22.04"
VEPFS = "vepfs-cnbj2c98dea54433"
VE = "/home/wangkuiyu/.local/bin/ve"


def portable(path: pathlib.Path) -> str:
    return str(path.resolve()).replace("/vepfs-mlp2/c20250508/400083", "/root/400083", 1)


def call(body: dict[str, object]) -> tuple[int, str]:
    command = [VE, "mlplatform20240701", "CreateJob", "--body", json.dumps(body)]
    if os.geteuid() == 0:
        command = ["runuser", "-u", "wangkuiyu", "--", *command]
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return result.returncode, result.stdout.strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", type=pathlib.Path, required=True)
    parser.add_argument("--job-number-start", type=int, required=True)
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args()
    manifest = args.experiment / "volc_manifest.tsv"
    rows = list(csv.DictReader(manifest.open(), delimiter="\t"))
    shards = sorted({int(row["shard"]) for row in rows})
    if shards != list(range(8)):
        raise SystemExit(f"expected shards 0-7, found {shards}")
    exp = portable(args.experiment)
    records = []
    failed = False
    for shard in shards:
        number = args.job_number_start + shard
        name = f"ky-{dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).strftime('%Y%m%d')}-{number:03d}"
        body = {
            "Name": name,
            "ProjectName": "default",
            "DryRun": not args.submit,
            "ResourceConfig": {
                "ResourceQueueId": QUEUE,
                "MaxRuntimeSeconds": 7200,
                "Roles": [{"Name": "worker", "Replicas": 1,
                           "Resource": {"Type": "Native", "InstanceTypeId": "ml.pni2.3xlarge"}}],
            },
            "RuntimeConfig": {
                "Framework": "Custom",
                "Image": {"Type": "VolcEngine", "Url": IMAGE},
                "Command": f"/bin/bash /root/400083/ova_p3_tetramer_design/af3_pipeline/run_protenix_volc_shard.sh {exp} {shard}",
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
        with (args.experiment / "volc_jobs.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, list(records[0]), delimiter="\t")
            writer.writeheader(); writer.writerows(records)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
