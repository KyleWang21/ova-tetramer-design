#!/usr/bin/env python3
"""Dry-run or submit up to eight OVA AF3 shards to Volcengine MLP."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import pathlib
import subprocess


# This project/account currently has CreateJob permission on the c20250508
# shared queue. The older queue012 ID in volcengine_mlp_queue_tasks.md is
# visible but returned 403 on 2026-08-29.
QUEUE = "q-20251224105744-qhrgl"
IMAGE = "vemlp-cn-beijing.cr.volces.com/preset-images/cuda:12.8.1-py3.12-ubuntu22.04"
VEPFS = "vepfs-cnbj2c98dea54433"
VE = "/home/wangkuiyu/.local/bin/ve"
ROOT = pathlib.Path(__file__).resolve().parents[1]


def portable(path: pathlib.Path) -> str:
    return str(path.resolve()).replace("/vepfs-mlp2/c20250508/400083", "/root/400083", 1)


def parse_shards(spec: str, available: set[int]) -> list[int]:
    result: set[int] = set()
    for part in spec.split(","):
        if "-" in part:
            lo, hi = map(int, part.split("-", 1))
            result.update(range(lo, hi + 1))
        else:
            result.add(int(part))
    unknown = result - available
    if unknown:
        raise SystemExit(f"shards have no manifest entries: {sorted(unknown)}")
    return sorted(result)


def payload(
    experiment: pathlib.Path, shard: int, dry_run: bool, job_number_start: int
) -> dict[str, object]:
    exp_portable = portable(experiment)
    # User-mandated queue naming convention.  Do not expose experiment names,
    # shard labels, or free-form tags in the Volcengine job name.
    beijing = dt.timezone(dt.timedelta(hours=8))
    date = dt.datetime.now(beijing).strftime("%Y%m%d")
    job_number = job_number_start + shard
    if job_number < 1:
        raise SystemExit(
            f"Volcengine job number must be positive; got {job_number}"
        )
    # The user-mandated form is ``ky-YYYYMMDD-number``.  Keep three-digit
    # zero-padding for the historical range, but do not reject the four-digit
    # serials used by the continuing experiment registry.
    name = f"ky-{date}-{job_number:03d}"
    return {
        "Name": name,
        "ProjectName": "default",
        "DryRun": dry_run,
        "ResourceConfig": {
            "ResourceQueueId": QUEUE,
            "MaxRuntimeSeconds": 43200,
            "Roles": [{
                "Name": "worker",
                "Replicas": 1,
                "Resource": {"Type": "Native", "InstanceTypeId": "ml.pni2.3xlarge"},
            }],
        },
        "RuntimeConfig": {
            "Framework": "Custom",
            "Image": {"Type": "VolcEngine", "Url": IMAGE},
            "Command": (
                "/bin/bash /root/400083/ova_p3_tetramer_design/af3_pipeline/run_worker.sh "
                f"{exp_portable} {shard}"
            ),
        },
        "StorageConfig": {"Storages": [{
            "Type": "Vepfs",
            "MountPath": "/root/400083",
            "ReadOnly": False,
            "Config": {"Vepfs": {
                "Id": VEPFS,
                "HostPath": f"/mnt/{VEPFS}",
                "SubPath": "c20250508/400083",
            }},
        }]},
    }


def call(body: dict[str, object]) -> tuple[int, str]:
    command = [VE, "mlplatform20240701", "CreateJob", "--body", json.dumps(body)]
    if os.geteuid() == 0:
        command = ["runuser", "-u", "wangkuiyu", "--", *command]
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return result.returncode, result.stdout.strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", type=pathlib.Path, required=True)
    parser.add_argument("--shards", default="0-7")
    parser.add_argument("--job-number-start", type=int, required=True)
    parser.add_argument("--record-file", default="jobs.tsv")
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args()
    superseded = args.experiment / "SUPERSEDED_DO_NOT_SUBMIT.md"
    if superseded.exists():
        raise SystemExit(f"refusing superseded experiment: {superseded}")
    manifest = args.experiment / "manifest.tsv"
    if not manifest.exists():
        raise SystemExit(f"missing {manifest}")
    rows = list(csv.DictReader(manifest.open(), delimiter="\t"))
    available = {int(row["shard"]) for row in rows}
    shards = parse_shards(args.shards, available)

    records: list[dict[str, object]] = []
    failed = False
    for shard in shards:
        body = payload(args.experiment, shard, not args.submit, args.job_number_start)
        rc, response = call(body)
        try:
            job_id = str((json.loads(response).get("Result") or {}).get("Id") or "")
        except json.JSONDecodeError:
            job_id = ""
        ok = (rc == 0) or (not args.submit and "DryRunOperation" in response)
        if args.submit and not job_id:
            ok = False
        failed |= not ok
        print(f"s{shard}\tname={body['Name']}\trc={rc}\tjob={job_id or '-'}\tok={ok}")
        if not ok:
            print(response)
        records.append({
            "submitted_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "shard": shard,
            "name": body["Name"],
            "job_id": job_id,
            "queue_id": QUEUE,
            "mode": "submit" if args.submit else "dry_run",
            "status": "ok" if ok else "error",
            "response": response.replace("\t", " ").replace("\n", " "),
        })
    if args.submit:
        jobs = args.experiment / args.record_file
        exists = jobs.exists() and jobs.stat().st_size > 0
        with jobs.open("a", newline="") as handle:
            writer = csv.DictWriter(handle, list(records[0]), delimiter="\t")
            if not exists:
                writer.writeheader()
            writer.writerows(records)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
