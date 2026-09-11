#!/usr/bin/env python3
"""Submit the first eight sparse C4 AF2 design tasks to c20250508."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import pathlib
import re

from submit_volc import QUEUE, call, payload


CUDA_IMAGE = "vemlp-cn-beijing.cr.volces.com/preset-images/cuda:12.9.0-py3.11-ubuntu22.04"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--experiment", type=pathlib.Path, required=True)
    ap.add_argument("--job-number-start", type=int, default=1)
    ap.add_argument("--worker-script", default="run_c4_lowmut_af2_worker.sh")
    ap.add_argument("--log-prefix", default="volc")
    ap.add_argument("--record-file", default="volc_jobs.tsv")
    ap.add_argument("--seed-offset", type=int, default=0)
    ap.add_argument("--run-prefix", default="run")
    ap.add_argument("--submit", action="store_true")
    args = ap.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_]+", args.run_prefix):
        raise SystemExit("--run-prefix must contain only letters, digits, and underscores")
    args.experiment.mkdir(parents=True, exist_ok=True)
    (args.experiment / "logs").mkdir(exist_ok=True)
    records = []
    failed = False
    experiment_portable = str(args.experiment.resolve()).replace(
        "/vepfs-mlp2/c20250508/400083", "/root/400083", 1
    )
    for task_index in range(8):
        body = payload(args.experiment, task_index, not args.submit, args.job_number_start)
        body["RuntimeConfig"]["Image"]["Url"] = CUDA_IMAGE
        body["ResourceConfig"]["MaxRuntimeSeconds"] = 43200
        body["RuntimeConfig"]["Command"] = (
            "/bin/bash -lc 'PYTHONUNBUFFERED=1 XLA_PYTHON_CLIENT_PREALLOCATE=false "
            f"OVA_SEED_OFFSET={args.seed_offset} OVA_IPTM90_RUN_PREFIX={args.run_prefix} "
            "/bin/bash /root/400083/ova_p3_tetramer_design/af3_pipeline/"
            f"{args.worker_script} {task_index} > "
            f"{experiment_portable}/logs/{args.log_prefix}_{task_index:02d}.log 2>&1'"
        )
        rc, response = call(body)
        try:
            job_id = str((json.loads(response).get("Result") or {}).get("Id") or "")
        except json.JSONDecodeError:
            job_id = ""
        ok = (rc == 0) or (not args.submit and "DryRunOperation" in response)
        if args.submit and not job_id:
            ok = False
        failed |= not ok
        print(f"task={task_index}\tname={body['Name']}\tjob={job_id or '-'}\tok={ok}")
        if not ok:
            print(response)
        records.append({
            "submitted_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "task_index": task_index,
            "name": body["Name"],
            "job_id": job_id,
            "queue_id": QUEUE,
            "mode": "submit" if args.submit else "dry_run",
            "status": "ok" if ok else "error",
            "response": response.replace("\t", " ").replace("\n", " "),
        })
    if args.submit:
        path = args.experiment / args.record_file
        exists = path.exists() and path.stat().st_size > 0
        with path.open("a", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
            if not exists:
                writer.writeheader()
            writer.writerows(records)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
