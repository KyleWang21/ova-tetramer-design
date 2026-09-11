#!/usr/bin/env python3
"""Submit eight tied-AF2 jobs only when c20250508 has >10 free GPUs."""

from __future__ import annotations

import argparse
import csv
import json
import os
import pathlib
import subprocess
import sys


QUEUE = "q-20251224105744-qhrgl"
VE = "/home/wangkuiyu/.local/bin/ve"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", type=pathlib.Path,
                        default=pathlib.Path("experiments/e155_c4_iptm90_crossmodel_design"))
    parser.add_argument("--job-number-start", type=int, default=376)
    parser.add_argument("--record-file", default="volc_design_jobs.tsv")
    parser.add_argument("--worker-script", default="run_iptm90_tied_design_worker.sh")
    parser.add_argument("--log-prefix", default="volc_iptm90")
    parser.add_argument("--seed-offset", type=int, default=1000)
    parser.add_argument("--run-prefix", default="run_volc")
    args = parser.parse_args()
    record = args.experiment / args.record_file
    if record.exists():
        rows = list(csv.DictReader(record.open(), delimiter="\t"))
        submitted = [row for row in rows if row.get("mode") == "submit" and row.get("status") == "ok"]
        if len(submitted) >= 8:
            print(f"already submitted {len(submitted)} successful Volc jobs; no action")
            return

    query = [VE, "mlplatform20240701", "GetResourceQueue", "--Id", QUEUE]
    if not pathlib.Path(VE).exists():
        raise SystemExit(f"missing Volcengine CLI: {VE}")
    if os.geteuid() == 0:
        query = ["runuser", "-u", "wangkuiyu", "--", *query]
    response = subprocess.run(query, check=True, text=True, stdout=subprocess.PIPE).stdout
    result = json.loads(response)["Result"]
    capability = int(result["QuotaCapability"]["GpuCount"])
    allocated = int(result["QuotaAllocated"]["GpuCount"])
    free = capability - allocated
    print(f"c20250508 capability={capability} allocated={allocated} free={free}")
    if free <= 10:
        print("capacity gate not met: require free GPU >10; no submission")
        return

    submitter = pathlib.Path(__file__).with_name("submit_volc_c4_lowmut_af2.py")
    command = [
        sys.executable, str(submitter), "--experiment", str(args.experiment),
        "--job-number-start", str(args.job_number_start),
        "--worker-script", args.worker_script,
        "--log-prefix", args.log_prefix, "--record-file", record.name, "--submit",
        "--seed-offset", str(args.seed_offset), "--run-prefix", args.run_prefix,
    ]
    subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
