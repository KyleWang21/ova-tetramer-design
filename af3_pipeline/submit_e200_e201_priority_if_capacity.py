#!/usr/bin/env python3
"""Fill eight Volc A100s with two E178 strict hits, four seeds each."""

from __future__ import annotations

import csv
import json
import os
import pathlib
import subprocess
import sys


QUEUE = "q-20251224105744-qhrgl"
VE = "/home/wangkuiyu/.local/bin/ve"
PROJECT = pathlib.Path(__file__).resolve().parents[1]


def queue_free_gpu() -> tuple[int, int, int]:
    command = [VE, "mlplatform20240701", "GetResourceQueue", "--Id", QUEUE]
    if os.geteuid() == 0:
        command = ["runuser", "-u", "wangkuiyu", "--", *command]
    response = subprocess.run(command, check=True, text=True, stdout=subprocess.PIPE).stdout
    result = json.loads(response)["Result"]
    capacity = int(result["QuotaCapability"]["GpuCount"])
    allocated = int(result["QuotaAllocated"]["GpuCount"])
    return capacity, allocated, capacity - allocated


def successful_rows(path: pathlib.Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as handle:
        return [
            row for row in csv.DictReader(handle, delimiter="\t")
            if row.get("mode") == "submit" and row.get("status") == "ok"
        ]


def submit(experiment: pathlib.Path, start: int) -> None:
    record = experiment / "volc_af3_jobs.tsv"
    existing = successful_rows(record)
    if existing:
        if len(existing) != 4:
            raise SystemExit(f"refusing partial submission: {record} ({len(existing)}/4)")
        print(f"already submitted {experiment}: 4/4")
        return
    subprocess.run([
        sys.executable, str(PROJECT / "af3_pipeline" / "submit_volc.py"),
        "--experiment", str(experiment), "--shards", "0-3",
        "--job-number-start", str(start), "--record-file", record.name, "--submit",
    ], check=True)
    if len(successful_rows(record)) != 4:
        raise SystemExit(f"submission audit failed: {record}")
    (experiment / "VOLC_AF3_SHARDS_SUBMITTED").touch()


def main() -> None:
    capacity, allocated, free = queue_free_gpu()
    print(f"c20250508 capability={capacity} allocated={allocated} free={free}")
    if free <= 10:
        print("priority remote gate not met: require free GPU >10")
        raise SystemExit(10)
    submit(PROJECT / "experiments" / "e200_c4_e178_04_multiseed", 667)
    submit(PROJECT / "experiments" / "e201_c4_e178_15_multiseed", 683)
    (PROJECT / "experiments" / "E200_E201_PRIORITY_VOLC_SUBMITTED").touch()
    print("submitted E200/E201 priority 8-GPU multiseed batch")


if __name__ == "__main__":
    main()
