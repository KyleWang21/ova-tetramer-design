#!/usr/bin/env python3
"""Submit independent AF2 homodimer hallucination seeds to the c20250508 queue."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import pathlib

from submit_volc import call, payload, portable, QUEUE


ROOT = pathlib.Path(__file__).resolve().parents[1]
ENV_PYTHON = "/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python"
PARAMS = "/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/params"
AF2_IMAGE = "vemlp-cn-beijing.cr.volces.com/preset-images/cuda:12.9.0-py3.11-ubuntu22.04"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--experiment", type=pathlib.Path, required=True)
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--job-number-start", type=int, required=True)
    ap.add_argument("--logits-iters", type=int, default=12)
    ap.add_argument("--soft-iters", type=int, default=6)
    ap.add_argument("--hard-iters", type=int, default=3)
    ap.add_argument("--submit", action="store_true")
    args = ap.parse_args()
    seeds = [int(value) for value in args.seeds.split(",") if value.strip()]
    if not 1 <= len(seeds) <= 8:
        raise SystemExit("one to eight seeds are required")
    args.experiment.mkdir(parents=True, exist_ok=True)
    records = []; failed = False
    for index, seed in enumerate(seeds):
        output = portable(args.experiment / f"volc_seed{seed:02d}")
        body = payload(args.experiment, index, not args.submit, args.job_number_start)
        body["RuntimeConfig"]["Image"]["Url"] = AF2_IMAGE
        body["ResourceConfig"]["MaxRuntimeSeconds"] = 43200
        inner = (
            f"{ENV_PYTHON} /root/400083/ova_p3_tetramer_design/af3_pipeline/af2_homodimer_hallucination.py "
            "--reference-fasta /root/400083/ova_p3_tetramer_design/OVA_P3-13R_四聚体候选_AA.fasta "
            f"--design-positions {portable(args.experiment / 'design_positions.txt')} "
            f"--params {PARAMS} --out {output} --seeds {seed} "
            f"--logits-iters {args.logits_iters} --soft-iters {args.soft_iters} "
            f"--hard-iters {args.hard_iters} --recycles 0"
        )
        body["RuntimeConfig"]["Command"] = (
            f"/bin/bash -lc 'mkdir -p {output} && PYTHONUNBUFFERED=1 "
            f"XLA_PYTHON_CLIENT_PREALLOCATE=false {inner} > {output}/run.log 2>&1'"
        )
        rc, response = call(body)
        try:
            job_id = str((json.loads(response).get("Result") or {}).get("Id") or "")
        except json.JSONDecodeError:
            job_id = ""
        ok = (rc == 0) or (not args.submit and "DryRunOperation" in response)
        if args.submit and not job_id: ok = False
        failed |= not ok
        print(f"seed={seed}\tname={body['Name']}\tjob={job_id or '-'}\tok={ok}")
        if not ok: print(response)
        records.append({
            "submitted_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "seed": seed,
            "name": body["Name"], "job_id": job_id, "queue_id": QUEUE,
            "mode": "submit" if args.submit else "dry_run", "status": "ok" if ok else "error",
            "response": response.replace("\t", " ").replace("\n", " "),
        })
    if args.submit:
        path = args.experiment / "volc_jobs.tsv"; exists = path.exists() and path.stat().st_size > 0
        with path.open("a", newline="") as handle:
            writer = csv.DictWriter(handle, list(records[0]), delimiter="\t")
            if not exists: writer.writeheader()
            writer.writerows(records)
    if failed: raise SystemExit(1)


if __name__ == "__main__":
    main()
