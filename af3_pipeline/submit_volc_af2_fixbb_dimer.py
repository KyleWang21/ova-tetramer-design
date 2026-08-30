#!/usr/bin/env python3
"""Submit independent AF2 fixed-backbone homodimer design seeds to Volcengine."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import pathlib

from submit_volc import call, payload, portable, QUEUE


ENV_PYTHON = "/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python"
PARAMS = "/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/params"
AF2_IMAGE = "vemlp-cn-beijing.cr.volces.com/preset-images/cuda:12.9.0-py3.11-ubuntu22.04"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--experiment", type=pathlib.Path, required=True)
    ap.add_argument("--target-pdb", type=pathlib.Path, required=True)
    ap.add_argument(
        "--reference-fasta",
        type=pathlib.Path,
        default=pathlib.Path("/root/400083/ova_p3_tetramer_design/OVA_P3-13R_四聚体候选_AA.fasta"),
    )
    ap.add_argument("--reference-name", default="D0-P3-13R")
    ap.add_argument("--chains", default="A,B")
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--job-number-start", type=int, required=True)
    ap.add_argument("--free-interface", action="store_true")
    ap.add_argument("--submit", action="store_true")
    args = ap.parse_args()
    seeds = [int(value) for value in args.seeds.split(",") if value.strip()]
    if not 1 <= len(seeds) <= 8: raise SystemExit("one to eight seeds are required")
    args.experiment.mkdir(parents=True, exist_ok=True)
    records = []; failed = False
    for index, seed in enumerate(seeds):
        output = portable(args.experiment / f"volc_seed{seed:02d}")
        body = payload(args.experiment, index, not args.submit, args.job_number_start)
        body["RuntimeConfig"]["Image"]["Url"] = AF2_IMAGE
        body["ResourceConfig"]["MaxRuntimeSeconds"] = 43200
        inner = (
            f"{ENV_PYTHON} /root/400083/ova_p3_tetramer_design/af3_pipeline/af2_fixbb_noncovalent_interface_design.py "
            f"--target-pdb {portable(args.target_pdb)} "
            f"--reference-fasta {portable(args.reference_fasta)} "
            f"--reference-name {args.reference_name} "
            f"--design-positions {portable(args.experiment / 'design_positions.txt')} "
            f"--params {PARAMS} --out {output} --chains A,B --seeds {seed} "
            "--logits-iters 10 --soft-iters 5 --hard-iters 3 --recycles 0 "
            "--init gumbel --target-contact-cutoff 12 --predicted-contact-cutoff 14 --w-target-interface 5"
        )
        inner = inner.replace("--chains A,B", f"--chains {args.chains}")
        if args.free_interface:
            inner += " --free-interface"
        body["RuntimeConfig"]["Command"] = (
            f"/bin/bash -lc 'mkdir -p {output} && PYTHONUNBUFFERED=1 "
            f"XLA_PYTHON_CLIENT_PREALLOCATE=false {inner} > {output}/run.log 2>&1'"
        )
        rc, response = call(body)
        try: job_id = str((json.loads(response).get("Result") or {}).get("Id") or "")
        except json.JSONDecodeError: job_id = ""
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
