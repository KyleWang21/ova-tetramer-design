#!/usr/bin/env python3
"""Flush stepwise diagnostics for ColabDesign/JAX queue startup crashes."""

from __future__ import annotations

import os
import platform
import sys


def mark(text: str) -> None:
    print(text, flush=True)


mark(f"python={sys.version}")
mark(f"platform={platform.platform()}")
mark(f"LD_LIBRARY_PATH={os.environ.get('LD_LIBRARY_PATH', '')}")
mark("before import jax")
import jax  # noqa: E402
import jaxlib  # noqa: E402
mark(f"jax={jax.__version__} jaxlib={jaxlib.__version__}")
mark("before jax.devices")
mark(f"devices={jax.devices()}")
mark("before import colabdesign")
import colabdesign  # noqa: E402
mark(f"colabdesign={colabdesign.__file__}")
mark("before import mk_afdesign_model")
from colabdesign import mk_afdesign_model  # noqa: E402
mark("before mk_afdesign_model")
model = mk_afdesign_model(
    protocol="hallucination", use_multimer=True, use_templates=False,
    data_dir="/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/params",
    num_recycles=0, model_names=["model_1_multimer_v3"],
)
mark("before prep_inputs")
model.prep_inputs(length=386, copies=2)
mark("after prep_inputs")
reference = ""
active = False
for line in open("/root/400083/ova_p3_tetramer_design/OVA_P3-13R_四聚体候选_AA.fasta"):
    if line.startswith(">"):
        if active: break
        active = line[1:].startswith("D0-P3-13R")
    elif active:
        reference += line.strip()
mark(f"reference_length={len(reference)} before predict")
model.predict(seq=reference, models=[0], num_recycles=0, verbose=False)
mark(f"after predict log={model.aux['log']}")
mark("before one design_logits step")
model.restart(seed=0, seq=reference, rm_aa="C")
model.design_logits(1, models=[0], num_models=1, sample_models=False,
                    dropout=False, save_best=True, verbose=1)
mark("after one design_logits step")
mark("diagnostic complete")
