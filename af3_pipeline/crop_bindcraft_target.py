#!/usr/bin/env python3
"""Crop a compact 3D target around OVA hotspot residues for BindCraft."""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
from Bio.PDB import PDBIO, PDBParser, Select


class ResidueSelect(Select):
    def __init__(self, chain: str, residues: set[int]):
        self.chain=chain;self.residues=residues
    def accept_chain(self, chain): return chain.id==self.chain
    def accept_residue(self, residue): return residue.id[0]==" " and int(residue.id[1]) in self.residues


def segments(numbers: list[int]) -> list[tuple[int,int]]:
    out=[]
    for n in numbers:
        if not out or n>out[-1][1]+1:out.append((n,n))
        else:out[-1]=(out[-1][0],n)
    return out


def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument("--pdb",type=pathlib.Path,required=True)
    ap.add_argument("--hotspots",default="99,155,157,182,184,334,336,338")
    ap.add_argument("--chain",default="A")
    ap.add_argument("--radius",type=float,default=16.0)
    ap.add_argument("--sequence-padding",type=int,default=2)
    ap.add_argument("--out",type=pathlib.Path,required=True)
    args=ap.parse_args();args.out.parent.mkdir(parents=True,exist_ok=True)
    structure=PDBParser(QUIET=True).get_structure("target",str(args.pdb));model=next(structure.get_models());chain=model[args.chain]
    residues={int(r.id[1]):r for r in chain if r.id[0]==" " and "CA" in r}
    hotspots=[int(x) for x in args.hotspots.split(",")]
    hot_xyz=np.stack([np.asarray(residues[x]["CB"].coord if "CB" in residues[x] else residues[x]["CA"].coord) for x in hotspots])
    keep=set()
    for rid,res in residues.items():
        xyz=np.asarray(res["CB"].coord if "CB" in res else res["CA"].coord)
        if float(np.min(np.linalg.norm(hot_xyz-xyz,axis=1)))<=args.radius:keep.add(rid)
    expanded=set(keep)
    for rid in keep:
        expanded.update(x for x in range(rid-args.sequence_padding,rid+args.sequence_padding+1) if x in residues)
    keep=expanded
    io=PDBIO();io.set_structure(structure);io.save(str(args.out),ResidueSelect(args.chain,keep))
    meta={"source":str(args.pdb),"chain":args.chain,"hotspots":hotspots,"radius":args.radius,
          "n_residues":len(keep),"segments":segments(sorted(keep)),"residues":sorted(keep)}
    args.out.with_suffix(".json").write_text(json.dumps(meta,indent=2)+"\n")
    print(json.dumps({k:v for k,v in meta.items() if k!="residues"},indent=2))


if __name__=="__main__":main()
