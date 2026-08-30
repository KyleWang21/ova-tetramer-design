#!/usr/bin/env python3
"""Extract selected protein chains from a PDB while preserving coordinates."""

from __future__ import annotations

import argparse
import pathlib

from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select


class SelectedChains(Select):
    def __init__(self, chains: set[str]) -> None:
        self.chains = chains

    def accept_chain(self, chain) -> bool:  # noqa: ANN001
        return str(chain.id) in self.chains

    def accept_residue(self, residue) -> bool:  # noqa: ANN001
        return residue.id[0] == " "


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", type=pathlib.Path, required=True)
    ap.add_argument("--chains", required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    args = ap.parse_args()
    selected = {item.strip() for item in args.chains.split(",") if item.strip()}
    parser = (
        MMCIFParser(QUIET=True)
        if args.pdb.suffix.lower() in {".cif", ".mmcif"}
        else PDBParser(QUIET=True)
    )
    structure = parser.get_structure("selected", str(args.pdb))
    present = {str(chain.id) for chain in structure.get_chains()}
    if not selected <= present:
        raise SystemExit(f"missing chains: {sorted(selected - present)}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO(); io.set_structure(structure); io.save(str(args.out), SelectedChains(selected))
    print(f"wrote chains {','.join(sorted(selected))} to {args.out}")


if __name__ == "__main__":
    main()
