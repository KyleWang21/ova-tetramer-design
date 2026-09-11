#!/usr/bin/env python3
"""Keep the first N complete records of an A3M alignment."""

from __future__ import annotations

import argparse
import pathlib

from prepare_inputs import parse_a3m_records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--records", type=int, required=True)
    args = parser.parse_args()
    records = parse_a3m_records(args.input.read_text())[: args.records]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("".join(f"{header}\n{sequence}\n" for header, sequence in records))
    print(f"wrote {len(records)} records to {args.out}")


if __name__ == "__main__":
    main()
