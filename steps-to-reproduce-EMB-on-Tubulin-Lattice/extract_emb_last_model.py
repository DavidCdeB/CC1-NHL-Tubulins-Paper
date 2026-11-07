#!/usr/bin/env python3
"""
extract_cc1_last_model.py

CC1_from_Prod_NPT_10ns.pdb contains many MODEL/ENDMDL frames.
For GROMACS setup you typically want ONE model.

This extracts the LAST model and writes it to a new PDB, setting a single chain ID (default Z).

Usage:
  python3 extract_cc1_last_model.py -i CC1_from_Prod_NPT_10ns.pdb -o cc1_lastmodel.pdb --chain Z
"""

from __future__ import annotations
import argparse

CHAIN_COL = 21  # 0-based index

def set_chain(line: str, chain: str) -> str:
    if len(line) <= CHAIN_COL:
        line = line.rstrip("\n").ljust(CHAIN_COL + 1) + "\n"
    return line[:CHAIN_COL] + chain + line[CHAIN_COL + 1:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--inp", required=True)
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--chain", default="Z")
    args = ap.parse_args()

    in_model = False
    any_model = False
    current = []
    chosen = []

    with open(args.inp, "r", encoding="utf-8", errors="replace") as fin:
        for line in fin:
            if line.startswith("MODEL"):
                any_model = True
                in_model = True
                current = []
                continue
            if line.startswith("ENDMDL"):
                if in_model:
                    chosen = current[:]  # keep the last one
                in_model = False
                continue

            if any_model:
                if in_model:
                    current.append(line)
            else:
                chosen.append(line)

    with open(args.out, "w", encoding="utf-8") as fout:
        for line in chosen:
            if line.startswith("END"):
                continue
            if line.startswith(("ATOM", "HETATM", "TER")):
                line = set_chain(line, args.chain)
            fout.write(line)
        fout.write("END\n")

    print(f"Wrote {args.out} (last model), chain={args.chain}")

if __name__ == "__main__":
    main()
