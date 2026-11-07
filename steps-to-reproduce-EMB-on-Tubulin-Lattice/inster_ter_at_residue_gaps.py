#!/usr/bin/env python3
"""
Insert TER records at internal residue-number gaps within chains.

Why:
- If a chain has residues like 37 then 47 (missing 38-46), pdb2gmx will
  otherwise connect them with a peptide bond, causing stretched bonds/exclusions
  and warnings like "largest distance between excluded atoms ... larger than cutoff".

This script preserves original residue numbers (important if you need resid 350-355).
It only splits chains by adding TER records.
"""

from __future__ import annotations
import argparse
from typing import Optional, Tuple, Dict, List

def is_atom_line(line: str) -> bool:
    return line.startswith("ATOM  ") or line.startswith("HETATM")

def get_chain_id(line: str) -> str:
    return line[21:22]

def get_resseq_icode(line: str) -> Tuple[int, str]:
    resseq_str = line[22:26]
    icode = line[26:27]
    try:
        resseq = int(resseq_str)
    except ValueError:
        resseq = 0
    return resseq, icode

def get_resname(line: str) -> str:
    return line[17:20]

def get_atom_serial(line: str) -> Optional[int]:
    s = line[6:11]
    try:
        return int(s)
    except ValueError:
        return None

def format_ter(serial: int, resname: str, chain: str, resseq: int, icode: str) -> str:
    # PDB TER format fields: "TER" + serial + resname + chain + resseq + icode
    # Columns:  1-6  7-11    18-20  22   23-26 27
    return f"TER   {serial:5d}      {resname:>3s} {chain}{resseq:4d}{icode}\n"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_pdb", required=True)
    ap.add_argument("--out_pdb", required=True)
    ap.add_argument(
        "--chains",
        default="",
        help="Optional: only process these chains (e.g. 'A' or 'AZ'). Default: all chains.",
    )
    ap.add_argument(
        "--min_jump",
        type=int,
        default=2,
        help="Insert TER if current_resseq - prev_resseq >= min_jump (default: 2, i.e. gaps > 1).",
    )
    args = ap.parse_args()

    target_chains = set(args.chains) if args.chains else None

    with open(args.in_pdb, "r") as f:
        lines = f.readlines()

    # Track last residue per chain
    last_res: Dict[str, Tuple[int, str, str]] = {}  # chain -> (resseq, icode, resname)
    last_atom_serial = 0
    inserted = 0
    report: List[str] = []

    out: List[str] = []

    for line in lines:
        if is_atom_line(line):
            serial = get_atom_serial(line)
            if serial is not None and serial > last_atom_serial:
                last_atom_serial = serial

            chain = get_chain_id(line)
            if target_chains is not None and chain not in target_chains:
                out.append(line)
                continue

            resseq, icode = get_resseq_icode(line)
            resname = get_resname(line)

            if chain in last_res:
                prev_resseq, prev_icode, prev_resname = last_res[chain]

                # New residue?
                if (resseq, icode, resname) != (prev_resseq, prev_icode, prev_resname):
                    jump = resseq - prev_resseq
                    if jump >= args.min_jump:
                        # Insert TER to split chain here
                        last_atom_serial += 1
                        out.append(format_ter(last_atom_serial, prev_resname, chain, prev_resseq, prev_icode))
                        inserted += 1
                        report.append(
                            f"Inserted TER in chain {chain}: {prev_resname.strip()} {prev_resseq}{prev_icode.strip() or ''}"
                            f"  ->  {resname.strip()} {resseq}{icode.strip() or ''} (jump {jump})"
                        )

            # Update last residue for this chain
            last_res[chain] = (resseq, icode, resname)
            out.append(line)

        else:
            # Pass through other records; keep existing TER/END/REMARK/etc.
            out.append(line)

    with open(args.out_pdb, "w") as f:
        f.writelines(out)

    print(f"Wrote: {args.out_pdb}")
    print(f"Inserted {inserted} TER record(s).")
    if report:
        print("Details:")
        for r in report[:50]:
            print("  " + r)
        if len(report) > 50:
            print(f"  ... {len(report)-50} more")
    else:
        print("No residue-number gaps detected (or chains filtered out).")

if __name__ == "__main__":
    main()

