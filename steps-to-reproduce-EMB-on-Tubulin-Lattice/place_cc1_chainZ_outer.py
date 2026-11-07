#!/usr/bin/env python3
"""
Translate ONLY chain Z (CC1) in a complex PDB so that:
  COM(chain A, res 350-355)  <->  COM(chain Z, res 24-30)
has a separation of sep_nm along the tubulin "radial" direction.

"Outer" radial direction is defined as:
  - Fit a principal axis to tubulin CA atoms (all chains except Z)
  - Compute radial vector from that axis line to the tubulin patch COM
  - Move CC1 patch COM to (tub_patch_COM + radial_dir * sep_nm)

PDB coordinates are assumed Å. sep_nm is in nm.
"""

from __future__ import annotations
import argparse
import math
from typing import List, Tuple

def dot(a, b): return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
def add(a, b): return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
def sub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def scale(a, s): return (a[0]*s, a[1]*s, a[2]*s)

def norm(a):
    n = math.sqrt(dot(a, a))
    if n < 1e-12:
        raise RuntimeError("Zero-length vector encountered.")
    return scale(a, 1.0/n)

def read_com(pdb_path: str, chain: str, r0: int, r1: int) -> Tuple[float,float,float]:
    xs = ys = zs = 0.0
    n = 0
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            ch = (line[21:22].strip() or " ")
            if ch != chain:
                continue
            rn = int(line[22:26])
            if rn < r0 or rn > r1:
                continue
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            xs += x; ys += y; zs += z
            n += 1
    if n == 0:
        raise RuntimeError(f"Empty selection for COM: chain {chain} resid {r0}-{r1}")
    return (xs/n, ys/n, zs/n)

def tubulin_ca_points(pdb_path: str, exclude_chain: str) -> List[Tuple[float,float,float]]:
    pts = []
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atom = line[12:16].strip()
            if atom != "CA":
                continue
            ch = (line[21:22].strip() or " ")
            if ch == exclude_chain:
                continue
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            pts.append((x,y,z))
    if len(pts) < 10:
        raise RuntimeError("Not enough tubulin CA atoms to estimate a principal axis.")
    return pts

def principal_axis(points: List[Tuple[float,float,float]]) -> Tuple[Tuple[float,float,float], Tuple[float,float,float]]:
    """Return (mean_point, dominant_eigenvector) of covariance via power iteration."""
    n = len(points)
    mx = sum(p[0] for p in points) / n
    my = sum(p[1] for p in points) / n
    mz = sum(p[2] for p in points) / n

    cxx = cyy = czz = cxy = cxz = cyz = 0.0
    for x,y,z in points:
        dx = x - mx; dy = y - my; dz = z - mz
        cxx += dx*dx; cyy += dy*dy; czz += dz*dz
        cxy += dx*dy; cxz += dx*dz; cyz += dy*dz

    inv = 1.0 / n
    cxx *= inv; cyy *= inv; czz *= inv
    cxy *= inv; cxz *= inv; cyz *= inv

    # Power iteration
    vx, vy, vz = 1.0, 0.3, 0.2
    for _ in range(80):
        nx = cxx*vx + cxy*vy + cxz*vz
        ny = cxy*vx + cyy*vy + cyz*vz
        nz = cxz*vx + cyz*vy + czz*vz
        l  = math.sqrt(nx*nx + ny*ny + nz*nz)
        if l < 1e-12:
            break
        vx, vy, vz = nx/l, ny/l, nz/l

    axis = norm((vx,vy,vz))
    mu = (mx,my,mz)
    return mu, axis

def apply_shift_to_chain(pdb_in: str, pdb_out: str, chain_to_shift: str, shift: Tuple[float,float,float]) -> None:
    sx, sy, sz = shift
    with open(pdb_in, "r", encoding="utf-8", errors="replace") as fin, \
         open(pdb_out, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith(("ATOM  ", "HETATM")):
                ch = (line[21:22].strip() or " ")
                if ch == chain_to_shift:
                    x = float(line[30:38]) + sx
                    y = float(line[38:46]) + sy
                    z = float(line[46:54]) + sz
                    # Keep fixed-width PDB columns
                    line = line[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[54:]
            fout.write(line)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_pdb", required=True)
    ap.add_argument("--out_pdb", required=True)
    ap.add_argument("--sep_nm", type=float, default=4.0)
    ap.add_argument("--side", choices=["outer","inner"], default="outer")

    ap.add_argument("--tub_chain", default="A")
    ap.add_argument("--tub_r0", type=int, default=350)
    ap.add_argument("--tub_r1", type=int, default=355)

    ap.add_argument("--cc_chain", default="Z")
    ap.add_argument("--cc_r0", type=int, default=24)
    ap.add_argument("--cc_r1", type=int, default=30)

    args = ap.parse_args()

    tub_patch = read_com(args.in_pdb, args.tub_chain, args.tub_r0, args.tub_r1)
    cc_patch  = read_com(args.in_pdb, args.cc_chain,  args.cc_r0,  args.cc_r1)

    ca = tubulin_ca_points(args.in_pdb, exclude_chain=args.cc_chain)
    mu, axis = principal_axis(ca)

    # Closest point on axis line (through mu) to the tubulin patch
    rel = sub(tub_patch, mu)
    axis_point = add(mu, scale(axis, dot(rel, axis)))

    radial_dir = norm(sub(tub_patch, axis_point))  # axis -> patch
    if args.side == "inner":
        radial_dir = scale(radial_dir, -1.0)

    desired_cc_patch = add(tub_patch, scale(radial_dir, args.sep_nm * 10.0))  # nm -> Å
    shift = sub(desired_cc_patch, cc_patch)

    print("Applying shift to chain", args.cc_chain, "(Å):", shift)
    apply_shift_to_chain(args.in_pdb, args.out_pdb, args.cc_chain, shift)
    print("Wrote", args.out_pdb)

if __name__ == "__main__":
    main()

