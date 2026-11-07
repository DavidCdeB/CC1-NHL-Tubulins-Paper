#!/usr/bin/env python3
import math
import argparse

def read_atoms(pdb, ca_only=False):
    atoms=[]
    with open(pdb) as f:
        for line in f:
            if not line.startswith(("ATOM  ","HETATM")):
                continue
            if ca_only and line[12:16].strip() != "CA":
                continue
            ch = line[21:22] or " "
            rn = int(line[22:26])
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            atoms.append((ch, rn, (x,y,z)))
    return atoms

def com(atoms, chain=None, r0=None, r1=None):
    xs=ys=zs=0.0; n=0
    for ch,rn,(x,y,z) in atoms:
        if chain is not None and ch != chain: continue
        if r0 is not None and not (r0 <= rn <= r1): continue
        xs+=x; ys+=y; zs+=z; n+=1
    if n==0:
        raise RuntimeError("Empty selection for COM (check chain/res range).")
    return (xs/n, ys/n, zs/n)

def dot(a,b): return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def sub(a,b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def add(a,b): return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
def scale(a,s): return (a[0]*s, a[1]*s, a[2]*s)
def norm(a):
    n=math.sqrt(dot(a,a))
    if n < 1e-12: raise RuntimeError("Zero vector.")
    return scale(a, 1.0/n)

def principal_axis(points):
    # covariance + power iteration
    n=len(points)
    mx=sum(p[0] for p in points)/n
    my=sum(p[1] for p in points)/n
    mz=sum(p[2] for p in points)/n
    cxx=cyy=czz=cxy=cxz=cyz=0.0
    for x,y,z in points:
        dx=x-mx; dy=y-my; dz=z-mz
        cxx+=dx*dx; cyy+=dy*dy; czz+=dz*dz
        cxy+=dx*dy; cxz+=dx*dz; cyz+=dy*dz
    inv=1.0/n
    cxx*=inv; cyy*=inv; czz*=inv; cxy*=inv; cxz*=inv; cyz*=inv

    vx,vy,vz = 1.0,0.3,0.2
    for _ in range(60):
        nx=cxx*vx + cxy*vy + cxz*vz
        ny=cxy*vx + cyy*vy + cyz*vz
        nz=cxz*vx + cyz*vy + czz*vz
        l=math.sqrt(nx*nx+ny*ny+nz*nz)
        if l < 1e-12: break
        vx,vy,vz = nx/l, ny/l, nz/l
    return (vx,vy,vz)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tub", default="tubulin_gmx.pdb")
    ap.add_argument("--cc1", default="cc1_gmx.pdb")
    ap.add_argument("--sep_nm", type=float, default=2.0)
    ap.add_argument("--side", choices=["outer","inner"], default="outer")
    args = ap.parse_args()

    tub_all = read_atoms(args.tub, ca_only=False)
    tub_ca  = read_atoms(args.tub, ca_only=True)
    cc1_all = read_atoms(args.cc1, ca_only=False)

    tub_com = com(tub_all)
    tub_patch_com = com(tub_all, chain="A", r0=350, r1=355)
    cc1_patch_com = com(cc1_all, chain="Z", r0=24, r1=30)

    axis = norm(principal_axis([p for _,_,p in tub_ca]))
    # closest point on axis (through tub_com) to the patch:
    rel = sub(tub_patch_com, tub_com)
    axis_point = add(tub_com, scale(axis, dot(rel, axis)))

    radial_dir = norm(sub(tub_patch_com, axis_point))  # axis -> patch
    if args.side == "inner":
        radial_dir = scale(radial_dir, -1.0)

    target_cc1_patch = add(tub_patch_com, scale(radial_dir, args.sep_nm))
    shift = sub(target_cc1_patch, cc1_patch_com)

    print(f"{shift[0]:.6f} {shift[1]:.6f} {shift[2]:.6f}")

if __name__ == "__main__":
    main()

