#!/usr/bin/env python3
# Compare a deterministic FDMSk solution to a Monte-Carlo reference (pure python).
#   MC file (.dat): lines "size  count  conc" (skip '#'); conc = number density n_size.
#   Det file: save_solution format — '#' lines, then "size time", then n[0..size-1],
#             where n[k] is the density of size k+1.
# Both are mass-normalised (sum_i i*n_i = 1), so they compare directly by size.
import sys, math

def load_mc(path):
    d = {}; mx = 0
    with open(path) as f:
        for ln in f:
            if ln.startswith('#') or not ln.strip():
                continue
            p = ln.split()
            s = int(p[0]); d[s] = float(p[2])
            if s > mx: mx = s
    return d, mx

def load_det(path):
    with open(path) as f:
        data = [l for l in f if not l.startswith('#')]
    return [float(x) for x in data[1:] if x.strip()]   # vals[k] = density of size k+1

mc, mc_max = load_mc(sys.argv[1])
det = load_det(sys.argv[2])                              # det density of size s = det[s-1]
S = max(len(det), mc_max)

mass_det = mass_mc = 0.0
num = den = 0.0            # density rel-L2
mnum = mden = 0.0         # mass-weighted rel-L2
for s in range(1, S + 1):
    a = det[s-1] if s-1 < len(det) else 0.0
    b = mc.get(s, 0.0)
    mass_det += s * a; mass_mc += s * b
    num += (a - b) ** 2; den += b * b
    mnum += (s * (a - b)) ** 2; mden += (s * b) ** 2

print(f"  det file : {sys.argv[2]}")
print(f"  mc  file : {sys.argv[1]}")
print(f"  mass: det={mass_det:.6f}  mc={mass_mc:.6f}")
print(f"  rel L2 (density)       = {math.sqrt(num/den):.4e}")
print(f"  rel L2 (mass-weighted) = {math.sqrt(mnum/mden):.4e}")
print(f"  {'size':>8} {'det n_i':>14} {'mc n_i':>14} {'ratio':>8}")
for s in [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000]:
    if s <= S:
        a = det[s-1] if s-1 < len(det) else 0.0
        b = mc.get(s, 0.0)
        r = (a / b) if b > 0 else float('nan')
        print(f"  {s:8d} {a:14.6e} {b:14.6e} {r:8.3f}")
