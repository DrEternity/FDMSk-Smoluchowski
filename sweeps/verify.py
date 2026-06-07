#!/usr/bin/env python3
"""Compare solver outputs (pure-python, no numpy).
Reports relative Frobenius vs a reference file and total mass = sum n[i]*(i+1).
Usage: verify.py <ref.txt> <file1.txt> [file2.txt ...]
File format: 4 comment lines, 1 'size time' line, then values."""
import sys, math

def load(f):
    with open(f) as fh:
        lines = fh.readlines()
    return [float(x) for x in lines[5:] if x.strip()]

def mass(v):
    return sum(v[i]*(i+1) for i in range(len(v)))

ref = load(sys.argv[1])
dref = math.sqrt(sum(x*x for x in ref))
print(f"{'file':<32} {'rel_frobenius':>14} {'mass':>14} {'mass_err':>12}")
print(f"{sys.argv[1].split('/')[-1]:<32} {0.0:>14.3e} {mass(ref):>14.10f} {'(ref)':>12}")
for f in sys.argv[2:]:
    v = load(f)
    n = min(len(v), len(ref))
    num = math.sqrt(sum((v[i]-ref[i])**2 for i in range(n)))
    rel = num/dref
    m = mass(v)
    print(f"{f.split('/')[-1]:<32} {rel:>14.3e} {m:>14.10f} {abs(m-1.0):>12.3e}")
