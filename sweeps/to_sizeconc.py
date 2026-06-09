# Convert solver output (det) to 2-column "size concentration" (nonzero only).
import sys
with open(sys.argv[1]) as f:
    lines = [l for l in f if not l.startswith('#')]
vals = [float(x) for x in lines[1:] if x.strip()]   # vals[k] = density of size k+1
out = sys.argv[2]
with open(out, 'w') as g:
    g.write("# size  concentration  (deterministic FDMSk, ballistic kernel)\n")
    for k, v in enumerate(vals):
        if v > 0.0:
            g.write(f"{k+1} {v:.16e}\n")
print(f"wrote {out}")
