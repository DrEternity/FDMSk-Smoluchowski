#!/bin/bash
#SBATCH --job-name=smol_1e7
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=c32m384
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --output=sweeps/run_t1e7_%j.out

# FINAL GOAL: solve to t=1e7 with mass conservation.
#
# Memory wall solved by TWO levers (ANALYSIS.md H8'/H9), both with mass identical and
# accuracy intact (Frobenius vs baseline ~1e-7):
#   - MSK_REL_TOL=1e-5: shrinks the build scratch (~rank).
#   - MSK_MIN_BLOCK=1024: fewer blocks => far fewer convolve buffers (3x less mem, 2x faster).
#   - MSK_NJOBS=2: kills the convolve buffer pile-up (nj=1 ->18GB vs nj=16 ->150GB @2^20).
# Combined, size 2^22 (the largest reached by t=1e7) fits in ~200GB << 376GB.
# Attempts #1 (rel_tol=1e-10) and #2 (rel_tol=1e-5, nj=16) both OOM'd at 2^21.
#
# max_size=2^23: front reaches only ~2M by t=1e7 (Part III'), so the run ends at size
# 2^22 with front << grid top => mass conserved. A MEMORY GUARD stops the run cleanly
# if usage approaches the node limit (so a wrong estimate fails safe, not via OOM-kill).
set -u
ROOT=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}
LIBS=${SMOL_LIBS:-/opt/ohpc/pub/compiler/gcc/14.2.0/lib64:/opt/ohpc/pub/libs/gnu14/openmpi5/fftw/3.3.10/lib:/opt/ohpc/pub/libs/gnu14/openblas/0.3.29/lib}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$ROOT"

MAX=8388608        # 2^23
T=10000000         # 1e7
TOL=1e-6           # ODE tol (mass conservation is set by grid size, not this)
RELTOL=${RELTOL:-1e-5}
# Optimum from the t=1e7 config sweep (ANALYSIS.md): mb=1024 minimises memory AND time;
# nj=4 is the fastest that fits (1976s, peak 280GB). For more headroom use nj=3 (2380s,
# 208GB) or nj=2 (3223s, 152GB). nj>=8 OOMs at 2^23 (convolve pile-up).
NJ=${NJ:-4}
MB=${MB:-1024}
GUARD_MB=${GUARD_MB:-350000}   # stop cleanly if used memory exceeds ~350GB
OUTFILE=$ROOT/build/new_reference_solution_atmos_T10000000.txt
RUNLOG=$ROOT/sweeps/out/t1e7_run.log

echo "=== t=1e7 run: host=$(hostname) date=$(date) ==="
echo "=== max_size=$MAX T=$T ode_tol=$TOL MSK_REL_TOL=$RELTOL MSK_NJOBS=$NJ MSK_MIN_BLOCK=$MB ==="
echo "=== node total mem: $(free -g | awk '/Mem:/{print $2}') GB, guard=${GUARD_MB}MB ==="

srun --cpus-per-task=64 --cpu-bind=none env LD_LIBRARY_PATH=$LIBS MSK_NJOBS=$NJ MSK_REL_TOL=$RELTOL MSK_MIN_BLOCK=$MB \
     "$ROOT/build_release/example" "$MAX" "$T" "$OUTFILE" 512 0.01 "$TOL" > "$RUNLOG" 2>&1 &
SP=$!
peak=0
while kill -0 $SP 2>/dev/null; do
    u=$(free -m | awk '/Mem:/{print $3}'); [ "$u" -gt "$peak" ] && peak=$u
    if [ "$u" -gt "$GUARD_MB" ]; then
        echo "!!! MEMORY GUARD: used ${u}MB > ${GUARD_MB}MB -> stopping run cleanly (downshift nj/raise mb)"
        kill $SP 2>/dev/null; sleep 5; kill -9 $SP 2>/dev/null
        break
    fi
    sz=$(grep -h "Building MSk for size" "$RUNLOG" 2>/dev/null | tail -1 | grep -oE "[0-9]+")
    mass=$(grep -hE "mass = " "$RUNLOG" 2>/dev/null | tail -1 | grep -oE "mass = [0-9.eE+-]+")
    printf "[%s] used=%6dMB peak=%6dMB size=%-8s %s\n" "$(date +%H:%M:%S)" "$u" "$peak" "${sz:-?}" "${mass:-}"
    sleep 15
done
wait $SP 2>/dev/null; rc=$?
echo "=== example exit rc=$rc, peak_used=${peak}MB ($(awk "BEGIN{printf \"%.1f\",$peak/1024}")GB) ==="
grep -hE "Integration complete|Final time|Final size|Final mass|Total mass|Negative comp" "$RUNLOG" | tail -10
echo "=== t=1e7 DONE: $(date) ==="
