#!/bin/bash
#SBATCH --job-name=smol_1e7
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=c32m384
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --output=sweeps/run_t1e7_%j.out
#
# THE production run: atmospheric kernel out to t=1e7 with mass conserved.
# Best configuration from ANALYSIS.md §3.4 (the t=1e7 config sweep). Submit from the
# repo root:  sbatch sweeps/run_t1e7.sbatch
#
# Why these parameters (all leave the answer unchanged — Frobenius vs baseline ~1e-7):
#   max_size = 2^23  the front reaches only ~2.3M by t=1e7 (ANALYSIS §2), so the run
#                    ends at size 2^22 with the front well below the grid top => mass
#                    has nowhere to leak and stays conserved.
#   min_block = 1024 fewer, larger MSk blocks => ~3x less convolve memory AND ~2x faster
#                    (ANALYSIS §3.3). 2048 is already slower.
#   rel_tol = 1e-5   looser kernel-approximation tolerance => ~16x smaller build scratch,
#                    no accuracy loss (ANALYSIS §3.3). NOT the same knob as ode_tol.
#   n_jobs = 4       fastest setting that still fits in 376 GB (~33 min, ~280 GB peak).
#                    The convolve buffer pile-up grows with n_jobs, so nj>=8 OOMs at 2^23.
#                    For more headroom: nj=3 (~40 min, ~208 GB) or nj=2 (~54 min, ~152 GB).
#   ode_tol = 1e-6   the atmospheric kernel is mild; mass conservation here is set by the
#                    grid size, not ode_tol (contrast the stiff ballistic kernel, §4).
#
# A MEMORY GUARD stops the run cleanly if usage approaches the node limit, so a wrong
# estimate fails safe (clean stop with a "downshift nj" hint) instead of an OOM-kill.
# Override any knob from the command line, e.g.:  sbatch --export=ALL,NJ=2 sweeps/run_t1e7.sbatch
set -u
ROOT=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}
# Solver binary: prefer build/ (the README default), fall back to build_release/.
BIN=${BIN:-$ROOT/build/example}
[ -x "$BIN" ] || BIN=$ROOT/build_release/example
# Cluster runtime libraries (gcc-14.2 libstdc++ + fftw + openblas); override SMOL_LIBS elsewhere.
LIBS=${SMOL_LIBS:-/opt/ohpc/pub/compiler/gcc/14.2.0/lib64:/opt/ohpc/pub/libs/gnu14/openmpi5/fftw/3.3.10/lib:/opt/ohpc/pub/libs/gnu14/openblas/0.3.29/lib}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$ROOT"

MAX=${MAX:-8388608}            # 2^23
T=${T:-10000000}               # 1e7
TOL=${TOL:-1e-6}               # ode_tol
RELTOL=${RELTOL:-1e-5}
NJ=${NJ:-4}
MB=${MB:-1024}
GUARD_MB=${GUARD_MB:-350000}   # stop cleanly if used memory exceeds ~350 GB
OUTFILE=${OUTFILE:-$ROOT/build/new_reference_solution_atmos_T10000000.txt}
RUNLOG=$ROOT/sweeps/out/t1e7_run.log
mkdir -p "$ROOT/sweeps/out"

echo "=== t=1e7 atmospheric run: host=$(hostname) date=$(date) ==="
echo "=== binary: $BIN ==="
echo "=== max_size=$MAX T=$T ode_tol=$TOL rel_tol=$RELTOL n_jobs=$NJ min_block=$MB ==="
echo "=== node total mem: $(free -g | awk '/Mem:/{print $2}') GB, guard=${GUARD_MB}MB ==="

srun --cpus-per-task=64 --cpu-bind=none env LD_LIBRARY_PATH=$LIBS \
     MSK_NJOBS=$NJ MSK_REL_TOL=$RELTOL MSK_MIN_BLOCK=$MB \
     "$BIN" "$MAX" "$T" "$OUTFILE" 512 0.01 "$TOL" > "$RUNLOG" 2>&1 &
SP=$!
peak=0
while kill -0 $SP 2>/dev/null; do
    u=$(free -m | awk '/Mem:/{print $3}'); [ "$u" -gt "$peak" ] && peak=$u
    if [ "$u" -gt "$GUARD_MB" ]; then
        echo "!!! MEMORY GUARD: used ${u}MB > ${GUARD_MB}MB -> stopping cleanly (downshift NJ or raise MB)"
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
