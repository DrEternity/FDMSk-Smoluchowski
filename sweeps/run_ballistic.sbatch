#!/bin/bash
#SBATCH --job-name=bal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=c32m384
#SBATCH --exclusive
#SBATCH --time=4:00:00

# Deterministic ballistic-kernel solve (for comparison against an external reference).
# Ballistic front grows very fast (~t^4.7) -> watch the size. Pass T and MAX via env:
#   sbatch --export=ALL,T=,MAX= --output=sweeps/ballistic/bal_t${T}_%j.out run_ballistic.sbatch
set -u
ROOT=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}
# Solver binary: prefer build/ (the README default), fall back to build_release/.
BIN=${BIN:-$ROOT/build/example}
[ -x "$BIN" ] || BIN=$ROOT/build_release/example
LIBS=${SMOL_LIBS:-/opt/ohpc/pub/compiler/gcc/14.2.0/lib64:/opt/ohpc/pub/libs/gnu14/openmpi5/fftw/3.3.10/lib:/opt/ohpc/pub/libs/gnu14/openblas/0.3.29/lib}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$ROOT"
T=${T:?set T} ; MAX=${MAX:?set MAX}
NJ=${NJ:-4} ; MB=${MB:-1024} ; RELTOL=${RELTOL:-1e-6} ; ODETOL=${ODETOL:-1e-9} ; GUARD_MB=${GUARD_MB:-340000}
TAG=${TAG:-}    # optional suffix to isolate parallel runs writing the same T
OUT=$ROOT/sweeps/ballistic/det_ballistic_t${T}${TAG}.txt
RUNLOG=$ROOT/sweeps/ballistic/run_t${T}${TAG}.log

echo "=== ballistic t=$T host=$(hostname) max_size=$MAX nj=$NJ mb=$MB rel_tol=$RELTOL ode_tol=$ODETOL $(date) ==="
srun --cpus-per-task=64 --cpu-bind=none env LD_LIBRARY_PATH=$LIBS \
     SMOL_KERNEL=ballistic SMOL_MOSAIC=monodiag MSK_NJOBS=$NJ MSK_REL_TOL=$RELTOL MSK_MIN_BLOCK=$MB \
     "$BIN" "$MAX" "$T" "$OUT" 256 0.01 "$ODETOL" > "$RUNLOG" 2>&1 &
SP=$!
peak=0
while kill -0 $SP 2>/dev/null; do
    u=$(free -m | awk '/Mem:/{print $3}'); [ "$u" -gt "$peak" ] && peak=$u
    if [ "$u" -gt "$GUARD_MB" ]; then echo "!!! MEMORY GUARD t=$T: ${u}MB>${GUARD_MB}MB -> stop"; kill $SP 2>/dev/null; sleep 3; kill -9 $SP 2>/dev/null; break; fi
    sz=$(grep -h "Building MSk for size" "$RUNLOG" 2>/dev/null | tail -1 | grep -oE "[0-9]+")
    printf "[%s] used=%6dMB peak=%6dMB size=%s\n" "$(date +%H:%M:%S)" "$u" "$peak" "${sz:-?}"
    sleep 10
done
wait $SP 2>/dev/null; rc=$?
gb=$(awk "BEGIN{printf \"%.1f\",$peak/1024}")
echo "#### BALLISTIC t=$T rc=$rc peak=${gb}GB $(grep -hE 'Final size|Total mass' "$RUNLOG" | tail -2 | tr '\n' ' ')"
echo "=== DONE $(date) ==="
