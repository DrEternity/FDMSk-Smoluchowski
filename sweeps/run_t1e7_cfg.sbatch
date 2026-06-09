#!/bin/bash
#SBATCH --job-name=c7
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=c32m384
#SBATCH --exclusive
#SBATCH --time=4:00:00

# One t=1e7 run for a given (MB, NJ), measuring wall + peak mem + mass, to find the
# optimal memory-acceptable config. Pass MB, NJ via: sbatch --export=ALL,MB=,NJ= \
#   --output=sweeps/cfg7_mb${MB}_nj${NJ}_%j.out run_t1e7_cfg.sbatch
# Memory guard stops cleanly before OOM (so a wrong estimate fails safe).
set -u
ROOT=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}
# Solver binary: prefer build/ (the README default), fall back to build_release/.
BIN=${BIN:-$ROOT/build/example}
[ -x "$BIN" ] || BIN=$ROOT/build_release/example
LIBS=${SMOL_LIBS:-/opt/ohpc/pub/compiler/gcc/14.2.0/lib64:/opt/ohpc/pub/libs/gnu14/openmpi5/fftw/3.3.10/lib:/opt/ohpc/pub/libs/gnu14/openblas/0.3.29/lib}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$ROOT"
MB=${MB:?set MB} ; NJ=${NJ:?set NJ}
MAX=8388608 ; T=10000000 ; RELTOL=1e-5 ; TOL=1e-6 ; GUARD_MB=340000
RUNLOG=$ROOT/sweeps/out/cfg7_mb${MB}_nj${NJ}.log

echo "=== t=1e7 cfg: host=$(hostname) MB=$MB NJ=$NJ rel_tol=$RELTOL date=$(date) ==="
srun --cpus-per-task=64 --cpu-bind=none env LD_LIBRARY_PATH=$LIBS MSK_NJOBS=$NJ MSK_REL_TOL=$RELTOL MSK_MIN_BLOCK=$MB \
     "$BIN" "$MAX" "$T" "/tmp/cfg7_out_${MB}_${NJ}.txt" 512 0.01 "$TOL" > "$RUNLOG" 2>&1 &
SP=$!
peak=0; guarded=0
while kill -0 $SP 2>/dev/null; do
    u=$(free -m | awk '/Mem:/{print $3}'); [ "$u" -gt "$peak" ] && peak=$u
    if [ "$u" -gt "$GUARD_MB" ]; then
        echo "!!! MEMORY GUARD MB=$MB NJ=$NJ: used ${u}MB > ${GUARD_MB}MB -> stop"; guarded=1
        kill $SP 2>/dev/null; sleep 3; kill -9 $SP 2>/dev/null; break
    fi
    sleep 3
done
wait $SP 2>/dev/null; rc=$?
wall=$(grep -hE "Total wall time" "$RUNLOG" | tail -1 | grep -oE "[0-9.]+")
mass=$(grep -hE "Total mass" "$RUNLOG" | tail -1 | grep -oE "[0-9.]+$")
fsz=$(grep -hE "Final size" "$RUNLOG" | tail -1 | grep -oE "[0-9]+$")
gb=$(awk "BEGIN{printf \"%.1f\", $peak/1024}")
echo "#### RESULT MB=$MB NJ=$NJ wall=${wall:-NA}s peak=${gb}GB mass=${mass:-NA} final_sz=${fsz:-NA} rc=$rc guarded=$guarded"
echo "=== DONE $(date) ==="
