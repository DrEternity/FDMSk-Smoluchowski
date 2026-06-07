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
# KEY FIX (ANALYSIS.md H8): rel_tol=1e-5 instead of 1e-10. The MSk build allocates
# per-thread scratch proportional to block RANK; at 1e-10 ranks are huge => ~470GB
# at size 2^21 => OOM (killed attempt #1, job 44950). At 1e-5 the same build is ~29GB,
# while accuracy is unchanged (Frobenius 1e-7 vs 1e-10, mass identical 0.9999995).
#
# max_size=2^23: front reaches only ~2M by t=1e7 (Part III'), so the run ends at
# size 2^22 with front << grid top => mass conserved. 2^23 is cheap headroom now and
# guards against faster-than-expected front growth. nj=16 (memory-independent; speed).
set -u
ROOT=/home/r.dyachenko/FDMSk-Smoluchowski
LIBS=/opt/ohpc/pub/compiler/gcc/14.2.0/lib64:/opt/ohpc/pub/libs/gnu14/openmpi5/fftw/3.3.10/lib:/opt/ohpc/pub/libs/gnu14/openblas/0.3.29/lib
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$ROOT"

MAX=8388608        # 2^23
T=10000000         # 1e7
TOL=1e-6           # ODE tol (mass conservation is set by grid size, not this)
RELTOL=${RELTOL:-1e-5}
NJ=${NJ:-16}
OUTFILE=$ROOT/build/new_reference_solution_atmos_T10000000.txt
RUNLOG=$ROOT/sweeps/out/t1e7_run.log

echo "=== t=1e7 run: host=$(hostname) date=$(date) ==="
echo "=== max_size=$MAX T=$T ode_tol=$TOL MSK_REL_TOL=$RELTOL MSK_NJOBS=$NJ ==="
echo "=== node total mem: $(free -g | awk '/Mem:/{print $2}') GB ==="

srun --cpus-per-task=64 --cpu-bind=none env LD_LIBRARY_PATH=$LIBS MSK_NJOBS=$NJ MSK_REL_TOL=$RELTOL \
     "$ROOT/build_release/example" "$MAX" "$T" "$OUTFILE" 512 0.01 "$TOL" > "$RUNLOG" 2>&1 &
SP=$!
peak=0
while kill -0 $SP 2>/dev/null; do
    u=$(free -m | awk '/Mem:/{print $3}'); [ "$u" -gt "$peak" ] && peak=$u
    sz=$(grep -h "Building MSk for size" "$RUNLOG" 2>/dev/null | tail -1 | grep -oE "[0-9]+")
    mass=$(grep -hE "mass = " "$RUNLOG" 2>/dev/null | tail -1 | grep -oE "mass = [0-9.eE+-]+")
    printf "[%s] used=%6dMB peak=%6dMB size=%-8s %s\n" "$(date +%H:%M:%S)" "$u" "$peak" "${sz:-?}" "${mass:-}"
    sleep 15
done
wait $SP; rc=$?
echo "=== example exit rc=$rc, peak_used=${peak}MB ($(awk "BEGIN{printf \"%.1f\",$peak/1024}")GB) ==="
grep -hE "Integration complete|Final time|Final size|Final mass|Total mass|Negative comp" "$RUNLOG" | tail -10
echo "=== t=1e7 DONE: $(date) ==="
