#!/bin/bash
#SBATCH --job-name=smol_1e7
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=c32m384
#SBATCH --exclusive
#SBATCH --time=48:00:00
#SBATCH --output=sweeps/run_t1e7_%j.out

# FINAL GOAL: solve to t=1e7 with mass conservation.
# max_size = 2^23 = 8388608: front growth (x2.7 per doubling, ANALYSIS.md Part III)
# reaches ~2.25M < size/2=4.19M at t=1e7, so no particles fall off the grid =>
# mass stays ~1. (max_size=2^20 would stall ~t=3e6 and bleed mass.)
# Launch pattern validated in sweep_large.sbatch (srun --cpus-per-task=64 needed
# or the batch step is pinned to 1 CPU). ode_tol=1e-6: mass conservation is set by
# the grid size, not tol; 1e-6 is ~4x faster than the 1e-9 reference quality.
set -u
ROOT=/home/r.dyachenko/FDMSk-Smoluchowski
LIBS=/opt/ohpc/pub/compiler/gcc/14.2.0/lib64:/opt/ohpc/pub/libs/gnu14/openmpi5/fftw/3.3.10/lib:/opt/ohpc/pub/libs/gnu14/openblas/0.3.29/lib
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
cd "$ROOT"

MAX=8388608        # 2^23
T=10000000         # 1e7
TOL=1e-6
NJ=${NJ:-16}    # sweep optimum at size 65536 AND 262144 (nj=32 regresses; HT+Amdahl)
INIT=512
MASS_THR=0.01
OUTFILE=$ROOT/build/new_reference_solution_atmos_T10000000.txt

echo "=== t=1e7 run: host=$(hostname) date=$(date) ==="
echo "=== max_size=$MAX  T=$T  ode_tol=$TOL  MSK_NJOBS=$NJ ==="
srun --cpus-per-task=64 --cpu-bind=none bash -c 'grep Cpus_allowed_list /proc/self/status'

srun --cpus-per-task=64 --cpu-bind=none env LD_LIBRARY_PATH=$LIBS MSK_NJOBS=$NJ \
     "$ROOT/build_release/example" "$MAX" "$T" "$OUTFILE" "$INIT" "$MASS_THR" "$TOL"

echo "=== t=1e7 DONE: $(date) ==="
