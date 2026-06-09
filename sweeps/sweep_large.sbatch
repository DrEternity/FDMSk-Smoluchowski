#!/bin/bash
#SBATCH --job-name=smol_lg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=c32m384
#SBATCH --exclusive
#SBATCH --time=40:00
#SBATCH --output=sweeps/sweep_large_%j.out

# Picks production n_jobs at a larger size (T=100000 -> size 262144) and validates
# the sbatch + "srun --cpus-per-task=64 --cpu-bind=none" launch pattern that the
# long t=1e7 job will reuse (plain sbatch pins the batch step to 1 CPU).
set -u
ROOT=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}
LIBS=${SMOL_LIBS:-/opt/ohpc/pub/compiler/gcc/14.2.0/lib64:/opt/ohpc/pub/libs/gnu14/openmpi5/fftw/3.3.10/lib:/opt/ohpc/pub/libs/gnu14/openblas/0.3.29/lib}
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
cd "$ROOT"
MAX=1048576 ; T=100000 ; TOL=1e-6
OUT=$ROOT/sweeps/out ; mkdir -p "$OUT"

echo "=== affinity check (should be 0-63) ==="
srun --cpus-per-task=64 --cpu-bind=none bash -c 'grep Cpus_allowed_list /proc/self/status'

LAUNCH="srun --cpus-per-task=64 --cpu-bind=none env LD_LIBRARY_PATH=$LIBS"
for NJ in 16 32; do
    echo "########## T=$T nj=$NJ  ($(date +%H:%M:%S)) ##########"
    $LAUNCH MSK_NJOBS=$NJ "$ROOT/build_release/example" "$MAX" "$T" "$OUT/out_lg_nj${NJ}.txt" 512 0.01 "$TOL" \
        > "$OUT/lg_nj${NJ}.log" 2>&1
    grep -E "n_jobs|Total wall time|Final size|Total mass|Negative comp" "$OUT/lg_nj${NJ}.log" | sed 's/^/    /'
done
echo "=== DONE $(date) ==="
