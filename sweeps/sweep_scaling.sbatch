#!/bin/bash
#SBATCH --job-name=smol_sweep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=c32m384
#SBATCH --exclusive
#SBATCH --output=sweeps/sweep_scaling_%j.out

# Scaling sweep: tests H1 (MSk n_jobs block-parallelism) and H2 (-O0 vs -O3/-march).
# Workload: T=10000 (reaches size 65536). Baseline O0 single-thread ~84s (ref table).
# All non-MSk thread sources pinned to 1 so the only parallelism axis is MSK_NJOBS.

set -u
ROOT=/home/r.dyachenko/FDMSk-Smoluchowski
export LD_LIBRARY_PATH=/opt/ohpc/pub/compiler/gcc/14.2.0/lib64:${LD_LIBRARY_PATH:-}
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

MAX=1048576
T=10000
TOL=1e-6
OUT=$ROOT/sweeps/out
mkdir -p "$OUT"
cd "$ROOT"

echo "=== HOST: $(hostname)  batch-nproc=$(nproc)  srun-nproc=$(srun --cpus-per-task=64 --cpu-bind=none nproc)  date=$(date) ==="
echo "=== workload: T=$T  ode_tol=$TOL  max_size=$MAX ==="
echo "=== OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS OMP_NUM_THREADS=$OMP_NUM_THREADS ==="

# Launch via srun --cpus-per-task=64 --cpu-bind=none so the process inherits the full 64-core cgroup
# (the batch step itself is affinity-pinned to 1 CPU by task/affinity).
run() {
    local tag=$1 ; shift
    local outfile=$1 ; shift
    local log=$OUT/${tag}.log
    echo "########## $tag ##########"
    srun --cpus-per-task=64 --cpu-bind=none "$@" "$MAX" "$T" "$outfile" 512 0.01 "$TOL" > "$log" 2>&1
    local rc=$?
    if [ $rc -ne 0 ]; then echo "    !!! exit code $rc"; tail -3 "$log" | sed 's/^/    ERR: /'; fi
    grep -E "n_jobs|Total wall time|Final size|Final time|Total mass|Negative comp" "$log" | sed 's/^/    /'
}

# H2: O0 baseline (hardcoded n_jobs=1, no MSK_NJOBS env support)
run "o0_baseline" "$OUT/out_o0.txt" "$ROOT/build/example_O0_baseline"

# H1 + H2: Release -O3 -march scaling over MSk worker threads
for NJ in 1 2 4 8 16 32; do
    export MSK_NJOBS=$NJ
    run "rel_nj${NJ}" "$OUT/out_rel_nj${NJ}.txt" "$ROOT/build_release/example"
done

echo "=== ALL DONE: $(date) ==="
