#!/bin/bash
# Sweep body — run via srun-as-allocator so the process gets full 64-core affinity:
#   srun --partition=c32m384 --exclusive --ntasks=1 --cpus-per-task=64 \
#        --time=40:00 --cpu-bind=none bash sweeps/sweep_body.sh > sweeps/<log> 2>&1 &
# Tests H1 (MSk n_jobs block-parallelism) and H2 (-O0 vs -O3/-march).

set -u
ROOT=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}
# Rebuilt binaries (build_release/build_o0) lost the rpath the original had, so we
# must add FFTW + OpenBLAS lib dirs explicitly (plus gcc-14 libstdc++).
export LD_LIBRARY_PATH=${SMOL_LIBS:-/opt/ohpc/pub/compiler/gcc/14.2.0/lib64:/opt/ohpc/pub/libs/gnu14/openmpi5/fftw/3.3.10/lib:/opt/ohpc/pub/libs/gnu14/openblas/0.3.29/lib}:${LD_LIBRARY_PATH:-}
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

MAX=1048576
T=${T:-10000}
TOL=${TOL:-1e-6}
OUT=$ROOT/sweeps/out
mkdir -p "$OUT"
cd "$ROOT"

echo "=== HOST: $(hostname)  nproc=$(nproc)  allowed=$(grep -o 'Cpus_allowed_list.*' /proc/self/status)  date=$(date) ==="
echo "=== workload: T=$T  ode_tol=$TOL  max_size=$MAX ==="
echo "=== OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS OMP_NUM_THREADS=$OMP_NUM_THREADS ==="

run() {
    local tag=$1 ; shift
    local outfile=$1 ; shift
    local log=$OUT/${tag}.log
    echo "########## $tag  ($(date +%H:%M:%S)) ##########"
    "$@" "$MAX" "$T" "$outfile" 512 0.01 "$TOL" > "$log" 2>&1
    local rc=$?
    if [ $rc -ne 0 ]; then echo "    !!! exit code $rc"; tail -3 "$log" | sed 's/^/    ERR: /'; fi
    grep -E "n_jobs|Total wall time|Final size|Final time|Total mass|Negative comp" "$log" | sed 's/^/    /'
}

# H2: O0 baseline (atmospheric, empty build type = no -O3; single MSk worker)
export MSK_NJOBS=1
run "o0_baseline" "$OUT/out_o0.txt" "$ROOT/build_o0/example"

# H1 + H2: Release -O3 -march scaling over MSk worker threads
for NJ in 1 2 4 8 16 32; do
    export MSK_NJOBS=$NJ
    run "rel_nj${NJ}" "$OUT/out_rel_nj${NJ}.txt" "$ROOT/build_release/example"
done

echo "=== ALL DONE: $(date) ==="
