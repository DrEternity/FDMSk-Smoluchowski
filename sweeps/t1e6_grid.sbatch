#!/bin/bash
#SBATCH --job-name=smol_g6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=c32m384
#SBATCH --exclusive
#SBATCH --time=3:00:00
#SBATCH --output=sweeps/t1e6_grid_%j.out

# Efficient-config search for t=1e6 (reaches size 1048576): the memory x threads x
# min_block trade-off. t=1e6 fits in RAM at all configs, so the target is WALL TIME;
# peak memory + mass are tracked too. rel_tol=1e-5 fixed (accuracy-neutral, faster).
set -u
ROOT=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}
# Solver binary: prefer build/ (the README default), fall back to build_release/.
BIN=${BIN:-$ROOT/build/example}
[ -x "$BIN" ] || BIN=$ROOT/build_release/example
LIBS=${SMOL_LIBS:-/opt/ohpc/pub/compiler/gcc/14.2.0/lib64:/opt/ohpc/pub/libs/gnu14/openmpi5/fftw/3.3.10/lib:/opt/ohpc/pub/libs/gnu14/openblas/0.3.29/lib}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$ROOT"
OUT=$ROOT/sweeps/out ; mkdir -p "$OUT"
MAX=1048576 ; T=1000000

echo "=== t=1e6 grid: host=$(hostname) date=$(date) (rel_tol=1e-5) ==="
printf "%-9s %-5s %-12s %-9s %-12s %-9s\n" "min_block" "nj" "wall(s)" "peak(GB)" "mass" "final_sz"

run() {
    local MB=$1 NJ=$2
    local log=$OUT/g6_mb${MB}_nj${NJ}.log
    srun --cpus-per-task=64 --cpu-bind=none env LD_LIBRARY_PATH=$LIBS MSK_NJOBS=$NJ MSK_REL_TOL=1e-5 MSK_MIN_BLOCK=$MB \
         "$BIN" "$MAX" "$T" "$OUT/g6_out_mb${MB}_nj${NJ}.txt" 512 0.01 1e-6 > "$log" 2>&1 &
    local SP=$!
    local peak=0 u
    while kill -0 $SP 2>/dev/null; do
        u=$(free -m | awk '/Mem:/{print $3}'); [ "$u" -gt "$peak" ] && peak=$u
        sleep 3
    done
    wait $SP
    local wall=$(grep -hE "Total wall time" "$log" | tail -1 | grep -oE "[0-9.]+")
    local mass=$(grep -hE "Total mass" "$log" | tail -1 | grep -oE "[0-9.]+$")
    local fsz=$(grep -hE "Final size" "$log" | tail -1 | grep -oE "[0-9]+$")
    local gb=$(awk "BEGIN{printf \"%.1f\", $peak/1024}")
    printf "%-9s %-5s %-12s %-9s %-12s %-9s\n" "$MB" "$NJ" "${wall:-ERR}" "$gb" "${mass:-?}" "${fsz:-?}"
}

for MB in 128 256 512 1024 2048; do
    for NJ in 4 8 16 32; do
        run $MB $NJ
    done
done
echo "=== DONE $(date) ==="
