#!/bin/bash
#SBATCH --job-name=smol_lev
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=c32m384
#SBATCH --exclusive
#SBATCH --time=40:00
#SBATCH --output=sweeps/mem_lever_%j.out

# Find a memory-feasible rel_tol for large max_size while conserving mass.
# Fixed size: max_size=524288 (baseline 47GB peak at rel_tol=1e-10). For each
# rel_tol, run to the can't-resize cap (fast), sample peak node memory, and read
# the final mass (the acceptance criterion). Exclusive node => node mem ~= job mem.
set -u
ROOT=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}
LIBS=${SMOL_LIBS:-/opt/ohpc/pub/compiler/gcc/14.2.0/lib64:/opt/ohpc/pub/libs/gnu14/openmpi5/fftw/3.3.10/lib:/opt/ohpc/pub/libs/gnu14/openblas/0.3.29/lib}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$ROOT"
MAX=524288 ; T=10000000 ; NJ=16
OUT=$ROOT/sweeps/out ; mkdir -p "$OUT"

for TOL in 1e-10 1e-7 1e-5 1e-3; do
    log=$OUT/lev_tol${TOL}.log
    srun --cpus-per-task=64 --cpu-bind=none env LD_LIBRARY_PATH=$LIBS MSK_NJOBS=$NJ MSK_REL_TOL=$TOL \
         "$ROOT/build_release/example" "$MAX" "$T" "$OUT/lev_out_${TOL}.txt" 512 0.01 1e-6 > "$log" 2>&1 &
    SP=$!
    peak=0
    while kill -0 $SP 2>/dev/null; do
        u=$(free -m | awk '/Mem:/{print $3}')
        [ "$u" -gt "$peak" ] && peak=$u
        sleep 3
    done
    wait $SP
    mass=$(grep -hE "Total mass" "$log" | tail -1)
    comp=$(grep -hE "Compression" "$log" | tail -1)
    wall=$(grep -hE "Total wall time" "$log" | tail -1)
    gb=$(awk "BEGIN{printf \"%.1f\", $peak/1024}")
    printf "### rel_tol=%-6s peak=%6d MB (%s GB) | %s | %s | %s\n" \
        "$TOL" "$peak" "$gb" "$mass" "$comp" "$wall"
done
echo "=== DONE $(date) ==="
