#!/bin/bash
#SBATCH --job-name=smol_mb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=c32m384
#SBATCH --exclusive
#SBATCH --time=45:00
#SBATCH --output=sweeps/mb_sweep_%j.out

# H9: does larger min_block cut the convolve memory peak? peak ~ 8B*sum_blocks(N-(i0+j0))
# and #blocks ~ 1/min_block. Fixed: max_size=1048576 (peak ~150GB at mb=128/nj16),
# nj=16, rel_tol=1e-5. Vary min_block; record peak mem, mass, wall, compression.
# Outputs saved per mb for an accuracy cross-check afterwards.
set -u
ROOT=/home/r.dyachenko/FDMSk-Smoluchowski
LIBS=/opt/ohpc/pub/compiler/gcc/14.2.0/lib64:/opt/ohpc/pub/libs/gnu14/openmpi5/fftw/3.3.10/lib:/opt/ohpc/pub/libs/gnu14/openblas/0.3.29/lib
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$ROOT"
OUT=$ROOT/sweeps/out ; mkdir -p "$OUT"

for MB in 128 256 512 1024; do
    log=$OUT/mb_${MB}.log
    srun --cpus-per-task=64 --cpu-bind=none env LD_LIBRARY_PATH=$LIBS MSK_NJOBS=16 MSK_REL_TOL=1e-5 MSK_MIN_BLOCK=$MB \
         "$ROOT/build_release/example" 1048576 10000000 "$OUT/mb_out_${MB}.txt" 512 0.01 1e-6 > "$log" 2>&1 &
    SP=$!
    peak=0
    while kill -0 $SP 2>/dev/null; do
        u=$(free -m | awk '/Mem:/{print $3}')
        [ "$u" -gt "$peak" ] && peak=$u
        sleep 2
    done
    wait $SP; rc=$?
    gb=$(awk "BEGIN{printf \"%.1f\", $peak/1024}")
    mass=$(grep -hE "Total mass" "$log" | tail -1 | grep -oE "[0-9.]+$")
    wall=$(grep -hE "Total wall time" "$log" | tail -1 | grep -oE "[0-9.]+ sec")
    comp=$(grep -hE "Compression" "$log" | tail -1 | grep -oE "[0-9.]+x")
    printf "### min_block=%-5s peak=%6d MB (%s GB) rc=%s mass=%s comp(last)=%s wall=%s\n" \
        "$MB" "$peak" "$gb" "$rc" "${mass:-?}" "${comp:-?}" "${wall:-?}"
done
echo "=== DONE $(date) ==="
