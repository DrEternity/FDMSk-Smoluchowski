#!/bin/bash
#SBATCH --job-name=smol_nji
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=c32m384
#SBATCH --exclusive
#SBATCH --time=40:00
#SBATCH --output=sweeps/nj_integ_%j.out

# Decisive: is the integration (convolve) memory spike per-thread or per-block?
# Run the resize path to max_size=1048576 (peak ~154GB at nj=16) for several n_jobs,
# sampling node memory every 2s to catch the transient convolve spike. If peak drops
# ~linearly with n_jobs => per-thread => low nj makes t=1e7 feasible.
set -u
ROOT=/home/r.dyachenko/FDMSk-Smoluchowski
LIBS=/opt/ohpc/pub/compiler/gcc/14.2.0/lib64:/opt/ohpc/pub/libs/gnu14/openmpi5/fftw/3.3.10/lib:/opt/ohpc/pub/libs/gnu14/openblas/0.3.29/lib
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$ROOT"
OUT=$ROOT/sweeps/out ; mkdir -p "$OUT"

for NJ in 16 4 1; do
    log=$OUT/nji_nj${NJ}.log
    srun --cpus-per-task=64 --cpu-bind=none env LD_LIBRARY_PATH=$LIBS MSK_NJOBS=$NJ MSK_REL_TOL=1e-5 \
         "$ROOT/build_release/example" 1048576 10000000 "$OUT/nji_out.txt" 512 0.01 1e-6 > "$log" 2>&1 &
    SP=$!
    peak=0
    while kill -0 $SP 2>/dev/null; do
        u=$(free -m | awk '/Mem:/{print $3}')
        [ "$u" -gt "$peak" ] && peak=$u
        sleep 2
    done
    wait $SP; rc=$?
    gb=$(awk "BEGIN{printf \"%.1f\", $peak/1024}")
    mass=$(grep -hE "Total mass" "$log" | tail -1)
    printf "### nj=%-3s peak=%6d MB (%s GB) rc=%s | %s\n" "$NJ" "$peak" "$gb" "$rc" "$mass"
done
echo "=== DONE $(date) ==="
