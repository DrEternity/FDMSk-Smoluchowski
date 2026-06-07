#!/bin/bash
#SBATCH --job-name=smol_pk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=c32m384
#SBATCH --exclusive
#SBATCH --time=25:00
#SBATCH --output=sweeps/peak_mem_%j.out

# Is the memory peak per-thread (n_jobs lever works) or per-N (won't help)?
# Build a size directly (initial=max=size) and watch memory. Record peak during the
# approximate() BUILD, then keep sampling ~30s into the first RHS (where the real run
# OOM'd), then KILL (the full integration to T would take forever at large size).
set -u
ROOT=/home/r.dyachenko/FDMSk-Smoluchowski
LIBS=/opt/ohpc/pub/compiler/gcc/14.2.0/lib64:/opt/ohpc/pub/libs/gnu14/openmpi5/fftw/3.3.10/lib:/opt/ohpc/pub/libs/gnu14/openblas/0.3.29/lib
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$ROOT"
OUT=$ROOT/sweeps/out ; mkdir -p "$OUT"
echo "node total: $(free -g | awk '/Mem:/{print $2}') GB"

probe() {
    local SIZE=$1 NJ=$2
    local log=$OUT/pk_${SIZE}_nj${NJ}.log
    srun --cpus-per-task=64 --cpu-bind=none env LD_LIBRARY_PATH=$LIBS MSK_NJOBS=$NJ MSK_REL_TOL=1e-5 \
         "$ROOT/build_release/example" "$SIZE" 1.0 "$OUT/pk_out.txt" "$SIZE" 0.01 1e-6 > "$log" 2>&1 &
    local SP=$! peak_build=0 peak=0 u phase=build rhs_t=0
    while kill -0 $SP 2>/dev/null; do
        u=$(free -m | awk '/Mem:/{print $3}'); [ "$u" -gt "$peak" ] && peak=$u
        if [ "$phase" = build ]; then
            peak_build=$peak
            grep -q "Integrating from" "$log" 2>/dev/null && phase=rhs
        else
            rhs_t=$((rhs_t+1)); [ "$rhs_t" -ge 15 ] && { kill $SP 2>/dev/null; break; }
        fi
        sleep 2
    done
    kill $SP 2>/dev/null; wait $SP 2>/dev/null
    local bt=$(grep -hoE "Approximation time: [0-9.]+" "$log" | head -1)
    printf "### size=%-8s nj=%-3s build_peak=%5d MB (%.0fGB)  build+rhs_peak=%5d MB (%.0fGB)  %s\n" \
        "$SIZE" "$NJ" "$peak_build" "$(awk "BEGIN{print $peak_build/1024}")" \
        "$peak" "$(awk "BEGIN{print $peak/1024}")" "$bt"
}

for SIZE in 1048576 2097152; do
    for NJ in 16 4 1; do probe $SIZE $NJ; done
done
echo "=== DONE $(date) ==="
