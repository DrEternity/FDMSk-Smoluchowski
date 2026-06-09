#!/bin/bash
#SBATCH --job-name=smol_mem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=c32m384
#SBATCH --exclusive
#SBATCH --time=30:00
#SBATCH --output=sweeps/mem_probe_%j.out

# Diagnose the t=1e7 OOM: run example until it reaches/dies at size 2^21-2^22,
# sampling node memory every 4s and tagging it with the current MSk size.
# Exclusive node => node memory used ~= this job's usage.
set -u
ROOT=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}
LIBS=${SMOL_LIBS:-/opt/ohpc/pub/compiler/gcc/14.2.0/lib64:/opt/ohpc/pub/libs/gnu14/openmpi5/fftw/3.3.10/lib:/opt/ohpc/pub/libs/gnu14/openblas/0.3.29/lib}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$ROOT"
RUNLOG=$ROOT/sweeps/out/mem_run.log

echo "=== node total mem ==="; free -g | awk '/Mem:/{print "total="$2"GB"}'
srun --cpus-per-task=64 --cpu-bind=none env LD_LIBRARY_PATH=$LIBS MSK_NJOBS=16 \
     "$ROOT/build_release/example" 8388608 4000000 "$ROOT/sweeps/out/mem_out.txt" 512 0.01 1e-6 \
     > "$RUNLOG" 2>&1 &
SP=$!
peak=0
while kill -0 $SP 2>/dev/null; do
    used=$(free -m | awk '/Mem:/{print $3}')
    [ "$used" -gt "$peak" ] && peak=$used
    sz=$(grep -h "Building MSk for size" "$RUNLOG" 2>/dev/null | tail -1 | grep -oE "[0-9]+")
    printf "[%s] size=%-8s used=%6d MB  peak=%6d MB\n" "$(date +%H:%M:%S)" "${sz:-?}" "$used" "$peak"
    sleep 4
done
wait $SP
echo "=== run ended, peak_used=${peak}MB ==="
grep -hE "Building MSk for size|Approximation time|Compression|Killed|Total mass|Front reached" "$RUNLOG" | tail -20
