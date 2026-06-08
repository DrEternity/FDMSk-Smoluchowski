#!/bin/bash
#SBATCH --job-name=v6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=c32m384
#SBATCH --exclusive
#SBATCH --time=30:00
#SBATCH --output=sweeps/verify_t1e6_%j.out
set -u
ROOT=/home/r.dyachenko/FDMSk-Smoluchowski
LIBS=/opt/ohpc/pub/compiler/gcc/14.2.0/lib64:/opt/ohpc/pub/libs/gnu14/openmpi5/fftw/3.3.10/lib:/opt/ohpc/pub/libs/gnu14/openblas/0.3.29/lib
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$ROOT"
# New (refactored) binary, same config as grid run g6_out_mb1024_nj16
srun --cpus-per-task=64 --cpu-bind=none env LD_LIBRARY_PATH=$LIBS MSK_NJOBS=16 MSK_REL_TOL=1e-5 MSK_MIN_BLOCK=1024 \
     "$ROOT/build_release/example" 1048576 1000000 "$ROOT/sweeps/out/verify_t1e6_new.txt" 512 0.01 1e-6
echo "=== compare new(refactored) vs grid(old binary) at t=1e6 ==="
python3 sweeps/verify.py sweeps/out/g6_out_mb1024_nj16.txt sweeps/out/verify_t1e6_new.txt
echo "=== DONE $(date) ==="
