#!/bin/bash
# ============================================================================
# Regression + performance benchmark for the FDMSk Smoluchowski solver.
#
# Runs a fixed set of SMALL cases on BOTH kernels (atmospheric + ballistic) and
# for each checks:
#   * mass conservation     (|Total mass - 1| < MASS_TOL)  — physics invariant (all cases)
#   * numerical regression  rel Frobenius vs a committed golden < FROB_TOL  (golden cases)
# and reports wall time + peak RAM (per-process VmRSS via /proc). Prints PASS/FAIL.
# Self-contained: no external data needed (goldens live in sweeps/bench/golden/).
#
# Use after ANY code change to confirm the method still works and is not slower.
#
# Usage:
#   bash sweeps/bench.sh                 # run all cases, compare to golden
#   bash sweeps/bench.sh --make-golden   # (re)generate golden files (after an intended change)
#
# Quick (login-runnable, ~2 min total). Memory here is the small-case footprint;
# for large-scale memory/stability regression use the dedicated slurm sweeps.
# ============================================================================
set -u
# Repo root derived from this script's location, so the benchmark is portable
# (works from any clone, not a hardcoded path).
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
# Solver binary: prefer build/ (the README default), fall back to build_release/.
BIN=${BIN:-$ROOT/build/example}
[ -x "$BIN" ] || BIN=$ROOT/build_release/example
# Cluster runtime libraries (gcc-14.2 libstdc++ + fftw + openblas). Override LD_LIBRARY_PATH
# in the environment if you build elsewhere.
LIBS=/opt/ohpc/pub/compiler/gcc/14.2.0/lib64:/opt/ohpc/pub/libs/gnu14/openmpi5/fftw/3.3.10/lib:/opt/ohpc/pub/libs/gnu14/openblas/0.3.29/lib
export LD_LIBRARY_PATH=$LIBS:${LD_LIBRARY_PATH:-} OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
GOLD=$ROOT/sweeps/bench/golden
WORK=$ROOT/sweeps/bench/work
mkdir -p "$GOLD" "$WORK"

MASS_TOL=${MASS_TOL:-1e-3}     # |mass-1| must be below this (catches blowups like 1.9)
FROB_TOL=${FROB_TOL:-1e-6}     # vs golden (parallel FP noise ~1e-8 << this)

# name        kernel     mosaic    t     max_size init  min_block rel_tol ode_tol nj  check
CASES=(
"atmos_t100   atmospheric tridiag  100   4096     512   128       1e-8    1e-7    4   golden"
"atmos_t1000  atmospheric tridiag  1000  16384    512   128       1e-8    1e-7    4   golden"
"bal_t10      ballistic   monodiag 10    8192     256   128       1e-8    1e-8    4   golden"
"bal_t20      ballistic   monodiag 20    262144   256   256       1e-7    1e-8    4   mass"
)

run_case() {  # fields -> runs solver, echoes "wall peakMB mass finalsize outfile"
    local name=$1 kernel=$2 mosaic=$3 t=$4 max=$5 init=$6 mb=$7 rel=$8 ode=$9 nj=${10}
    local out=$WORK/${name}.txt log=$WORK/${name}.log
    SMOL_KERNEL=$kernel SMOL_MOSAIC=$mosaic MSK_MIN_BLOCK=$mb MSK_REL_TOL=$rel MSK_NJOBS=$nj \
    SMOL_QUIET=1 \
        "$BIN" "$max" "$t" "$out" "$init" 0.01 "$ode" > "$log" 2>&1 &
    local pid=$! peak=0 r
    while kill -0 $pid 2>/dev/null; do
        r=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null)
        [ -n "$r" ] && [ "$r" -gt "$peak" ] && peak=$r
        sleep 0.3
    done
    wait $pid
    local wall mass fsz
    wall=$(grep -hE "Total wall time" "$log" | tail -1 | grep -oE "[0-9.]+")
    mass=$(grep -hE "Total mass" "$log" | tail -1 | grep -oE "[0-9.]+$")
    fsz=$(grep -hE "Final size" "$log" | tail -1 | grep -oE "[0-9]+$")
    echo "${wall:-NA} $((peak/1024)) ${mass:-NA} ${fsz:-NA} $out"
}

if [ "${1:-}" = "--make-golden" ]; then
    echo "Generating golden references..."
    for c in "${CASES[@]}"; do
        set -- $c; [ "${11}" = "golden" ] || continue
        echo "  $1 ..."; read -r w p m f o < <(run_case $c)
        cp "$o" "$GOLD/$1.txt"; echo "    golden saved (mass=$m, wall=${w}s)"
    done
    echo "Done. Commit sweeps/bench/golden/*.txt"
    exit 0
fi

echo "=== FDMSk benchmark  (bin: $BIN)  $(date) ==="
printf "%-13s %-11s %8s %9s %-13s %-14s %s\n" case kernel "wall(s)" "peakRAM" mass check result
printf '%.0s-' {1..86}; echo
pass=0; total=0
for c in "${CASES[@]}"; do
    set -- $c
    name=$1; kernel=$2; check=${11}
    read -r wall peakMB mass fsz outf < <(run_case $c)
    total=$((total+1))
    # mass check (physics invariant)
    massok=$(awk -v m="$mass" -v tol="$MASS_TOL" 'BEGIN{print (m!="NA" && (m-1<tol && 1-m<tol))?1:0}')
    # numerical-regression check
    detail=""; corrok=1
    if [ "$check" = "golden" ]; then
        if [ -f "$GOLD/$name.txt" ]; then
            fr=$(python3 "$ROOT/sweeps/verify.py" "$GOLD/$name.txt" "$outf" 2>/dev/null | awk 'NR==3{print $2}')
            detail="frob=${fr:-NA}"
            corrok=$(awk -v f="${fr:-9}" -v tol="$FROB_TOL" 'BEGIN{print (f+0<tol)?1:0}')
        else detail="no-golden"; corrok=0; fi
    else
        detail="(mass-only)"   # mass check above is the correctness criterion
    fi
    if [ "$massok" = 1 ] && [ "$corrok" = 1 ]; then res="PASS"; pass=$((pass+1)); else res="FAIL"; fi
    [ "$massok" = 1 ] || detail="$detail MASS!"
    printf "%-13s %-11s %8s %7sMB  %-13s %-14s %s\n" "$name" "$kernel" "$wall" "$peakMB" "$mass" "$detail" "$res"
done
printf '%.0s-' {1..86}; echo
echo "OVERALL: $pass/$total passed"
[ "$pass" = "$total" ]
