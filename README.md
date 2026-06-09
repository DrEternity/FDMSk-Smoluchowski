# FDMSk-Smoluchowski

Deterministic solver for the **Smoluchowski coagulation equation** using a
**mosaic-skeleton (MSk)** low-rank compression of the dense coagulation-kernel matrix.
Companion code for the paper *"Mosaic-skeleton approximation is all you need for
Smoluchowski equations"* ([arXiv:2501.10206](https://arxiv.org/abs/2501.10206)).

The equation is integrated as a system of ODEs for the size concentrations `n_k(t)`
(monodisperse start: `n[0]=1`). Each right-hand side evaluation is a **gain**
(`½·Σ_{i+j=k} K(i,j) n_i n_j`, a convolution) minus a **loss** (`n_k·Σ_j K(k,j) n_j`,
a matvec); both use the MSk-compressed kernel `K`. An adaptive RK4 stepper grows the
grid (`size`) on the fly as the particle "front" advances.

---

## Repository layout

```
src/                C++ sources
  main.cpp          the `example` program — solver / reference-solution generator
  smoluchowski.hpp  the solver: modeling() + adaptive RK4 + RHS operator
  config.hpp        SolverConfig — ALL hyperparameters in one place
  bench.cpp         the `bench` program — MSk matvec/convolve micro-benchmark
CMakeLists.txt      build (fetches the MSk library `zaimsk` automatically)
README.md           this file
ANALYSIS.md         optimization journal + Part VII = guide for a NEW kernel (read this!)
OPTIMIZATION_STORY_MONTECARLO.md   the sister Monte-Carlo project's story
sweeps/             slurm experiment scripts, the benchmark, and python helpers
  bench.sh          regression+perf benchmark (run after any change)
  run_t1e7.sbatch   large production run (atmospheric, t=1e7)
  run_ballistic.sbatch   ballistic-kernel runs
  *.sbatch          parameter sweeps (memory, threads, min_block)
  *.py              verify.py / mc_compare.py / to_sizeconc.py (comparison helpers)
external/  build/   fetched library / build trees (git-ignored)
```

---

## Built on

The numerical heavy lifting is **not** in this repo — it lives in the mosaic-skeleton
library, which CMake fetches and builds automatically. This project is the thin
Smoluchowski layer on top of it.

| dependency | role here |
|---|---|
| **zaimsk** — the Mosaic-Skeleton (MSk) library ([gitlab.com/bulatral/mosaic-skeleton](https://gitlab.com/bulatral/mosaic-skeleton), branch `dev`) | the core: low-rank/block compression of the dense `N×N` kernel matrix, the FFT-based block **convolution** (gain term), the **matvec** (loss term), and the worker thread pool (`MSK_NJOBS`). Fetched via CMake `FetchContent` on the first configure. We also use its helper headers `smart_array` (RAII buffers) and `cpp_blas` (the `BLAS::axpy/copy/nrm2` used by the RK4 stepper). **We do not modify this library.** |
| **FFTW3** | the FFTs behind MSk's block convolution (built with `MSK_USE_FFT=ON`) |
| **BLAS + LAPACK** | dense / low-rank block linear algebra and the matvec inside MSk |
| **pthreads** | MSk's block-parallel thread pool |

So the only first-class code in this repository is `src/` (the solver + config + the
`example`/`bench` programs); everything mathematically heavy is delegated to `zaimsk`,
which in turn stands on FFTW3 / BLAS / LAPACK. MPI is disabled (`MSK_USE_MPI=OFF`).

## Requirements

* C, C++ **and Fortran** compilers (the MSk library has a Fortran component)
* CMake ≥ 3.20
* **BLAS, LAPACK, FFTW3, pthreads** (system libraries linked by `zaimsk`)
* `git` + network access on the **first** configure (CMake clones `zaimsk` into `external/`)

---

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target example -j
```

This produces `build/example` (and `build/bench`). The build defaults to an optimized
`Release` and tunes for the host CPU (`-march=native`); pass `-DSMOL_NATIVE=OFF` for a
portable binary.

<details>
<summary><b>Notes for the `nieve` cluster</b> (where this work was done)</summary>

The system CMake/GCC are too old; use the modules and set the runtime library path:

```bash
export PATH=/opt/ohpc/pub/utils/cmake/4.1.2/bin:/opt/ohpc/pub/compiler/gcc/14.2.0/bin:$PATH
export LD_LIBRARY_PATH=/opt/ohpc/pub/compiler/gcc/14.2.0/lib64:\
/opt/ohpc/pub/libs/gnu14/openmpi5/fftw/3.3.10/lib:\
/opt/ohpc/pub/libs/gnu14/openblas/0.3.29/lib:$LD_LIBRARY_PATH

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target example -j 16
```

`LD_LIBRARY_PATH` is also needed at **run** time. Long/large runs go through SLURM —
see the ready-made scripts in `sweeps/` (they pin threads and add a memory guard).
</details>

---

## Run

```bash
cd build
./example                     # default: atmospheric kernel, small demo run
```

All parameters are read by `SolverConfig` (see [`src/config.hpp`](src/config.hpp)) with
the precedence **built-in default → positional CLI → environment variable** (env wins,
which makes parameter sweeps trivial). The program prints the full configuration at
startup.

### Hyperparameters (environment variables — the recommended interface)

| group | variable | meaning | default |
|---|---|---|---|
| problem | `SMOL_MAX_SIZE` | grid cap (power of two) | 65536 |
| | `SMOL_TIME` | final time `T` | 1000 |
| | `SMOL_INIT_SIZE` | starting grid size | 512 |
| | `SMOL_KERNEL` | `atmospheric` or `ballistic` | atmospheric |
| | `SMOL_MOSAIC` | `tridiag` (ρ=2) or `monodiag` (ρ=1) | tridiag |
| | `SMOL_OUTPUT` | output file path | reference_solution_atmos.txt |
| integrator | `SMOL_DT` | initial RK4 step | 1e-4 |
| | `SMOL_ODE_TOL` | step-error tolerance (↓ for stiff kernels!) | 1e-6 |
| | `SMOL_MASS_GUARD` | abort if mass exceeds initial × this (instability guard; ≤1 disables) | 1.02 |
| mosaic-MSk | `MSK_REL_TOL` | block approximation tolerance | 1e-10 |
| | `MSK_MIN_BLOCK` | min mosaic block (memory/speed lever) | 128 |
| | `MSK_MAX_RANK` | block rank cap (0 = unlimited) | 0 |
| | `MSK_NJOBS` | MSk worker threads | 1 |

The legacy positional form is still accepted (same order as before):
`./example  max_size  time  output  initial_size  <unused>  ode_tol`.

### Examples

```bash
# Atmospheric kernel to t=1000 on a 16384 grid, 8 threads, output to a file
SMOL_KERNEL=atmospheric SMOL_MAX_SIZE=16384 SMOL_TIME=1000 \
MSK_NJOBS=8 SMOL_OUTPUT=atmos_t1000.txt ./example

# Ballistic kernel to t=20  (this kernel is STIFF — use a small ode_tol, see below)
SMOL_KERNEL=ballistic SMOL_MOSAIC=monodiag SMOL_MAX_SIZE=262144 SMOL_TIME=20 \
SMOL_ODE_TOL=1e-8 MSK_NJOBS=4 SMOL_OUTPUT=ball_t20.txt ./example
```

The output file lists the concentrations `n[k]` (size `k+1`) one per line, after a short
header. `sweeps/to_sizeconc.py` converts it to a two-column `size concentration` table.

---

## Verify a change (benchmark)

After any code change, run the regression+performance benchmark:

```bash
bash sweeps/bench.sh
```

It runs both kernels on small cases (~1.5 min) and checks **mass conservation** plus
**correctness** (atmospheric vs a committed golden, ballistic vs a Monte-Carlo
reference), reporting wall time and peak RAM with an overall `PASS/FAIL` and exit code.

---

## ⚠️ Before a large run — read this

Memory grows **explosively** (~×3 per grid doubling) and the explicit RK4 can go
**unstable** for fast-growing ("stiff") kernels — both bit us hard. The hard-won
playbook for choosing `max_size`, `ode_tol`, `MSK_NJOBS`, `MSK_MIN_BLOCK` and for
running a **new kernel** safely is in **[ANALYSIS.md](ANALYSIS.md), Part VII** (with a
dedicated memory-caution section and a during-run control checklist). In short:
classify the kernel by its homogeneity `λ`; probe the front/memory on a *small* run
first; then launch via the `sweeps/` SLURM scripts with the built-in memory guard.

---

## License

See [LICENSE](LICENSE).
