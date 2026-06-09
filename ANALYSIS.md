# ANALYSIS — distilled experiments

Summary of the work on speeding up the deterministic Smoluchowski solver (FDM +
mosaic-skeleton compression of the kernel matrix, library `zaimsk`). Goals: reach
**t = 10⁷ with mass conservation** and find efficient configurations. Below is the useful
distillation only: key mechanisms, numbers and tables. **The practical guide for running
a NEW kernel is in §6** (read it before a large run).

Two kernels were studied, differing in homogeneity `λ` (`K(ai,aj)=a^λ·K(i,j)`):
**atmospheric** `λ=−5/9` and **ballistic** `λ=5/6`. `λ` predicts both the front-growth
speed (memory) and the stiffness (integrator stability) — see §6.

---

## 1. Problem setup (brief)

- Coagulation as a system of ODEs for the size concentrations `n_k(t)`; monodisperse
  start `n[0]=1`.
- RHS = **gain** (convolution `Σ_{i+j=k} K(i,j)n_i n_j`, via FFT — `convolve`) − **loss**
  (matvec `n_k·Σ_j K(k,j)n_j` — `matvec`).
- The `N×N` kernel matrix is dense but smooth off-diagonal → compressed by
  mosaic-skeleton: Toeplitz blocks (FFT) + low-rank/dense blocks (BLAS).
- Integrator — **adaptive RK4** with step control on the local error (`ode_tol`).
  The grid is **doubled** when the "front" (where 99.9% of the mass accumulates)
  reaches `size/2`.
- **Mass `Σ n_k·k` is the key invariant and the acceptance criterion.** For `λ≤1` it is
  conserved (≈1); any GROWTH above 1 is a numerical artifact (instability), see §4.

---

## 2. Reference solutions and cost profile

Two reference series (in `build/`, up to 1048576 points), differing ONLY in `ode_tol`:

| series | ode_tol | purpose |
|---|---|---|
| `reference_solution_atmos_T*` | 1e-6 | fast run |
| `new_reference_solution_atmos_T*` | 1e-9 | accurate reference (~4× slower) |

Timings of the original (−O0, single-thread) run:

| T | reference (1e-6) | new_reference (1e-9) | final size |
|---:|---:|---:|---:|
| 10 | 4.9 s | 171 s | 512 |
| 100 | 6.9 s | 216 s | 2048 |
| 1000 | 14.6 s | 318 s | 8192 |
| 10000 | 84 s | 795 s | 65536 |
| 100000 | 351 s | 1968 s | 262144 |
| 1000000 | 1835 s | 7165 s | 1048576 |

**Cost profile (t=10⁶):** time splits **≈ evenly** between building the MSk approximation
(approx) and integration (the `convolve`/`matvec` calls) — both phases need optimizing.
Inside the RHS, `convolve` (gain) ≫ `matvec` (loss) (~18× at size 16384). Large epochs
dominate: the last two (524288 + 1048576) ≈ **78%** of total time.

**Front growth (atmospheric kernel) → choice of `max_size`.** The time `t` at which the
front reaches `size/2` and triggers a resize grows **×~2.7 per size doubling**
(characteristic size ∝ t^0.7):

| size | 2¹⁹ | 2²⁰ | 2²¹ | 2²² | 2²³ |
|---|---|---|---|---|---|
| front→size/2 at t≈ | 4.6e5 | 1.2e6 | 3.3e6 | 9e6 | 2.4e7 |

⇒ by t=10⁷ the front reaches ~2.3M; final size = 2²² (a resize into 2²³ would happen only
at t≈2.4e7). Mass is conserved as long as the front stays below the top of the grid.

---

## 3. Optimization (atmospheric kernel)

### 3.1 MSk parallelism + compiler flags (T=10⁴, size→65536)

The baseline repository built with `-O0` and ran in a **single thread**
(`approximate(...,1)`). Clean measurement on an exclusive slurm node (mass 0.99999902 in
every row, Frobenius between configs ≲ 5e-13 — the physics does not change):

| config | wall, s | × vs −O0 |
|---|---:|---:|
| −O0, 1 thread | 75.4 | 1.00 |
| Release `-O3 -march=native`, 1 thread | 43.8 | 1.72 |
| + nj=4 | 17.3 | 4.36 |
| + nj=8 | 13.3 | 5.67 |
| **+ nj=16** | **12.9** | **5.85** |
| + nj=32 | 14.3 | 5.27 (regression) |

**−O3/−march = 1.72×; MSk parallelism = a further ~3.4×.** At this size the optimum is
nj≈16; nj=32 regresses (HT contention + Amdahl on the serial BLAS-1 operations of RK4).
At larger sizes the parallel fraction grows (already 7.06× vs −O0 at size 262144).

### 3.2 Efficient config for t=10⁶ (min_block × n_jobs grid)

Memory is not the constraint here (everything fits) ⇒ the target is wall time. Mass is
0.999938 everywhere, Frobenius between configs ≤4e-7 (they do NOT change the answer).
**wall, s / peak, GB:**

| min_block \ nj | 4 | 8 | 16 | 32 |
|---|---|---|---|---|
| 128 | 367/109 | 312/136 | **297**/150 | 325/157 |
| 256 | 281/60 | 229/77 | **204**/85 | 217/89 |
| 512 | 318/31 | 225/49 | 187/58 | **175**/60 |
| 1024 | 429/36 | 275/39 | 200/48 | **174/50** ⭐ |
| 2048 | 673/59 | 404/59 | 265/60 | 206/66 |

**Optimum for t=10⁶: `min_block=1024, nj=32` → 174 s, 50 GB** (×1.7 vs the default
mb=128/nj=16; ~×10 vs the original −O0/1-thread). `min_block` is a lever for both memory
and speed (but 2048 is already slower: dense diagonal blocks too large). The optimal nj
GROWS with min_block (fewer blocks → no convolve pile-up → threads scale).

### 3.3 Memory — the wall and the levers (critical for size ≥ 2²¹)

Peak memory is the **`convolve` transient**: each block writes into its own output buffer
of length `N−(i0+j0)`, and the buffers pile up "in flight" (many producer threads, one
accumulator). Peak ≈ (#blocks)×(pile-up factor of n_jobs), and it grows **×~3 per size
doubling** — this, not the compressed matrix (a few GB), is the memory killer. At
`rel_tol=1e-10` a huge build scratch (∝ rank) adds on top. Measurements at size 1048576
(mass identical, accuracy unaffected):

| lever | values → peak memory |
|---|---|
| **n_jobs** (pile-up) | 16 → 150 GB; 4 → 111 GB; **1 → 18 GB** |
| **min_block** (#blocks) | 128 → 150 GB; 512 → 61; **1024 → 51 GB** (+ ×2 speed) |
| **rel_tol** (build scratch) | 1e-10 → ~470 GB@2²¹ (OOM); **1e-5 → ×16 smaller**, Frobenius vs 1e-10 = 1e-7 |

Conclusion: **large size → LOW nj** (2–4), **`min_block=1024`**, **`rel_tol=1e-5`**.
(Mirror image of t=10⁶, where the system fits and a high nj is best.)

### 3.4 t=10⁷ — configs and reaching the goal

Configs at size 2²³ (mass 0.99977 everywhere; nj≥8 → OOM from convolve pile-up):

| min_block | nj | wall | peak, GB | verdict |
|---:|---:|---:|---:|---|
| **1024** | **4** | **33 min** | **280** | ⭐ fastest that fits |
| 1024 | 3 | 40 min | 208 | conservative |
| 1024 | 2 | 54 min | 152 | max memory headroom |
| 512 | 4 | — | 333 | ✗ OOM |

🎯 **t=10⁷ achieved** (`max_size=2²³, rel_tol=1e-5, min_block=1024, nj=2`, memory guard):
final time 10⁷, final size 2²³, **Total mass = 0.99977** (drift 2.3e-4 — mass conserved),
0 negatives, peak 152 GB, ~54 min. The mass drift accumulates over 10⁷ of time
(asymmetry of the approximated kernel + clamping of negatives); on the t≤10⁶ references
it was ~1e-6.

> Memory GREW during the run for two reasons: (1) step-ups at each resize (convolve
> buffer ∝ size); (2) glibc-arena ratchet (freed multi-MB buffers are not returned to the
> OS). Mitigation if needed:
> `MALLOC_MMAP_THRESHOLD_=131072 MALLOC_TRIM_THRESHOLD_=131072`.

---

## 4. Ballistic kernel (t = 10, 20, 40, 50)

The front grows **very fast (~t^4.7)**: t=10 → ~1.2K, t=20 → ~65K, t=40 → ~0.4M (2²⁰),
t=50 → ~1.3M (2²²). Config: `mosaic=monodiag (rho=1)`, otherwise like the atmospheric one.

**🔴 Main lesson (stiffness): the ballistic kernel needs a SMALL `ode_tol`.** The first
t=40/50 run at `ode_tol=1e-7` gave a **MASS BLOW-UP** (1.90 and 2.44 instead of 1.0). This
is not physics (λ=5/6<1, no gelation) but an instability of explicit RK4: the controller
watches ACCURACY, not STABILITY → the step grew past the stability limit → oscillations
with negatives → the `y<0→0` clamp adds mass → runaway. **Fix: `ode_tol=1e-10`** kept the
step bounded ⇒ mass 0.99997. (The atmospheric kernel was fine even at 1e-6 — ballistic is
stiffer.) This lesson is built into the solver as a **mass guard** (`SMOL_MASS_GUARD`, §6).

**Results against an external reference solution** (mass conserved):

| t | rel L2 (density) | mass | final size | peak GB | ode_tol |
|---:|---:|---:|---:|---:|---:|
| 10 | 0.34% | 0.99999 | 4096 | <1 | 1e-7 |
| 20 | 0.57% | 0.99998 | 65536 | ~1 | 1e-7 |
| 40 | 1.70% | 0.99997 | 2²⁰ | 16 | 1e-10 |
| 50 | 2.28% | 0.99997 | 2²² | 95 | 1e-10 |

Bulk agreement is good (det/reference ratios ~0.95–1.06). The large mass-weighted L2
(t=50) is the statistical noise of the reference in the tail (single particles at size
~10⁶), not a solver error.

---

## 5. Methodology (brief)

- **Hardware:** Intel Xeon Gold 6226R (Cascade Lake), 32 physical cores (64 HT), AVX-512,
  376 GB RAM, partition `c32m384`. Build host = run host ⇒ `-march=native` is safe.
- **Where to measure:** only on a dedicated slurm node (the login node is noisy). Each
  measurement is a separate, exclusive sbatch.
- **Correctness:** Frobenius `||n−n_ref||/||n_ref||` ≲ 1e-6 (the ODE-tolerance scale) plus
  mass conservation. Parallel reduction / `-march` perturb FP by ~1e-13, not the physics.
- **Cluster gotchas (for reproducibility):** an `srun` step inside `sbatch` needs an
  EXPLICIT `--cpus-per-task=64 --cpu-bind=none` (otherwise it is pinned to 1 core);
  `/usr/bin/time` is absent on the nodes — measure from the program's output; rebuilt
  binaries need `LD_LIBRARY_PATH` (gcc-14.2 libstdc++ + openblas + fftw); sample peak
  memory frequently (near the peak it jumps >16 GB between samples).

---

## 6. GUIDE: running a simulation for a NEW kernel

The unifying parameter is the **homogeneity `λ`** (`K(ai,aj)=a^λ K(i,j)`). Compute it
FIRST: it predicts both the front-growth speed (memory) and the stiffness (stability).

> ⚠️ **The most dangerous thing in practice is MEMORY.** It grows explosively (~×3 per
> size doubling) and "blows up in an instant" — it caused all our failures (OOM, runaway).
> Rule: **measure first (a cheap probe of the front and peak memory), then launch — and
> always under a memory guard.**

**Step 0 — classify by `λ`:**
- `λ > 1` ⇒ **gelation**: mass physically decreases after t_gel — mass-as-invariant is no
  longer the criterion (a separate gel-fraction accounting is needed). We had none such.
- `λ ≤ 1` ⇒ mass is conserved ⇒ `Total mass ≈ 1` is the hard criterion; GROWTH >1 = artifact.
- The larger `λ`, the faster the front grows AND the stiffer the system.

**Step 1 — stability (`ode_tol`):** explicit RK4 is stable only for a small enough step,
but the controller watches accuracy, not stability. The instability symptom is **`mass`
growing above 1** (accelerating) + the front overshooting the real one. Rule: the larger
`λ`/size, the smaller `ode_tol` (atmospheric was fine at 1e-6; ballistic needed 1e-10).
Start at `ode_tol≈1e-9…1e-10`. The solver has a built-in **🛡 mass guard**
(`SMOL_MASS_GUARD`, default 1.02): the run aborts with a "reduce ode_tol" message before a
full runaway.

**Step 2 — `max_size` (front ⇄ mass):** mass is conserved only if the front (99.9% of the
mass) stays below the top of the grid. Estimate the front at the target `t` with a cheap
small run / extrapolation, and take `max_size` with margin. Memory grows ×~3 per doubling
— `max_size` is memory lever #1.

**Step 3 — memory (`n_jobs`, `min_block`, `rel_tol`):**
- `min_block` ↑ ⇒ fewer blocks ⇒ less memory and often faster. **Default 1024** (2048 is
  already slower). Does not change the answer.
- `n_jobs` is regime-dependent: system fits in RAM (size ≤ ~2²⁰) → HIGH nj (16–32); large
  (≥2²¹) → LOW nj (2–4), otherwise the convolve pile-up causes OOM.
- `rel_tol` looser (≤1e-5) ⇒ smaller build scratch, no accuracy loss. **Do not confuse
  with `ode_tol`** (rel_tol — kernel-approximation accuracy; ode_tol — time-stepping
  stability).
- Always set a **memory guard** (stop at ~340 of 376 GB).

**Step 4 — mosaic (`rho`):** chosen for the kernel's singularity near the diagonal
(atmospheric — `tridiag` ρ=2, ballistic — `monodiag` ρ=1). At a fixed `rel_tol` it has
little effect on accuracy.

**During-run control checklist:**

| metric | normal | alarm |
|---|---|---|
| `Total mass` | ≈1, stable | GROWTH >1 ⇒ instability (↓ode_tol); strong drop ⇒ leakage (↑max_size) |
| `Negative components` | 0 / tiny | many ⇒ oscillations (↓ode_tol) |
| `Final size` vs expected front | matches | overshoot ⇒ fake mass / instability |
| peak RAM | < guard | approaching guard ⇒ ↓nj or ↑min_block |
| `step` over time | grows moderately, plateaus | unbounded growth ⇒ runaway imminent |

**Recipe:** (1) compute `λ`; (2) cheap recon run at small `t` — tune `ode_tol` until mass
stops growing, measure the front and peak memory; (3) choose `max_size` with margin;
(4) `min_block=1024`, `rel_tol=1e-5`, `n_jobs` per the memory budget; (5) launch via slurm
with a memory guard + live monitoring of mass/size/RAM; (6) compare against a reference
(`sweeps/ref_compare.py`). After any code change — `bash sweeps/bench.sh` (regression).

**Calibration from the two kernels:**

| | atmospheric | ballistic |
|---|---|---|
| `λ` | −5/9 (decreasing) | 5/6 (increasing) |
| front growth | slow (t=10⁷ → 2.3M) | fast (t=50 → 1.3M, ~t^4.7) |
| stiffness / `ode_tol` | mild, 1e-6 OK | stiff, needs 1e-10 |
| mosaic | tridiag (ρ=2) | monodiag (ρ=1) |
| mass result | 0.99977 @ t=10⁷ | 0.99997 @ t≤50 |
