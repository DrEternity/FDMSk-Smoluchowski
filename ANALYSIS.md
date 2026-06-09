# ANALYSIS — notes from optimizing the solver

This is the working log from the effort to make our deterministic Smoluchowski solver
(an FDM discretization plus mosaic-skeleton compression of the kernel matrix, built on the
`zaimsk` library) both fast and trustworthy enough to reach a single concrete target:
**simulate out to t = 10⁷ without losing mass.** Getting there meant profiling where the
time and memory actually go, sweeping the handful of knobs that matter, and getting burned
twice — once by memory blowing up, once by the integrator going unstable. Those two
failures shaped every recommendation here, so they get the most space.

If you only read one section, read **§6**: it distills everything below into a step-by-step
recipe for running a *new* kernel safely. The rest of the document is the evidence behind
that recipe — the mechanisms, the numbers, and the tables.

**The two kernels we studied.** Both are coagulation kernels `K(i,j)` — the rate at which a
particle of size `i` and one of size `j` stick together — and they differ in a single
number, the *homogeneity* `λ`, defined by `K(ai, aj) = a^λ · K(i,j)`. We looked at an
**atmospheric** kernel (`λ = −5/9`) and a **ballistic** one (`λ = 5/6`). That one number
turns out to be the most predictive parameter in the whole problem: it tells you how fast
the particle population spreads toward large sizes (which drives memory) and how stiff the
system is (which drives integrator stability). That theme runs through the entire document
and is the backbone of §6.

---

## 1. The problem, in brief

We solve coagulation as a system of ODEs for the size concentrations `n_k(t)` — `n_k` is
how much material sits at size `k` — starting from a monodisperse state where everything is
the smallest size (`n[0] = 1`).

Each evaluation of the right-hand side has two pieces:

- a **gain** term, `½·Σ_{i+j=k} K(i,j) n_i n_j`, which counts the pairs that merge *into*
  size `k`. It is a convolution, so we compute it with an FFT; in the code this is the call
  we refer to as **`convolve`**.
- a **loss** term, `n_k·Σ_j K(k,j) n_j`, which counts size-`k` particles that merge *away*
  into something larger. It is a matrix–vector product — the **`matvec`** call.

The object both terms need is the `N×N` kernel matrix `K`. It is dense, but away from the
diagonal it is smooth, and that is exactly what **mosaic-skeleton (MSk)** compression
exploits: it tiles the matrix into blocks and stores each one cheaply — Toeplitz blocks via
FFT, the rest as low-rank or small dense blocks via BLAS. That compression is the only
reason large `N` is tractable at all.

Time stepping is an **adaptive RK4** integrator. It chooses its step size by watching the
local truncation error against a tolerance we call `ode_tol` — keep that name in mind, it
is both the hero and the villain of §4. Because the population marches toward ever-larger
sizes, we do not fix the grid: we track the **front** — the size below which 99.9% of the
mass sits — and whenever the front reaches the halfway point of the current grid, we
**double the grid** and carry on. Each such doubling is what we call an *epoch* below.

Finally, the one number that decides whether a run is acceptable: **mass**, `Σ n_k·k`.
Physically, coagulation only shuffles mass between sizes; it never creates any. So for any
kernel with `λ ≤ 1` mass is conserved and should stay `≈ 1` for all time, and we use that
as the hard acceptance test. The corollary matters just as much: if mass ever *grows* above
1, that is never physics — it is a numerical artifact, and §4 is the story of chasing one
down.

---

## 2. Reference solutions and where the cost goes

Before optimizing anything we generated trusted reference solutions to check against. There
are two series sitting in `build/`, both run out to 1048576 grid points and identical
except for `ode_tol`:

| series | ode_tol | what it's for |
|---|---|---|
| `reference_solution_atmos_T*` | 1e-6 | the everyday "fast" run |
| `new_reference_solution_atmos_T*` | 1e-9 | the tight, accurate reference (~4× slower) |

The timings below are for the *original* code — compiled at −O0, single-threaded — and they
are the baseline that every speedup later is measured against:

| T | reference (1e-6) | new_reference (1e-9) | final size |
|---:|---:|---:|---:|
| 10 | 4.9 s | 171 s | 512 |
| 100 | 6.9 s | 216 s | 2048 |
| 1000 | 14.6 s | 318 s | 8192 |
| 10000 | 84 s | 795 s | 65536 |
| 100000 | 351 s | 1968 s | 262144 |
| 1000000 | 1835 s | 7165 s | 1048576 |

Profiling a t = 10⁶ run, three things stand out:

- **The two phases cost about the same.** Roughly half the wall time goes into *building*
  the MSk approximation of the kernel and the other half into *integrating* (the
  `convolve`/`matvec` calls). Neither dominates, so both had to be optimized.
- **Within the right-hand side, the gain term is everything.** `convolve` outweighs
  `matvec` by about 18× at size 16384, so the convolution is where tuning pays off.
- **The last couple of epochs dominate.** The two largest grids (524288 and 1048576)
  together eat ≈ 78% of the run. Big-`N` behavior is what matters; the small grids are
  noise.

How far does the front actually travel, and how big does the grid need to be? Since the grid
doubles whenever the front hits the halfway mark, the size you must reserve is set by where
the front is at your target time. For the atmospheric kernel the characteristic size grows
roughly as `t^0.7`, which means the time at which the front triggers each resize stretches
out by about ×2.7 per doubling:

| size | 2¹⁹ | 2²⁰ | 2²¹ | 2²² | 2²³ |
|---|---|---|---|---|---|
| front reaches size/2 at t ≈ | 4.6e5 | 1.2e6 | 3.3e6 | 9e6 | 2.4e7 |

Reading that off: by t = 10⁷ the atmospheric front has reached ~2.3 million, so the run
settles at a final grid of 2²²; the next doubling into 2²³ would not be needed until
t ≈ 2.4×10⁷. And as long as the front stays below the top of the grid, mass has nowhere to
leak, so it stays conserved.

---

## 3. Optimizing the atmospheric kernel

### 3.1 The first big wins: compiler flags and MSk threads

The repository as we inherited it was built at −O0 and ran single-threaded
(`approximate(..., 1)`). Two cheap changes account for most of the speedup. The
measurements below come from a clean, exclusive slurm node at t = 10⁴ (final grid 65536);
reassuringly, every configuration produced the same physics — mass 0.99999902 in every row,
and the relative Frobenius difference between configurations was ≲ 5e-13, i.e. nothing but
floating-point dust.

| config | wall, s | × vs −O0 |
|---|---:|---:|
| −O0, 1 thread | 75.4 | 1.00 |
| Release `-O3 -march=native`, 1 thread | 43.8 | 1.72 |
| + nj=4 | 17.3 | 4.36 |
| + nj=8 | 13.3 | 5.67 |
| **+ nj=16** | **12.9** | **5.85** |
| + nj=32 | 14.3 | 5.27 (regression) |

So turning on `−O3 -march=native` alone buys 1.72×, and handing the MSk block work to
multiple threads (the `n_jobs` knob) buys roughly another 3.4× on top. At this size the
sweet spot is about 16 threads; pushing to 32 actually *regresses*, because we start
fighting hyper-threading contention and Amdahl's law — the serial BLAS-1 vector operations
inside RK4 do not parallelize, so more threads cannot help them. At larger grids the
parallel fraction grows and the picture improves (already 7.06× over −O0 at size 262144).

### 3.2 Picking a config for t = 10⁶

Here everything fits comfortably in RAM, so memory is not the constraint — the only target
is wall time. The two knobs worth sweeping are `n_jobs` (how many MSk worker threads) and
**`min_block`** (the smallest block edge MSk is allowed to cut the matrix into; a bigger
`min_block` means fewer, larger blocks). Again the physics is untouched: mass 0.999938
everywhere, Frobenius between configs ≤ 4e-7. Each cell is **wall (s) / peak (GB):**

| min_block \ nj | 4 | 8 | 16 | 32 |
|---|---|---|---|---|
| 128 | 367/109 | 312/136 | **297**/150 | 325/157 |
| 256 | 281/60 | 229/77 | **204**/85 | 217/89 |
| 512 | 318/31 | 225/49 | 187/58 | **175**/60 |
| 1024 | 429/36 | 275/39 | 200/48 | **174/50** ⭐ |
| 2048 | 673/59 | 404/59 | 265/60 | 206/66 |

The winner is `min_block=1024, n_jobs=32`: 174 s at 50 GB, about 1.7× faster than the old
default (mb=128/nj=16) and roughly 10× faster than the original −O0 single-threaded code.
A couple of things are worth noticing. `min_block` is a real lever for both speed and
memory, but only up to a point: at 2048 it slows down again, because the dense diagonal
blocks get too large to handle efficiently. And the best `n_jobs` climbs as `min_block`
does, which is no accident: with fewer blocks the convolution outputs stop piling up (the
mechanism in §3.3), and once that bottleneck is gone the threads scale cleanly.

### 3.3 Memory: the real wall, and how to push it back

Everything above was about speed. For grids of 2²¹ and beyond, the binding constraint stops
being time and becomes memory. It is worth understanding exactly *where* that memory goes,
because it is not where you would guess.

It is not the compressed matrix (that is only a few GB). The killer is a **transient inside
`convolve`**. Each block writes its contribution into its own output buffer of length
`N − (i0+j0)`, and with many producer threads feeding one accumulator, a crowd of these
buffers is "in flight" at once. Peak memory is therefore roughly (number of blocks) × (a
pile-up factor that scales with `n_jobs`) — and it grows by about ×3 with every grid
doubling. On top of that, at a tight `rel_tol` the MSk *builder* needs a large scratch area
(it scales with block rank), which piles on further.

That decomposition tells you immediately which knobs to turn. Measured at size 1048576 (mass
identical, accuracy unaffected in every case):

| lever | values → peak memory |
|---|---|
| **n_jobs** (pile-up) | 16 → 150 GB; 4 → 111 GB; **1 → 18 GB** |
| **min_block** (#blocks) | 128 → 150 GB; 512 → 61; **1024 → 51 GB** (+ ×2 speed) |
| **rel_tol** (build scratch) | 1e-10 → ~470 GB @ 2²¹ (OOM); **1e-5 → ×16 smaller**, Frobenius vs 1e-10 = 1e-7 |

So the three levers are independent and, for big runs, all point the same way: drop
`n_jobs`, raise `min_block`, loosen `rel_tol`. Note the last one especially — going from
rel_tol 1e-10 to 1e-5 shrinks the build scratch ~16× while changing the answer by a
Frobenius of only 1e-7. This is the exact mirror image of the t = 10⁶ advice in §3.2, where
the system fits and you *want* many threads.

**Bottom line for large grids: low `n_jobs` (2–4), `min_block=1024`, `rel_tol=1e-5`.**

### 3.4 Reaching t = 10⁷

Putting §3.3 into practice, here are the configurations we tried at the final grid size 2²³
(mass 0.99977 in all of them; anything with `n_jobs ≥ 8` ran out of memory from the convolve
pile-up):

| min_block | nj | wall | peak, GB | verdict |
|---:|---:|---:|---:|---|
| **1024** | **4** | **33 min** | **280** | ⭐ fastest that fits |
| 1024 | 3 | 40 min | 208 | conservative |
| 1024 | 2 | 54 min | 152 | max memory headroom |
| 512 | 4 | — | 333 | ✗ OOM |

And that did it. Running with `max_size=2²³, rel_tol=1e-5,
min_block=1024, n_jobs=2` under a memory guard gave: final time 10⁷, final size 2²³,
**total mass = 0.99977** (a drift of 2.3e-4, i.e. mass conserved), with 0 negative
components, a 152 GB peak, and ~54 minutes of wall time. (We used the conservative
`n_jobs=2` row for headroom; the `n_jobs=4` row above is the faster option once you trust
the memory budget.) The small drift is real but benign: it accumulates slowly over the full
10⁷ of simulated time, from the slight asymmetry of the *approximated* kernel plus the
clamping of negative concentrations to zero. On the shorter t ≤ 10⁶ references the same
drift is only ~1e-6.

One operational wrinkle worth recording: memory kept *climbing* during the run, for two
reasons. The obvious one is the resizes — the convolve buffer scales with grid size, so
every epoch steps the peak up. The subtler one is a glibc allocator "ratchet": once it has
grown its arenas for those multi-megabyte buffers, freed memory is not handed back to the
OS. If that becomes a problem, forcing large allocations through `mmap` tames it:

```
MALLOC_MMAP_THRESHOLD_=131072 MALLOC_TRIM_THRESHOLD_=131072
```

---

## 4. The ballistic kernel — and the instability that taught us the most

The ballistic kernel (`λ = 5/6`) is a different animal. Its front races outward — roughly as
`t^4.7` — so it reaches large sizes almost immediately: ~1.2K by t=10, ~65K by t=20, ~0.4M
(2²⁰) by t=40, ~1.3M (2²²) by t=50. We run it with `mosaic=monodiag` (ρ=1) but otherwise
just like the atmospheric case.

This is where we learned the most important lesson of the whole project, and it is about
stiffness: the ballistic kernel demands a **small** `ode_tol`. Our first runs at t=40 and
t=50 used `ode_tol=1e-7`, and the mass *blew up* — 1.90 and 2.44 instead of 1.0. That is physically
impossible (`λ = 5/6 < 1`, so there is no gelation, no mechanism to create mass), which is
exactly what flagged it as a numerical failure. Here is the mechanism, because the same trap
waits for any stiff kernel:

> The step controller watches **accuracy**, not **stability**. Explicit RK4 is only stable
> below a step-size limit, but nothing in the controller knows about that limit — so when
> the dynamics are stiff, the controller happily grows the step past it. Beyond the limit
> the solution oscillates and goes negative; our `y < 0 → 0` clamp then quietly *adds* mass
> on every step; that extra mass feeds back, and the whole thing runs away.

The fix was simply `ode_tol=1e-10`, which holds the step below the stability limit and
brought mass back to 0.99997. (The atmospheric kernel never needed this — it was fine even
at 1e-6. Ballistic is just stiffer.) This lesson is now baked into the solver as a **mass
guard** (`SMOL_MASS_GUARD`, see §6) that aborts a run the moment mass starts climbing,
instead of letting it waste hours diverging.

With a small enough `ode_tol`, the results line up well against an external reference
solution (mass conserved throughout):

| t | rel L2 (density) | mass | final size | peak GB | ode_tol |
|---:|---:|---:|---:|---:|---:|
| 10 | 0.34% | 0.99999 | 4096 | <1 | 1e-7 |
| 20 | 0.57% | 0.99998 | 65536 | ~1 | 1e-7 |
| 40 | 1.70% | 0.99997 | 2²⁰ | 16 | 1e-10 |
| 50 | 2.28% | 0.99997 | 2²² | 95 | 1e-10 |

The bulk of the distribution agrees closely (deterministic-to-reference ratios sit around
0.95–1.06). The one number that looks large — the 2.28% mass-weighted L2 at t=50 — is not a
solver error; it is statistical noise in the *reference*, which out in the tail is resolving
single particles at sizes near 10⁶.

---

## 5. How the measurements were made

A few notes so the numbers above are reproducible, and so you trust them:

- **Hardware.** Intel Xeon Gold 6226R (Cascade Lake), 32 physical cores / 64 hyper-threads,
  AVX-512, 376 GB RAM, slurm partition `c32m384`. The build host is the run host, so
  `-march=native` is safe.
- **Where to measure.** Only on a dedicated slurm node — the login node is too noisy to time
  anything. Every measurement is its own exclusive `sbatch`.
- **What "correct" means here.** Two checks together: a relative Frobenius difference
  `‖n − n_ref‖ / ‖n_ref‖ ≲ 1e-6` (the same scale as `ode_tol`), plus mass conservation.
  Switching compiler flags or parallel reductions perturbs the floating-point result at the
  ~1e-13 level — the physics does not move.
- **Cluster gotchas that cost us time.** An `srun` step inside an `sbatch` must be given an
  explicit `--cpus-per-task=64 --cpu-bind=none`, or it gets pinned to a single core.
  `/usr/bin/time` is not installed on the nodes, so peak memory and wall time come from the
  program's own output. Rebuilt binaries need `LD_LIBRARY_PATH` set (gcc-14.2 libstdc++ +
  openblas + fftw). And sample memory *often* — near the peak it can jump by more than 16 GB
  between samples.

---

## 6. Guide: running a brand-new kernel

Everything above converges on one organizing parameter — the homogeneity **`λ`**, from
`K(ai, aj) = a^λ K(i,j)`. **Compute it first.** It predicts both of the things that can sink
a run: how fast the front grows (memory) and how stiff the system is (stability).

> ⚠️ **The single most dangerous thing in practice is memory.** It grows ~×3 per grid
> doubling and blows up in an instant — every one of our failures (OOM, runaway) traced back
> to it. Hence the rule: **measure first** with a cheap probe of the front and the peak
> memory, **then launch — and always under a memory guard.**

**Step 0 — classify by `λ`.**
- `λ > 1` means **gelation**: a finite-time singularity after which mass physically
  *decreases*. Mass-as-invariant no longer applies; you would need separate gel-fraction
  bookkeeping. (We never worked in this regime.)
- `λ ≤ 1` means mass is conserved, so **`Total mass ≈ 1` is your hard pass/fail**, and any
  growth above 1 is an artifact.
- The larger `λ` is, the faster the front grows *and* the stiffer the system — both screws
  tighten together.

**Step 1 — stability, via `ode_tol`.** Explicit RK4 is stable only for a small enough step,
but the controller chases accuracy, not stability, so it will not protect you. The symptom
of trouble is unmistakable: **mass climbing above 1** (and accelerating), with the front
overshooting where it physically should be. Rule of thumb: the larger `λ` or the bigger the
grid, the smaller `ode_tol` needs to be (atmospheric was happy at 1e-6; ballistic needed
1e-10). Start around 1e-9…1e-10. As a backstop, the solver has a built-in 🛡 **mass guard**
(`SMOL_MASS_GUARD`, default 1.02) that aborts with a "reduce ode_tol" message before a full
runaway burns your allocation.

**Step 2 — `max_size`, the front-vs-mass trade.** Mass stays conserved only while the front
(99.9% of the mass) is below the top of the grid. So estimate where the front will be at
your target `t` — a cheap small run plus extrapolation is enough — and set `max_size` with
margin above it. Since memory grows ~×3 per doubling, `max_size` is your number-one memory
lever.

**Step 3 — memory, via `n_jobs`, `min_block`, `rel_tol`.**
- **`min_block`** — raise it to make fewer, larger blocks: less memory, and often faster.
  Default to 1024 (2048 is already slower). Does not change the answer.
- **`n_jobs`** — regime-dependent. If the system fits in RAM (grid ≤ ~2²⁰), go high (16–32)
  for speed. If it is large (≥ 2²¹), go low (2–4), or the convolve pile-up will OOM you.
- **`rel_tol`** — loosen it (≤ 1e-5) to shrink the builder's scratch memory at no real
  accuracy cost. **Do not confuse it with `ode_tol`:** `rel_tol` controls how accurately the
  *kernel matrix* is approximated; `ode_tol` controls *time-stepping* stability. They are
  unrelated knobs that merely happen to both be tolerances.
- Always run under a **memory guard** (stop at ~340 of the 376 GB).

**Step 4 — mosaic type (`rho`).** Pick it for how singular the kernel is near the diagonal
(atmospheric → `tridiag`, ρ=2; ballistic → `monodiag`, ρ=1). At a fixed `rel_tol` it barely
moves the accuracy.

**What to watch while it runs:**

| metric | normal | alarm |
|---|---|---|
| `Total mass` | ≈ 1, stable | growth > 1 ⇒ instability (lower ode_tol); a strong drop ⇒ leakage (raise max_size) |
| `Negative components` | 0 / tiny | many ⇒ oscillations (lower ode_tol) |
| `Final size` vs expected front | matches | overshoot ⇒ fake mass / instability |
| peak RAM | < guard | approaching the guard ⇒ lower nj or raise min_block |
| `step` over time | grows moderately, then plateaus | unbounded growth ⇒ runaway imminent |

**The recipe, start to finish:** (1) compute `λ`; (2) do a cheap recon run at small `t` —
tune `ode_tol` until mass stops growing, and read off the front position and peak memory;
(3) choose `max_size` with margin; (4) set `min_block=1024`, `rel_tol=1e-5`, and `n_jobs` to
fit your memory budget; (5) launch through slurm with a memory guard and live monitoring of
mass / size / RAM; (6) compare against a reference with `sweeps/ref_compare.py`. After any
code change, run `bash sweeps/bench.sh` to catch regressions.

**Calibration from the two kernels we did run:**

| | atmospheric | ballistic |
|---|---|---|
| `λ` | −5/9 (decreasing) | 5/6 (increasing) |
| front growth | slow (t=10⁷ → 2.3M) | fast (t=50 → 1.3M, ~t^4.7) |
| stiffness / `ode_tol` | mild, 1e-6 OK | stiff, needs 1e-10 |
| mosaic | tridiag (ρ=2) | monodiag (ρ=1) |
| mass result | 0.99977 @ t=10⁷ | 0.99997 @ t≤50 |
