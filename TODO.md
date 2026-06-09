# TODO

Ideas not yet implemented (see ANALYSIS.md for context).

## Benchmark: `MODE=full` for large-scale memory/stability regression
`sweeps/bench.sh` currently covers **correctness** on small login-runnable cases. It does
NOT catch the failures that actually bit us — OOM and mass blow-up at large size. Add a
`MODE=full bash sweeps/bench.sh` that submits a couple of SLURM cases (e.g. atmospheric
at size 2^21 and ballistic at t=40) and checks **peak memory** (against the node limit)
and **mass conservation** there, with PASS/FAIL. Reuse `run_t1e7_cfg.sbatch` /
`run_ballistic.sbatch` patterns (memory guard + `/proc` sampling) and the SolverConfig
env interface.

## Visualization: `sweeps/plot_compare.py` (solution vs external reference)
Small matplotlib script (the `smol_research` conda env has matplotlib) plotting the
size distribution `n(size)` on log-log for the deterministic solution vs an external
reference, for both kernels and several times. A visual check is far more intuitive than
the rel-L2 number and is useful for the paper/slides. Inputs: the 2-column
`sweeps/ballistic/conc_det_ballistic_t*.dat` (via `to_sizeconc.py`) and the external
reference `.dat` files.
