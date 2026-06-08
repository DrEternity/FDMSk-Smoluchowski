#ifndef SMOLUCHOWSKI_CONFIG_HPP
#define SMOLUCHOWSKI_CONFIG_HPP
//
// Single place for every hyperparameter of the solver.
//
// Precedence: built-in default  ->  legacy positional CLI  ->  environment var.
// (env wins, so parameter sweeps just export variables.)
//
// Environment variables
//   Problem:    SMOL_MAX_SIZE  SMOL_TIME  SMOL_INIT_SIZE  SMOL_KERNEL(atmos|ballistic)
//               SMOL_MOSAIC(tridiag|monodiag)  SMOL_OUTPUT
//   Integrator: SMOL_DT  SMOL_ODE_TOL
//   Mosaic-MSk: MSK_REL_TOL  MSK_MIN_BLOCK  MSK_MAX_RANK  MSK_NJOBS
//
// Legacy positional CLI (kept for backward compatibility, same order as before):
//   ./example  max_size  time  output  initial_size  [unused]  ode_tol
//   (the 5th slot was `mass_threshold`, now unused; kept so arg order is unchanged.)
//
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>
#include <iomanip>

struct SolverConfig {
    // --- Problem definition ---
    int         max_size     = 65536;   // grid cap (power of two)
    double      time         = 1000.0;  // final time T
    int         initial_size = 512;     // starting grid size
    int         kernel       = 0;       // 0 = atmospheric, 1 = ballistic
    int         mosaic       = 1;       // 0 = monodiag (rho=1), 1 = tridiag (rho=2)
    std::string output_file  = "reference_solution_atmos.txt";
    // --- Time integrator (adaptive RK4 with mass/front monitoring) ---
    double      dt           = 1e-4;    // initial step size
    double      ode_tol      = 1e-6;    // step-error tolerance
    // --- Mosaic-skeleton approximation ---
    double      rel_tol      = 1e-10;   // block approximation tolerance
    uint64_t    min_block    = 128;     // min mosaic block (min_rows = min_cols)
    uint64_t    max_rank     = 0;       // block rank cap (0 = unlimited)
    uint64_t    n_jobs       = 1;       // MSk worker threads (approx + matvec + convolve)
};

namespace cfg_detail {
    inline bool get_d(const char* k, double& v) {
        if (const char* e = std::getenv(k)) { v = std::atof(e); return true; } return false; }
    inline bool get_i(const char* k, long& v) {
        if (const char* e = std::getenv(k)) { v = std::atol(e); return true; } return false; }
}

// Legacy positional arguments (order preserved; slot 5 is intentionally ignored).
inline void config_from_cli(SolverConfig& c, int argc, char* argv[]) {
    if (argc > 1) c.max_size     = std::stoi(argv[1]);
    if (argc > 2) c.time         = std::stod(argv[2]);
    if (argc > 3) c.output_file  = argv[3];
    if (argc > 4) c.initial_size = std::stoi(argv[4]);
    // argv[5] was mass_threshold — unused now (resize is driven by the front), kept for arg-order compat
    if (argc > 6) c.ode_tol      = std::stod(argv[6]);
}

// Environment overrides (highest precedence). Same validation as before for the
// three established knobs; new knobs added.
inline void config_from_env(SolverConfig& c) {
    using namespace cfg_detail;
    double d; long i;
    // problem
    if (get_i("SMOL_MAX_SIZE",  i) && i > 0) c.max_size     = (int)i;
    if (get_d("SMOL_TIME",      d) && d > 0) c.time         = d;
    if (get_i("SMOL_INIT_SIZE", i) && i > 0) c.initial_size = (int)i;
    if (const char* e = std::getenv("SMOL_KERNEL"))
        c.kernel = (std::strcmp(e, "ballistic") == 0 || std::strcmp(e, "1") == 0) ? 1 : 0;
    if (const char* e = std::getenv("SMOL_MOSAIC"))
        c.mosaic = (std::strcmp(e, "monodiag") == 0 || std::strcmp(e, "0") == 0) ? 0 : 1;
    if (const char* e = std::getenv("SMOL_OUTPUT")) c.output_file = e;
    // integrator
    if (get_d("SMOL_DT",      d) && d > 0) c.dt      = d;
    if (get_d("SMOL_ODE_TOL", d) && d > 0) c.ode_tol = d;
    // mosaic-skeleton
    if (get_d("MSK_REL_TOL",   d) && d > 0)  c.rel_tol   = d;
    if (get_i("MSK_MIN_BLOCK", i) && i >= 1) c.min_block = (uint64_t)i;
    if (get_i("MSK_MAX_RANK",  i) && i >= 0) c.max_rank  = (uint64_t)i;
    if (get_i("MSK_NJOBS",     i) && i >= 1) c.n_jobs    = (uint64_t)i;
}

inline SolverConfig load_config(int argc, char* argv[]) {
    SolverConfig c;
    config_from_cli(c, argc, argv);   // legacy positional
    config_from_env(c);               // env overrides win
    return c;
}

inline void print_config(const SolverConfig& c) {
    std::cout << "=== Solver configuration ===\n";
    std::cout << "  [problem]    kernel=" << (c.kernel ? "ballistic" : "atmospheric")
              << "  mosaic=" << (c.mosaic ? "tridiag(rho=2)" : "monodiag(rho=1)")
              << "  max_size=" << c.max_size
              << "  initial_size=" << c.initial_size
              << "  time=" << c.time << "\n";
    std::cout << "  [integrator] dt=" << c.dt << "  ode_tol=" << c.ode_tol << "\n";
    std::cout << "  [mosaic-MSk] rel_tol=" << c.rel_tol
              << "  min_block=" << c.min_block
              << "  max_rank=" << c.max_rank << (c.max_rank ? "" : " (unlimited)")
              << "  n_jobs=" << c.n_jobs << "\n";
    std::cout << "  [output]     " << c.output_file << "\n";
    std::cout << "============================" << std::endl;
}

#endif // SMOLUCHOWSKI_CONFIG_HPP
