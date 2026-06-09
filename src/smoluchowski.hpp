#ifndef MOSAIC_SKELETON_MATRIX_PARTICLES_HPP
#define MOSAIC_SKELETON_MATRIX_PARTICLES_HPP

#include <cstdint>
#include <cstdlib>
#include <functional>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "msk/matrix/basic.hpp"
#include "msk/oracle/elementary.hpp"
#include "msk/mosaic_skeleton.hpp"
#include "msk/utility/testing.hpp"
#include "smart_array/smart_array.hpp"
#include "cpp_blas/BLAS.hpp"

// ---------------------------------------------------------------------------
// Smoluchowski RHS operator
// ---------------------------------------------------------------------------
template<typename DataType, class Matrix>
class Smoluch {
    const uint64_t _size;
    const Matrix &_matrix;
    const DataType _src{0};
    DataType *_tmp_y{nullptr};
    DataType *_tmp_x{nullptr};
    mutable uint64_t _rhs_count{0};
    mutable double _rhs_total_time{0.0};
public:
    Smoluch(const uint64_t size, const Matrix &matrix, DataType src = 0)
        : _size{size}, _matrix{matrix}, _src{src},
          _tmp_y{new DataType[_size]}, _tmp_x{new DataType[_size]} {}
    ~Smoluch() {
        delete[] _tmp_y;
        delete[] _tmp_x;
    }
    void apply(const DataType *x, DataType *y) const {
        auto t_start = std::chrono::steady_clock::now();

        for (uint64_t i = 0; i < _size; ++i) {
            _tmp_x[i] = (x[i] > 0.0) ? x[i] : 0.0;
        }
        _tmp_y[0] = 0;
        _matrix.convolve(_tmp_x, _tmp_x, _tmp_y + 1, _size - 1);
        _matrix.matvec(_tmp_x, y);
        for (uint64_t i = 0; i < _size; ++i) {
            y[i] = 0.5 * _tmp_y[i] - _tmp_x[i] * y[i];
        }
        y[0] += _src;

        auto t_end = std::chrono::steady_clock::now();
        _rhs_total_time += std::chrono::duration<double>(t_end - t_start).count();
        _rhs_count++;
    }
    uint64_t rhs_count() const { return _rhs_count; }
    double rhs_total_time() const { return _rhs_total_time; }
    double rhs_avg_time() const { return _rhs_count > 0 ? _rhs_total_time / _rhs_count : 0.0; }
};

// ---------------------------------------------------------------------------
// Utilities
// ---------------------------------------------------------------------------
inline double compute_mass(const double* n, uint64_t size) {
    double mass = 0.0;
    for (uint64_t i = 0; i < size; ++i) {
        mass += n[i] * static_cast<double>(i + 1);
    }
    return mass;
}

// ---------------------------------------------------------------------------
// Modified adaptive RK4 (based on Bulat's runge_kutta4_adaptive template)
//
// Additions vs original:
//   - Starts from t_start (not 0)
//   - After each accepted step, checks mass; if relative drop > mass_threshold,
//     returns with need_resize = true
//   - Returns current step size (for reuse after resize)
//   - Logs mass instead of norm
// ---------------------------------------------------------------------------
enum class StepperResult {
    finished,       // reached T
    mass_dropped    // mass dropped below threshold => need resize
};

template<typename RealType, class Operator>
StepperResult runge_kutta4_adaptive_mass(
        uint64_t size, Operator &A, RealType *y,
        RealType T_total,          // total end time (absolute)
        RealType &t_current,       // current time (in/out)
        RealType &step,
	RealType tol)            // current step size (in/out)
{
    SmartArray<RealType> y_step(size), A_y_step(size), dy(size), dy_check(size);

    RealType T = T_total - t_current;  // remaining time
    RealType t_local{0};
    RealType next_t_verbose{T / 100};
    uint64_t steps_total{0}, steps_hit{0};

    while (t_local < T) {
        if (step > T - t_local) {
            step = T - t_local;
        }
        if (t_local >= next_t_verbose) {
            RealType mass = compute_mass(y, size);
            std::cout << "t = " << (t_current + t_local) << " :: step = " << step
                      << ", mass = " << mass
                      << "; steps: total = " << steps_total
                      << ", successful = " << steps_hit << std::endl;
            if (std::isnan(mass)) {
                break;
            }
            next_t_verbose += T / 100;
        }

        // step 1
        A.apply(y, A_y_step.data());
        BLAS::copy(size, A_y_step.data(), 1, dy.data(), 1);
        BLAS::copy(size, A_y_step.data(), 1, dy_check.data(), 1);

        // step 2
        BLAS::copy(size, y, 1, y_step.data(), 1);
        BLAS::axpy(size, step / 2, A_y_step.data(), 1, y_step.data(), 1);
        A.apply(y_step.data(), A_y_step.data());
        BLAS::axpy(size, RealType{2}, A_y_step.data(), 1, dy.data(), 1);

        // step 3
        BLAS::copy(size, y, 1, y_step.data(), 1);
        BLAS::axpy(size, step / 2, A_y_step.data(), 1, y_step.data(), 1);
        A.apply(y_step.data(), A_y_step.data());
        BLAS::axpy(size, RealType{2}, A_y_step.data(), 1, dy.data(), 1);

        // step 4
        BLAS::copy(size, y, 1, y_step.data(), 1);
        BLAS::axpy(size, step, A_y_step.data(), 1, y_step.data(), 1);
        A.apply(y_step.data(), A_y_step.data());
        BLAS::axpy(size, RealType{1}, A_y_step.data(), 1, dy.data(), 1);

        // check step tolerance
        BLAS::axpy(size, RealType{-1}, A_y_step.data(), 1, dy_check.data(), 1);
        RealType dy_err = step / 6 * BLAS::nrm2(size, dy_check.data(), 1);

        if (dy_err < tol) {
            t_local += step;
            BLAS::axpy(size, step / 6, dy.data(), 1, y, 1);
	    for (uint64_t i = 0; i < size; i++) {
    		if (y[i] < 0.0) y[i] = 0.0;
	    }
            ++steps_hit;

	    // --- Front check after accepted step ---
            double cumulative = 0.0;
            double total = compute_mass(y, size);
            uint64_t i_front = 0;
            for (uint64_t i = 0; i < size; i++) {
                cumulative += y[i] * static_cast<double>(i + 1);
                if (cumulative >= (1.0 - 1e-3) * total) {
                    i_front = i;
                    break;
                }
            }
            if (i_front > size / 2) {
                t_current += t_local;
                std::cout << "  *** Front reached " << i_front
                          << " / " << size
                          << " at t = " << t_current
                          << ", triggering resize ***" << std::endl;
                std::cout << "Total steps: " << steps_total
                          << ", successful steps: " << steps_hit << std::endl;
                return StepperResult::mass_dropped;
            }
	}

        // update step size
        step *= std::min(RealType{4}, std::max(RealType{0.25},
                std::pow(tol / (2 * dy_err), 0.25)));
        if (t_local + step > T) {
            step = T - t_local;
        }
        ++steps_total;
    }

    t_current += t_local;
    std::cout << "Total steps: " << steps_total
              << ", successful steps: " << steps_hit << std::endl;
    return StepperResult::finished;
}

// ---------------------------------------------------------------------------
// MSk approximation builder
// ---------------------------------------------------------------------------
enum class MosaicType {
    monodiag,
    tridiag
};

// ---------------------------------------------------------------------------
// Main modeling function with adaptive matrix sizing
// ---------------------------------------------------------------------------
// Hyperparameters are passed explicitly (see config.hpp / SolverConfig). The previous
// scattered getenv() reads were moved to config.hpp; modeling() is now pure plumbing.
double* modeling(unsigned int max_size,
                 std::function<double(double, double)> kernel,
                 double rel_tol, double *n_0, double time, double first_step,
                 MosaicType mosaic_type = MosaicType::tridiag,
                 unsigned int initial_size = 512,
                 double ode_tol = 1e-6,
                 uint64_t n_jobs = 1,
                 uint64_t min_block = 128,
                 uint64_t max_rank = 0) {

    auto new_kernel = [kernel](uint64_t i, uint64_t j) {
        return kernel(i + 1, j + 1);
    };

    uint64_t current_size = initial_size;
    double t_current = 0.0;
    double step = first_step;
    double initial_mass = compute_mass(n_0, max_size);

    auto t_global_start = std::chrono::steady_clock::now();

    std::cout << "\n=== Adaptive Integration ===" << std::endl;
    std::cout << "Initial size: " << current_size << ", max size: " << max_size << std::endl;
    std::cout << "Initial mass: " << initial_mass << std::endl;
    std::cout << "MSk n_jobs: " << n_jobs << ", min_block: " << min_block
              << ", max_rank: " << max_rank << std::endl;
    std::cout << "ODE tol: " << ode_tol << ", MSk rel_tol: " << rel_tol << std::endl;

    while (t_current < time && current_size <= static_cast<uint64_t>(max_size)) {
        // --- Build MSk approximation for current_size ---
        std::cout << "\n--- Building MSk for size " << current_size << " ---" << std::endl;

        MSk::oracle::Parameters params;
        params.min_rows = min_block;
        params.min_cols = min_block;
        switch (mosaic_type) {
        case MosaicType::monodiag: params.rho = 1.0; break;
        case MosaicType::tridiag:  params.rho = 2.0; break;
        default: params.rho = 2.0;
        }
        params.block_comp = {0.0, rel_tol, max_rank};

        MSk::oracle::Elementary<double> oracle{params};
        MSk::matrix::Elementary<double> matrix(current_size, current_size, new_kernel);

        auto t_approx_start = std::chrono::steady_clock::now();
        auto msk = MSk::MosaicSkeleton<double>::approximate(matrix, oracle, n_jobs);
        auto t_approx_end = std::chrono::steady_clock::now();
        double t_approx = std::chrono::duration<double>(t_approx_end - t_approx_start).count();

        double compression = MSk::check_compression<double>(msk, nullptr);
        std::cout << "  Approximation time: " << t_approx << " sec" << std::endl;
        std::cout << "  Compression: " << 1.0 / compression
                  << "x (" << compression * 100.0 << "% of full)" << std::endl;

        Smoluch<double, MSk::MosaicSkeleton<double>> smoluch(current_size, msk);

        std::cout << "--- Integrating from t = " << t_current
                  << ", size = " << current_size << " ---" << std::endl;

        // --- Run adaptive RK4 with mass monitoring ---
        // Note: stepper operates on first current_size elements of n_0
        StepperResult result = runge_kutta4_adaptive_mass<double>(
            current_size, smoluch, n_0,
            time, t_current, step,
            ode_tol);

        std::cout << "  RHS evaluations: " << smoluch.rhs_count()
                  << ", avg time: " << smoluch.rhs_avg_time() << " sec" << std::endl;
        std::cout << "  Current mass: " << compute_mass(n_0, current_size)
                  << ", t = " << t_current << std::endl;

        if (result == StepperResult::finished) {
            break;
        }

        // --- Resize: double the matrix ---
        if (current_size * 2 <= static_cast<uint64_t>(max_size)) {
            uint64_t new_size = current_size * 2;
            std::cout << "\n*** Resizing: " << current_size << " -> " << new_size << " ***" << std::endl;
            // New elements are already zero (n_0 was zero-initialized for max_size)
            current_size = new_size;
        } else {
            std::cout << "\n*** Cannot resize further (max_size = " << max_size
                      << "), continuing with current size ***" << std::endl;
            // Continue with current size until time ends
            StepperResult final_result = runge_kutta4_adaptive_mass<double>(
                current_size, smoluch, n_0,
                time, t_current, step,
                ode_tol);  // relax threshold
            break;
        }
    }

    auto t_global_end = std::chrono::steady_clock::now();
    double total_time = std::chrono::duration<double>(t_global_end - t_global_start).count();
    std::cout << "\n=== Integration complete ===" << std::endl;
    std::cout << "  Final time: " << t_current << std::endl;
    std::cout << "  Final size: " << current_size << std::endl;
    std::cout << "  Final mass: " << compute_mass(n_0, current_size) << std::endl;
    std::cout << "  Total wall time: " << total_time << " sec" << std::endl;

    return n_0;
}

#endif
