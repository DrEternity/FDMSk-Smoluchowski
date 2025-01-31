#ifndef MOSAIC_SKELETON_MATRIX_PARTICLES_HPP
#define MOSAIC_SKELETON_MATRIX_PARTICLES_HPP


#include <cstdint>
#include <functional>
#include <chrono>

#include "msk/matrix/basic.hpp"
#include "msk/oracle/elementary.hpp"
#include "msk/mosaic_skeleton.hpp"
#include "msk/solvers/ode/runge_kutta.hpp"
#include "msk/utility/testing.hpp"



template<typename DataType, class Matrix>
class Smoluch {
    const uint64_t _size;
    const Matrix &_matrix;
    const DataType _src{0};
    DataType *_tmp_y{nullptr};
public:
    Smoluch(const uint64_t size, const Matrix &matrix, DataType src = 0) : _size{size}, _matrix{matrix}, _src{src}, _tmp_y{new DataType[_size]} {}
    ~Smoluch() {
        delete [] _tmp_y;
    }

    void apply(const DataType *x, DataType *y) const {
        _tmp_y[0] = 0;
        _matrix.convolve(x, x, _tmp_y + 1, _size - 1);
        _matrix.matvec(x, y);
        for (uint64_t i = 0; i < _size; ++i) {
            y[i] = 0.5 * _tmp_y[i] - x[i] * y[i];
        }
        y[0] += _src;
    }
};

enum class MosaicType {
    monodiag,
    tridiag
};

double* modeling(unsigned int size, std::function<double(double, double)> kernel, double rel_tol, double *n_0, double time, double first_step, MosaicType mosaic_type = MosaicType::monodiag) {
    auto new_kernel = [kernel](uint64_t i, uint64_t j) { // because Bulat MSk approx have zero indexing
        return kernel(i + 1, j + 1);
    };

    double abs_tol = 0.0;
    uint64_t max_rank = 0;

    MSk::oracle::Parameters params;
    params.min_rows = 32;
    params.min_cols = 32;
    switch (mosaic_type) {
    case MosaicType::monodiag:
        params.rho = 1.0;
        break;
    case MosaicType::tridiag:
        params.rho = 2.0;
        break;
    default:
        params.rho = 1.0;
    }
    params.block_comp = {0.0, rel_tol, max_rank};

    MSk::oracle::Elementary<double> oracle{params};

    MSk::matrix::Elementary<double> matrix(size, size, new_kernel);

    auto t_start = std::chrono::steady_clock::now();
    auto msk = MSk::MosaicSkeleton<double>::approximate(matrix, oracle, 1);
    auto t_end = std::chrono::steady_clock::now();
    double t_diff = std::chrono::duration<double>(t_end - t_start).count();
    std::cout << "Kernel approximated, time = " << t_diff  << " sec" << std::endl;
    const bool check_compression{true}, check_error{true};
    if (check_compression) {
        double compression = MSk::check_compression<double>(msk, nullptr);
        std::cout << "Compression rate = " << 1.0 / compression << " (" << compression * 100.0 << "\% of full matrix)" << std::endl;
    }
    if (check_error) {
        double error = MSk::check_error<double, MSk::matrix::Elementary<double>>(msk, matrix, nullptr);
        std::cout << "Approximation error = " << error << std::endl;
    }

    Smoluch<double, MSk::MosaicSkeleton<double>> smoluch(size, msk);

    t_start = std::chrono::steady_clock::now();
    MSk::solvers::ode::runge_kutta4_adaptive(size, smoluch, n_0, time, first_step);
    t_end = std::chrono::steady_clock::now();
    t_diff = std::chrono::duration<double>(t_end - t_start).count();
    std::cout << "ODE solved, time = " << t_diff  << " sec" << std::endl;

    return n_0;
}


#endif