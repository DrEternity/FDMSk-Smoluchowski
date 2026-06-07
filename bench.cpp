/*
 * Benchmark: MSk approximation + matvec/convolve speed
 * for atmospheric kernel at different min_block sizes
 *
 * Usage: ./bench [size]
 *   size - matrix size (default 65536)
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <vector>
#include <cstdint>
#include <functional>

#include "msk/matrix/basic.hpp"
#include "msk/oracle/elementary.hpp"
#include "msk/mosaic_skeleton.hpp"
#include "msk/utility/testing.hpp"

// Atmospheric kernel (1-indexed)
double kernel_atmos(uint64_t i, uint64_t j) {
    if (i == j) return 4.0;
    double num = (i + j) * std::pow(std::pow(i, 1.0/3.0) + std::pow(j, 1.0/3.0), 2.0/3.0);
    double denom = std::pow(i * j, 5.0/9.0) * std::abs(std::pow(i, 2.0/3.0) - std::pow(j, 2.0/3.0));
    return num / denom;
}

int main(int argc, char* argv[]) {
    uint64_t size = 65536;
    if (argc > 1) size = std::stoull(argv[1]);

    int n_runs = 10;
    if (argc > 2) n_runs = std::stoi(argv[2]);

    double rel_tol = 1e-10;

    auto new_kernel = [](uint64_t i, uint64_t j) {
        return kernel_atmos(i + 1, j + 1);  // 0-indexed -> 1-indexed
    };

    std::vector<uint64_t> block_sizes = {32, 64, 128, 256, 512};

    // Prepare test vectors
    double* x = new double[size];
    double* y = new double[size];
    double* conv_out = new double[size];
    for (uint64_t i = 0; i < size; i++) {
        x[i] = (i == 0) ? 1.0 : std::exp(-0.01 * i);  // some non-trivial input
    }

    std::cout << "=== MSk Benchmark ===" << std::endl;
    std::cout << "Matrix size: " << size << "x" << size << std::endl;
    std::cout << "Kernel: Atmospheric" << std::endl;
    std::cout << "MSk rel_tol: " << rel_tol << std::endl;
    std::cout << "Averaging over " << n_runs << " runs" << std::endl;
    std::cout << "=====================\n" << std::endl;

    std::cout << std::left
              << std::setw(12) << "min_block"
              << std::setw(16) << "approx_time(s)"
              << std::setw(16) << "compression(%)"
              << std::setw(16) << "matvec_avg(s)"
              << std::setw(16) << "matvec_std(s)"
              << std::setw(16) << "conv_avg(s)"
              << std::setw(16) << "conv_std(s)"
              << std::endl;
    std::cout << std::string(108, '-') << std::endl;

    for (uint64_t block_sz : block_sizes) {
        // Build approximation
        MSk::oracle::Parameters params;
        params.min_rows = block_sz;
        params.min_cols = block_sz;
        params.rho = 2.0;
        uint64_t max_rank = 0;
        params.block_comp = {0.0, rel_tol, max_rank};

        MSk::oracle::Elementary<double> oracle{params};
        MSk::matrix::Elementary<double> matrix(size, size, new_kernel);

        auto t0 = std::chrono::steady_clock::now();
        auto msk = MSk::MosaicSkeleton<double>::approximate(matrix, oracle, 1);
        auto t1 = std::chrono::steady_clock::now();
        double approx_time = std::chrono::duration<double>(t1 - t0).count();

        double compression = MSk::check_compression<double>(msk, nullptr);

        // Warmup
        msk.matvec(x, y);
        msk.convolve(x, x, conv_out + 1, size - 1);

        // Benchmark matvec
        std::vector<double> matvec_times(n_runs);
        for (int r = 0; r < n_runs; r++) {
            auto ts = std::chrono::steady_clock::now();
            msk.matvec(x, y);
            auto te = std::chrono::steady_clock::now();
            matvec_times[r] = std::chrono::duration<double>(te - ts).count();
        }

        // Benchmark convolve
        std::vector<double> conv_times(n_runs);
        for (int r = 0; r < n_runs; r++) {
            conv_out[0] = 0;
            auto ts = std::chrono::steady_clock::now();
            msk.convolve(x, x, conv_out + 1, size - 1);
            auto te = std::chrono::steady_clock::now();
            conv_times[r] = std::chrono::duration<double>(te - ts).count();
        }

        // Compute mean and std
        auto mean_std = [](const std::vector<double>& v) -> std::pair<double, double> {
            double sum = 0, sum2 = 0;
            for (double x : v) { sum += x; sum2 += x * x; }
            double mean = sum / v.size();
            double var = sum2 / v.size() - mean * mean;
            return {mean, std::sqrt(std::max(0.0, var))};
        };

        auto [mv_mean, mv_std] = mean_std(matvec_times);
        auto [cv_mean, cv_std] = mean_std(conv_times);

        std::cout << std::left
                  << std::setw(12) << block_sz
                  << std::setw(16) << std::fixed << std::setprecision(4) << approx_time
                  << std::setw(16) << std::fixed << std::setprecision(4) << compression * 100.0
                  << std::setw(16) << std::scientific << std::setprecision(4) << mv_mean
                  << std::setw(16) << std::scientific << std::setprecision(4) << mv_std
                  << std::setw(16) << std::scientific << std::setprecision(4) << cv_mean
                  << std::setw(16) << std::scientific << std::setprecision(4) << cv_std
                  << std::endl;
    }

    delete[] x;
    delete[] y;
    delete[] conv_out;

    std::cout << "\nDone!" << std::endl;
    return 0;
}
