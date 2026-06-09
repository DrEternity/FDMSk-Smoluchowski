/*
 * Generate Reference Solution for Atmospheric Kernel
 */

#include "smoluchowski.hpp"
#include "config.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

double kernel_atmos(uint64_t i, uint64_t j) {
    uint64_t ip1 = i;
    uint64_t jp1 = j;

    if (ip1 == jp1) {
        return 4.0;
    }

    double num = (ip1 + jp1) * std::pow(std::pow(ip1, 1.0 / 3.0) + std::pow(jp1, 1.0 / 3.0), 2.0 / 3.0);
    double denom = std::pow(ip1 * jp1, 5.0 / 9.0) * std::abs(std::pow(ip1, 2.0 / 3.0) - std::pow(jp1, 2.0 / 3.0));
    return num / denom;
}

double kernel_ballistic(uint64_t i, uint64_t j) {
    if (i == j) return 4.0;
    double si = std::pow(i, 1.0 / 3.0);
    double sj = std::pow(j, 1.0 / 3.0);
    double sum = si + sj;
    double diff = std::abs(std::pow(i, 1.0 / 6.0) - std::pow(j, 1.0 / 6.0));
    return sum * sum * diff;
}

void save_solution(const std::string& filename, double* n, int size, double time) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "ERROR: Cannot open file " << filename << std::endl;
        return;
    }
    file << "# Reference solution for atmospheric kernel\n";
    file << "# Size: " << size << "\n";
    file << "# Time: " << time << "\n";
    file << "# Format: index concentration\n";
    file << size << " " << std::scientific << std::setprecision(16) << time << "\n";
    for (int i = 0; i < size; i++) {
        file << std::scientific << std::setprecision(16) << n[i] << "\n";
    }
    file.close();
    std::cout << "Solution saved to: " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    // All hyperparameters live in one place now: SolverConfig (config.hpp).
    // Defaults -> legacy positional CLI -> environment variables (env wins).
    SolverConfig cfg = load_config(argc, argv);
    print_config(cfg);

    // Convenience locals reused by the summary/save code below.
    const int max_size = cfg.max_size;
    const double time = cfg.time;
    const std::string output_file = cfg.output_file;

    // Initial condition: monodisperse (allocate max_size, zero-initialized)
    double* n_0 = new double[max_size]();
    n_0[0] = 1.0;

    std::cout << "Initial condition: n[0] = 1, others = 0" << std::endl;
    std::cout << "Starting ODE integration..." << std::endl;

    std::function<double(double, double)> kernel =
        cfg.kernel ? std::function<double(double, double)>(kernel_ballistic)
                   : std::function<double(double, double)>(kernel_atmos);
    MosaicType mosaic_type = cfg.mosaic ? MosaicType::tridiag : MosaicType::monodiag;

    double* n_solution = modeling(max_size, kernel, cfg.rel_tol, n_0, time, cfg.dt,
                                  mosaic_type, cfg.initial_size, cfg.ode_tol,
                                  cfg.n_jobs, cfg.min_block, cfg.max_rank, cfg.mass_guard);


    // Print summary
    double total_mass = 0.0;
    double total_concentration = 0.0;
    int nonzero_count = 0;
    int negative_count = 0;
    double max_concentration = 0.0;
    int max_index = 0;

    for (int i = 0; i < max_size; i++) {
        total_mass += n_solution[i] * (i + 1);
        total_concentration += n_solution[i];
        if (n_solution[i] > 1e-15) nonzero_count++;
        if (n_solution[i] < 0.0) negative_count++;
        if (n_solution[i] > max_concentration) {
            max_concentration = n_solution[i];
            max_index = i;
        }
    }

    std::cout << "\nSolution summary:" << std::endl;
    std::cout << "  Total mass (sum n[i]*(i+1)): " << total_mass << std::endl;
    std::cout << "  Total concentration (sum n[i]): " << total_concentration << std::endl;
    std::cout << "  Non-zero components (>1e-15): " << nonzero_count << std::endl;
    std::cout << "  Negative components: " << negative_count << std::endl;
    std::cout << "  Max concentration: " << max_concentration << " at index " << max_index << std::endl;

    std::cout << "\nFirst 10 values:" << std::endl;
    for (int i = 0; i < std::min(10, max_size); i++) {
        std::cout << "  n[" << i << "] = " << std::scientific << std::setprecision(6) << n_solution[i] << std::endl;
    }

    save_solution(output_file, n_solution, max_size, time);
    delete[] n_0;

    std::cout << "\nDone!" << std::endl;
    return 0;
}
