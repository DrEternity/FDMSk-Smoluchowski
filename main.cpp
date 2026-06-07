/*
 * Generate Reference Solution for Atmospheric Kernel
 */

#include "python_var.hpp"
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
    int max_size = 65536;
    double time = 1000.0;
    double dt = 0.0001;
    double rel_tol = 1e-10;
    double ode_tol = 1e-6;
    int initial_size = 512;
    double mass_threshold = 1e-2;
    std::string output_file = "reference_solution_atmos.txt";

    // Parse: max_size time output_file initial_size mass_threshold ode_tol
    if (argc > 1) max_size = std::stoi(argv[1]);
    if (argc > 2) time = std::stod(argv[2]);
    if (argc > 3) output_file = argv[3];
    if (argc > 4) initial_size = std::stoi(argv[4]);
    if (argc > 5) mass_threshold = std::stod(argv[5]);
    if (argc > 6) ode_tol = std::stod(argv[6]);

    std::cout << "=== Reference Solution Generator ===" << std::endl;
    std::cout << "Kernel: Atmospheric" << std::endl;
    std::cout << "Max size: " << max_size << std::endl;
    std::cout << "Initial size: " << initial_size << std::endl;
    std::cout << "Time: " << time << std::endl;
    std::cout << "MSk tolerance: " << rel_tol << std::endl;
    std::cout << "ODE tolerance: " << ode_tol << std::endl;
    std::cout << "Mass threshold: " << mass_threshold << std::endl;
    std::cout << "dt: " << dt << std::endl;
    std::cout << "Output: " << output_file << std::endl;
    std::cout << "====================================\n" << std::endl;

    // Initial condition: monodisperse (allocate max_size, zero-initialized)
    double* n_0 = new double[max_size]();
    n_0[0] = 1.0;

    std::cout << "Initial condition: n[0] = " << n_0[0] << ", others = 0" << std::endl;
    std::cout << "Starting ODE integration..." << std::endl;

    // Atmospheric kernel + tridiag (rho=2.0): the configuration used to generate
    // the reference_solution_atmos / new_reference_solution_atmos series.
    MosaicType mosaic_type{MosaicType::tridiag};
    double* n_solution = modeling(max_size, kernel_atmos, rel_tol, n_0, time, dt,
                                  mosaic_type, initial_size, mass_threshold, ode_tol);
    //MosaicType mosaic_type{MosaicType::monodiag};
    //double* n_solution = modeling(max_size, kernel_ballistic, rel_tol, n_0, time, dt,
    //                              mosaic_type, initial_size, mass_threshold, ode_tol);


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
