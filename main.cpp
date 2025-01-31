#include "python_var.hpp"
#include <iostream>


double kernel_atmos(uint64_t i, uint64_t j) {
    if (i == j) {
        return 4;
    }
    double num = (i + j) * std::pow(std::pow(i, 1.0 / 3.0) + std::pow(j, 1.0 / 3.0), 2.0 / 3.0);
    double denom = std::pow(i * j, 5.0 / 9.0) * std::abs(std::pow(i, 2.0 / 3.0) - std::pow(j, 2.0 / 3.0));
    return num / denom;
}

int main()
{
    int size = 256;
    double rel_tol = 1e-6;
    double* n_0 = new double[size];
    n_0[0] = 1;
    double time = 1.0;
    double dt = 0.001;
    MosaicType mosaic_type{MosaicType::monodiag};

    n_0 = modeling(size, kernel_atmos, rel_tol, n_0, time, dt, mosaic_type);

    // for (int i = 0; i < size; i++) {
    //     std::cout << n_0[i] << " ";
    // }
    // std::cout << std::endl;

    delete[] n_0;

    return 0;
}
