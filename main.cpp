#include "python_var.hpp"
#include <iostream>

double kernel(uint64_t i, uint64_t j)
{
    return 1.0 / (i + j);
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

    n_0 = modeling(size, kernel, rel_tol, n_0, time, dt, mosaic_type);

    // for (int i = 0; i < size; i++) {
    //     std::cout << n_0[i] << " ";
    // }
    // std::cout << std::endl;

    delete[] n_0;

    return 0;
}
