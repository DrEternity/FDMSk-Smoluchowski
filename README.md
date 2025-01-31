# FDMSk-Smoluchowski
This repository contains the code for the paper titled "Mosaic-skeleton approximation is all you need for Smoluchowski equations", see at [arXiv.org](https://arxiv.org/abs/2501.10206).


## Usage guide
The C++ example is built using `CMake`.
Requirements:
* C++, C and Fortran compilers;
* CMake, BLAS, LAPACK, fftw3, pthreads libraries.


Build
```bash
cmake -S . -B build
cmake --build build --target example
```
and run
```bash
./example
```
