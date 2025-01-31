cmake_minimum_required(VERSION 3.20)
project(mosaic_skeleton_examples
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)


if(NOT TARGET zaimsk)
        find_package(zaimsk)
endif()
if(NOT zaimsk_FOUND)
        include(FetchContent)

        set(MSK_USE_FFT ON CACHE BOOL "USE FFT")
        set(MSK_USE_MPI OFF CACHE BOOL "Do not use MPI")
        set(MSK_BUILD_DOCS OFF CACHE BOOL "Do not build docs")
        set(MSK_BUILD_C_INTERFACE OFF CACHE BOOL "Do not build C interface")
        set(MSK_BUILD_FORTRAN_INTERFACE OFF CACHE BOOL "Do not build Fortran interface")

        FetchContent_Declare(
                zaimsk
                GIT_REPOSITORY https://gitlab.com/bulatral/mosaic-skeleton.git
                GIT_TAG        main
                PREFIX ${PROJECT_SOURCE_DIR}/external/mosaic-skeleton
                SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/mosaic-skeleton
                BINARY_DIR ${PROJECT_BINARY_DIR}/external/mosaic-skeleton
                TMP_DIR ${PROJECT_BINARY_DIR}/external/mosaic-skeleton/tmp
                STAMP_DIR ${PROJECT_BINARY_DIR}/external/mosaic-skeleton/stamp
                INSTALL_DIR ${PROJECT_SOURCE_DIR}/external/mosaic-skeleton
                LOG_DIR ${PROJECT_BINARY_DIR}/external/mosaic-skeleton/log
                CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
        )
        FetchContent_MakeAvailable(zaimsk)
        set(zaimsk_FOUND zaimsk_POPULATED)

        include(FindPackageHandleStandardArgs)
        find_package_handle_standard_args(zaimsk DEFAULT_MSG zaimsk_FOUND)
endif()


set(EXAMPLE_SOURCE
        main.cpp)
add_executable(example ${EXAMPLE_SOURCE})
target_link_libraries(example PRIVATE zaimsk::zaimsk)
