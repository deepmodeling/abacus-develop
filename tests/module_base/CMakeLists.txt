########################################
# CMake build system
# This file is part of ABACUS
cmake_minimum_required(VERSION 3.18)
########################################
# project name
project(module_base_test)
# enable cmake test
enable_testing()
# link the libraries anyway
set(CMAKE_LINK_WHAT_YOU_USE TRUE)
# set extra modules
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../../modules)

# load the package needed by the current UT module
find_package(MPI)
find_package(OpenMP)
find_package(Threads)
find_package(MPI REQUIRED)
find_package(ELPA REQUIRED)
find_package(GTest REQUIRED)
find_package(Cereal REQUIRED)

# compile the base module
set(libname "module_base")
set(LIB_BASE_DIR ${CMAKE_SOURCE_DIR}/../../source/module_base)
include_directories(${LIB_BASE_DIR})
file(GLOB LIB_SRC 
     ${LIB_BASE_DIR}/complexarray.cpp 
     ${LIB_BASE_DIR}/global_variable.cpp
     ${LIB_BASE_DIR}/math_polyint.cpp
     ${LIB_BASE_DIR}/matrix.cpp 
     ${LIB_BASE_DIR}/sph_bessel.cpp
     ${LIB_BASE_DIR}/tool_title.cpp 
     ${LIB_BASE_DIR}/complexmatrix.cpp
     ${LIB_BASE_DIR}/intarray.cpp 
     ${LIB_BASE_DIR}/math_sphbes.cpp 
     ${LIB_BASE_DIR}/memory.cpp
     ${LIB_BASE_DIR}/sph_bessel_recursive-d1.cpp
     ${LIB_BASE_DIR}/ylm.cpp
     ${LIB_BASE_DIR}/element_basis_index.cpp
     ${LIB_BASE_DIR}/integral.cpp
     ${LIB_BASE_DIR}/math_ylmreal.cpp
     ${LIB_BASE_DIR}/sph_bessel_recursive-d2.cpp
     ${LIB_BASE_DIR}/export.cpp
     ${LIB_BASE_DIR}/inverse_matrix.cpp
     ${LIB_BASE_DIR}/timer.cpp
     ${LIB_BASE_DIR}/global_file.cpp
     ${LIB_BASE_DIR}/main.cpp
     ${LIB_BASE_DIR}/mathzone.cpp
     ${LIB_BASE_DIR}/polint.cpp
     ${LIB_BASE_DIR}/tool_check.cpp
     ${LIB_BASE_DIR}/global_function.cpp
     ${LIB_BASE_DIR}/math_integral.cpp
     ${LIB_BASE_DIR}/matrix3.cpp
     ${LIB_BASE_DIR}/realarray.cpp
     ${LIB_BASE_DIR}/tool_quit.cpp
     ${LIB_BASE_DIR}/opt_CG.cpp
     ${LIB_BASE_DIR}/opt_DCsrch.cpp
)
add_library(${libname} ${LIB_SRC})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


# compile the test executable
file(GLOB TEST_SRC 
     test_blas_connector.cpp
     test_matrix3.cpp
)

add_executable(runUnitTests ${TEST_SRC})


target_link_libraries(
    runUnitTests 
    gtest 
    gtest_main 
    ${libname} 
    pthread 
    coverage_config
)

# link the third party library
if(DEFINED MKL_DIR)
    set(MKL_DIR_CACHE "${MKL_DIR}")
endif()

if(DEFINED MKL_DIR_CACHE)
    find_package(IntelMKL REQUIRED)
    add_definitions(-D__MKL -DMKL_ILP64)
    target_link_libraries(
        runUnitTests
        -lifcore
        MPI::MPI_CXX
        IntelMKL::MKL
        OpenMP::OpenMP_CXX
    )
else()
    find_package(LAPACK REQUIRED)
    find_package(ScaLAPACK REQUIRED)
    target_link_libraries(
        runUnitTests
        FFTW3::FFTW3
        MPI::MPI_CXX
        ScaLAPACK::ScaLAPACK
        LAPACK::LAPACK
        OpenMP::OpenMP_CXX
    )
endif()

add_test(runUnitTests runUnitTests)

# coverage setting
add_library(coverage_config INTERFACE)
target_compile_options(coverage_config INTERFACE
    -O0        # no optimization
    -g         # generate debug info
    --coverage # sets all required flags
)
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
target_link_options(coverage_config INTERFACE --coverage)
else()
target_link_libraries(coverage_config INTERFACE --coverage)
endif()

# setting gtest environ
if(NOT GTEST_LIBRARIES)
  configure_file(../../cmake/googletest.cmake.in googletest-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
  if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
  if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
  endif()
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src ${CMAKE_CURRENT_BINARY_DIR}/googletest-build EXCLUDE_FROM_ALL)
else ()
  include_directories(${GTEST_INCLUDE_DIRS})
endif ()
