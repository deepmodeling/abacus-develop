#include "module_base/complexmatrix.h"
#include "module_base/lapack_connector.h"
#include "module_hsolver/kernels/dngvd_op.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_psi/kernels/memory_op.h"

#include <algorithm>
#include <complex>
#include <fstream>
#include <gtest/gtest.h>
#include <iostream>
#include <vector>

class TestModuleHsolverMathDngvd : public ::testing::Test
{
  protected:
    using resize_memory_op_Z = psi::memory::resize_memory_op<std::complex<double>, psi::DEVICE_GPU>;
    using delete_memory_op_Z = psi::memory::delete_memory_op<std::complex<double>, psi::DEVICE_GPU>;
    using resize_memory_op_D = psi::memory::resize_memory_op<double, psi::DEVICE_GPU>;
    using delete_memory_op_D = psi::memory::delete_memory_op<double, psi::DEVICE_GPU>;
    // from CPU to GPU
    using synchronize_memory_op_C2G_Z
        = psi::memory::synchronize_memory_op<std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
    using synchronize_memory_op_C2G_D = psi::memory::synchronize_memory_op<double, psi::DEVICE_GPU, psi::DEVICE_CPU>;
    using synchronize_memory_op_G2C_Z
        = psi::memory::synchronize_memory_op<std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_GPU>;
    using synchronize_memory_op_G2C_D = psi::memory::synchronize_memory_op<double, psi::DEVICE_CPU, psi::DEVICE_GPU>;

    const psi::DEVICE_CPU* cpu_ctx = {};
    const psi::DEVICE_GPU* gpu_ctx = {};

    // prepare A & B in CPU
    std::vector<std::complex<double>> matrix_A = {
        {-0.351417, -1.73472},
        {-8.32667, 2.3744},
        {4.16334, 3.64292},
        {5.20417, -3.85976},
        {-8.32667, -2.3744},
        {0.551651, -2.60209},
        {2.08167, 1.9082},
        {-6.93889, 1.04083},
        {4.16334, -3.64292},
        {2.08167, -1.9082},
        {0.551651, -2.25514},
        {-1.31839, 5.20417},
        {5.20417, 3.85976},
        {-6.93889, -1.04083},
        {-1.31839, -5.20417},
        {0.551651, -3.64292}

        // {-4.280e-01,3.084e-17}, {-1.288e-16,9.021e-17}, {5.204e-18,-3.990e-17}, {-5.204e-17,-2.776e-17},
        // {-1.288e-16,-9.021e-17}, {4.574e-01,-8.687e-17}, {-6.939e-18,1.908e-17}, {8.327e-17,-1.041e-17},
        // {5.204e-18,3.990e-17}, {-6.939e-18,-1.908e-17}, {4.574e-01,-5.783e-17},{-6.939e-18,4.163e-17},
        // {-5.204e-17,2.776e-17}, {8.327e-17,1.041e-17}, {-6.939e-18,-4.163e-17}, {4.574e-01,-3.451e-17}

        // {-4.280e-01,3.084e-17}, {-1.288e-16,9.021e-17}, {5.204e-18,-3.990e-17}, {-1.145e-16,2.255e-17},
        // {3.946e-17,-2.949e-17}, {4.574e-01,-8.687e-17}, {2.082e-17,1.908e-17}, {4.163e-17,-4.163e-17},
        // {3.296e-17,-9.454e-17}, {-6.939e-18,-1.908e-17}, {4.574e-01,-5.783e-17},{-1.388e-17,6.939e-18},
        // {-5.204e-17,2.776e-17}, {8.327e-17,1.041e-17}, {-6.939e-18,-4.163e-17}, {4.574e-01,-3.451e-17}
    };
    std::vector<std::complex<double>> matrix_B = {
        {1, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {1, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {1, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {1, 0}
        // {1.000e+00,0.000e+00}, {-7.069e-17,-3.123e-17}, {-5.204e-17,0.000e+00}, {5.551e-17,-2.082e-17},
        // {-7.069e-17,3.123e-17}, {1.000e+00,0.000e+00}, {1.110e-16,-7.286e-17}, {8.327e-17,-1.110e-16},
        // {-5.204e-17,0.000e+00}, {1.110e-16,7.286e-17}, {1.000e+00,0.000e+00}, {-9.714e-17,2.776e-17},
        // {5.551e-17,2.082e-17}, {8.327e-17,1.110e-16}, {-9.714e-17,-2.776e-17}, {1.000e+00,0.000e+00}
    };
    const int matrix_size = 16;

    // prepare W & V in CPU in dngv_op
    std::vector<double> W_dngv_op = {0.0, 0.0, 0.0, 0.0};
    std::vector<std::complex<double>> matrix_V_dngv_op = {{0.0, 0.0},
                                                     {0.0, 0.0},
                                                     {0.0, 0.0},
                                                     {0.0, 0.0},
                                                     {0.0, 0.0},
                                                     {0.0, 0.0},
                                                     {0.0, 0.0},
                                                     {0.0, 0.0},
                                                     {0.0, 0.0},
                                                     {0.0, 0.0},
                                                     {0.0, 0.0},
                                                     {0.0, 0.0},
                                                     {0.0, 0.0},
                                                     {0.0, 0.0},
                                                     {0.0, 0.0},
                                                     {0.0, 0.0}};

    // prepare W & V in CPU in dngvx_op
    std::vector<double> W_DNGVX = {0.0, 0.0};
    std::vector<std::complex<double>> matrix_V_DNGVX = {{0.0, 0.0},
                                                   {0.0, 0.0},
                                                   {0.0, 0.0},
                                                   {0.0, 0.0},
                                                   {0.0, 0.0},
                                                   {0.0, 0.0},
                                                   {0.0, 0.0},
                                                   {0.0, 0.0},
                                                   {0.0, 0.0},
                                                   {0.0, 0.0},
                                                   {0.0, 0.0},
                                                   {0.0, 0.0},
                                                   {0.0, 0.0},
                                                   {0.0, 0.0},
                                                   {0.0, 0.0},
                                                   {0.0, 0.0}};
};



#if __UT_USE_CUDA || __UT_USE_ROCM

TEST_F(TestModuleHsolverMathDngvd, transpose_gpu)
{
    // prepare transpose in GPU
    std::vector<std::complex<double>> transpose = {
        {-0.351417, -1.73472},
        {-8.32667, 2.3744},
        {4.16334, 3.64292},
        {-0.351417, -1.73472},
        {-8.32667, 2.3744},
        {4.16334, 3.64292},
        // {-0.351417,-1.73472}, {-8.32667,2.3744}, {4.16334,3.64292}
    };
    std::complex<double>* device_transpose = nullptr;
    resize_memory_op_Z()(gpu_ctx, device_transpose, matrix_size);
    synchronize_memory_op_C2G_Z()(gpu_ctx, cpu_ctx, device_transpose, transpose.data(), transpose.size());

    // run
    hsolver::createGpuBlasHandle();
    hsolver::matrixTranspose_op<double, psi::DEVICE_GPU>()(gpu_ctx, 2, 3, device_transpose, device_transpose);
    hsolver::destoryBLAShandle();

    // copy transpose data from GPU to CPU
    std::vector<std::complex<double>> transpose_result = {
        {0.0, 0.0},
        {0.0, 0.0},
        {0.0, 0.0},
        {0.0, 0.0},
        {0.0, 0.0},
        {0.0, 0.0},
        // {0.0,  0.0}, {0.0,  0.0}, {0.0,  0.0}
    };
    synchronize_memory_op_G2C_Z()(cpu_ctx, gpu_ctx, transpose_result.data(), device_transpose, transpose.size());

    // std::vector<std::complex<double> > test_result = {
    //     {-0.351417,-1.73472}, {-0.351417,-1.73472}, {-0.351417,-1.73472},
    //     {-8.32667,2.3744}, {-8.32667,2.3744}, {-8.32667,2.3744},
    //     {4.16334,3.64292}, {4.16334,3.64292}, {4.16334,3.64292},
    // };
    std::vector<std::complex<double>> test_result = {
        {-0.351417, -1.73472},
        {-0.351417, -1.73472},
        {-8.32667, 2.3744},
        {-8.32667, 2.3744},
        {4.16334, 3.64292},
        {4.16334, 3.64292},
    };

    // for (int i = 0; i < 3; i++)
    // {
    //    for (int j = 0; j < 2; j++)
    //     {
    //         std::cout << transpose_result[i * 2 + j];
    //     }
    //     std::cout << std::endl;
    // }

    for (int i = 0; i < transpose_result.size(); i++)
    {
        EXPECT_LT(fabs(test_result[i].imag()) - fabs(transpose_result[i].imag()), 1e-8);
        EXPECT_LT(fabs(test_result[i].real()) - fabs(transpose_result[i].real()), 1e-8);
    }
}

#endif // __UT_USE_CUDA || __UT_USE_ROCM
