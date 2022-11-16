#include "module_base/complexmatrix.h"
#include "module_base/lapack_connector.h"
#include "module_hsolver/include/dngvd_op.h"
#include "module_hsolver/include/math_kernel.h"
#include "module_psi/include/memory.h"

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
    std::vector<complex<double>> matrix_A = {
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
    std::vector<complex<double>> matrix_B = {
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
    std::vector<complex<double>> matrix_V_dngv_op = {{0.0, 0.0},
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
    std::vector<complex<double>> matrix_V_DNGVX = {{0.0, 0.0},
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

TEST_F(TestModuleHsolverMathDngvd, dngv_cpu)
{
    // (1)
    ModuleBase::ComplexMatrix B(4, 4);
    ModuleBase::ComplexMatrix V(4, 4);
    // ==================================
    int col = 4;
    int row = 4;
    int lwork = 0;
    int nb = LapackConnector::ilaenv(1, "ZHETRD", "L", col, -1, -1, -1);
    if (nb < 1)
    {
        nb = std::max(1, col);
    }
    if (nb == 1 || nb >= col)
    {
        lwork = 2 * col; // mohan modify 2009-08-02
    }
    else
    {
        lwork = (nb + 1) * col;
    }
    std::complex<double>* work = new std::complex<double>[lwork];
    ModuleBase::GlobalFunc::ZEROS(work, lwork);
    int info = 0;
    int rwork_dim;
    rwork_dim = 3 * col - 2;
    double* rwork = new double[rwork_dim];
    ModuleBase::GlobalFunc::ZEROS(rwork, rwork_dim);
    // V.c = matrix_A.data();
    // B.c = matrix_B.data();
    for (int i = 0; i < matrix_size; i++)
    {
        V.c[i] = matrix_A[i];
    }
    for (int i = 0; i < matrix_size; i++)
    {
        B.c[i] = matrix_B[i];
    }
    LapackConnector::zhegv(1, 'V', 'U', col, V, col, B, col, W_dngv_op.data(), work, lwork, rwork, info);
    // ==================================

    // (2)
    std::vector<double> W_result = {0.0, 0.0, 0.0, 0.0};
    std::vector<complex<double>> V_result = {{0.0, 0.0},
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
    hsolver::dngv_op<double, psi::DEVICE_CPU>()(cpu_ctx,
                                                4,
                                                4,
                                                matrix_A.data(),
                                                matrix_B.data(),
                                                W_result.data(),
                                                V_result.data());

    // output
    // std::cout << W_dngv_op[0] << "\t" <<  W_dngv_op[1] <<  W_dngv_op[2] <<  W_dngv_op[3] << std::endl;
    // std::cout << W_result[0] << "\t" <<  W_result[1] <<  W_result[2] <<  W_result[3] << std::endl;
    // std::cout << "CPU vector" << std::endl;
    // for (int i = 0; i < 4; i++)
    // {
    //     for (int j = 0; j < 4; j++)
    //     {
    //         std::cout << V_result[i * 4 + j] << ", ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "GPU vector" << std::endl;
    // for (int i = 0; i < 4; i++)
    // {
    //     for (int j = 0; j < 4; j++)
    //     {
    //         std::cout << V(i, j) << ", ";
    //     }
    //     std::cout << std::endl;
    // }

    // test
    for (int i = 0; i < W_dngv_op.size(); i++)
    {
        EXPECT_LT(fabs(W_dngv_op[i] - W_result[i]), 1e-8);
    }
    for (int i = 0; i < V_result.size(); i++)
    {
        EXPECT_LT(fabs(V.c[i].imag() - V_result[i].imag()), 1e-8);
        EXPECT_LT(fabs(V.c[i].real() - V_result[i].real()), 1e-8);
    }
    delete[] rwork;
    delete[] work;
}

TEST_F(TestModuleHsolverMathDngvd, dngvx_cpu)
{
    // (1)
    ModuleBase::ComplexMatrix A(4, 4);
    ModuleBase::ComplexMatrix B(4, 4);
    ModuleBase::ComplexMatrix V(4, 4);
    // ==================================
    int lwork;
    int info = 0;
    std::string name1 = "ZHETRD";
    std::string name2 = "L";
    int col = 4;
    int nb = LapackConnector::ilaenv(1, name1.c_str(), name2.c_str(), col, -1, -1, -1);
    if (nb < 1)
    {
        nb = std::max(1, col);
    }
    if (nb == 1 || nb >= col)
    {
        lwork = 2 * col; // qianrui fix a bug 2021-7-25 : lwork should be at least max(1,2*n)
    }
    else
    {
        lwork = (nb + 1) * col;
    }
    std::complex<double>* work = new std::complex<double>[2 * lwork];
    assert(work != 0);
    double* rwork = new double[7 * col];
    assert(rwork != 0);
    int* iwork = new int[5 * col];
    assert(iwork != 0);
    int* ifail = new int[col];
    assert(ifail != 0);
    ModuleBase::GlobalFunc::ZEROS(work, lwork); // qianrui change it, only first lwork numbers are used in zhegvx
    ModuleBase::GlobalFunc::ZEROS(rwork, 7 * col);
    ModuleBase::GlobalFunc::ZEROS(iwork, 5 * col);
    ModuleBase::GlobalFunc::ZEROS(ifail, col);
    for (int i = 0; i < matrix_size; i++)
    {
        A.c[i] = matrix_A[i];
    }
    for (int i = 0; i < matrix_size; i++)
    {
        B.c[i] = matrix_B[i];
    }
    for (int i = 0; i < matrix_size; i++)
    {
        V.c[i] = 0.0;
    }
    int m = 2;
    LapackConnector::zhegvx(
        1, // ITYPE = 1:  A*x = (lambda)*B*x
        'V', // JOBZ = 'V':  Compute eigenvalues and eigenvectors.
        'I', // RANGE = 'I': the IL-th through IU-th eigenvalues will be found.
        'L', // UPLO = 'L':  Lower triangles of A and B are stored.
        col, // N = base
        A, // A is COMPLEX*16 array  dimension (LDA, N)
        col, // LDA = base
        B, // B is COMPLEX*16 array, dimension (LDB, N)
        col, // LDB = base
        0.0, // Not referenced if RANGE = 'A' or 'I'.
        0.0, // Not referenced if RANGE = 'A' or 'I'.
        1, // IL: If RANGE='I', the index of the smallest eigenvalue to be returned. 1 <= IL <= IU <= N,
        m, // IU: If RANGE='I', the index of the largest eigenvalue to be returned. 1 <= IL <= IU <= N,
        0.0, // ABSTOL
        m, // M: The total number of eigenvalues found.  0 <= M <= N. if RANGE = 'I', M = IU-IL+1.
        W_DNGVX.data(), // W store eigenvalues
        V, // store eigenvector
        col, // LDZ: The leading dimension of the array Z.
        work,
        lwork,
        rwork,
        iwork,
        ifail,
        info);

    // ==================================

    // (2)
    std::vector<double> W_result = {0.0, 0.0};
    std::vector<complex<double>> V_result = {{0.0, 0.0},
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
    hsolver::dngvx_op<double, psi::DEVICE_CPU>()(cpu_ctx,
                                                 4,
                                                 4,
                                                 matrix_A.data(),
                                                 matrix_B.data(),
                                                 2,
                                                 W_result.data(),
                                                 V_result.data());

    // output
    // std::cout << W_result[0] << "\t" <<  W_result[1]  << std::endl;
    // std::cout << W_DNGVX[0] << "\t" <<  W_DNGVX[1]<< std::endl;
    // std::cout << "CPU vector" << std::endl;
    // for (int i = 0; i < 4; i++)
    // {
    //     for (int j = 0; j < 4; j++)
    //     {
    //         std::cout << V_result[i * 4 + j] << ", ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "GPU vector" << std::endl;
    // for (int i = 0; i < 4; i++)
    // {
    //     for (int j = 0; j < 4; j++)
    //     {
    //         std::cout << V(i, j) << ", ";
    //     }
    //     std::cout << std::endl;
    // }

    // test
    for (int i = 0; i < W_result.size(); i++)
    {
        EXPECT_LT(fabs(W_DNGVX[i] - W_result[i]), 1e-8);
    }
    for (int i = 0; i < V_result.size(); i++)
    {
        EXPECT_LT(fabs(V.c[i].imag() - V_result[i].imag()), 1e-8);
        EXPECT_LT(fabs(V.c[i].real() - V_result[i].real()), 1e-8);
    }
    delete[] work;
    delete[] rwork;
    delete[] iwork;
    delete[] ifail;
}

#if __UT_USE_CUDA || __UT_USE_ROCM

// computes all the eigenvalues and eigenvectors
TEST_F(TestModuleHsolverMathDngvd, dngv_gpu)
{
    // prepare A & B in GPU
    std::complex<double>* device_matrix_A = nullptr;
    std::complex<double>* device_matrix_B = nullptr;
    resize_memory_op_Z()(gpu_ctx, device_matrix_A, matrix_size);
    resize_memory_op_Z()(gpu_ctx, device_matrix_B, matrix_size);
    synchronize_memory_op_C2G_Z()(gpu_ctx, cpu_ctx, device_matrix_A, matrix_A.data(), matrix_size);
    synchronize_memory_op_C2G_Z()(gpu_ctx, cpu_ctx, device_matrix_B, matrix_B.data(), matrix_size);

    // prepare W & V in GPU in dngv_op
    double* device_W_dngv_op = nullptr;
    resize_memory_op_D()(gpu_ctx, device_W_dngv_op, W_dngv_op.size());
    psi::memory::set_memory_op<double, psi::DEVICE_GPU>()(gpu_ctx, device_W_dngv_op, 0, W_dngv_op.size());
    std::complex<double>* device_matrix_V_dngv_op = nullptr;
    resize_memory_op_Z()(gpu_ctx, device_matrix_V_dngv_op, matrix_V_dngv_op.size());
    psi::memory::set_memory_op<complex<double>, psi::DEVICE_GPU>()(gpu_ctx,
                                                                   device_matrix_V_dngv_op,
                                                                   0,
                                                                   matrix_V_dngv_op.size());

    // run in GPU
    hsolver::dngv_op<double, psi::DEVICE_GPU>()(gpu_ctx,
                                                4,
                                                4,
                                                device_matrix_A,
                                                device_matrix_B,
                                                device_W_dngv_op,
                                                device_matrix_V_dngv_op);
    // copy W data from GPU to CPU
    std::vector<double> W_result = {0.0, 0.0, 0.0, 0.0};
    synchronize_memory_op_G2C_D()(cpu_ctx, gpu_ctx, W_result.data(), device_W_dngv_op, W_result.size());
    // copy V data from GPU to CPU
    std::vector<complex<double>> V_result = {{0.0, 0.0},
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
    synchronize_memory_op_G2C_Z()(cpu_ctx, gpu_ctx, V_result.data(), device_matrix_V_dngv_op, V_result.size());

    // run in CPU
    hsolver::dngv_op<double, psi::DEVICE_CPU>()(cpu_ctx,
                                                4,
                                                4,
                                                matrix_A.data(),
                                                matrix_B.data(),
                                                W_dngv_op.data(),
                                                matrix_V_dngv_op.data());

    // // CPU
    // std::cout << W_dngv_op[0] << "\t" <<  W_dngv_op[1] << "\t"  <<  W_dngv_op[2] << "\t"  <<  W_dngv_op[3] <<
    // std::endl;
    // // GPU
    // std::cout << W_result[0] << "\t" <<  W_result[1] << "\t"  <<  W_result[2] << "\t"  <<  W_result[3] << std::endl;
    // std::cout << "CPU vector" << std::endl;
    // for (int i = 0; i < 4; i++)
    // {
    //     for (int j = 0; j < 4; j++)
    //     {
    //         std::cout << matrix_V_dngv_op[i * 4 + j] << ", ";

    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "GPU vector" << std::endl;
    // for (int i = 0; i < 4; i++)
    // {
    //     for (int j = 0; j < 4; j++)
    //     {
    //         std::cout << V_result[i * 4 + j] << ", ";
    //     }
    //     std::cout << std::endl;
    // }

    // we need to compare
    //          1. W with W_result
    //          2. matrix_V with V_result
    for (int i = 0; i < W_dngv_op.size(); i++)
    {
        EXPECT_LT(fabs(W_dngv_op[i] - W_result[i]), 1e-8);
    }
    for (int i = 0; i < matrix_V_dngv_op.size(); i++)
    {
        EXPECT_LT(fabs(matrix_V_dngv_op[i].imag()) - fabs(V_result[i].imag()), 1e-8);
        EXPECT_LT(fabs(matrix_V_dngv_op[i].real()) - fabs(V_result[i].real()), 1e-8);
    }

    // delete values in GPU
    delete_memory_op_Z()(gpu_ctx, device_matrix_A);
    delete_memory_op_Z()(gpu_ctx, device_matrix_B);
    delete_memory_op_Z()(gpu_ctx, device_matrix_V_dngv_op);
    delete_memory_op_D()(gpu_ctx, device_W_dngv_op);
}

// computes the first m eigenvalues ​​and their corresponding eigenvectors
TEST_F(TestModuleHsolverMathDngvd, dngvx_gpu)
{
    // prepare A & B in GPU
    std::complex<double>* device_matrix_A = nullptr;
    std::complex<double>* device_matrix_B = nullptr;
    resize_memory_op_Z()(gpu_ctx, device_matrix_A, matrix_size);
    resize_memory_op_Z()(gpu_ctx, device_matrix_B, matrix_size);
    synchronize_memory_op_C2G_Z()(gpu_ctx, cpu_ctx, device_matrix_A, matrix_A.data(), matrix_size);
    synchronize_memory_op_C2G_Z()(gpu_ctx, cpu_ctx, device_matrix_B, matrix_B.data(), matrix_size);
    // prepare W & V in GPU in dngvx_op
    double* device_W_DNGVX = nullptr;
    resize_memory_op_D()(gpu_ctx, device_W_DNGVX, W_DNGVX.size());
    synchronize_memory_op_C2G_D()(gpu_ctx, cpu_ctx, device_W_DNGVX, W_DNGVX.data(), W_DNGVX.size());
    std::complex<double>* device_matrix_V_DNGVX = nullptr;
    resize_memory_op_Z()(gpu_ctx, device_matrix_V_DNGVX, matrix_V_DNGVX.size());
    synchronize_memory_op_C2G_Z()(gpu_ctx,
                                  cpu_ctx,
                                  device_matrix_V_DNGVX,
                                  matrix_V_DNGVX.data(),
                                  matrix_V_DNGVX.size());

    // run in GPU
    hsolver::dngvx_op<double, psi::DEVICE_GPU>()(gpu_ctx,
                                                 4,
                                                 4,
                                                 device_matrix_A,
                                                 device_matrix_B,
                                                 2,
                                                 device_W_DNGVX,
                                                 device_matrix_V_DNGVX);
    // copy W data from GPU to CPU
    std::vector<double> W_result = {0.0, 0.0};
    synchronize_memory_op_G2C_D()(cpu_ctx, gpu_ctx, W_result.data(), device_W_DNGVX, W_result.size());
    // copy V data from GPU to CPU
    std::vector<complex<double>> V_result = {{0.0, 0.0},
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
    synchronize_memory_op_G2C_Z()(cpu_ctx, gpu_ctx, V_result.data(), device_matrix_V_DNGVX, V_result.size());

    // run in CPU
    hsolver::dngvx_op<double, psi::DEVICE_CPU>()(cpu_ctx,
                                                 4,
                                                 4,
                                                 matrix_A.data(),
                                                 matrix_B.data(),
                                                 2,
                                                 W_DNGVX.data(),
                                                 matrix_V_DNGVX.data());

    // std::cout << "GPU::" << std::endl;
    // std::cout << W_result[0] << "\t" <<  W_result[1] << std::endl;
    // std::cout << "CPU::" << std::endl;
    // std::cout << W_DNGVX[0] << "\t" <<  W_DNGVX[1] << std::endl;
    // std::cout << "GPU vector" << std::endl;
    // for (int i = 0; i < 4; i++)
    // {
    //     for (int j = 0; j < 4; j++)
    //     {
    //         std::cout << V_result[i * 4 + j] << ", ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "CPU vector" << std::endl;
    // for (int i = 0; i < 4; i++)
    // {
    //     for (int j = 0; j < 4; j++)
    //     {
    //         std::cout << matrix_V_DNGVX[i * 4 + j] << ", ";
    //     }
    //     std::cout << std::endl;
    // }

    // we need to compare
    //          1. W with W_result
    //          2. matrix_V with V_result
    for (int i = 0; i < W_DNGVX.size(); i++)
    {
        EXPECT_LT(fabs(W_DNGVX[i] - W_result[i]), 1e-8);
    }
    for (int i = 0; i < matrix_V_dngv_op.size(); i++)
    {
        EXPECT_LT(fabs(matrix_V_dngv_op[i].imag()) - fabs(V_result[i].imag()), 1e-8);
        EXPECT_LT(fabs(matrix_V_dngv_op[i].real()) - fabs(V_result[i].real()), 1e-8);
    }

    // delete values in GPU
    delete_memory_op_Z()(gpu_ctx, device_matrix_A);
    delete_memory_op_Z()(gpu_ctx, device_matrix_B);
    delete_memory_op_Z()(gpu_ctx, device_matrix_V_DNGVX);
    delete_memory_op_D()(gpu_ctx, device_W_DNGVX);
}

TEST_F(TestModuleHsolverMathDngvd, transpose_gpu)
{
    // prepare transpose in GPU
    std::vector<complex<double>> transpose = {
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
    hsolver::createBLAShandle();
    hsolver::matrixTranspose_op<double, psi::DEVICE_GPU>()(gpu_ctx, 2, 3, device_transpose, device_transpose);
    hsolver::destoryBLAShandle();

    // copy transpose data from GPU to CPU
    std::vector<complex<double>> transpose_result = {
        {0.0, 0.0},
        {0.0, 0.0},
        {0.0, 0.0},
        {0.0, 0.0},
        {0.0, 0.0},
        {0.0, 0.0},
        // {0.0,  0.0}, {0.0,  0.0}, {0.0,  0.0}
    };
    synchronize_memory_op_G2C_Z()(cpu_ctx, gpu_ctx, transpose_result.data(), device_transpose, transpose.size());

    // std::vector<complex<double> > test_result = {
    //     {-0.351417,-1.73472}, {-0.351417,-1.73472}, {-0.351417,-1.73472},
    //     {-8.32667,2.3744}, {-8.32667,2.3744}, {-8.32667,2.3744},
    //     {4.16334,3.64292}, {4.16334,3.64292}, {4.16334,3.64292},
    // };
    std::vector<complex<double>> test_result = {
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
