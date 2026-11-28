#include "source_hsolver/kernels/cuda/diag_cusolver.cuh"
#include "gtest/gtest.h"
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <memory>
#include <string>

namespace
{
bool reject_allocations = false;
int rejected_allocations = 0;
cudaMemcpyKind failing_copy = cudaMemcpyDefault;
bool fail_gemm = false;
}

extern "C" cudaError_t __real_cudaMalloc(void** pointer, std::size_t bytes);
extern "C" cudaError_t __wrap_cudaMalloc(void** pointer, std::size_t bytes)
{
    if (reject_allocations)
    {
        ++rejected_allocations;
        return cudaErrorMemoryAllocation;
    }
    return __real_cudaMalloc(pointer, bytes);
}

extern "C" cudaError_t __real_cudaMemcpy(void* destination, const void* source, std::size_t bytes, cudaMemcpyKind kind);
extern "C" cudaError_t __wrap_cudaMemcpy(void* destination, const void* source, std::size_t bytes, cudaMemcpyKind kind)
{
    if (kind == failing_copy)
    {
        failing_copy = cudaMemcpyDefault;
        return cudaErrorInvalidValue;
    }
    return __real_cudaMemcpy(destination, source, bytes, kind);
}

extern "C" cublasStatus_t __real_cublasDgemm_v2(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb,
                                               int m, int n, int k, const double* alpha, const double* a, int lda,
                                               const double* b, int ldb, const double* beta, double* c, int ldc);
extern "C" cublasStatus_t __wrap_cublasDgemm_v2(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb,
                                               int m, int n, int k, const double* alpha, const double* a, int lda,
                                               const double* b, int ldb, const double* beta, double* c, int ldc)
{
    if (fail_gemm)
    {
        fail_gemm = false;
        return CUBLAS_STATUS_EXECUTION_FAILED;
    }
    return __real_cublasDgemm_v2(handle, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

class CusolverDensityFailure : public testing::Test
{
protected:
    void SetUp() override
    {
        solver.reset(new Diag_Cusolver_gvd());
        const double h[4] = {1.0, 0.0, 0.0, 2.0};
        const double s[4] = {1.0, 0.0, 0.0, 1.0};
        double eigenvalues[2] = {};
        double wave[4] = {};
        solver->solve_full(2, 2, h, s, eigenvalues, wave, true);
    }

    void TearDown() override
    {
        reject_allocations = false;
        rejected_allocations = 0;
        failing_copy = cudaMemcpyDefault;
        fail_gemm = false;
        // The production caller terminates; consume injected errors between tests.
        cudaGetLastError();
        solver.reset();
    }

    void expect_failure(const std::string& stage, const std::string& detail)
    {
        const double weights[2] = {2.0, 0.125};
        double density[4] = {};
        std::string error;
        const bool success = solver->density_from_retained(2, weights, density, error);
        EXPECT_FALSE(success);
        EXPECT_NE(error.find(stage), std::string::npos) << error;
        EXPECT_NE(error.find(detail), std::string::npos) << error;
        EXPECT_NE(error.find("n=2, bands=2"), std::string::npos) << error;
        EXPECT_NE(error.find("device cpu"), std::string::npos) << error;
    }

    std::unique_ptr<Diag_Cusolver_gvd> solver;
};

TEST_F(CusolverDensityFailure, ReusesSolverArraysWithoutDeviceAllocations)
{
    // The DM stage must not allocate any new device arrays in the CUDA kernel code.
    reject_allocations = true;
    const double weights[2] = {2.0, 0.125};
    double density[4] = {};
    std::string error;
    const bool success = solver->density_from_retained(2, weights, density, error);
    EXPECT_TRUE(success) << error;
    EXPECT_EQ(rejected_allocations, 0);
    EXPECT_DOUBLE_EQ(density[0], 2.0);
    EXPECT_DOUBLE_EQ(density[3], 0.125);
}

TEST_F(CusolverDensityFailure, OccupationUploadReportsStage)
{
    failing_copy = cudaMemcpyHostToDevice;
    expect_failure("upload occupations", "invalid argument");
}

TEST_F(CusolverDensityFailure, DensityDownloadReportsStage)
{
    failing_copy = cudaMemcpyDeviceToHost;
    expect_failure("download density matrix", "invalid argument");
}

TEST_F(CusolverDensityFailure, GemmReportsLibraryError)
{
    fail_gemm = true;
    expect_failure("form density matrix (cublasDgemm)", "CUBLAS_STATUS_EXECUTION_FAILED");
}
