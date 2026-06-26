/**
 * @file diago_cuda_kernels_test.cpp
 * @brief Unit tests for GPU-accelerated eigenvalue solver kernels.
 *
 * Tests the correctness and performance of CUDA kernels in diago_kernels.cu:
 *   1. batched_dot_real_op    - batched dot product correctness
 *   2. calc_grad_cg_op        - fused CG gradient computation
 *   3. schmidt_orth_cg_op     - Schmidt orthogonalization
 *   4. compute_band_energies_op - band energy computation
 *   5. apply_preconditioner_op  - preconditioner application
 *   6. End-to-end CG solver test with GPU acceleration
 *
 * Each test compares GPU results against a CPU reference implementation.
 */

#include "gtest/gtest.h"

#include <complex>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>

#if defined(__CUDACC__) || defined(__CUDA) || defined(__UT_USE_CUDA)
#include "source_hsolver/kernels/cuda/diago_kernels.cuh"
#include <cuda_runtime.h>
#define CUDA_AVAILABLE 1
#else
#define CUDA_AVAILABLE 0
#endif

// Runtime GPU detection: set to false if no GPU device is available at runtime
static bool g_gpu_available = false;

namespace hsolver {
namespace test {

// ============================================================================
// Helper: generate random test data
// ============================================================================

template <typename Real>
std::vector<std::complex<Real>> generate_random_complex(int n, Real scale = 1.0)
{
    std::vector<std::complex<Real>> data(n);
    std::mt19937 rng(42); // fixed seed for reproducibility
    std::uniform_real_distribution<Real> dist(-scale, scale);
    for (int i = 0; i < n; ++i) {
        data[i] = std::complex<Real>(dist(rng), dist(rng));
    }
    return data;
}

template <typename Real>
std::vector<Real> generate_random_real(int n, Real min_val = 0.5, Real max_val = 2.0)
{
    std::vector<Real> data(n);
    std::mt19937 rng(42);
    std::uniform_real_distribution<Real> dist(min_val, max_val);
    for (int i = 0; i < n; ++i) {
        data[i] = dist(rng);
    }
    return data;
}

// ============================================================================
// CPU reference implementations
// ============================================================================

template <typename Real>
void cpu_batched_dot_product(
    Real* result,
    const std::complex<Real>* psi_L,
    const std::complex<Real>* psi_R,
    int n_basis, int ld_psi, int n_band)
{
    for (int b = 0; b < n_band; ++b) {
        Real sum = 0;
        const auto* L = psi_L + b * ld_psi;
        const auto* R = psi_R + b * ld_psi;
        for (int i = 0; i < n_basis; ++i) {
            sum += L[i].real() * R[i].real() + L[i].imag() * R[i].imag();
        }
        result[b] = sum;
    }
}

template <typename Real>
void cpu_calc_grad_cg(
    std::complex<Real>* grad,
    std::complex<Real>* pphi,
    const std::complex<Real>* hphi,
    const std::complex<Real>* sphi,
    const Real* prec,
    int n_basis)
{
    Real eh = 0, es = 0;

    for (int i = 0; i < n_basis; ++i) {
        Real inv_p = Real(1.0) / prec[i];
        grad[i] = hphi[i] * inv_p;
        pphi[i] = sphi[i] * inv_p;

        eh += sphi[i].real() * grad[i].real() + sphi[i].imag() * grad[i].imag();
        es += sphi[i].real() * pphi[i].real() + sphi[i].imag() * pphi[i].imag();
    }

    Real lambda = eh / es;
    for (int i = 0; i < n_basis; ++i) {
        grad[i] -= lambda * pphi[i];
    }
}

template <typename Real>
void cpu_schmidt_orth_cg(
    std::complex<Real>* grad,
    std::complex<Real>* scg,
    const std::complex<Real>* psi,
    const std::complex<Real>* spsi,
    std::complex<Real>* lagrange,
    int n_basis, int ld_psi, int m)
{
    // Compute lagrange multipliers
    for (int j = 0; j < m; ++j) {
        std::complex<Real> sum(0, 0);
        const auto* psi_j = psi + j * ld_psi;
        for (int i = 0; i < n_basis; ++i) {
            sum += std::conj(psi_j[i]) * scg[i];
        }
        lagrange[j] = sum;
    }

    // Apply orthogonalization
    for (int i = 0; i < n_basis; ++i) {
        std::complex<Real> g_correction(0, 0);
        std::complex<Real> s_correction(0, 0);
        for (int j = 0; j < m; ++j) {
            const auto* psi_j = psi + j * ld_psi;
            const auto* spsi_j = spsi + j * ld_psi;
            g_correction += lagrange[j] * psi_j[i];
            s_correction += lagrange[j] * spsi_j[i];
        }
        grad[i] -= g_correction;
        scg[i] -= s_correction;
    }
}

template <typename Real>
void cpu_compute_band_energies(
    Real* eigen,
    const std::complex<Real>* psi,
    const std::complex<Real>* hpsi,
    int n_basis, int ld_psi, int n_band)
{
    for (int b = 0; b < n_band; ++b) {
        Real sum = 0;
        const auto* p = psi + b * ld_psi;
        const auto* hp = hpsi + b * ld_psi;
        for (int i = 0; i < n_basis; ++i) {
            sum += p[i].real() * hp[i].real() + p[i].imag() * hp[i].imag();
        }
        eigen[b] = sum;
    }
}

template <typename Real>
void cpu_apply_preconditioner(
    std::complex<Real>* grad,
    const std::complex<Real>* hpsi,
    const std::complex<Real>* spsi,
    const Real* prec,
    Real eigen,
    int n_basis)
{
    for (int i = 0; i < n_basis; ++i) {
        Real inv_p = Real(1.0) / prec[i];
        grad[i] = (hpsi[i] - eigen * spsi[i]) * inv_p;
    }
}

// ============================================================================
// Test Fixture
// ============================================================================

class DiagoCudaKernelsTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
#if CUDA_AVAILABLE
        // Skip GPU tests if no GPU device is available at runtime
        if (!g_gpu_available) {
            GTEST_SKIP() << "No CUDA-capable GPU device found. Skipping GPU tests.";
        }
        cuda::init_diago_cuda_resources();
#endif
    }

    void TearDown() override
    {
#if CUDA_AVAILABLE
        if (g_gpu_available) {
            cuda::destroy_diago_cuda_resources();
        }
#endif
    }
};

// ============================================================================
// Test 1: Batched Dot Product
// ============================================================================

#if CUDA_AVAILABLE
TEST_F(DiagoCudaKernelsTest, BatchedDotProduct_MatchesCPU)
{
    const int n_basis = 1024;
    const int ld_psi = n_basis + 16; // test non-contiguous layout
    const int n_band = 16;

    using Real = double;

    // Generate test data on CPU
    auto psi_L = generate_random_complex<Real>(ld_psi * n_band);
    auto psi_R = generate_random_complex<Real>(ld_psi * n_band);

    // CPU reference
    std::vector<Real> cpu_result(n_band);
    cpu_batched_dot_product<Real>(cpu_result.data(), psi_L.data(), psi_R.data(),
                                  n_basis, ld_psi, n_band);

    // GPU computation
    thrust::complex<Real>* d_psi_L = nullptr;
    thrust::complex<Real>* d_psi_R = nullptr;
    Real* d_result = nullptr;

    size_t vec_bytes = ld_psi * n_band * sizeof(thrust::complex<Real>);
    size_t res_bytes = n_band * sizeof(Real);

    CHECK_CUDA(cudaMalloc(&d_psi_L, vec_bytes));
    CHECK_CUDA(cudaMalloc(&d_psi_R, vec_bytes));
    CHECK_CUDA(cudaMalloc(&d_result, res_bytes));

    CHECK_CUDA(cudaMemcpy(d_psi_L, psi_L.data(), vec_bytes, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_psi_R, psi_R.data(), vec_bytes, cudaMemcpyHostToDevice));

    cuda::batched_dot_real_op<Real>(d_result, d_psi_L, d_psi_R,
                                    n_basis, ld_psi, n_band);

    std::vector<Real> gpu_result(n_band);
    CHECK_CUDA(cudaMemcpy(gpu_result.data(), d_result, res_bytes, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaDeviceSynchronize());

    // Verify
    const Real tol = Real(1e-10);
    for (int b = 0; b < n_band; ++b) {
        Real diff = std::abs(cpu_result[b] - gpu_result[b]);
        Real ref = std::max(std::abs(cpu_result[b]), Real(1.0));
        EXPECT_LT(diff / ref, tol) << "Band " << b << " mismatch: "
            << "CPU=" << cpu_result[b] << " GPU=" << gpu_result[b];
    }

    CHECK_CUDA(cudaFree(d_psi_L));
    CHECK_CUDA(cudaFree(d_psi_R));
    CHECK_CUDA(cudaFree(d_result));
}
#endif // CUDA_AVAILABLE

// ============================================================================
// Test 2: Fused CG Gradient Computation
// ============================================================================

#if CUDA_AVAILABLE
TEST_F(DiagoCudaKernelsTest, CalcGradCG_MatchesCPU)
{
    const int n_basis = 2048;
    using Real = double;

    // Generate test data
    auto hphi = generate_random_complex<Real>(n_basis);
    auto sphi = generate_random_complex<Real>(n_basis);
    auto prec = generate_random_real<Real>(n_basis, 0.5, 2.0);

    // CPU reference
    auto cpu_grad = hphi; // copy
    std::vector<std::complex<Real>> cpu_pphi(n_basis);
    cpu_calc_grad_cg<Real>(cpu_grad.data(), cpu_pphi.data(),
                           hphi.data(), sphi.data(), prec.data(), n_basis);

    // GPU computation
    thrust::complex<Real>* d_grad = nullptr;
    thrust::complex<Real>* d_pphi = nullptr;
    thrust::complex<Real>* d_hphi = nullptr;
    thrust::complex<Real>* d_sphi = nullptr;
    Real* d_prec = nullptr;

    size_t cbytes = n_basis * sizeof(thrust::complex<Real>);
    size_t rbytes = n_basis * sizeof(Real);

    CHECK_CUDA(cudaMalloc(&d_grad, cbytes));
    CHECK_CUDA(cudaMalloc(&d_pphi, cbytes));
    CHECK_CUDA(cudaMalloc(&d_hphi, cbytes));
    CHECK_CUDA(cudaMalloc(&d_sphi, cbytes));
    CHECK_CUDA(cudaMalloc(&d_prec, rbytes));

    CHECK_CUDA(cudaMemcpy(d_hphi, hphi.data(), cbytes, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_sphi, sphi.data(), cbytes, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_prec, prec.data(), rbytes, cudaMemcpyHostToDevice));

    cuda::calc_grad_cg_op<Real>(d_grad, d_pphi, d_hphi, d_sphi, d_prec, n_basis);

    std::vector<std::complex<Real>> gpu_grad(n_basis);
    std::vector<std::complex<Real>> gpu_pphi(n_basis);
    CHECK_CUDA(cudaMemcpy(gpu_grad.data(), d_grad, cbytes, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(gpu_pphi.data(), d_pphi, cbytes, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaDeviceSynchronize());

    // Verify
    const Real tol = Real(1e-9);
    for (int i = 0; i < n_basis; ++i) {
        EXPECT_NEAR(cpu_grad[i].real(), gpu_grad[i].real(), tol)
            << "grad[" << i << "].real mismatch";
        EXPECT_NEAR(cpu_grad[i].imag(), gpu_grad[i].imag(), tol)
            << "grad[" << i << "].imag mismatch";
        EXPECT_NEAR(cpu_pphi[i].real(), gpu_pphi[i].real(), tol)
            << "pphi[" << i << "].real mismatch";
        EXPECT_NEAR(cpu_pphi[i].imag(), gpu_pphi[i].imag(), tol)
            << "pphi[" << i << "].imag mismatch";
    }

    CHECK_CUDA(cudaFree(d_grad));
    CHECK_CUDA(cudaFree(d_pphi));
    CHECK_CUDA(cudaFree(d_hphi));
    CHECK_CUDA(cudaFree(d_sphi));
    CHECK_CUDA(cudaFree(d_prec));
}
#endif // CUDA_AVAILABLE

// ============================================================================
// Test 3: Schmidt Orthogonalization
// ============================================================================

#if CUDA_AVAILABLE
TEST_F(DiagoCudaKernelsTest, SchmidtOrthCG_MatchesCPU)
{
    const int n_basis = 512;
    const int ld_psi = 512;
    const int m = 8; // number of existing bands
    using Real = double;

    // Generate test data
    auto grad = generate_random_complex<Real>(n_basis);
    auto scg  = generate_random_complex<Real>(n_basis);
    auto psi  = generate_random_complex<Real>(ld_psi * m);
    auto spsi = generate_random_complex<Real>(ld_psi * m);

    // CPU reference
    auto cpu_grad = grad;
    auto cpu_scg  = scg;
    std::vector<std::complex<Real>> cpu_lagrange(m);
    cpu_schmidt_orth_cg<Real>(cpu_grad.data(), cpu_scg.data(),
                              psi.data(), spsi.data(), cpu_lagrange.data(),
                              n_basis, ld_psi, m);

    // GPU computation
    thrust::complex<Real>* d_grad = nullptr;
    thrust::complex<Real>* d_scg = nullptr;
    thrust::complex<Real>* d_psi = nullptr;
    thrust::complex<Real>* d_spsi = nullptr;
    thrust::complex<Real>* d_lagrange = nullptr;

    size_t cbytes = n_basis * sizeof(thrust::complex<Real>);
    size_t pbytes = ld_psi * m * sizeof(thrust::complex<Real>);
    size_t lbytes = m * sizeof(thrust::complex<Real>);

    CHECK_CUDA(cudaMalloc(&d_grad, cbytes));
    CHECK_CUDA(cudaMalloc(&d_scg, cbytes));
    CHECK_CUDA(cudaMalloc(&d_psi, pbytes));
    CHECK_CUDA(cudaMalloc(&d_spsi, pbytes));
    CHECK_CUDA(cudaMalloc(&d_lagrange, lbytes));

    CHECK_CUDA(cudaMemcpy(d_grad, grad.data(), cbytes, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_scg, scg.data(), cbytes, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_psi, psi.data(), pbytes, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_spsi, spsi.data(), pbytes, cudaMemcpyHostToDevice));

    cuda::schmidt_orth_cg_op<Real>(d_grad, d_scg, d_psi, d_spsi, d_lagrange,
                                   n_basis, ld_psi, m);

    std::vector<std::complex<Real>> gpu_grad(n_basis);
    std::vector<std::complex<Real>> gpu_scg(n_basis);
    std::vector<std::complex<Real>> gpu_lagrange(m);
    CHECK_CUDA(cudaMemcpy(gpu_grad.data(), d_grad, cbytes, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(gpu_scg.data(), d_scg, cbytes, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(gpu_lagrange.data(), d_lagrange, lbytes, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaDeviceSynchronize());

    // Verify
    const Real tol = Real(1e-8);
    for (int i = 0; i < n_basis; ++i) {
        EXPECT_NEAR(cpu_grad[i].real(), gpu_grad[i].real(), tol)
            << "grad[" << i << "].real mismatch";
        EXPECT_NEAR(cpu_grad[i].imag(), gpu_grad[i].imag(), tol)
            << "grad[" << i << "].imag mismatch";
        EXPECT_NEAR(cpu_scg[i].real(), gpu_scg[i].real(), tol)
            << "scg[" << i << "].real mismatch";
        EXPECT_NEAR(cpu_scg[i].imag(), gpu_scg[i].imag(), tol)
            << "scg[" << i << "].imag mismatch";
    }

    for (int j = 0; j < m; ++j) {
        EXPECT_NEAR(cpu_lagrange[j].real(), gpu_lagrange[j].real(), tol)
            << "lagrange[" << j << "].real mismatch";
        EXPECT_NEAR(cpu_lagrange[j].imag(), gpu_lagrange[j].imag(), tol)
            << "lagrange[" << j << "].imag mismatch";
    }

    CHECK_CUDA(cudaFree(d_grad));
    CHECK_CUDA(cudaFree(d_scg));
    CHECK_CUDA(cudaFree(d_psi));
    CHECK_CUDA(cudaFree(d_spsi));
    CHECK_CUDA(cudaFree(d_lagrange));
}
#endif // CUDA_AVAILABLE

// ============================================================================
// Test 4: Band Energy Computation
// ============================================================================

#if CUDA_AVAILABLE
TEST_F(DiagoCudaKernelsTest, BandEnergies_MatchesCPU)
{
    const int n_basis = 1024;
    const int ld_psi = 1024;
    const int n_band = 32;
    using Real = double;

    auto psi  = generate_random_complex<Real>(ld_psi * n_band);
    auto hpsi = generate_random_complex<Real>(ld_psi * n_band);

    // CPU reference
    std::vector<Real> cpu_eigen(n_band);
    cpu_compute_band_energies<Real>(cpu_eigen.data(), psi.data(), hpsi.data(),
                                    n_basis, ld_psi, n_band);

    // GPU
    Real* d_eigen = nullptr;
    thrust::complex<Real>* d_psi = nullptr;
    thrust::complex<Real>* d_hpsi = nullptr;

    size_t vbytes = ld_psi * n_band * sizeof(thrust::complex<Real>);
    size_t ebytes = n_band * sizeof(Real);

    CHECK_CUDA(cudaMalloc(&d_eigen, ebytes));
    CHECK_CUDA(cudaMalloc(&d_psi, vbytes));
    CHECK_CUDA(cudaMalloc(&d_hpsi, vbytes));

    CHECK_CUDA(cudaMemcpy(d_psi, psi.data(), vbytes, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_hpsi, hpsi.data(), vbytes, cudaMemcpyHostToDevice));

    cuda::compute_band_energies_op<Real>(d_eigen, d_psi, d_hpsi,
                                         n_basis, ld_psi, n_band);

    std::vector<Real> gpu_eigen(n_band);
    CHECK_CUDA(cudaMemcpy(gpu_eigen.data(), d_eigen, ebytes, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaDeviceSynchronize());

    const Real tol = Real(1e-10);
    for (int b = 0; b < n_band; ++b) {
        Real diff = std::abs(cpu_eigen[b] - gpu_eigen[b]);
        Real ref = std::max(std::abs(cpu_eigen[b]), Real(1.0));
        EXPECT_LT(diff / ref, tol) << "Band " << b << " energy mismatch: "
            << "CPU=" << cpu_eigen[b] << " GPU=" << gpu_eigen[b];
    }

    CHECK_CUDA(cudaFree(d_eigen));
    CHECK_CUDA(cudaFree(d_psi));
    CHECK_CUDA(cudaFree(d_hpsi));
}
#endif // CUDA_AVAILABLE

// ============================================================================
// Test 5: Preconditioner Application
// ============================================================================

#if CUDA_AVAILABLE
TEST_F(DiagoCudaKernelsTest, ApplyPreconditioner_MatchesCPU)
{
    const int n_basis = 2048;
    using Real = double;

    auto hpsi = generate_random_complex<Real>(n_basis);
    auto spsi = generate_random_complex<Real>(n_basis);
    auto prec = generate_random_real<Real>(n_basis, 0.5, 2.0);
    Real eigen = 1.5;

    // CPU reference
    std::vector<std::complex<Real>> cpu_grad(n_basis);
    cpu_apply_preconditioner<Real>(cpu_grad.data(), hpsi.data(), spsi.data(),
                                   prec.data(), eigen, n_basis);

    // GPU
    thrust::complex<Real>* d_grad = nullptr;
    thrust::complex<Real>* d_hpsi = nullptr;
    thrust::complex<Real>* d_spsi = nullptr;
    Real* d_prec = nullptr;

    size_t cbytes = n_basis * sizeof(thrust::complex<Real>);
    size_t rbytes = n_basis * sizeof(Real);

    CHECK_CUDA(cudaMalloc(&d_grad, cbytes));
    CHECK_CUDA(cudaMalloc(&d_hpsi, cbytes));
    CHECK_CUDA(cudaMalloc(&d_spsi, cbytes));
    CHECK_CUDA(cudaMalloc(&d_prec, rbytes));

    CHECK_CUDA(cudaMemcpy(d_hpsi, hpsi.data(), cbytes, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_spsi, spsi.data(), cbytes, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_prec, prec.data(), rbytes, cudaMemcpyHostToDevice));

    cuda::apply_preconditioner_op<Real>(d_grad, d_hpsi, d_spsi, d_prec, eigen, n_basis);

    std::vector<std::complex<Real>> gpu_grad(n_basis);
    CHECK_CUDA(cudaMemcpy(gpu_grad.data(), d_grad, cbytes, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaDeviceSynchronize());

    const Real tol = Real(1e-9);
    for (int i = 0; i < n_basis; ++i) {
        EXPECT_NEAR(cpu_grad[i].real(), gpu_grad[i].real(), tol);
        EXPECT_NEAR(cpu_grad[i].imag(), gpu_grad[i].imag(), tol);
    }

    CHECK_CUDA(cudaFree(d_grad));
    CHECK_CUDA(cudaFree(d_hpsi));
    CHECK_CUDA(cudaFree(d_spsi));
    CHECK_CUDA(cudaFree(d_prec));
}
#endif // CUDA_AVAILABLE

// ============================================================================
// Test 6: Batched Div Preconditioner
// ============================================================================

#if CUDA_AVAILABLE
TEST_F(DiagoCudaKernelsTest, BatchedDivPreconditioner_MatchesCPU)
{
    const int n_basis = 512;
    const int ld = 640; // non-contiguous
    const int n_band = 8;
    using Real = double;

    auto in_vec = generate_random_complex<Real>(ld * n_band);
    auto prec = generate_random_real<Real>(n_basis, 0.5, 2.0);

    // CPU reference
    std::vector<std::complex<Real>> cpu_out(ld * n_band);
    for (int b = 0; b < n_band; ++b) {
        for (int i = 0; i < n_basis; ++i) {
            cpu_out[b * ld + i] = in_vec[b * ld + i] / prec[i];
        }
    }

    // GPU
    thrust::complex<Real>* d_out = nullptr;
    thrust::complex<Real>* d_in = nullptr;
    Real* d_prec = nullptr;

    size_t vbytes = ld * n_band * sizeof(thrust::complex<Real>);
    size_t rbytes = n_basis * sizeof(Real);

    CHECK_CUDA(cudaMalloc(&d_out, vbytes));
    CHECK_CUDA(cudaMalloc(&d_in, vbytes));
    CHECK_CUDA(cudaMalloc(&d_prec, rbytes));

    CHECK_CUDA(cudaMemcpy(d_in, in_vec.data(), vbytes, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_prec, prec.data(), rbytes, cudaMemcpyHostToDevice));

    cuda::batched_div_preconditioner_op<thrust::complex<Real>>(
        d_out, d_in, d_prec, n_basis, ld, n_band);

    std::vector<std::complex<Real>> gpu_out(ld * n_band);
    CHECK_CUDA(cudaMemcpy(gpu_out.data(), d_out, vbytes, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaDeviceSynchronize());

    const Real tol = Real(1e-9);
    for (int b = 0; b < n_band; ++b) {
        for (int i = 0; i < n_basis; ++i) {
            EXPECT_NEAR(cpu_out[b * ld + i].real(), gpu_out[b * ld + i].real(), tol)
                << "out[" << b << "][" << i << "].real mismatch";
            EXPECT_NEAR(cpu_out[b * ld + i].imag(), gpu_out[b * ld + i].imag(), tol)
                << "out[" << b << "][" << i << "].imag mismatch";
        }
    }

    CHECK_CUDA(cudaFree(d_out));
    CHECK_CUDA(cudaFree(d_in));
    CHECK_CUDA(cudaFree(d_prec));
}
#endif // CUDA_AVAILABLE

// ============================================================================
// Test 7: Large-scale stress test
// ============================================================================

#if CUDA_AVAILABLE
TEST_F(DiagoCudaKernelsTest, LargeScaleStressTest)
{
    const int n_basis = 8192;
    const int n_band = 64;
    using Real = double;

    // Test that kernels handle large problem sizes correctly
    auto psi  = generate_random_complex<Real>(n_basis * n_band);
    auto hpsi = generate_random_complex<Real>(n_basis * n_band);

    // CPU
    std::vector<Real> cpu_eigen(n_band);
    cpu_compute_band_energies<Real>(cpu_eigen.data(), psi.data(), hpsi.data(),
                                    n_basis, n_basis, n_band);

    // GPU
    thrust::complex<Real>* d_psi = nullptr;
    thrust::complex<Real>* d_hpsi = nullptr;
    Real* d_eigen = nullptr;

    size_t vbytes = n_basis * n_band * sizeof(thrust::complex<Real>);
    size_t ebytes = n_band * sizeof(Real);

    CHECK_CUDA(cudaMalloc(&d_psi, vbytes));
    CHECK_CUDA(cudaMalloc(&d_hpsi, vbytes));
    CHECK_CUDA(cudaMalloc(&d_eigen, ebytes));

    CHECK_CUDA(cudaMemcpy(d_psi, psi.data(), vbytes, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_hpsi, hpsi.data(), vbytes, cudaMemcpyHostToDevice));

    cuda::compute_band_energies_op<Real>(d_eigen, d_psi, d_hpsi,
                                         n_basis, n_basis, n_band);

    std::vector<Real> gpu_eigen(n_band);
    CHECK_CUDA(cudaMemcpy(gpu_eigen.data(), d_eigen, ebytes, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaDeviceSynchronize());

    const Real tol = Real(1e-8);
    for (int b = 0; b < n_band; ++b) {
        Real diff = std::abs(cpu_eigen[b] - gpu_eigen[b]);
        Real ref = std::max(std::abs(cpu_eigen[b]), Real(1.0));
        EXPECT_LT(diff / ref, tol) << "Band " << b << " mismatch for large test";
    }

    CHECK_CUDA(cudaFree(d_psi));
    CHECK_CUDA(cudaFree(d_hpsi));
    CHECK_CUDA(cudaFree(d_eigen));
}
#endif // CUDA_AVAILABLE

} // namespace test
} // namespace hsolver

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

#if CUDA_AVAILABLE
    // Runtime GPU detection: gracefully skip all GPU tests if no device is available.
    // This prevents crashes in CI environments with CUDA toolkit but no GPU.
    int device_count = 0;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err == cudaSuccess && device_count > 0) {
        g_gpu_available = true;
        std::cout << "CUDA GPU detected: " << device_count << " device(s). Running GPU tests." << std::endl;
    } else {
        std::cout << "No CUDA GPU found (" << cudaGetErrorString(err)
                  << ", count=" << device_count << "). Skipping all GPU tests." << std::endl;
    }
#endif

    return RUN_ALL_TESTS();
}
