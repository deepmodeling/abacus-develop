/**
 * diago_ppcg_float_test.cpp — single-precision unit test for DiagoPPCG.
 *
 * Exercises the std::complex<float> instantiation of the BLOCK_SUBSPACE and
 * CONJUGATE_GRADIENT strategies on dense matrices with analytical reference
 * eigenvalues.  Tolerances are looser than the double-precision suite because
 * single precision has roughly 7 significant digits.
 */

#include "../diago_ppcg.h"

#include <gtest/gtest.h>
#include <cmath>
#include <complex>
#include <limits>
#include <random>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using T    = std::complex<float>;
using Real = float;

// -----------------------------------------------------------------------------
// Helper: dense H-matrix times a set of column vectors (column-major H).
// -----------------------------------------------------------------------------
static void dense_h_multiply(const T* H_data, int n_dim,
                             const T* in, T* out, int ld, int ncol)
{
    for (int j = 0; j < ncol; ++j)
    {
        for (int i = 0; i < n_dim; ++i)
        {
            T sum = T(0.0f, 0.0f);
            for (int k = 0; k < n_dim; ++k)
            {
                sum += H_data[i + k * n_dim] * in[k + j * ld];
            }
            out[i + j * ld] = sum;
        }
    }
}

// Orthonormalize columns of psi in-place (S = I).
static void gram_schmidt(std::vector<T>& psi, int ld, int n_dim, int nband)
{
    for (int j = 0; j < nband; ++j)
    {
        for (int k = 0; k < j; ++k)
        {
            T dot = T(0.0f, 0.0f);
            for (int i = 0; i < n_dim; ++i)
            {
                dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
            }
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] -= dot * psi[i + k * ld];
            }
        }
        Real nrm = 0.0f;
        for (int i = 0; i < n_dim; ++i)
        {
            nrm += std::norm(psi[i + j * ld]);
        }
        nrm = std::sqrt(nrm);
        for (int i = 0; i < n_dim; ++i)
        {
            psi[i + j * ld] /= nrm;
        }
    }
}

// -----------------------------------------------------------------------------
// Diagonal matrix: H = diag(1, 2, 3).  All eigenvalues are computed (nband ==
// n_dim), so there is no ambiguity about which end of the spectrum to converge
// to; single-precision Rayleigh-Ritz can otherwise drift toward the upper
// eigenvalues on some platforms.
// -----------------------------------------------------------------------------
TEST(DiagoPPCGFloatTest, DiagonalBlockSubspace)
{
    const int n_dim = 3;
    const int nband = 3;
    const int ld = n_dim;

    std::vector<T> H_mat(n_dim * n_dim, T(0.0f, 0.0f));
    for (int i = 0; i < n_dim; ++i)
    {
        H_mat[i + i * n_dim] = T(Real(i + 1), 0.0f);
    }

    std::vector<Real> prec(n_dim);
    for (int i = 0; i < n_dim; ++i)
    {
        prec[i] = Real(i + 1);
    }

    const Real exact[3] = {1.0f, 2.0f, 3.0f};
    std::vector<double> ethr(nband, 1e-4);

    std::mt19937 rng(42);
    std::uniform_real_distribution<Real> dist(-1.0f, 1.0f);
    std::vector<T> psi(ld * nband, T(0.0f, 0.0f));
    for (int j = 0; j < nband; ++j)
    {
        for (int i = 0; i < n_dim; ++i)
        {
            psi[i + j * ld] = T(dist(rng), 0.0f);
        }
    }
    gram_schmidt(psi, ld, n_dim, nband);

    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0f);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-5f,
        /* max_iter = */ 100,
        /* sbsize   = */ 3,
        /* rr_step  = */ 3,
        /* gamma_g0 = */ false,
        hsolver::PpcgStrategy::BLOCK_SUBSPACE);

    auto h_op = [&](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim,
                                  psi_run.data(), eval.data(), ethr, prec.data());

    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(double(eval[i]), double(exact[i]), 1e-4)
            << "Diagonal float BLOCK: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, 100.0) << "Diagonal float BLOCK: too many iterations";
}

// -----------------------------------------------------------------------------
// Tridiagonal Laplacian: H[i,i]=2, H[i,i±1]=-1, exact λ_k = 2 - 2cos(kπ/(n+1))
// -----------------------------------------------------------------------------
TEST(DiagoPPCGFloatTest, TridiagonalBlockSubspace)
{
    const int n_dim = 10;
    const int nband = 3;
    const int ld = n_dim;

    std::vector<T> H_mat(n_dim * n_dim, T(0.0f, 0.0f));
    for (int i = 0; i < n_dim; ++i)
    {
        H_mat[i + i * n_dim] = T(2.0f, 0.0f);
        if (i > 0)
        {
            H_mat[i + (i - 1) * n_dim] = T(-1.0f, 0.0f);
        }
        if (i < n_dim - 1)
        {
            H_mat[i + (i + 1) * n_dim] = T(-1.0f, 0.0f);
        }
    }

    std::vector<Real> prec(n_dim, 2.0f);
    std::vector<Real> exact(nband);
    for (int k = 0; k < nband; ++k)
    {
        exact[k] = 2.0f - 2.0f * std::cos(Real(k + 1) * M_PI
                                           / Real(n_dim + 1));
    }
    std::vector<double> ethr(nband, 1e-4);

    std::mt19937 rng(42);
    std::uniform_real_distribution<Real> dist(-1.0f, 1.0f);
    std::vector<T> psi(ld * nband, T(0.0f, 0.0f));
    for (int j = 0; j < nband; ++j)
    {
        for (int i = 0; i < n_dim; ++i)
        {
            psi[i + j * ld] = T(dist(rng), 0.0f);
        }
    }
    gram_schmidt(psi, ld, n_dim, nband);

    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0f);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-5f,
        /* max_iter = */ 100,
        /* sbsize   = */ 4,
        /* rr_step  = */ 4,
        /* gamma_g0 = */ false,
        hsolver::PpcgStrategy::BLOCK_SUBSPACE);

    auto h_op = [&](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim,
                                  psi_run.data(), eval.data(), ethr, prec.data());

    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(double(eval[i]), double(exact[i]), 1e-4)
            << "Tridiagonal float BLOCK: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, 100.0) << "Tridiagonal float BLOCK: too many iterations";
}

// -----------------------------------------------------------------------------
// CONJUGATE_GRADIENT fallback strategy on the diagonal matrix.
// -----------------------------------------------------------------------------
TEST(DiagoPPCGFloatTest, ConjugateGradientFallback)
{
    const int n_dim = 5;
    const int nband = 3;
    const int ld = n_dim;

    std::vector<T> H_mat(n_dim * n_dim, T(0.0f, 0.0f));
    for (int i = 0; i < n_dim; ++i)
    {
        H_mat[i + i * n_dim] = T(Real(i + 1), 0.0f);
    }

    std::vector<Real> prec(n_dim);
    for (int i = 0; i < n_dim; ++i)
    {
        prec[i] = Real(i + 1);
    }

    const Real exact[3] = {1.0f, 2.0f, 3.0f};
    std::vector<double> ethr(nband, 1e-4);

    std::mt19937 rng(42);
    std::uniform_real_distribution<Real> dist(-1.0f, 1.0f);
    std::vector<T> psi(ld * nband, T(0.0f, 0.0f));
    for (int j = 0; j < nband; ++j)
    {
        for (int i = 0; i < n_dim; ++i)
        {
            psi[i + j * ld] = T(dist(rng), 0.0f);
        }
    }
    gram_schmidt(psi, ld, n_dim, nband);

    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0f);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-5f,
        /* max_iter = */ 200,
        /* sbsize   = */ 3,
        /* rr_step  = */ 3,
        /* gamma_g0 = */ false,
        hsolver::PpcgStrategy::CONJUGATE_GRADIENT);

    auto h_op = [&](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim,
                                  psi_run.data(), eval.data(), ethr, prec.data());

    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(double(eval[i]), double(exact[i]), 1e-4)
            << "Diagonal float CG: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, 200.0) << "Diagonal float CG: too many iterations";
}

// -----------------------------------------------------------------------------
// Non-finite input validation (throws).
// -----------------------------------------------------------------------------
TEST(DiagoPPCGFloatTest, NonFiniteInputThrows)
{
    const int n_dim = 5;
    const int nband = 3;
    const int ld = n_dim;

    std::vector<T> H_mat(n_dim * n_dim, T(0.0f, 0.0f));
    for (int i = 0; i < n_dim; ++i)
    {
        H_mat[i + i * n_dim] = T(Real(i + 1), 0.0f);
    }

    std::vector<Real> prec(n_dim, 1.0f);
    std::vector<T> psi(ld * nband, T(1.0f, 0.0f));
    std::vector<Real> eval(nband, 0.0f);
    std::vector<double> ethr(nband, 1e-4);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-5f, 100, 3, 3, false,
        hsolver::PpcgStrategy::BLOCK_SUBSPACE);

    auto h_op = [&](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    std::vector<double> bad_ethr = ethr;
    bad_ethr[0] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(solver.diag(h_op, nullptr, ld, nband, n_dim,
                             psi.data(), eval.data(), bad_ethr, prec.data()),
                 std::invalid_argument);

    std::vector<Real> bad_prec = prec;
    bad_prec[0] = std::numeric_limits<Real>::quiet_NaN();
    EXPECT_THROW(solver.diag(h_op, nullptr, ld, nband, n_dim,
                             psi.data(), eval.data(), ethr, bad_prec.data()),
                 std::invalid_argument);
}
