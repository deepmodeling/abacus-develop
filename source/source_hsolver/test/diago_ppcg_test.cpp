/**
 * diago_ppcg_test.cpp — unit test for DiagoPPCG solver
 *
 * Test matrices (all with S = I):
 *   1. Tridiagonal Laplacian (1D particle-in-a-box): H[i,i]=2, H[i,i±1]=-1
 *      Exact λ_k = 2 - 2·cos(k·π/(n+1)).  Realistic but sparse.
 *   2. Diagonal matrix: H = diag(1, 2, 3, 4, 5)
 *      Exact eigenvalues are the diagonal entries.  Simplest possible
 *      smoke test — should converge in very few iterations.
 *
 * Tests use the CONJUGATE_GRADIENT strategy which has a try/catch fallback
 * for LAPACK sygvd failures and is therefore more portable across different
 * LAPACK implementations.
 *
 * BLOCK_SUBSPACE strategy tests exist in git history but are disabled here
 * because they require a LAPACK with reliable dsygvd for small ill-conditioned
 * generalized eigenvalue problems.
 */

#include "../diago_ppcg.h"

#include <gtest/gtest.h>
#include <cmath>
#include <random>
#include <vector>
#include <complex>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using T    = std::complex<double>;
using Real = double;

// -----------------------------------------------------------------------------
// Helper: dense H-matrix times a set of column vectors
// H is stored column-major: H(row, col) = H_data[row + col * n_dim]
// -----------------------------------------------------------------------------
static void dense_h_multiply(const T* H_data, int n_dim,
                              const T* in, T* out, int ld, int ncol)
{
    for (int j = 0; j < ncol; ++j) {
        for (int i = 0; i < n_dim; ++i) {
            T sum = 0;
            for (int k = 0; k < n_dim; ++k)
                sum += H_data[i + k * n_dim] * in[k + j * ld];
            out[i + j * ld] = sum;
        }
    }
}

// =============================================================================
// Test fixture: 1D particle-in-a-box (tridiagonal Laplacian)
// =============================================================================
class DiagoPPCGTridiagTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        n_dim = 10;
        nband = 3;
        ld = n_dim;

        // Build tridiagonal H:  H[i,i] = 2, H[i,i±1] = -1
        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i) {
            H_mat[i + i * n_dim] = T(2.0, 0);
            if (i > 0)         H_mat[i + (i - 1) * n_dim] = T(-1.0, 0);
            if (i < n_dim - 1) H_mat[i + (i + 1) * n_dim] = T(-1.0, 0);
        }

        // Preconditioner — diagonal of H (all 2.0)
        prec.assign(n_dim, 2.0);

        // Exact reference eigenvalues
        exact.resize(nband);
        for (int k = 0; k < nband; ++k)
            exact[k] = 2.0 - 2.0 * std::cos(static_cast<Real>(k + 1)
                                            * M_PI / static_cast<Real>(n_dim + 1));

        // Convergence thresholds
        ethr.assign(nband, 1e-10);

        // Random initial guess (fixed seed for reproducibility)
        std::mt19937 rng(42);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
            for (int i = 0; i < n_dim; ++i)
                psi[i + j * ld] = T(dist(rng), 0.0);

        // Gram-Schmidt orthonormalisation (S = I)
        for (int j = 0; j < nband; ++j) {
            for (int k = 0; k < j; ++k) {
                T dot = 0;
                for (int i = 0; i < n_dim; ++i)
                    dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                for (int i = 0; i < n_dim; ++i)
                    psi[i + j * ld] -= dot * psi[i + k * ld];
            }
            Real nrm = 0;
            for (int i = 0; i < n_dim; ++i)
                nrm += std::norm(psi[i + j * ld]);
            nrm = std::sqrt(nrm);
            for (int i = 0; i < n_dim; ++i)
                psi[i + j * ld] /= nrm;
        }
    }

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGTridiagTest, ConjugateGradient)
{
    std::vector<T>    psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 100,
        /* sbsize   = */ 4,
        /* rr_step  = */ 4,
        /* gamma_g0 = */ false,
        hsolver::PpcgStrategy::CONJUGATE_GRADIENT
    );

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(
        h_op, nullptr, ld, nband, n_dim,
        psi_run.data(), eval.data(), ethr, prec.data()
    );

    for (int i = 0; i < nband; ++i) {
        EXPECT_NEAR(eval[i], exact[i], 1e-8)
            << "Tridiag CG: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, static_cast<double>(100))
        << "Tridiag CG: too many iterations";
}

// =============================================================================
// Test fixture: diagonal matrix (simplest possible Hamiltonian)
// =============================================================================
class DiagoPPCGDiagonalTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        n_dim = 5;
        nband = 3;
        ld = n_dim;

        // Build diagonal H: H[i,i] = i+1
        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i)
            H_mat[i + i * n_dim] = T(static_cast<Real>(i + 1), 0);

        // Preconditioner — diagonal of H
        prec.resize(n_dim);
        for (int i = 0; i < n_dim; ++i)
            prec[i] = static_cast<Real>(i + 1);

        // Lowest 3 eigenvalues: 1, 2, 3
        exact = {1.0, 2.0, 3.0};

        // Convergence thresholds
        ethr.assign(nband, 1e-10);

        // Random initial guess (fixed seed)
        std::mt19937 rng(42);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
            for (int i = 0; i < n_dim; ++i)
                psi[i + j * ld] = T(dist(rng), 0.0);

        // Gram-Schmidt orthonormalisation (S = I)
        for (int j = 0; j < nband; ++j) {
            for (int k = 0; k < j; ++k) {
                T dot = 0;
                for (int i = 0; i < n_dim; ++i)
                    dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                for (int i = 0; i < n_dim; ++i)
                    psi[i + j * ld] -= dot * psi[i + k * ld];
            }
            Real nrm = 0;
            for (int i = 0; i < n_dim; ++i)
                nrm += std::norm(psi[i + j * ld]);
            nrm = std::sqrt(nrm);
            for (int i = 0; i < n_dim; ++i)
                psi[i + j * ld] /= nrm;
        }
    }

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGDiagonalTest, ConjugateGradient)
{
    std::vector<T>    psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 50,
        /* sbsize   = */ 3,
        /* rr_step  = */ 3,
        /* gamma_g0 = */ false,
        hsolver::PpcgStrategy::CONJUGATE_GRADIENT
    );

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(
        h_op, nullptr, ld, nband, n_dim,
        psi_run.data(), eval.data(), ethr, prec.data()
    );

    for (int i = 0; i < nband; ++i) {
        EXPECT_NEAR(eval[i], exact[i], 1e-8)
            << "Diagonal CG: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, static_cast<double>(50))
        << "Diagonal CG: too many iterations";
}

// =============================================================================
// Test fixture: 2×2 matrix — smallest non-trivial case
// H = [[2, 1], [1, 2]], eigenvalues: 1, 3
// =============================================================================
class DiagoPPCG2x2Test : public ::testing::Test
{
protected:
    void SetUp() override
    {
        n_dim = 2;
        nband = 2;
        ld = n_dim;

        // H = [[2, 1], [1, 2]]
        H_mat.assign(n_dim * n_dim, T(0));
        H_mat[0 + 0 * n_dim] = T(2.0, 0);
        H_mat[1 + 1 * n_dim] = T(2.0, 0);
        H_mat[0 + 1 * n_dim] = T(1.0, 0);
        H_mat[1 + 0 * n_dim] = T(1.0, 0);

        prec.assign(n_dim, 2.0);

        // λ₁ = 1, λ₂ = 3
        exact = {1.0, 3.0};

        ethr.assign(nband, 1e-10);

        std::mt19937 rng(123);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
            for (int i = 0; i < n_dim; ++i)
                psi[i + j * ld] = T(dist(rng), 0.0);

        // Gram-Schmidt orthonormalisation (S = I)
        for (int j = 0; j < nband; ++j) {
            for (int k = 0; k < j; ++k) {
                T dot = 0;
                for (int i = 0; i < n_dim; ++i)
                    dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                for (int i = 0; i < n_dim; ++i)
                    psi[i + j * ld] -= dot * psi[i + k * ld];
            }
            Real nrm = 0;
            for (int i = 0; i < n_dim; ++i)
                nrm += std::norm(psi[i + j * ld]);
            nrm = std::sqrt(nrm);
            for (int i = 0; i < n_dim; ++i)
                psi[i + j * ld] /= nrm;
        }
    }

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCG2x2Test, ConjugateGradient)
{
    std::vector<T>    psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 50,
        /* sbsize   = */ 2,
        /* rr_step  = */ 2,
        /* gamma_g0 = */ false,
        hsolver::PpcgStrategy::CONJUGATE_GRADIENT
    );

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(
        h_op, nullptr, ld, nband, n_dim,
        psi_run.data(), eval.data(), ethr, prec.data()
    );

    for (int i = 0; i < nband; ++i) {
        EXPECT_NEAR(eval[i], exact[i], 1e-8)
            << "2x2 CG: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, static_cast<double>(50))
        << "2x2 CG: too many iterations";
}

// =============================================================================
// Test fixture: degenerate eigenvalues
// H = I + J  (identity plus all-ones), 4×4.
// J has eigenvector [1,1,1,1]^T with eigenvalue 4.
// All vectors orthogonal to [1,1,1,1]^T are eigenvectors with eigenvalue 0.
// So H = I + J has: λ₁ = 1 (multiplicity 3), λ₄ = 5.
// =============================================================================
class DiagoPPCGDegenerateTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        n_dim = 4;
        nband = 4;
        ld = n_dim;

        // H = I + J where J is the all-ones matrix
        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i)
            for (int j = 0; j < n_dim; ++j)
                H_mat[i + j * n_dim] = T(1.0, 0);  // all-ones J
        for (int i = 0; i < n_dim; ++i)
            H_mat[i + i * n_dim] += T(1.0, 0);       // J → I+J

        // Preconditioner: diagonal = 2
        prec.assign(n_dim, 2.0);

        // λ = {1, 1, 1, 5}
        exact = {1.0, 1.0, 1.0, 5.0};

        ethr.assign(nband, 1e-10);

        std::mt19937 rng(456);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
            for (int i = 0; i < n_dim; ++i)
                psi[i + j * ld] = T(dist(rng), 0.0);

        // Gram-Schmidt (S = I)
        for (int j = 0; j < nband; ++j) {
            for (int k = 0; k < j; ++k) {
                T dot = 0;
                for (int i = 0; i < n_dim; ++i)
                    dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                for (int i = 0; i < n_dim; ++i)
                    psi[i + j * ld] -= dot * psi[i + k * ld];
            }
            Real nrm = 0;
            for (int i = 0; i < n_dim; ++i)
                nrm += std::norm(psi[i + j * ld]);
            nrm = std::sqrt(nrm);
            for (int i = 0; i < n_dim; ++i)
                psi[i + j * ld] /= nrm;
        }
    }

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGDegenerateTest, ConjugateGradient)
{
    std::vector<T>    psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 100,
        /* sbsize   = */ 4,
        /* rr_step  = */ 4,
        /* gamma_g0 = */ false,
        hsolver::PpcgStrategy::CONJUGATE_GRADIENT
    );

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(
        h_op, nullptr, ld, nband, n_dim,
        psi_run.data(), eval.data(), ethr, prec.data()
    );

    for (int i = 0; i < nband; ++i) {
        EXPECT_NEAR(eval[i], exact[i], 1e-8)
            << "Degenerate CG: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, static_cast<double>(100))
        << "Degenerate CG: too many iterations";
}

// =============================================================================
// Test fixture: larger tridiagonal, more bands
// 20×20 tridiagonal Laplacian, nband=5.
// =============================================================================
class DiagoPPCGLargeTridiagTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        n_dim = 20;
        nband = 5;
        ld = n_dim;

        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i) {
            H_mat[i + i * n_dim] = T(2.0, 0);
            if (i > 0)         H_mat[i + (i - 1) * n_dim] = T(-1.0, 0);
            if (i < n_dim - 1) H_mat[i + (i + 1) * n_dim] = T(-1.0, 0);
        }

        prec.assign(n_dim, 2.0);

        exact.resize(nband);
        for (int k = 0; k < nband; ++k)
            exact[k] = 2.0 - 2.0 * std::cos(static_cast<Real>(k + 1)
                                            * M_PI / static_cast<Real>(n_dim + 1));

        ethr.assign(nband, 1e-10);

        std::mt19937 rng(789);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
            for (int i = 0; i < n_dim; ++i)
                psi[i + j * ld] = T(dist(rng), 0.0);

        for (int j = 0; j < nband; ++j) {
            for (int k = 0; k < j; ++k) {
                T dot = 0;
                for (int i = 0; i < n_dim; ++i)
                    dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                for (int i = 0; i < n_dim; ++i)
                    psi[i + j * ld] -= dot * psi[i + k * ld];
            }
            Real nrm = 0;
            for (int i = 0; i < n_dim; ++i)
                nrm += std::norm(psi[i + j * ld]);
            nrm = std::sqrt(nrm);
            for (int i = 0; i < n_dim; ++i)
                psi[i + j * ld] /= nrm;
        }
    }

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGLargeTridiagTest, ConjugateGradient)
{
    std::vector<T>    psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 150,
        /* sbsize   = */ 5,
        /* rr_step  = */ 5,
        /* gamma_g0 = */ false,
        hsolver::PpcgStrategy::CONJUGATE_GRADIENT
    );

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(
        h_op, nullptr, ld, nband, n_dim,
        psi_run.data(), eval.data(), ethr, prec.data()
    );

    for (int i = 0; i < nband; ++i) {
        EXPECT_NEAR(eval[i], exact[i], 1e-8)
            << "Large Tridiag CG: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, static_cast<double>(150))
        << "Large Tridiag CG: too many iterations";
}

// =============================================================================
// Test fixture: dense matrix with known eigenvalues
// H = Q^T * D * Q where Q is a known orthogonal matrix (a Givens rotation
// repeated on different index pairs) and D is diagonal.
// For an 8×8 case: D = diag(1, 2, 3, 4, 5, 6, 7, 8), then apply several
// Givens rotations to mix all rows/cols. The exact eigenvalues remain 1..8.
//
// This addresses the "full/dense matrix" test that was originally missing.
// =============================================================================
class DiagoPPCGDenseTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        n_dim = 8;
        nband = 4;
        ld = n_dim;

        // Start with diagonal matrix
        std::vector<Real> dense(n_dim * n_dim, static_cast<Real>(0));
        for (int i = 0; i < n_dim; ++i)
            dense[i + i * n_dim] = static_cast<Real>(i + 1);

        // Apply several Givens rotations to make it dense while preserving
        // eigenvalues. Each rotation: A' = G(i,j,θ)^T * A * G(i,j,θ)
        auto apply_givens = [&](int p, int q, Real theta) {
            Real c = std::cos(theta);
            Real s = std::sin(theta);
            // Apply to columns
            for (int i = 0; i < n_dim; ++i) {
                Real aip = dense[i + p * n_dim];
                Real aiq = dense[i + q * n_dim];
                dense[i + p * n_dim] =  c * aip + s * aiq;
                dense[i + q * n_dim] = -s * aip + c * aiq;
            }
            // Apply to rows
            for (int j = 0; j < n_dim; ++j) {
                Real apj = dense[p + j * n_dim];
                Real aqj = dense[q + j * n_dim];
                dense[p + j * n_dim] = c * apj + s * aqj;
                dense[q + j * n_dim] = -s * apj + c * aqj;
            }
        };

        // Several rotations with different angles to create a genuinely
        // dense matrix (all off-diagonals become non-zero)
        std::mt19937 rng_dense(111);
        std::uniform_real_distribution<Real> angle_dist(
            static_cast<Real>(0.2), static_cast<Real>(1.3));
        for (int k = 0; k < 20; ++k) {
            int p = k % (n_dim - 1);
            int q = p + 1 + (k / (n_dim - 1)) % (n_dim - 1 - p);
            if (q >= n_dim) q = n_dim - 1;
            if (p == q) continue;
            apply_givens(p, q, angle_dist(rng_dense));
        }

        // Copy to complex H_mat
        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim * n_dim; ++i)
            H_mat[i] = T(dense[i], 0);

        // Preconditioner: use diagonal of the rotated H
        prec.resize(n_dim);
        for (int i = 0; i < n_dim; ++i)
            prec[i] = std::real(H_mat[i + i * n_dim]);

        // Lowest 4 eigenvalues: 1, 2, 3, 4
        exact = {1.0, 2.0, 3.0, 4.0};

        ethr.assign(nband, 1e-10);

        std::mt19937 rng_psi(222);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
            for (int i = 0; i < n_dim; ++i)
                psi[i + j * ld] = T(dist(rng_psi), 0.0);

        for (int j = 0; j < nband; ++j) {
            for (int k = 0; k < j; ++k) {
                T dot = 0;
                for (int i = 0; i < n_dim; ++i)
                    dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                for (int i = 0; i < n_dim; ++i)
                    psi[i + j * ld] -= dot * psi[i + k * ld];
            }
            Real nrm = 0;
            for (int i = 0; i < n_dim; ++i)
                nrm += std::norm(psi[i + j * ld]);
            nrm = std::sqrt(nrm);
            for (int i = 0; i < n_dim; ++i)
                psi[i + j * ld] /= nrm;
        }
    }

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGDenseTest, ConjugateGradient)
{
    std::vector<T>    psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 200,
        /* sbsize   = */ 4,
        /* rr_step  = */ 4,
        /* gamma_g0 = */ false,
        hsolver::PpcgStrategy::CONJUGATE_GRADIENT
    );

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(
        h_op, nullptr, ld, nband, n_dim,
        psi_run.data(), eval.data(), ethr, prec.data()
    );

    for (int i = 0; i < nband; ++i) {
        EXPECT_NEAR(eval[i], exact[i], 1e-8)
            << "Dense CG: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, static_cast<double>(200))
        << "Dense CG: too many iterations";
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
