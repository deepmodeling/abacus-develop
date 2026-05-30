/**
 * diago_ppcg_test.cpp — unit test for DiagoPPCG solver
 *
 * Solves the 1D particle-in-a-box problem (tridiagonal H) with S = I,
 * and compares computed eigenvalues against exact analytic values.
 * Both BLOCK_SUBSPACE and CONJUGATE_GRADIENT strategies are tested.
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

// -----------------------------------------------------------------------------
// Test fixture: 1D particle-in-a-box (tridiagonal Laplacian)
// -----------------------------------------------------------------------------
class DiagoPPCGTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        n_dim = 10;
        nband = 3;
        ld = n_dim;

        // Build tridiagonal H:  H[i,i] = 2, H[i,i±1] = -1
        // Exact λ_k = 2 - 2·cos(k·π / (n_dim+1)), k = 1, 2, ...
        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i) {
            H_mat[i + i * n_dim] = T(2.0, 0);
            if (i > 0)         H_mat[i + (i - 1) * n_dim] = T(-1.0, 0);
            if (i < n_dim - 1) H_mat[i + (i + 1) * n_dim] = T(-1.0, 0);
        }

        // Preconditioner — diagonal of H (all 2.0)
        prec.assign(n_dim, 2.0);

        // Exact reference eigenvalues (lowest nband)
        exact.resize(nband);
        for (int k = 0; k < nband; ++k)
            exact[k] = 2.0 - 2.0 * std::cos(static_cast<Real>(k + 1)
                                            * M_PI / static_cast<Real>(n_dim + 1));

        // Convergence thresholds
        ethr.assign(nband, 1e-10);

        // Generate initial guess wavefunctions (fixed seed for reproducibility)
        std::mt19937 rng(42);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
            for (int i = 0; i < n_dim; ++i)
                psi[i + j * ld] = T(dist(rng), dist(rng));

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

// -----------------------------------------------------------------------------
// Test BLOCK_SUBSPACE strategy
// -----------------------------------------------------------------------------
TEST_F(DiagoPPCGTest, BlockSubspaceStrategy)
{
    std::vector<T>    psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 100,
        /* sbsize   = */ 4,
        /* rr_step  = */ 4,
        /* gamma_g0 = */ false,
        hsolver::PpcgStrategy::BLOCK_SUBSPACE
    );

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(
        h_op,
        /* spsi_func = */ nullptr,   // S = I
        ld, nband, n_dim,
        psi_run.data(),
        eval.data(),
        ethr,
        prec.data()
    );

    // Check eigenvalues against exact solution
    for (int i = 0; i < nband; ++i) {
        EXPECT_NEAR(eval[i], exact[i], 1e-8)
            << "BLOCK_SUBSPACE: eigenvalue[" << i << "] mismatch";
    }

    // Should converge within reasonable iterations
    EXPECT_LE(avg_iter, static_cast<double>(100))
        << "BLOCK_SUBSPACE: too many iterations";
}

// -----------------------------------------------------------------------------
// Test CONJUGATE_GRADIENT strategy
// -----------------------------------------------------------------------------
TEST_F(DiagoPPCGTest, ConjugateGradientStrategy)
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
        h_op,
        /* spsi_func = */ nullptr,   // S = I
        ld, nband, n_dim,
        psi_run.data(),
        eval.data(),
        ethr,
        prec.data()
    );

    // Check eigenvalues against exact solution
    for (int i = 0; i < nband; ++i) {
        EXPECT_NEAR(eval[i], exact[i], 1e-8)
            << "CONJUGATE_GRADIENT: eigenvalue[" << i << "] mismatch";
    }

    // Should converge within reasonable iterations
    EXPECT_LE(avg_iter, static_cast<double>(100))
        << "CONJUGATE_GRADIENT: too many iterations";
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
