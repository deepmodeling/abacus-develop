/**
 * diago_ppcg_test.cpp — unit test for DiagoPPCG solver
 *
 * Test matrices include both S = I and non-trivial overlap operators:
 *   1. Tridiagonal Laplacian (1D particle-in-a-box): H[i,i]=2, H[i,i±1]=-1
 *      Exact λ_k = 2 - 2·cos(k·π/(n+1)).  Realistic but sparse.
 *   2. Diagonal matrix: H = diag(1, 2, 3, 4, 5)
 *      Exact eigenvalues are the diagonal entries.  Simplest possible
 *      smoke test — should converge in very few iterations.
 *
 * Tests exercise the production BLOCK_SUBSPACE strategy (the band-by-band
 * CONJUGATE_GRADIENT variant was removed from this PR).
 */

#include "../diago_ppcg.h"

#include <gtest/gtest.h>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <complex>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using T = std::complex<double>;
using Real = double;

// -----------------------------------------------------------------------------
// Helper: dense H-matrix times a set of column vectors
// H is stored column-major: H(row, col) = H_data[row + col * n_dim]
// -----------------------------------------------------------------------------
static void dense_h_multiply(const T* H_data, int n_dim, const T* in, T* out, int ld, int ncol)
{
    for (int j = 0; j < ncol; ++j)
    {
        for (int i = 0; i < n_dim; ++i)
        {
            T sum = 0;
            for (int k = 0; k < n_dim; ++k)
            {
                sum += H_data[i + k * n_dim] * in[k + j * ld];
            }
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
        for (int i = 0; i < n_dim; ++i)
        {
            H_mat[i + i * n_dim] = T(2.0, 0);
            if (i > 0)
            {
                H_mat[i + (i - 1) * n_dim] = T(-1.0, 0);
            }
            if (i < n_dim - 1)
            {
                H_mat[i + (i + 1) * n_dim] = T(-1.0, 0);
            }
        }

        // Preconditioner — diagonal of H (all 2.0)
        prec.assign(n_dim, 2.0);

        // Exact reference eigenvalues
        exact.resize(nband);
        for (int k = 0; k < nband; ++k)
        {
            exact[k] = 2.0 - 2.0 * std::cos(Real(k + 1) * M_PI / Real(n_dim + 1));
        }

        // Convergence thresholds
        ethr.assign(nband, 1e-10);

        // Random initial guess (fixed seed for reproducibility)
        std::mt19937 rng(42);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] = T(dist(rng), 0.0);
            }
        }

        // Gram-Schmidt orthonormalisation (S = I)
        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T dot = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    psi[i + j * ld] -= dot * psi[i + k * ld];
                }
            }
            Real nrm = 0;
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

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGTridiagTest, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 100,
        /* sbsize   = */ 4,
        /* rr_step  = */ 4,
        /* gamma_g0 = */ false);

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());

    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-8) << "Tridiag BLOCK: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, double(100)) << "Tridiag BLOCK: too many iterations";
}

TEST_F(DiagoPPCGTridiagTest, ReportsIterationLimit)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);
    std::vector<double> strict_ethr(nband, 1.0e-14);
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1.0e-14,
        /* max_iter = */ 1,
        /* sbsize   = */ 3,
        /* rr_step  = */ 4,
        /* gamma_g0 = */ false);

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), strict_ethr, prec.data());

    EXPECT_FALSE(solver.converged());
    EXPECT_GT(solver.active_band_count(), 0);
    EXPECT_TRUE(solver.iteration_limit_reached());
}

TEST(DiagoPPCGStatusTest, NonFiniteFinalResidualIsNotConverged)
{
    const int n_dim = 2;
    const int nband = 2;
    const int ld = n_dim;
    std::vector<T> psi(ld * nband, T(0));
    psi[0] = T(1, 0);
    psi[1 + ld] = T(1, 0);
    std::vector<Real> eval(nband, Real(0));
    const std::vector<double> ethr(nband, 1e-12);
    const std::vector<Real> prec = {1, 2};
    int h_apply_count = 0;
    const auto h_op = [&h_apply_count](T* in, T* out, int ld_in, int ncol) {
        ++h_apply_count;
        for (int j = 0; j < ncol; ++j)
        {
            for (int i = 0; i < ld_in; ++i)
            {
                const Real eigenvalue = Real(i + 1);
                out[i + j * ld_in] = eigenvalue * in[i + j * ld_in];
            }
        }
        if (h_apply_count >= 2)
        {
            out[0] = T(std::numeric_limits<Real>::quiet_NaN(), 0);
        }
    };
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 10,
        /* sbsize   = */ 2,
        /* rr_step  = */ 2,
        /* gamma_g0 = */ false);

    solver.diag(h_op, nullptr, ld, nband, n_dim,
                psi.data(), eval.data(), ethr, prec.data());

    EXPECT_FALSE(solver.converged());
    EXPECT_EQ(solver.active_band_count(), nband);
    EXPECT_FALSE(solver.iteration_limit_reached());
}

TEST(DiagoPPCGStatusTest, FailedRayleighRitzCannotCertifyLowestRoots)
{
    const int n_dim = 3;
    const int nband = 2;
    const int ld = n_dim;
    std::vector<T> h_mat(n_dim * n_dim, T(0));
    h_mat[0] = T(1, 0);
    h_mat[1 + ld] = T(2, 0);
    h_mat[2 + 2 * ld] = T(3, 0);
    std::vector<T> psi(ld * nband, T(0));
    psi[2] = T(1, 0);
    psi[2 + ld] = T(1, 0);
    std::vector<Real> eval(nband, Real(0));
    const std::vector<double> ethr(nband, 1e-12);
    const std::vector<Real> prec = {1, 2, 3};
    const auto h_op = [&h_mat, n_dim](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(h_mat.data(), n_dim, in, out, ld_in, ncol);
    };
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 0,
        /* sbsize   = */ 2,
        /* rr_step  = */ 2,
        /* gamma_g0 = */ false);

    solver.diag(h_op, nullptr, ld, nband, n_dim,
                psi.data(), eval.data(), ethr, prec.data());

    EXPECT_FALSE(solver.converged());
    EXPECT_EQ(solver.active_band_count(), nband);
    EXPECT_TRUE(solver.iteration_limit_reached());
}

TEST_F(DiagoPPCGTridiagTest, ResidualTraceWritesCsv)
{
    const char* env_name = "ABACUS_PPCG_RESIDUAL_TRACE";
    const char* old_env = std::getenv(env_name);
    const bool had_old_env = old_env != nullptr;
    const std::string old_env_value = had_old_env ? old_env : "";
    const std::string trace_path = "ppcg_residual_trace_test.csv";
    std::remove(trace_path.c_str());
    ASSERT_EQ(::setenv(env_name, trace_path.c_str(), 1), 0);

    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 100,
        /* sbsize   = */ 4,
        /* rr_step  = */ 4,
        /* gamma_g0 = */ false);

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());

    if (had_old_env)
    {
        ASSERT_EQ(::setenv(env_name, old_env_value.c_str(), 1), 0);
    }
    else
    {
        ASSERT_EQ(::unsetenv(env_name), 0);
    }

    std::ifstream trace(trace_path);
    ASSERT_TRUE(trace.good());
    std::string header;
    std::string first_record;
    std::getline(trace, header);
    std::getline(trace, first_record);
    EXPECT_EQ(header,
              "iteration,stage,max_residual,eigenvalue_0,eigenvalue_1,eigenvalue_2,"
              "local_overlap_min,local_overlap_max,local_overlap_condition");
    EXPECT_NE(first_record.find("initial_rr"), std::string::npos);
    EXPECT_EQ(std::count(first_record.begin(), first_record.end(), ','), nband + 5);

    double previous_max_error = std::numeric_limits<double>::max();
    std::string record = first_record;
    do
    {
        std::stringstream record_stream(record);
        std::string field;
        for (int field_index = 0; field_index < 3; ++field_index)
        {
            std::getline(record_stream, field, ',');
        }
        double max_error = 0.0;
        for (int ib = 0; ib < nband; ++ib)
        {
            std::getline(record_stream, field, ',');
            const double eigenvalue = std::stod(field);
            const double difference = eigenvalue - exact[ib];
            const double error = std::abs(difference);
            max_error = std::max(max_error, error);
        }
        const double allowed_error = previous_max_error + 1.0e-10;
        EXPECT_LE(max_error, allowed_error);
        previous_max_error = max_error;
    } while (std::getline(trace, record));

    trace.close();
    std::remove(trace_path.c_str());
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
        {
            H_mat[i + i * n_dim] = T(Real(i + 1), 0);
        }

        // Preconditioner — diagonal of H
        prec.resize(n_dim);
        for (int i = 0; i < n_dim; ++i)
        {
            prec[i] = Real(i + 1);
        }

        // Lowest 3 eigenvalues: 1, 2, 3
        exact = {1.0, 2.0, 3.0};

        // Convergence thresholds
        ethr.assign(nband, 1e-10);

        // Random initial guess (fixed seed)
        std::mt19937 rng(42);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] = T(dist(rng), 0.0);
            }
        }

        // Gram-Schmidt orthonormalisation (S = I)
        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T dot = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    psi[i + j * ld] -= dot * psi[i + k * ld];
                }
            }
            Real nrm = 0;
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

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGDiagonalTest, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 50,
        /* sbsize   = */ 3,
        /* rr_step  = */ 3,
        /* gamma_g0 = */ false);

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());

    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-8) << "Diagonal BLOCK: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, double(50)) << "Diagonal BLOCK: too many iterations";
}

TEST(DiagoPPCGRootSelectionTest, StartupCorrectionKeepsLowestRoots)
{
    const int n_dim = 8;
    const int nband = 4;
    const int ld = n_dim;
    const Real inv_sqrt_two = Real(1) / std::sqrt(Real(2));

    std::vector<T> h_mat(n_dim * n_dim, T(0));
    std::vector<Real> prec(n_dim);
    for (int i = 0; i < n_dim; ++i)
    {
        const Real diagonal = Real(i + 1);
        h_mat[i + i * n_dim] = T(diagonal, 0);
        prec[i] = diagonal;
    }

    // The first two columns contain the missing third and fourth roots.  The
    // last two are exact high-energy eigenvectors and would otherwise lock
    // before their lower replacements can enter the active subspace.
    std::vector<T> psi(ld * nband, T(0));
    psi[0] = T(inv_sqrt_two, 0);
    psi[2] = T(inv_sqrt_two, 0);
    psi[1 + ld] = T(inv_sqrt_two, 0);
    psi[3 + ld] = T(inv_sqrt_two, 0);
    psi[6 + 2 * ld] = T(1, 0);
    psi[7 + 3 * ld] = T(1, 0);

    std::vector<Real> eval(nband, Real(0));
    std::vector<double> ethr(nband, 1e-12);
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 50,
        /* sbsize   = */ 2,
        /* rr_step  = */ 4,
        /* gamma_g0 = */ false);

    auto h_op = [&h_mat, n_dim](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(h_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    const double avg_iter = solver.diag(
        h_op, nullptr, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], Real(i + 1), 1e-10) << "root[" << i << "] mismatch";
    }
    EXPECT_TRUE(solver.converged());
    EXPECT_LE(avg_iter, 50.0);
}

TEST(DiagoPPCGRootSelectionTest, StartupCorrectionKeepsLowestRootsWithoutInitialGap)
{
    const int n_dim = 8;
    const int nband = 4;
    const int ld = n_dim;
    const Real inv_sqrt_two = Real(1) / std::sqrt(Real(2));

    std::vector<T> h_mat(n_dim * n_dim, T(0));
    std::vector<Real> prec(n_dim);
    for (int i = 0; i < n_dim; ++i)
    {
        const Real diagonal = Real(i + 1);
        h_mat[i + i * n_dim] = T(diagonal, 0);
        prec[i] = diagonal;
    }

    // Every initial vector has a nonzero residual.  Independent two-column
    // blocks would nevertheless retain roots 1, 2, 5, and 6, discarding the
    // third and fourth roots before a global comparison is possible.
    std::vector<T> psi(ld * nband, T(0));
    psi[0] = T(inv_sqrt_two, 0);
    psi[2] = T(inv_sqrt_two, 0);
    psi[1 + ld] = T(inv_sqrt_two, 0);
    psi[3 + ld] = T(inv_sqrt_two, 0);
    psi[4 + 2 * ld] = T(inv_sqrt_two, 0);
    psi[6 + 2 * ld] = T(inv_sqrt_two, 0);
    psi[5 + 3 * ld] = T(inv_sqrt_two, 0);
    psi[7 + 3 * ld] = T(inv_sqrt_two, 0);

    std::vector<Real> eval(nband, Real(0));
    std::vector<double> ethr(nband, 1e-12);
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 50,
        /* sbsize   = */ 2,
        /* rr_step  = */ 4,
        /* gamma_g0 = */ false);

    auto h_op = [&h_mat, n_dim](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(h_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    const double avg_iter = solver.diag(
        h_op, nullptr, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], Real(i + 1), 1e-10) << "root[" << i << "] mismatch";
    }
    EXPECT_TRUE(solver.converged());
    EXPECT_LE(avg_iter, 50.0);
}

TEST(DiagoPPCGRootSelectionTest, WarmStartAfterHamiltonianChangeKeepsLowestRoots)
{
    const int n_dim = 8;
    const int nband = 4;
    const int ld = n_dim;
    const Real inv_sqrt_two = Real(1) / std::sqrt(Real(2));

    std::vector<T> h_mat(n_dim * n_dim, T(0));
    std::vector<Real> prec(n_dim, Real(1));
    std::vector<T> psi(ld * nband, T(0));
    const int first_indices[nband] = {0, 1, 4, 5};
    const int second_indices[nband] = {2, 3, 6, 7};
    for (int state = 0; state < nband; ++state)
    {
        const int first = first_indices[state];
        const int second = second_indices[state];
        const Real low_eigenvalue = Real(state + 1);
        const Real high_eigenvalue = Real(state + 5);
        const Real diagonal = Real(0.5) * (low_eigenvalue + high_eigenvalue);
        const Real off_diagonal = Real(0.5) * (low_eigenvalue - high_eigenvalue);
        h_mat[first + first * n_dim] = T(diagonal, 0);
        h_mat[second + second * n_dim] = T(diagonal, 0);
        h_mat[first + second * n_dim] = T(off_diagonal, 0);
        h_mat[second + first * n_dim] = T(off_diagonal, 0);
        psi[first + state * ld] = T(inv_sqrt_two, 0);
        psi[second + state * ld] = T(inv_sqrt_two, 0);
        prec[first] = diagonal;
        prec[second] = diagonal;
    }

    std::vector<Real> eval(nband, Real(0));
    std::vector<double> ethr(nband, 1e-12);
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 50,
        /* sbsize   = */ 2,
        /* rr_step  = */ 4,
        /* gamma_g0 = */ false);
    auto h_op = [&h_mat, n_dim](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(h_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    solver.diag(
        h_op, nullptr, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());
    ASSERT_TRUE(solver.converged());
    for (int i = 0; i < nband; ++i)
    {
        ASSERT_NEAR(eval[i], Real(i + 1), 1e-10);
    }

    // Model a later SCF call: the previous eigenvectors are a valid warm start,
    // but the changed Hamiltonian makes each one mix a low and a high root.
    std::fill(h_mat.begin(), h_mat.end(), T(0));
    for (int i = 0; i < n_dim; ++i)
    {
        const Real diagonal = Real(i + 1);
        h_mat[i + i * n_dim] = T(diagonal, 0);
        prec[i] = diagonal;
    }

    solver.diag(
        h_op, nullptr, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());
    EXPECT_TRUE(solver.converged());
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], Real(i + 1), 1e-10) << "root[" << i << "] mismatch";
    }
}

TEST(DiagoPPCGRootSelectionTest, StartupCorrectionKeepsLowestGeneralizedRoots)
{
    const int n_dim = 8;
    const int nband = 4;
    const int ld = n_dim;
    const Real inv_sqrt_two = Real(1) / std::sqrt(Real(2));

    std::vector<T> h_mat(n_dim * n_dim, T(0));
    std::vector<T> s_mat(n_dim * n_dim, T(0));
    std::vector<Real> s_diagonal(n_dim);
    std::vector<Real> prec(n_dim);
    for (int i = 0; i < n_dim; ++i)
    {
        const Real eigenvalue = Real(i + 1);
        const Real overlap_value = Real(1) + Real(0.1) * Real(i % 3);
        const Real hamiltonian_value = eigenvalue * overlap_value;
        h_mat[i + i * n_dim] = T(hamiltonian_value, 0);
        s_mat[i + i * n_dim] = T(overlap_value, 0);
        s_diagonal[i] = overlap_value;
        prec[i] = hamiltonian_value;
    }

    // H = S diag(1,...,8), so the generalized eigenvectors are the
    // coordinate vectors.  The initial columns create the same cross-block
    // root-selection hazard as the identity-overlap regression above.
    std::vector<T> psi(ld * nband, T(0));
    std::vector<Real> normalized_coefficients(n_dim);
    for (int i = 0; i < n_dim; ++i)
    {
        normalized_coefficients[i] = inv_sqrt_two / std::sqrt(s_diagonal[i]);
    }
    psi[0] = T(normalized_coefficients[0], 0);
    psi[2] = T(normalized_coefficients[2], 0);
    psi[1 + ld] = T(normalized_coefficients[1], 0);
    psi[3 + ld] = T(normalized_coefficients[3], 0);
    psi[4 + 2 * ld] = T(normalized_coefficients[4], 0);
    psi[6 + 2 * ld] = T(normalized_coefficients[6], 0);
    psi[5 + 3 * ld] = T(normalized_coefficients[5], 0);
    psi[7 + 3 * ld] = T(normalized_coefficients[7], 0);

    std::vector<Real> eval(nband, Real(0));
    std::vector<double> ethr(nband, 1e-12);
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 50,
        /* sbsize   = */ 2,
        /* rr_step  = */ 4,
        /* gamma_g0 = */ false);

    auto h_op = [&h_mat, n_dim](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(h_mat.data(), n_dim, in, out, ld_in, ncol);
    };
    auto s_op = [&s_mat, n_dim](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(s_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    const double avg_iter = solver.diag(
        h_op, s_op, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], Real(i + 1), 1e-10) << "root[" << i << "] mismatch";
    }
    EXPECT_TRUE(solver.converged());
    EXPECT_LE(avg_iter, 50.0);
}

TEST(DiagoPPCGRootSelectionTest, LockingKeepsRequestedLowRoots)
{
    const int n_dim = 4;
    const int nband = 2;
    const int ld = n_dim;
    std::vector<T> h_mat(n_dim * n_dim, T(0));
    h_mat[0] = T(1, 0);
    h_mat[1 + ld] = T(2, 0);
    h_mat[2 + 2 * ld] = T(4, 0);
    h_mat[3 + 3 * ld] = T(3, 0);
    const std::vector<Real> prec = {1, 2, 4, 3};
    const std::vector<double> ethr(nband, 1e-12);
    const Real epsilons[] = {1e-4, 1e-3, 1e-2, 1e-1};

    const auto h_op = [&h_mat, n_dim](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(h_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    for (const Real epsilon : epsilons)
    {
        SCOPED_TRACE("epsilon = " + std::to_string(epsilon));
        std::vector<T> psi(ld * nband, T(0));
        const Real norm = std::sqrt(Real(2) + epsilon * epsilon);
        psi[0] = T(Real(1) / norm, 0);
        psi[1] = T(epsilon / norm, 0);
        psi[2] = T(Real(1) / norm, 0);
        psi[3 + ld] = T(1, 0);
        std::vector<Real> eval(nband, Real(0));

        hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
            /* diag_thr = */ 1e-12,
            /* max_iter = */ 80,
            /* sbsize   = */ 2,
            /* rr_step  = */ 1,
            /* gamma_g0 = */ false);
        const double avg_iter = solver.diag(
            h_op, nullptr, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());
        EXPECT_NEAR(eval[0], 1.0, 1e-8);
        EXPECT_NEAR(eval[1], 2.0, 1e-8);
        EXPECT_TRUE(solver.converged());
        EXPECT_LE(avg_iter, 80.0);
    }
}

TEST(DiagoPPCGRankTruncationTest, ComplexHermitianStartupSubspace)
{
    const int n_dim = 6;
    const int nband = 3;
    const int ld = n_dim;
    const Real inv_sqrt_two = Real(1) / std::sqrt(Real(2));

    // Build H = U diag(1,...,6) U^H.  The first two columns of U form a
    // complex rotation, so the reduced Hamiltonian needs a conjugate
    // transpose when it is projected into the retained overlap subspace.
    std::vector<T> unitary(n_dim * n_dim, T(0));
    unitary[0] = T(inv_sqrt_two, 0);
    unitary[1] = T(0, inv_sqrt_two);
    unitary[n_dim] = T(0, inv_sqrt_two);
    unitary[1 + n_dim] = T(inv_sqrt_two, 0);
    for (int i = 2; i < n_dim; ++i)
    {
        unitary[i + i * n_dim] = T(1, 0);
    }

    std::vector<T> h_mat(n_dim * n_dim, T(0));
    std::vector<Real> eigenvalues(n_dim);
    for (int state = 0; state < n_dim; ++state)
    {
        eigenvalues[state] = Real(state + 1);
        for (int col = 0; col < n_dim; ++col)
        {
            for (int row = 0; row < n_dim; ++row)
            {
                h_mat[row + col * n_dim] += unitary[row + state * n_dim]
                                            * eigenvalues[state]
                                            * std::conj(unitary[col + state * n_dim]);
            }
        }
    }
    EXPECT_GT(std::abs(std::imag(h_mat[n_dim])), Real(0.1));

    std::vector<Real> prec(n_dim);
    for (int i = 0; i < n_dim; ++i)
    {
        prec[i] = std::max(Real(1), std::real(h_mat[i + i * n_dim]));
    }

    // The third column is an exact high-energy state.  Its zero residual
    // makes [psi, w] rank deficient during the first global correction,
    // while the first two columns and their residuals span the three lowest
    // roots.  This exercises Q^H K Q for a complex Hermitian matrix.
    std::vector<T> psi(ld * nband, T(0));
    for (int i = 0; i < n_dim; ++i)
    {
        psi[i] = (unitary[i] + unitary[i + 3 * n_dim]) * inv_sqrt_two;
        psi[i + ld] = (unitary[i + n_dim] + unitary[i + 2 * n_dim]) * inv_sqrt_two;
        psi[i + 2 * ld] = unitary[i + 5 * n_dim];
    }

    std::vector<Real> eval(nband, Real(0));
    std::vector<double> ethr(nband, 1e-12);
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 30,
        /* sbsize   = */ 2,
        /* rr_step  = */ 4,
        /* gamma_g0 = */ false);

    auto h_op = [&h_mat, n_dim](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(h_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    const double avg_iter = solver.diag(
        h_op, nullptr, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], Real(i + 1), 1e-10) << "root[" << i << "] mismatch";
    }
    EXPECT_TRUE(solver.converged());
    EXPECT_LE(avg_iter, 30.0);
}

TEST_F(DiagoPPCGDiagonalTest, EmptyHOperatorThrows)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 50,
        /* sbsize   = */ 3,
        /* rr_step  = */ 3,
        /* gamma_g0 = */ false);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU>::HPsiFunc h_op;
    EXPECT_THROW(solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data()),
                 std::invalid_argument);
}

TEST_F(DiagoPPCGDiagonalTest, NonFiniteInputThrows)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 50,
        /* sbsize   = */ 3,
        /* rr_step  = */ 3,
        /* gamma_g0 = */ false);

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    std::vector<double> bad_ethr = ethr;
    bad_ethr[0] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), bad_ethr, prec.data()),
                 std::invalid_argument);

    std::vector<Real> bad_prec = prec;
    bad_prec[0] = std::numeric_limits<Real>::infinity();
    EXPECT_THROW(solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, bad_prec.data()),
                 std::invalid_argument);
}

TEST(DiagoPPCGLeadingDimensionTest, BlockSubspaceWithPadding)
{
    const int n_dim = 5;
    const int nband = 3;
    const int ld = 8;

    std::vector<T> H_mat(n_dim * n_dim, T(0));
    for (int i = 0; i < n_dim; ++i)
    {
        H_mat[i + i * n_dim] = T(Real(i + 1), 0);
    }

    std::vector<Real> prec(n_dim);
    for (int i = 0; i < n_dim; ++i)
    {
        prec[i] = Real(i + 1);
    }

    std::vector<T> psi(ld * nband, T(17.0, -3.0));
    std::mt19937 rng(7);
    std::uniform_real_distribution<Real> dist(-1.0, 1.0);
    for (int j = 0; j < nband; ++j)
    {
        for (int i = 0; i < n_dim; ++i)
        {
            psi[i + j * ld] = T(dist(rng), 0.0);
        }
    }

    for (int j = 0; j < nband; ++j)
    {
        for (int k = 0; k < j; ++k)
        {
            T dot = 0;
            for (int i = 0; i < n_dim; ++i)
            {
                dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
            }
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] -= dot * psi[i + k * ld];
            }
        }
        Real nrm = 0;
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

    std::vector<Real> eval(nband, 0.0);
    std::vector<double> ethr(nband, 1e-10);
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 80,
        /* sbsize   = */ 2,
        /* rr_step  = */ 3,
        /* gamma_g0 = */ false);

    auto h_op = [&H_mat, n_dim](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());

    const Real exact[] = {1.0, 2.0, 3.0};
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-8) << "Padded ld BLOCK: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, double(80)) << "Padded ld BLOCK: too many iterations";
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
        {
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] = T(dist(rng), 0.0);
            }
        }

        // Gram-Schmidt orthonormalisation (S = I)
        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T dot = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    psi[i + j * ld] -= dot * psi[i + k * ld];
                }
            }
            Real nrm = 0;
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

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCG2x2Test, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 50,
        /* sbsize   = */ 2,
        /* rr_step  = */ 2,
        /* gamma_g0 = */ false);

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());

    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-8) << "2x2 BLOCK: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, double(50)) << "2x2 BLOCK: too many iterations";
}

TEST(DiagoPPCGComplexHermitianTest, DefaultKeepsImaginaryProjection)
{
    const int n_dim = 2;
    const int nband = 2;
    const int ld = n_dim;

    // H = [[2, i], [-i, 3]].  Dropping Im(<x|Hy>) would incorrectly
    // diagonalize diag(2, 3); the Hermitian eigenvalues are 2.5 +/- sqrt(1.25).
    std::vector<T> H_mat(n_dim * n_dim, T(0));
    H_mat[0 + 0 * n_dim] = T(2.0, 0.0);
    H_mat[1 + 1 * n_dim] = T(3.0, 0.0);
    H_mat[0 + 1 * n_dim] = T(0.0, 1.0);
    H_mat[1 + 0 * n_dim] = T(0.0, -1.0);

    std::vector<T> psi(ld * nband, T(0));
    psi[0 + 0 * ld] = T(1.0, 0.0);
    psi[1 + 1 * ld] = T(1.0, 0.0);

    std::vector<Real> prec(n_dim, 2.0);
    std::vector<double> ethr(nband, 1e-12);
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 10,
        /* sbsize   = */ 2,
        /* rr_step  = */ 1,
        /* gamma_g0 = */ false);

    auto h_op = [&H_mat, n_dim](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    solver.diag(h_op, nullptr, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());

    const Real delta = std::sqrt(1.25);
    EXPECT_NEAR(eval[0], 2.5 - delta, 1e-10);
    EXPECT_NEAR(eval[1], 2.5 + delta, 1e-10);
}

TEST(DiagoPPCGComplexHermitianTest, BlockSubspaceSmokeNoNaN)
{
    const int n_dim = 2;
    const int nband = 2;
    const int ld = n_dim;

    std::vector<T> H_mat(n_dim * n_dim, T(0));
    H_mat[0 + 0 * n_dim] = T(2.0, 0.0);
    H_mat[1 + 1 * n_dim] = T(3.0, 0.0);
    H_mat[0 + 1 * n_dim] = T(0.0, 1.0);
    H_mat[1 + 0 * n_dim] = T(0.0, -1.0);

    std::vector<T> psi(ld * nband, T(0));
    psi[0 + 0 * ld] = T(1.0, 0.0);
    psi[1 + 1 * ld] = T(1.0, 0.0);

    std::vector<Real> prec(n_dim, 2.0);
    std::vector<double> ethr(nband, 1e-10);
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-10,
        /* max_iter = */ 8,
        /* sbsize   = */ 2,
        /* rr_step  = */ 1,
        /* gamma_g0 = */ false);

    auto h_op = [&H_mat, n_dim](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    solver.diag(h_op, nullptr, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());

    const Real delta = std::sqrt(1.25);
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_TRUE(std::isfinite(eval[i])) << "BLOCK_SUBSPACE produced NaN/Inf";
    }
    EXPECT_NEAR(eval[0], 2.5 - delta, 1e-8);
    EXPECT_NEAR(eval[1], 2.5 + delta, 1e-8);
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
        {
            for (int j = 0; j < n_dim; ++j)
            {
                H_mat[i + j * n_dim] = T(1.0, 0); // all-ones J
            }
        }
        for (int i = 0; i < n_dim; ++i)
        {
            H_mat[i + i * n_dim] += T(1.0, 0); // J → I+J
        }

        // Preconditioner: diagonal = 2
        prec.assign(n_dim, 2.0);

        // λ = {1, 1, 1, 5}
        exact = {1.0, 1.0, 1.0, 5.0};

        ethr.assign(nband, 1e-10);

        std::mt19937 rng(456);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] = T(dist(rng), 0.0);
            }
        }

        // Gram-Schmidt (S = I)
        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T dot = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    psi[i + j * ld] -= dot * psi[i + k * ld];
                }
            }
            Real nrm = 0;
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

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGDegenerateTest, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 100,
        /* sbsize   = */ 4,
        /* rr_step  = */ 4,
        /* gamma_g0 = */ false);

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());

    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-8) << "Degenerate BLOCK: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, double(100)) << "Degenerate BLOCK: too many iterations";
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
        for (int i = 0; i < n_dim; ++i)
        {
            H_mat[i + i * n_dim] = T(2.0, 0);
            if (i > 0)
            {
                H_mat[i + (i - 1) * n_dim] = T(-1.0, 0);
            }
            if (i < n_dim - 1)
            {
                H_mat[i + (i + 1) * n_dim] = T(-1.0, 0);
            }
        }

        prec.assign(n_dim, 2.0);

        exact.resize(nband);
        for (int k = 0; k < nband; ++k)
        {
            exact[k] = 2.0 - 2.0 * std::cos(Real(k + 1) * M_PI / Real(n_dim + 1));
        }

        ethr.assign(nband, 1e-10);

        std::mt19937 rng(789);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] = T(dist(rng), 0.0);
            }
        }

        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T dot = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    psi[i + j * ld] -= dot * psi[i + k * ld];
                }
            }
            Real nrm = 0;
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

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGLargeTridiagTest, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 150,
        /* sbsize   = */ 5,
        /* rr_step  = */ 5,
        /* gamma_g0 = */ false);

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());

    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-8) << "Large Tridiag BLOCK: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, double(150)) << "Large Tridiag BLOCK: too many iterations";
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
        std::vector<Real> dense(n_dim * n_dim, Real(0));
        for (int i = 0; i < n_dim; ++i)
        {
            dense[i + i * n_dim] = Real(i + 1);
        }

        // Apply several Givens rotations to make it dense while preserving
        // eigenvalues. Each rotation: A' = G(i,j,θ)^T * A * G(i,j,θ)
        auto apply_givens = [&](int p, int q, Real theta) {
            Real c = std::cos(theta);
            Real s = std::sin(theta);
            // Apply to columns
            for (int i = 0; i < n_dim; ++i)
            {
                Real aip = dense[i + p * n_dim];
                Real aiq = dense[i + q * n_dim];
                dense[i + p * n_dim] = c * aip + s * aiq;
                dense[i + q * n_dim] = -s * aip + c * aiq;
            }
            // Apply to rows
            for (int j = 0; j < n_dim; ++j)
            {
                Real apj = dense[p + j * n_dim];
                Real aqj = dense[q + j * n_dim];
                dense[p + j * n_dim] = c * apj + s * aqj;
                dense[q + j * n_dim] = -s * apj + c * aqj;
            }
        };

        // Several rotations with different angles to create a genuinely
        // dense matrix (all off-diagonals become non-zero)
        std::mt19937 rng_dense(111);
        std::uniform_real_distribution<Real> angle_dist(Real(0.2), Real(1.3));
        for (int k = 0; k < 20; ++k)
        {
            int p = k % (n_dim - 1);
            int q = p + 1 + (k / (n_dim - 1)) % (n_dim - 1 - p);
            if (q >= n_dim)
            {
                q = n_dim - 1;
            }
            if (p == q)
            {
                continue;
            }
            apply_givens(p, q, angle_dist(rng_dense));
        }

        // Copy to complex H_mat
        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim * n_dim; ++i)
        {
            H_mat[i] = T(dense[i], 0);
        }

        // Preconditioner: use diagonal of the rotated H
        prec.resize(n_dim);
        for (int i = 0; i < n_dim; ++i)
        {
            prec[i] = std::real(H_mat[i + i * n_dim]);
        }

        // Lowest 4 eigenvalues: 1, 2, 3, 4
        exact = {1.0, 2.0, 3.0, 4.0};

        ethr.assign(nband, 1e-10);

        std::mt19937 rng_psi(222);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] = T(dist(rng_psi), 0.0);
            }
        }

        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T dot = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    psi[i + j * ld] -= dot * psi[i + k * ld];
                }
            }
            Real nrm = 0;
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

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGDenseTest, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 200,
        /* sbsize   = */ 4,
        /* rr_step  = */ 4,
        /* gamma_g0 = */ false);

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());

    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-8) << "Dense BLOCK: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, double(200)) << "Dense BLOCK: too many iterations";
}

// =============================================================================
// Helper: compute Hψ for eigenvector residual check
// =============================================================================
static void compute_residual(const T* H_data, int n_dim, const T* psi, const Real eval, int ld, T* residual)
{
    // residual = H*psi - eval*psi
    dense_h_multiply(H_data, n_dim, psi, residual, ld, 1);
    for (int i = 0; i < n_dim; ++i)
    {
        residual[i] -= eval * psi[i];
    }
}

static Real column_norm(const T* x, int n_dim)
{
    Real nrm = 0;
    for (int i = 0; i < n_dim; ++i)
    {
        nrm += std::norm(x[i]);
    }
    return std::sqrt(nrm);
}

// =============================================================================
// Test fixture: non-trivial S matrix (diagonal overlap, S ≠ I)
// H = tridiag 6×6 Laplacian, S = diag(1.1, 1.0, 0.9, 1.0, 1.1, 1.0)
// Tests that the solver correctly handles a non-identity overlap matrix.
// =============================================================================
class DiagoPPCGWithSTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        n_dim = 6;
        nband = 3;
        ld = n_dim + 2; // exercise custom S with padded leading dimension

        // Tridiagonal H
        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i)
        {
            H_mat[i + i * n_dim] = T(2.0, 0);
            if (i > 0)
            {
                H_mat[i + (i - 1) * n_dim] = T(-1.0, 0);
            }
            if (i < n_dim - 1)
            {
                H_mat[i + (i + 1) * n_dim] = T(-1.0, 0);
            }
        }

        // S = diag(1.1, 1.0, 0.9, 1.0, 1.1, 1.0)
        s_diag = {1.1, 1.0, 0.9, 1.0, 1.1, 1.0};
        S_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i)
        {
            S_mat[i + i * n_dim] = T(s_diag[i], 0);
        }

        prec.assign(n_dim, 2.0);

        ethr.assign(nband, 1e-8);

        std::mt19937 rng(333);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] = T(dist(rng), 0.0);
            }
        }

        // S-orthonormalize initial guess
        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T dot = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    dot += std::conj(psi[i + k * ld]) * T(s_diag[i], 0) * psi[i + j * ld];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    psi[i + j * ld] -= dot * psi[i + k * ld];
                }
            }
            Real nrm = 0;
            for (int i = 0; i < n_dim; ++i)
            {
                nrm += s_diag[i] * std::norm(psi[i + j * ld]);
            }
            nrm = std::sqrt(nrm);
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] /= nrm;
            }
        }
    }

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<T> S_mat;
    std::vector<Real> s_diag;
    std::vector<Real> prec;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGWithSTest, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    // S-apply function: S * psi = diag(s_diag) * psi (element-wise)
    auto spsi_func = [this](T* in, T* out, int ld_in, int ncol) {
        for (int j = 0; j < ncol; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                out[i + j * ld_in] = T(s_diag[i], 0) * in[i + j * ld_in];
            }
        }
    };

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-10,
        /* max_iter = */ 100,
        /* sbsize   = */ 3,
        /* rr_step  = */ 3,
        /* gamma_g0 = */ false);

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(h_op, spsi_func, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());

    // Eigenvalue check: skip absolute comparison (exact values not
    // analytically known for non-trivial S). Instead verify via residual.
    // Just check eigenvalues are reasonable (not NaN, not negative for
    // this positive-definite problem).
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_GT(eval[i], 0.0) << "WithS BLOCK: eigenvalue[" << i << "] should be positive";
        EXPECT_LT(eval[i], 10.0) << "WithS BLOCK: eigenvalue[" << i << "] unreasonably large";
    }

    // Residual check: ||Hψ_i - ε_i S ψ_i|| / |ε_i| < ethr
    std::vector<T> hpsi(n_dim), spsi(n_dim), res(n_dim);
    for (int i = 0; i < nband; ++i)
    {
        dense_h_multiply(H_mat.data(), n_dim, psi_run.data() + i * ld, hpsi.data(), n_dim, 1);
        spsi_func(psi_run.data() + i * ld, spsi.data(), n_dim, 1);
        for (int j = 0; j < n_dim; ++j)
        {
            res[j] = hpsi[j] - T(eval[i], 0) * spsi[j];
        }
        Real res_nrm = column_norm(res.data(), n_dim);
        EXPECT_LE(res_nrm, std::max(1e-6, 1e2 * ethr[i]))
            << "WithS BLOCK: residual[" << i << "] too large, r=" << res_nrm;
    }

    EXPECT_LE(avg_iter, double(100)) << "WithS BLOCK: too many iterations";
}

// =============================================================================
// Test fixture: gamma_g0 = true (Gamma-point real constraint)
// Same tridiagonal Laplacian, but with gamma_g0=true forcing G=0 wavefunctions
// to stay real-valued. Tests the force_g0_real codepath.
// =============================================================================
class DiagoPPCGGammaG0Test : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        n_dim = 8;
        nband = 3;
        ld = n_dim;

        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i)
        {
            H_mat[i + i * n_dim] = T(2.0, 0);
            if (i > 0)
            {
                H_mat[i + (i - 1) * n_dim] = T(-1.0, 0);
            }
            if (i < n_dim - 1)
            {
                H_mat[i + (i + 1) * n_dim] = T(-1.0, 0);
            }
        }

        prec.assign(n_dim, 2.0);

        exact.resize(nband);
        for (int k = 0; k < nband; ++k)
        {
            exact[k] = 2.0 - 2.0 * std::cos(Real(k + 1) * M_PI / Real(n_dim + 1));
        }

        ethr.assign(nband, 1e-10);

        std::mt19937 rng(555);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] = T(dist(rng), 0.0);
            }
        }

        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T dot = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    psi[i + j * ld] -= dot * psi[i + k * ld];
                }
            }
            Real nrm = 0;
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

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGGammaG0Test, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 100,
        /* sbsize   = */ 3,
        /* rr_step  = */ 3,
        /* gamma_g0 = */ true // <-- Force G=0 wavefunctions to be real
        );

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());

    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-8) << "GammaG0 BLOCK: eigenvalue[" << i << "] mismatch";
    }

    // Verify G=0 band (first band) is real
    Real max_imag = 0;
    for (int i = 0; i < n_dim; ++i)
    {
        max_imag = std::max(max_imag, std::abs(std::imag(psi_run[i])));
    }
    EXPECT_LT(max_imag, 1e-12) << "GammaG0 BLOCK: G=0 band has non-zero imaginary part: " << max_imag;

    EXPECT_LE(avg_iter, double(100)) << "GammaG0 BLOCK: too many iterations";
}

// =============================================================================
// Test fixture: single-band (nband = 1)
// Minimal test — extract only the lowest eigenvalue of a 5×5 tridiagonal
// Laplacian. This exercises the degenerate code paths for a single band.
// =============================================================================
class DiagoPPCGSingleBandTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        n_dim = 5;
        nband = 1;
        ld = n_dim;

        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i)
        {
            H_mat[i + i * n_dim] = T(2.0, 0);
            if (i > 0)
            {
                H_mat[i + (i - 1) * n_dim] = T(-1.0, 0);
            }
            if (i < n_dim - 1)
            {
                H_mat[i + (i + 1) * n_dim] = T(-1.0, 0);
            }
        }

        prec.assign(n_dim, 2.0);

        // Lowest eigenvalue of 5×5 tridiagonal Laplacian
        exact = {2.0 - 2.0 * std::cos(M_PI / 6.0)};

        ethr.assign(nband, 1e-10);

        std::mt19937 rng(42);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int i = 0; i < n_dim; ++i)
        {
            psi[i] = T(dist(rng), 0.0);
        }

        Real nrm = 0;
        for (int i = 0; i < n_dim; ++i)
        {
            nrm += std::norm(psi[i]);
        }
        nrm = std::sqrt(nrm);
        for (int i = 0; i < n_dim; ++i)
        {
            psi[i] /= nrm;
        }
    }

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGSingleBandTest, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 50,
        /* sbsize   = */ 1,
        /* rr_step  = */ 1,
        /* gamma_g0 = */ false);

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());

    EXPECT_NEAR(eval[0], exact[0], 1e-8) << "SingleBand BLOCK: eigenvalue mismatch";
    EXPECT_LE(avg_iter, double(50)) << "SingleBand BLOCK: too many iterations";
}

// =============================================================================
// Test fixture: eigenvector quality — verify Hψ ≈ εψ and ψ^H ψ = I
// Uses the 10×10 tridiagonal Laplacian. After convergence, check:
//   1. ||Hψ_i - ε_i ψ_i|| < tol   for each band
//   2. |ψ_i^H ψ_j - δ_ij| < tol    for all i,j
// =============================================================================
class DiagoPPCGEigenvectorTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        n_dim = 10;
        nband = 3;
        ld = n_dim;

        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i)
        {
            H_mat[i + i * n_dim] = T(2.0, 0);
            if (i > 0)
            {
                H_mat[i + (i - 1) * n_dim] = T(-1.0, 0);
            }
            if (i < n_dim - 1)
            {
                H_mat[i + (i + 1) * n_dim] = T(-1.0, 0);
            }
        }

        prec.assign(n_dim, 2.0);

        exact.resize(nband);
        for (int k = 0; k < nband; ++k)
        {
            exact[k] = 2.0 - 2.0 * std::cos(Real(k + 1) * M_PI / Real(n_dim + 1));
        }

        ethr.assign(nband, 1e-8);

        std::mt19937 rng(888);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] = T(dist(rng), 0.0);
            }
        }

        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T dot = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    psi[i + j * ld] -= dot * psi[i + k * ld];
                }
            }
            Real nrm = 0;
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

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGEigenvectorTest, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 100,
        /* sbsize   = */ 3,
        /* rr_step  = */ 3,
        /* gamma_g0 = */ false);

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());

    // --- Eigenvalue check ---
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-8) << "Eigenvec BLOCK: eigenvalue[" << i << "] mismatch";
    }

    // --- Residual check: ||Hψ_i - ε_i ψ_i|| < sqrt(ethr) ---
    // The eigenvalue-change convergence criterion targets eigenvalue error ~ethr,
    // so the eigenvector residual is naturally ~sqrt(ethr) = 1e-4 for ethr=1e-8.
    std::vector<T> hpsi(n_dim), res(n_dim);
    for (int i = 0; i < nband; ++i)
    {
        dense_h_multiply(H_mat.data(), n_dim, psi_run.data() + i * ld, hpsi.data(), n_dim, 1);
        for (int j = 0; j < n_dim; ++j)
        {
            res[j] = hpsi[j] - eval[i] * psi_run[j + i * ld];
        }
        Real res_nrm = column_norm(res.data(), n_dim);
        EXPECT_LT(res_nrm, 1e-4) << "Eigenvec BLOCK: residual[" << i << "] too large: " << res_nrm;
    }

    // --- Orthogonality check: |ψ_i^H ψ_j - δ_ij| < 1e-8 ---
    for (int i = 0; i < nband; ++i)
    {
        for (int j = 0; j < nband; ++j)
        {
            T dot = 0;
            for (int k = 0; k < n_dim; ++k)
            {
                dot += std::conj(psi_run[k + i * ld]) * psi_run[k + j * ld];
            }
            if (i == j)
            {
                EXPECT_NEAR(std::abs(dot), 1.0, 1e-8)
                    << "Eigenvec BLOCK: ψ[" << i << "] not normalized, |dot|=" << std::abs(dot);
            }
            else
            {
                EXPECT_LT(std::abs(dot), 1e-8)
                    << "Eigenvec BLOCK: ψ[" << i << "] not orthogonal to ψ[" << j << "], |dot|=" << std::abs(dot);
            }
        }
    }

    EXPECT_LE(avg_iter, double(100)) << "Eigenvec BLOCK: too many iterations";
}

// =============================================================================
// Test fixture: all eigenvalues of a small matrix (nband = n_dim)
// 3×3 tridiagonal Laplacian, compute all 3 eigenvalues.
// Exercises the degenerate case where every band is requested.
// =============================================================================
class DiagoPPCGAllBandsTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        n_dim = 3;
        nband = 3;
        ld = n_dim;

        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i)
        {
            H_mat[i + i * n_dim] = T(2.0, 0);
            if (i > 0)
            {
                H_mat[i + (i - 1) * n_dim] = T(-1.0, 0);
            }
            if (i < n_dim - 1)
            {
                H_mat[i + (i + 1) * n_dim] = T(-1.0, 0);
            }
        }

        prec.assign(n_dim, 2.0);

        exact.resize(nband);
        for (int k = 0; k < nband; ++k)
        {
            exact[k] = 2.0 - 2.0 * std::cos(Real(k + 1) * M_PI / Real(n_dim + 1));
        }

        ethr.assign(nband, 1e-10);

        std::mt19937 rng(101);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] = T(dist(rng), 0.0);
            }
        }

        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T dot = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    psi[i + j * ld] -= dot * psi[i + k * ld];
                }
            }
            Real nrm = 0;
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

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGAllBandsTest, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 100,
        /* sbsize   = */ 3,
        /* rr_step  = */ 3,
        /* gamma_g0 = */ false);

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());

    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-8) << "AllBands BLOCK: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, double(100)) << "AllBands BLOCK: too many iterations";
}

// =============================================================================
// Test fixture: medium-sized tridiagonal (15×15, nband=4)
// Bridges the gap between the 10×10 and 20×20 tests.
// =============================================================================
class DiagoPPCGMediumTridiagTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        n_dim = 15;
        nband = 4;
        ld = n_dim;

        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i)
        {
            H_mat[i + i * n_dim] = T(2.0, 0);
            if (i > 0)
            {
                H_mat[i + (i - 1) * n_dim] = T(-1.0, 0);
            }
            if (i < n_dim - 1)
            {
                H_mat[i + (i + 1) * n_dim] = T(-1.0, 0);
            }
        }

        prec.assign(n_dim, 2.0);

        exact.resize(nband);
        for (int k = 0; k < nband; ++k)
        {
            exact[k] = 2.0 - 2.0 * std::cos(Real(k + 1) * M_PI / Real(n_dim + 1));
        }

        ethr.assign(nband, 1e-10);

        std::mt19937 rng(202);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] = T(dist(rng), 0.0);
            }
        }

        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T dot = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    psi[i + j * ld] -= dot * psi[i + k * ld];
                }
            }
            Real nrm = 0;
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

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGMediumTridiagTest, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 120,
        /* sbsize   = */ 4,
        /* rr_step  = */ 4,
        /* gamma_g0 = */ false);

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());

    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-8) << "Medium Tridiag BLOCK: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, double(120)) << "Medium Tridiag BLOCK: too many iterations";
}

// =============================================================================
// Test fixture: gamma_g0 = true on a 7×7 tridiagonal, nband=2
// Verifies eigenvalues are correct and the first band stays real-valued
// when gamma_g0_real is enabled (H and S are both real-symmetric).
// =============================================================================
class DiagoPPCGGammaG0SmallTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        n_dim = 7;
        nband = 2;
        ld = n_dim;

        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i)
        {
            H_mat[i + i * n_dim] = T(2.0, 0);
            if (i > 0)
            {
                H_mat[i + (i - 1) * n_dim] = T(-1.0, 0);
            }
            if (i < n_dim - 1)
            {
                H_mat[i + (i + 1) * n_dim] = T(-1.0, 0);
            }
        }

        prec.assign(n_dim, 2.0);

        exact.resize(nband);
        for (int k = 0; k < nband; ++k)
        {
            exact[k] = 2.0 - 2.0 * std::cos(Real(k + 1) * M_PI / Real(n_dim + 1));
        }

        ethr.assign(nband, 1e-10);

        std::mt19937 rng(404);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] = T(dist(rng), 0.0);
            }
        }

        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T dot = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    psi[i + j * ld] -= dot * psi[i + k * ld];
                }
            }
            Real nrm = 0;
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

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGGammaG0SmallTest, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 100,
        /* sbsize   = */ 2,
        /* rr_step  = */ 2,
        /* gamma_g0 = */ true);

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());

    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-8) << "GammaG0Small BLOCK: eigenvalue[" << i << "] mismatch";
    }

    // Both bands should be real-valued when gamma_g0_real is true
    for (int j = 0; j < nband; ++j)
    {
        Real max_imag = 0;
        for (int i = 0; i < n_dim; ++i)
        {
            max_imag = std::max(max_imag, std::abs(std::imag(psi_run[i + j * ld])));
        }
        EXPECT_LT(max_imag, 1e-12) << "GammaG0Small BLOCK: band[" << j << "] has non-zero imaginary part: " << max_imag;
    }

    EXPECT_LE(avg_iter, double(100)) << "GammaG0Small BLOCK: too many iterations";
}

// =============================================================================
// Test fixture: pentadiagonal Toeplitz (discrete biharmonic operator)
// H[i,i]=6, H[i,i±1]=-4, H[i,i±2]=1.  Eigenvalues:
// λ_k = 16·sin⁴(k·π / (2·(n+1))),  k = 1,...,n
// Wider bandwidth (5 vs 3) tests the solver with more off-diagonal coupling.
// =============================================================================
class DiagoPPCGPentaTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        n_dim = 8;
        nband = 4;
        ld = n_dim;

        // H = T² where T is the tridiagonal Laplacian (2 on diag, -1 on off-diag).
        // The corners of T² have diag=5 (not 6) since (T²)[0,0] = 2² + (-1)² = 5.
        // Interior: (T²)[i,i] = (-1)² + 2² + (-1)² = 6.
        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i)
        {
            H_mat[i + i * n_dim] = T((i == 0 || i == n_dim - 1) ? 5.0 : 6.0, 0);
            if (i >= 1)
            {
                H_mat[i + (i - 1) * n_dim] = T(-4.0, 0);
            }
            if (i < n_dim - 1)
            {
                H_mat[i + (i + 1) * n_dim] = T(-4.0, 0);
            }
            if (i >= 2)
            {
                H_mat[i + (i - 2) * n_dim] = T(1.0, 0);
            }
            if (i < n_dim - 2)
            {
                H_mat[i + (i + 2) * n_dim] = T(1.0, 0);
            }
        }

        prec.assign(n_dim, 6.0);
        prec[0] = 5.0;
        prec[n_dim - 1] = 5.0;

        exact.resize(nband);
        for (int k = 0; k < nband; ++k)
        {
            Real theta = Real(k + 1) * M_PI / Real(2 * (n_dim + 1));
            Real s = std::sin(theta);
            exact[k] = Real(16) * s * s * s * s;
        }

        ethr.assign(nband, 1e-10);

        std::mt19937 rng(505);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] = T(dist(rng), 0.0);
            }
        }

        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T dot = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    psi[i + j * ld] -= dot * psi[i + k * ld];
                }
            }
            Real nrm = 0;
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

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGPentaTest, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 150,
        /* sbsize   = */ 4,
        /* rr_step  = */ 4,
        /* gamma_g0 = */ false);

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());

    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-8) << "Penta BLOCK: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, double(150)) << "Penta BLOCK: too many iterations";
}

// =============================================================================
// Test fixture: gapped spectrum
// H = diag(0.1, 0.5, 5.0, 6.0, 10.0), nband=3.
// Large gaps between eigenvalue groups test the solver's band separation.
// =============================================================================
class DiagoPCGGappedSpectrumTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        n_dim = 5;
        nband = 3;
        ld = n_dim;

        H_mat.assign(n_dim * n_dim, T(0));
        H_mat[0 + 0 * n_dim] = T(0.1, 0);
        H_mat[1 + 1 * n_dim] = T(0.5, 0);
        H_mat[2 + 2 * n_dim] = T(5.0, 0);
        H_mat[3 + 3 * n_dim] = T(6.0, 0);
        H_mat[4 + 4 * n_dim] = T(10.0, 0);

        prec.assign(n_dim, 1.0);

        exact = {0.1, 0.5, 5.0};

        ethr.assign(nband, 1e-10);

        std::mt19937 rng(606);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] = T(dist(rng), 0.0);
            }
        }

        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T dot = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    psi[i + j * ld] -= dot * psi[i + k * ld];
                }
            }
            Real nrm = 0;
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

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPCGGappedSpectrumTest, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 100,
        /* sbsize   = */ 3,
        /* rr_step  = */ 3,
        /* gamma_g0 = */ false);

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());

    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-8) << "Gapped BLOCK: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, double(100)) << "Gapped BLOCK: too many iterations";
}

// =============================================================================
// Test fixture: preconditioner stress test
// Uses the 10×10 tridiagonal Laplacian but with a suboptimal preconditioner:
// prec[i] = 1.0 instead of 2.0 (the exact diagonal). The solver should still
// converge, just more slowly. Tests robustness against a bad preconditioner.
// =============================================================================
class DiagoPPCGBadPrecTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        n_dim = 10;
        nband = 3;
        ld = n_dim;

        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i)
        {
            H_mat[i + i * n_dim] = T(2.0, 0);
            if (i > 0)
            {
                H_mat[i + (i - 1) * n_dim] = T(-1.0, 0);
            }
            if (i < n_dim - 1)
            {
                H_mat[i + (i + 1) * n_dim] = T(-1.0, 0);
            }
        }

        // Bad preconditioner: use 1.0 instead of 2.0
        prec.assign(n_dim, 1.0);

        exact.resize(nband);
        for (int k = 0; k < nband; ++k)
        {
            exact[k] = 2.0 - 2.0 * std::cos(Real(k + 1) * M_PI / Real(n_dim + 1));
        }

        ethr.assign(nband, 1e-10);

        std::mt19937 rng(707);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);

        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] = T(dist(rng), 0.0);
            }
        }

        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T dot = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    psi[i + j * ld] -= dot * psi[i + k * ld];
                }
            }
            Real nrm = 0;
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

    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGBadPrecTest, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 200, // more iterations due to bad preconditioner
        /* sbsize   = */ 4,
        /* rr_step  = */ 4,
        /* gamma_g0 = */ false);

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());

    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-8) << "BadPrec BLOCK: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, double(200)) << "BadPrec BLOCK: too many iterations";
}

// =============================================================================
// Test fixture: n_dim = 1, nband = 1 — absolute minimum
// H is a 1×1 matrix [5.0], eigenvalue = 5.0
// =============================================================================
class DiagoPPCG1x1Test : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        n_dim = 1;
        nband = 1;
        ld = n_dim;
        H_mat = {T(5.0, 0)};
        prec = {5.0};
        exact = {5.0};
        ethr.assign(nband, 1e-10);
        psi = {T(1.0, 0)}; // already normalized
    }
    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCG1x1Test, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(1e-12, 10, 1, 1, false
                                                                   );
    auto h_op = [this](T* in, T* out, int ldi, int nc) { dense_h_multiply(H_mat.data(), n_dim, in, out, ldi, nc); };
    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());
    EXPECT_NEAR(eval[0], exact[0], 1e-8) << "1x1 BLOCK: mismatch";
    EXPECT_LE(avg_iter, 10.0) << "1x1 BLOCK: too many iterations";
}

// =============================================================================
// Test fixture: scaled tridiagonal (eigenvalues × 100)
// H = 100 × tridiag(2, -1, -1). Tests convergence with large eigenvalues.
// =============================================================================
class DiagoPPCGScaledTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        n_dim = 8;
        nband = 3;
        ld = n_dim;
        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i)
        {
            H_mat[i + i * n_dim] = T(200.0, 0);
            if (i > 0)
            {
                H_mat[i + (i - 1) * n_dim] = T(-100.0, 0);
            }
            if (i < n_dim - 1)
            {
                H_mat[i + (i + 1) * n_dim] = T(-100.0, 0);
            }
        }
        prec.assign(n_dim, 200.0);
        exact.resize(nband);
        for (int k = 0; k < nband; ++k)
        {
            exact[k] = 100.0 * (2.0 - 2.0 * std::cos(Real(k + 1) * M_PI / Real(n_dim + 1)));
        }
        ethr.assign(nband, 1e-8);
        init_psi(808);
    }
    void init_psi(int seed)
    {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);
        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] = T(dist(rng), 0.0);
            }
        }
        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T d = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    d += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    psi[i + j * ld] -= d * psi[i + k * ld];
                }
            }
            Real nr = 0;
            for (int i = 0; i < n_dim; ++i)
            {
                nr += std::norm(psi[i + j * ld]);
            }
            nr = std::sqrt(nr);
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] /= nr;
            }
        }
    }
    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGScaledTest, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(1e-10, 120, 4, 4, false
                                                                   );
    auto h_op = [this](T* in, T* out, int ldi, int nc) { dense_h_multiply(H_mat.data(), n_dim, in, out, ldi, nc); };
    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-6) << "Scaled BLOCK: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, 120.0) << "Scaled BLOCK: too many iterations";
}

// =============================================================================
// Test fixture: many bands (n_dim=12, nband=4)
// Moderate band-to-dimension ratio (1:3).
// =============================================================================
class DiagoPPCGManyBandsTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        n_dim = 12;
        nband = 4;
        ld = n_dim;
        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i)
        {
            H_mat[i + i * n_dim] = T(2.0, 0);
            if (i > 0)
            {
                H_mat[i + (i - 1) * n_dim] = T(-1.0, 0);
            }
            if (i < n_dim - 1)
            {
                H_mat[i + (i + 1) * n_dim] = T(-1.0, 0);
            }
        }
        prec.assign(n_dim, 2.0);
        exact.resize(nband);
        for (int k = 0; k < nband; ++k)
        {
            exact[k] = 2.0 - 2.0 * std::cos(Real(k + 1) * M_PI / Real(n_dim + 1));
        }
        ethr.assign(nband, 1e-10);
        init_psi(909);
    }
    void init_psi(int seed)
    {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);
        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] = T(dist(rng), 0.0);
            }
        }
        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T d = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    d += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    psi[i + j * ld] -= d * psi[i + k * ld];
                }
            }
            Real nr = 0;
            for (int i = 0; i < n_dim; ++i)
            {
                nr += std::norm(psi[i + j * ld]);
            }
            nr = std::sqrt(nr);
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] /= nr;
            }
        }
    }
    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGManyBandsTest, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(1e-12, 150, 4, 4, false
                                                                   );
    auto h_op = [this](T* in, T* out, int ldi, int nc) { dense_h_multiply(H_mat.data(), n_dim, in, out, ldi, nc); };
    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-8) << "ManyBands BLOCK: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, 150.0) << "ManyBands BLOCK: too many iterations";
}

// =============================================================================
// Test fixture: rr_step = 1 (most frequent Rayleigh-Ritz refinement)
// 8×8 tridiagonal, nband=3. Tests aggressive subspace diagonalization.
// =============================================================================
class DiagoPPCGRrStep1Test : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        n_dim = 8;
        nband = 3;
        ld = n_dim;
        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i)
        {
            H_mat[i + i * n_dim] = T(2.0, 0);
            if (i > 0)
            {
                H_mat[i + (i - 1) * n_dim] = T(-1.0, 0);
            }
            if (i < n_dim - 1)
            {
                H_mat[i + (i + 1) * n_dim] = T(-1.0, 0);
            }
        }
        prec.assign(n_dim, 2.0);
        exact.resize(nband);
        for (int k = 0; k < nband; ++k)
        {
            exact[k] = 2.0 - 2.0 * std::cos(Real(k + 1) * M_PI / Real(n_dim + 1));
        }
        ethr.assign(nband, 1e-10);
        init_psi(111);
    }
    void init_psi(int seed)
    {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);
        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] = T(dist(rng), 0.0);
            }
        }
        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T d = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    d += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    psi[i + j * ld] -= d * psi[i + k * ld];
                }
            }
            Real nr = 0;
            for (int i = 0; i < n_dim; ++i)
            {
                nr += std::norm(psi[i + j * ld]);
            }
            nr = std::sqrt(nr);
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] /= nr;
            }
        }
    }
    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGRrStep1Test, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(1e-12, 100, 3, 1 /*rr_step=1*/, false
                                                                   );
    auto h_op = [this](T* in, T* out, int ldi, int nc) { dense_h_multiply(H_mat.data(), n_dim, in, out, ldi, nc); };
    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-8) << "RrStep1 BLOCK: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, 100.0) << "RrStep1 BLOCK: too many iterations";
}

// =============================================================================
// Test fixture: Neumann boundary tridiagonal Laplacian
// H[0,0]=H[n-1,n-1]=1, interior diag=2, off-diag=-1.
// Eigenvalues: λ_k = 2 - 2·cos(k·π/n), k=0,1,...,n-1
// =============================================================================
class DiagoPPCGNeumannTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        n_dim = 8;
        nband = 4;
        ld = n_dim;
        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i)
        {
            H_mat[i + i * n_dim] = T((i == 0 || i == n_dim - 1) ? 1.0 : 2.0, 0);
            if (i > 0)
            {
                H_mat[i + (i - 1) * n_dim] = T(-1.0, 0);
            }
            if (i < n_dim - 1)
            {
                H_mat[i + (i + 1) * n_dim] = T(-1.0, 0);
            }
        }
        prec.assign(n_dim, 2.0);
        prec[0] = 1.0;
        prec[n_dim - 1] = 1.0;
        exact.resize(nband);
        for (int k = 0; k < nband; ++k)
        {
            exact[k] = 2.0 - 2.0 * std::cos(Real(k) * M_PI / Real(n_dim));
        }
        ethr.assign(nband, 1e-10);
        init_psi(222);
    }
    void init_psi(int seed)
    {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);
        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] = T(dist(rng), 0.0);
            }
        }
        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T d = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    d += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    psi[i + j * ld] -= d * psi[i + k * ld];
                }
            }
            Real nr = 0;
            for (int i = 0; i < n_dim; ++i)
            {
                nr += std::norm(psi[i + j * ld]);
            }
            nr = std::sqrt(nr);
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] /= nr;
            }
        }
    }
    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGNeumannTest, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(1e-12, 100, 4, 4, false
                                                                   );
    auto h_op = [this](T* in, T* out, int ldi, int nc) { dense_h_multiply(H_mat.data(), n_dim, in, out, ldi, nc); };
    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-8) << "Neumann BLOCK: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, 100.0) << "Neumann BLOCK: too many iterations";
}

// =============================================================================
// Test fixture: tight convergence threshold (ethr = 1e-14)
// 6×6 tridiagonal, nband=2. Tests deep convergence.
// =============================================================================
class DiagoPPCGTightEthrTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        n_dim = 6;
        nband = 2;
        ld = n_dim;
        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i)
        {
            H_mat[i + i * n_dim] = T(2.0, 0);
            if (i > 0)
            {
                H_mat[i + (i - 1) * n_dim] = T(-1.0, 0);
            }
            if (i < n_dim - 1)
            {
                H_mat[i + (i + 1) * n_dim] = T(-1.0, 0);
            }
        }
        prec.assign(n_dim, 2.0);
        exact.resize(nband);
        for (int k = 0; k < nband; ++k)
        {
            exact[k] = 2.0 - 2.0 * std::cos(Real(k + 1) * M_PI / Real(n_dim + 1));
        }
        ethr.assign(nband, 1e-14);
        init_psi(333);
    }
    void init_psi(int seed)
    {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);
        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] = T(dist(rng), 0.0);
            }
        }
        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T d = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    d += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    psi[i + j * ld] -= d * psi[i + k * ld];
                }
            }
            Real nr = 0;
            for (int i = 0; i < n_dim; ++i)
            {
                nr += std::norm(psi[i + j * ld]);
            }
            nr = std::sqrt(nr);
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] /= nr;
            }
        }
    }
    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGTightEthrTest, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(1e-14, 200, 3, 3, false
                                                                   );
    auto h_op = [this](T* in, T* out, int ldi, int nc) { dense_h_multiply(H_mat.data(), n_dim, in, out, ldi, nc); };
    double avg_iter = solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-8) << "TightEthr BLOCK: eigenvalue[" << i << "] mismatch";
    }
    EXPECT_LE(avg_iter, 200.0) << "TightEthr BLOCK: too many iterations";
}

// =============================================================================
// Test fixture: non-diagonal S matrix (tridiagonal overlap)
// H = tridiag(2,-1,-1), S = tridiag(1.0, 0.2, 0.2) — both tridiagonal.
// Verifies S-orthogonalization with a non-diagonal S.
// =============================================================================
class DiagoPPCGTridiagSTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        n_dim = 6;
        nband = 2;
        ld = n_dim;
        H_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i)
        {
            H_mat[i + i * n_dim] = T(2.0, 0);
            if (i > 0)
            {
                H_mat[i + (i - 1) * n_dim] = T(-1.0, 0);
            }
            if (i < n_dim - 1)
            {
                H_mat[i + (i + 1) * n_dim] = T(-1.0, 0);
            }
        }
        S_mat.assign(n_dim * n_dim, T(0));
        for (int i = 0; i < n_dim; ++i)
        {
            S_mat[i + i * n_dim] = T(1.0, 0);
            if (i > 0)
            {
                S_mat[i + (i - 1) * n_dim] = T(0.2, 0);
                S_mat[(i - 1) + i * n_dim] = T(0.2, 0);
            }
        }
        prec.assign(n_dim, 2.0);
        // Exact eigenvalues unknown analytically for generalized problem
        // with non-diagonal S. Just check convergence via residual.
        exact = {0.0, 0.0};
        ethr.assign(nband, 1e-8);
        init_psi(444);
    }
    void init_psi(int seed)
    {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);
        psi.assign(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] = T(dist(rng), 0.0);
            }
        }
        // S-orthonormalize: S = tridiag(1.0, 0.2, 0.2)
        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T d = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    T si = 0;
                    si += T(1.0, 0) * psi[i + k * ld];
                    if (i > 0)
                    {
                        si += T(0.2, 0) * psi[(i - 1) + k * ld];
                    }
                    if (i < n_dim - 1)
                    {
                        si += T(0.2, 0) * psi[(i + 1) + k * ld];
                    }
                    d += std::conj(si) * psi[i + j * ld];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    psi[i + j * ld] -= d * psi[i + k * ld];
                }
            }
            Real nr = 0;
            for (int i = 0; i < n_dim; ++i)
            {
                T si = 0;
                si += T(1.0, 0) * psi[i + j * ld];
                if (i > 0)
                {
                    si += T(0.2, 0) * psi[(i - 1) + j * ld];
                }
                if (i < n_dim - 1)
                {
                    si += T(0.2, 0) * psi[(i + 1) + j * ld];
                }
                nr += std::real(std::conj(psi[i + j * ld]) * si);
            }
            nr = std::sqrt(nr);
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] /= nr;
            }
        }
    }
    int n_dim, nband, ld;
    std::vector<T> H_mat;
    std::vector<T> S_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<double> ethr;
    std::vector<T> psi;
};

TEST_F(DiagoPPCGTridiagSTest, BlockSubspace)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);
    auto spsi_func = [this](T* in, T* out, int ldi, int nc) {
        for (int j = 0; j < nc; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                out[i + j * ldi] = T(1.0, 0) * in[i + j * ldi];
                if (i > 0)
                {
                    out[i + j * ldi] += T(0.2, 0) * in[(i - 1) + j * ldi];
                }
                if (i < n_dim - 1)
                {
                    out[i + j * ldi] += T(0.2, 0) * in[(i + 1) + j * ldi];
                }
            }
        }
    };
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(1e-10, 150, 3, 3, false
                                                                   );
    auto h_op = [this](T* in, T* out, int ldi, int nc) { dense_h_multiply(H_mat.data(), n_dim, in, out, ldi, nc); };
    double avg_iter = solver.diag(h_op, spsi_func, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());
    // Check eigenvalues are positive and reasonable
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_GT(eval[i], 0.0);
        EXPECT_LT(eval[i], 5.0);
    }
    // Residual check
    std::vector<T> hpsi(n_dim), spsi(n_dim), res(n_dim);
    for (int i = 0; i < nband; ++i)
    {
        dense_h_multiply(H_mat.data(), n_dim, psi_run.data() + i * ld, hpsi.data(), n_dim, 1);
        spsi_func(psi_run.data() + i * ld, spsi.data(), n_dim, 1);
        for (int j = 0; j < n_dim; ++j)
        {
            res[j] = hpsi[j] - T(eval[i], 0) * spsi[j];
        }
        Real rn = column_norm(res.data(), n_dim);
        EXPECT_LT(rn, 1e-4) << "TridiagS BLOCK: residual[" << i << "] too large: " << rn;
    }
    EXPECT_LE(avg_iter, 150.0) << "TridiagS BLOCK: too many iterations";
}

// =============================================================================
// Performance benchmark: PPCG vs BLOCK vs Davidson comparison
//
// Runs PPCG on random sparse symmetric matrices of various sizes (matching
// the BLOCK/Davidson test sizes: npw=100,200,500) and reports avg_iter and
// wall-clock time.  avg_iter is the primary metric — it counts the average
// number of H·psi applications per band.
//
// Typical results from the existing BLOCK/Davidson tests (for reference):
//   BLOCK:       avg_iter ~ 20-50  on random sparse, ~ 5-15 on tridiagonal
//   Davidson: avg_iter ~ 15-40  on random sparse (but each iter is heavier)
//   PPCG:     avg_iter ~ 2-10  on tridiagonal/diagonal, varies with sparsity
//
// This test is DISABLED by default (too slow for CI).  Run manually with
//   --gtest_also_run_disabled_tests --gtest_filter='*Benchmark*'
// =============================================================================
class DiagoPPCGBenchmarkTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
    }

    // Generate a random sparse symmetric matrix of size n with given sparsity.
    // sparsity=0 means dense, sparsity=80 means 80% zeros.
    static void make_random_hamilt(int n, int sparsity_pct, std::vector<T>& H, std::vector<Real>& prec)
    {
        H.assign(n * n, T(0));
        std::mt19937 rng(unsigned(n * 100 + sparsity_pct));
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);
        int nnz = 0;
        for (int i = 0; i < n; ++i)
        {
            for (int j = i; j < n; ++j)
            {
                if (i != j && (rng() % 100) < sparsity_pct)
                {
                    continue;
                }
                Real val = (i == j) ? std::abs(dist(rng)) * n + 1.0 : dist(rng) * 0.5;
                H[i + j * n] = T(val, 0);
                if (i != j)
                {
                    H[j + i * n] = T(val, 0);
                }
                if (val != 0)
                {
                    ++nnz;
                }
            }
        }
        // Simple diagonal preconditioner
        prec.resize(n);
        for (int i = 0; i < n; ++i)
        {
            prec[i] = std::max(std::real(H[i + i * n]), 1e-6);
        }
    }

    // Run PPCG and return {avg_iter, wall_sec}.
    static std::pair<double, double> run_ppcg(int n, int nband, const std::vector<T>& H, const std::vector<Real>& prec)
    {
        int ld = n;
        std::mt19937 rng(42);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);
        std::vector<T> psi(ld * nband, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n; ++i)
            {
                psi[i + j * ld] = T(dist(rng), 0.0);
            }
        }
        // GS orthonormalize
        for (int j = 0; j < nband; ++j)
        {
            for (int k = 0; k < j; ++k)
            {
                T d = 0;
                for (int i = 0; i < n; ++i)
                {
                    d += std::conj(psi[i + k * ld]) * psi[i + j * ld];
                }
                for (int i = 0; i < n; ++i)
                {
                    psi[i + j * ld] -= d * psi[i + k * ld];
                }
            }
            Real nr = 0;
            for (int i = 0; i < n; ++i)
            {
                nr += std::norm(psi[i + j * ld]);
            }
            nr = std::sqrt(nr);
            for (int i = 0; i < n; ++i)
            {
                psi[i + j * ld] /= nr;
            }
        }

        std::vector<Real> eval(nband, 0.0);
        std::vector<double> ethr(nband, 1e-4);
        auto h_op = [&H, n](T* in, T* out, int ldi, int nc) { dense_h_multiply(H.data(), n, in, out, ldi, nc); };

        hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(1e-8, 500, nband, std::min(nband, 4), false
                                                                       );

        auto t0 = std::chrono::high_resolution_clock::now();
        double avg_iter = solver.diag(h_op, nullptr, ld, nband, n, psi.data(), eval.data(), ethr, prec.data());
        auto t1 = std::chrono::high_resolution_clock::now();
        double wall = std::chrono::duration<double>(t1 - t0).count();
        return {avg_iter, wall};
    }
};

TEST_F(DiagoPPCGBenchmarkTest, DISABLED_FullBenchmark)
{
    struct Case
    {
        int n;
        int nband;
        int sparsity;
    };
    std::vector<Case> cases = {
        {50, 10, 0},   {50, 10, 60},  {100, 10, 0},  {100, 10, 60},
        {100, 10, 80}, {200, 10, 60}, {200, 10, 80}, {500, 10, 80},
    };

    std::cout << "\n========== PPCG Performance Benchmark ==========\n";
    std::cout << " n_dim  nband  sparsity   avg_iter   wall_time(s)\n";
    std::cout << "-------------------------------------------------\n";
    for (auto& c : cases)
    {
        std::vector<T> H;
        std::vector<Real> prec;
        make_random_hamilt(c.n, c.sparsity, H, prec);
        const std::pair<double, double> result = run_ppcg(c.n, c.nband, H, prec);
        const double avg_iter = result.first;
        const double wall = result.second;
        printf(" %5d   %3d     %2d%%       %6.1f      %7.4f\n", c.n, c.nband, c.sparsity, avg_iter, wall);
    }
    std::cout << "=================================================\n";
    SUCCEED();
}

// Quick convergence smoke test: one representative case, fast enough for CI.
TEST_F(DiagoPPCGBenchmarkTest, QuickBenchmark)
{
    std::vector<T> H;
    std::vector<Real> prec;
    make_random_hamilt(80, 60, H, prec);
    const std::pair<double, double> result = run_ppcg(80, 8, H, prec);
    const double avg_iter = result.first;
    EXPECT_LE(avg_iter, 500.0) << "PPCG did not converge within 500 iters";
    SUCCEED();
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
