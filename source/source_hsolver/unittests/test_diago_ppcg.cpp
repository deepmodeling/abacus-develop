/**
 * test_diago_ppcg.cpp — correctness and regression tests for DiagoPPCG
 *
 * Test matrices include both S = I and non-trivial overlap operators:
 *   1. Tridiagonal Laplacian (1D particle-in-a-box): H[i,i]=2, H[i,i±1]=-1
 *      Exact λ_k = 2 - 2·cos(k·π/(n+1)).  Realistic but sparse.
 *   2. Diagonal matrix: H = diag(1, 2, 3, 4, 5)
 *      Exact eigenvalues are the diagonal entries.  Simplest possible
 *      smoke test — should converge in very few iterations.
 *
 * Analytical matrix cases share setup below; numerical failure regressions keep
 * their explicit initial vectors and assertions.
 */

#include "diago_ppcg_test_support.h"
#include "source_base/timer.h"
#ifdef __MPI
#include "source_basis/module_pw/test/test_tool.h"
#endif

#include <gtest/gtest.h>
#include <algorithm>
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

using ppcg_test::DiagoPPCGMatrixTest;
using ppcg_test::PPCG_PI;
using ppcg_test::Real;
using ppcg_test::T;
using ppcg_test::dense_h_multiply;
using ppcg_test::initialize_vectors;

// Like the CG test's preparation table, each row describes one complete solve.
struct PPCGLaplacianCase
{
    const char* name;
    int dimension;
    int bands;
    unsigned seed;
    double band_threshold;
    double matrix_scale;
    double diag_threshold;
    int max_iterations;
    int block_size;
    int ritz_interval;
    double eigenvalue_tolerance;
};

static void PrintTo(const PPCGLaplacianCase& config, std::ostream* out)
{
    *out << config.name;
}

class DiagoPPCGLaplacianTest : public DiagoPPCGMatrixTest,
                               public ::testing::WithParamInterface<PPCGLaplacianCase>
{
};

TEST_P(DiagoPPCGLaplacianTest, MatchesAnalyticalEigenvalues)
{
    const PPCGLaplacianCase& config = GetParam();
    initialize(config.dimension, config.bands, config.seed, config.band_threshold);
    set_laplacian(config.matrix_scale);
    check_eigenvalues(config.diag_threshold, config.max_iterations, config.block_size,
                      config.ritz_interval, config.eigenvalue_tolerance);
}

static const PPCGLaplacianCase laplacian_cases[] = {
    // name, dimension, bands, seed, band threshold, scale, diag threshold, budget, block, Ritz interval, tolerance
    {"Tridiagonal10", 10, 3, 42, 1e-10, 1.0, 1e-12, 100, 4, 4, 1e-8},
    {"SingleBand", 5, 1, 42, 1e-10, 1.0, 1e-12, 50, 1, 1, 1e-8},
    {"AllBands", 3, 3, 101, 1e-10, 1.0, 1e-12, 100, 3, 3, 1e-8},
    {"Laplacian15", 15, 4, 202, 1e-10, 1.0, 1e-12, 120, 4, 4, 1e-8},
    {"Laplacian20", 20, 5, 789, 1e-10, 1.0, 1e-12, 150, 5, 5, 1e-8},
    {"ScaledLaplacian", 8, 3, 808, 1e-8, 100.0, 1e-10, 120, 4, 4, 1e-6},
    {"ManyBands", 12, 4, 909, 1e-10, 1.0, 1e-12, 150, 4, 4, 1e-8},
    {"RitzEveryIteration", 8, 3, 111, 1e-10, 1.0, 1e-12, 100, 3, 1, 1e-8},
    {"TightThreshold", 6, 2, 333, 1e-14, 1.0, 1e-14, 200, 3, 3, 1e-8}};

INSTANTIATE_TEST_SUITE_P(
    Laplacians,
    DiagoPPCGLaplacianTest,
    ::testing::ValuesIn(laplacian_cases),
    [](const ::testing::TestParamInfo<PPCGLaplacianCase>& info) { return info.param.name; });

// Test fixture: 1D particle-in-a-box (tridiagonal Laplacian)
class DiagoPPCGTridiagTest : public DiagoPPCGMatrixTest
{
  protected:
    void SetUp() override
    {
        initialize(10, 3, 42, 1e-10);
        set_laplacian(1.0);
    }
};

TEST_F(DiagoPPCGTridiagTest, ReportsIterationLimit)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);
    std::vector<double> strict_ethr(nband, 1.0e-14);
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1.0e-14,
        /* max_iter = */ 1,
        /* sbsize   = */ 3,
        /* rr_step  = */ 4);

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), strict_ethr, prec.data());

    EXPECT_FALSE(solver.converged());
    EXPECT_GT(solver.active_band_count(), 0);
    EXPECT_TRUE(solver.iteration_limit_reached());
    EXPECT_GE(solver.first_active_band(), 0);
    EXPECT_GT(solver.max_residual(), Real(1.0e-14));
}

TEST_F(DiagoPPCGTridiagTest, ReusesProjectedOperatorImages)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);
    const std::vector<double> strict_ethr(nband, 1.0e-14);
    int h_calls = 0;
    int s_calls = 0;
    const auto h_op = [this, &h_calls](T* in, T* out, int ld_in, int ncol) {
        ++h_calls;
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };
    const auto s_op = [this, &s_calls](T* in, T* out, int ld_in, int ncol) {
        ++s_calls;
        for (int j = 0; j < ncol; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                out[i + j * ld_in] = (1.0 + 0.1 * i) * in[i + j * ld_in];
            }
        }
    };
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(1.0e-14, 2, 3, 16);
    solver.diag(h_op, s_op, ld, nband, n_dim,
                psi_run.data(), eval.data(), strict_ethr, prec.data());

    ASSERT_TRUE(solver.iteration_limit_reached());
    // Initial/final X applications, one W application per sweep, and one
    // cancellation-triggered P refresh. No unconditional HP or projected SW
    // recomputation is needed, including when S is not I.
    EXPECT_EQ(h_calls, 5);
    EXPECT_EQ(s_calls, 5);
    double max_residual = 0.0;
    for (int j = 0; j < nband; ++j)
    {
        double residual_squared = 0.0;
        for (int i = 0; i < n_dim; ++i)
        {
            T hpsi = T(0);
            for (int k = 0; k < n_dim; ++k)
            {
                hpsi += H_mat[i + k * n_dim] * psi_run[k + j * ld];
            }
            const T spsi = (1.0 + 0.1 * i) * psi_run[i + j * ld];
            residual_squared += std::norm(hpsi - eval[j] * spsi);
        }
        const Real residual_norm = std::sqrt(residual_squared);
        max_residual = std::max(max_residual, residual_norm);
    }
    EXPECT_NEAR(solver.max_residual(), max_residual, 1.0e-13);
}

TEST(DiagoPPCGReliableUpdateTest, RefreshesBeforeUsingNearPrecisionResidual)
{
    const int dimension = 2;
    const std::vector<Real> prec(dimension, 1.0);
    const std::vector<double> thresholds(1, 1.0e-14);
    const double perturbations[] = {1.0e-5, 1.0e-10};
    for (const double perturbation : perturbations)
    {
        SCOPED_TRACE(perturbation);
        const double normalization = std::sqrt(1.0 + perturbation * perturbation);
        std::vector<T> psi = {T(1.0 / normalization), T(perturbation / normalization)};
        std::vector<Real> eigenvalues(1, 0.0);
        int h_calls = 0;
        const auto h_op = [&h_calls](T* in, T* out, int ld, int ncol) {
            ++h_calls;
            for (int j = 0; j < ncol; ++j)
            {
                out[j * ld] = in[j * ld];
                out[1 + j * ld] = 3.0 * in[1 + j * ld];
            }
        };
        hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
            1.0e-14, 1, 1, 1);
        solver.diag(h_op, nullptr, dimension, 1, dimension,
                    psi.data(), eigenvalues.data(), thresholds, prec.data());
        EXPECT_TRUE(solver.converged());
        EXPECT_LT(solver.max_residual(), 1.0e-14);
        EXPECT_NEAR(eigenvalues[0], 1.0, 1.0e-14);
        const T residual0 = (1.0 - eigenvalues[0]) * psi[0];
        const T residual1 = (3.0 - eigenvalues[0]) * psi[1];
        const double residual = std::sqrt(std::norm(residual0) + std::norm(residual1));
        EXPECT_LT(residual, 1.0e-14);
        // Initial X, W, tentative-convergence validation, and final validation.
        // Only the near-precision solve adds a reliable X/P refresh before W.
        const int expected_calls = perturbation < 1.0e-8 ? 6 : 4;
        EXPECT_EQ(h_calls, expected_calls);
    }
}

TEST(DiagoPPCGPeriodicRitzTest, HonorsIntervalAndChecksFinalRitzResidual)
{
    const int dimension = 48;
    const int bands = 8;
    std::mt19937 random(184);
    std::uniform_real_distribution<double> uniform(-1.0, 1.0);
    std::vector<T> initial(dimension * bands);
    for (T& value : initial)
    {
        const double real = uniform(random);
        const double imaginary = uniform(random);
        value = T(real, imaginary);
    }
    const std::vector<Real> prec(dimension, 2.0);
    const std::vector<double> thresholds(bands, 1.0e-14);
    struct Configuration
    {
        int sweeps;
        int interval;
    };
    // Ten sweeps need no periodic refresh; eighteen exercise the one at 16.
    // Each configuration additionally encounters one projection cancellation.
    // Each run encounters one competing discarded direction and
    // therefore add one conditional expanded Ritz solve, without another H call.
    const Configuration configurations[] = {{10, 1}, {10, 4}, {10, 16},
                                             {18, 1}, {18, 4}, {18, 16}};
    for (const Configuration& configuration : configurations)
    {
        const int sweeps = configuration.sweeps;
        const int interval = configuration.interval;
        SCOPED_TRACE(sweeps);
        SCOPED_TRACE(interval);
        std::vector<T> psi = initial;
        std::vector<Real> eigenvalues(bands, 0.0);
        int applications = 0;
        const auto h_op = [&applications, dimension](T* in, T* out, int ld, int ncol) {
            ++applications;
            for (int j = 0; j < ncol; ++j)
            {
                for (int i = 0; i < dimension; ++i)
                {
                    T value = 2.0 * in[i + j * ld];
                    if (i > 0)
                    {
                        value -= in[i - 1 + j * ld];
                    }
                    if (i + 1 < dimension)
                    {
                        value -= in[i + 1 + j * ld];
                    }
                    out[i + j * ld] = value;
                }
            }
        };
        hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
            1.0e-14, sweeps, 2, interval);
        solver.diag(h_op, nullptr, dimension, bands, dimension,
                    psi.data(), eigenvalues.data(), thresholds, prec.data());
        ASSERT_TRUE(solver.iteration_limit_reached());
        const int scheduled = sweeps / interval;
        const int final_extra = sweeps % interval == 0 ? 0 : 1;
        const int startup_extra = interval == 1 ? 0 : 1;
        const int competing_extra = 1;
        const int expected_ritz = 1 + scheduled + final_extra + startup_extra + competing_extra;
        const int operator_refreshes = sweeps / 16;
        const int cancellation_refreshes = 1;
        const int expected_applications = 2 + sweeps + 2 * operator_refreshes
                                          + cancellation_refreshes;
        EXPECT_EQ(solver.root_recovery_count(), competing_extra);
        EXPECT_EQ(solver.rayleigh_ritz_count(), expected_ritz);
        EXPECT_EQ(applications, expected_applications);
        EXPECT_TRUE(std::is_sorted(eigenvalues.begin(), eigenvalues.end()));
        std::vector<T> hpsi(psi.size());
        h_op(psi.data(), hpsi.data(), dimension, bands);
        double residual = 0.0;
        for (int j = 0; j < bands; ++j)
        {
            double norm_squared = 0.0;
            for (int i = 0; i < dimension; ++i)
            {
                const T value = hpsi[i + j * dimension] - eigenvalues[j] * psi[i + j * dimension];
                norm_squared += std::norm(value);
            }
            const Real residual_norm = std::sqrt(norm_squared);
            residual = std::max(residual, residual_norm);
            for (int k = 0; k < bands; ++k)
            {
                T overlap = T(0);
                for (int i = 0; i < dimension; ++i)
                {
                    overlap += std::conj(psi[i + k * dimension]) * psi[i + j * dimension];
                }
                const T expected = j == k ? T(1) : T(0);
                const double error = std::abs(overlap - expected);
                EXPECT_LT(error, 1.0e-11);
            }
        }
        EXPECT_NEAR(solver.max_residual(), residual, 1.0e-13);
    }
}

TEST(DiagoPPCGPeriodicRitzTest, StartupReclassifiesDisplacedLockedPrefix)
{
    const int dimension = 5;
    const int bands = 2;
    std::vector<T> psi(dimension * bands, T(0));
    // The initial exact eigenvalue 2 appears converged, but it is not the
    // lowest root. The startup correction discovers a lower approximate
    // root, which must replace and unlock this column before the next sweep.
    psi[1] = T(1);
    const double scale = 1.0 / std::sqrt(3.0);
    psi[dimension] = T(scale);
    psi[2 + dimension] = T(scale);
    psi[4 + dimension] = T(scale);
    std::vector<Real> eigenvalues(bands, 0.0);
    const std::vector<Real> prec(dimension, 1.0);
    const std::vector<double> thresholds(bands, 1.0e-10);
    const auto h_op = [dimension](T* in, T* out, int ld, int ncol) {
        for (int j = 0; j < ncol; ++j)
        {
            for (int i = 0; i < dimension; ++i)
            {
                out[i + j * ld] = double(i + 1) * in[i + j * ld];
            }
        }
    };
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(1.0e-12, 5, 1, 16);
    solver.diag(h_op, nullptr, dimension, bands, dimension,
                psi.data(), eigenvalues.data(), thresholds, prec.data());
    EXPECT_TRUE(solver.converged());
    EXPECT_NEAR(eigenvalues[0], 1.0, 1.0e-10);
    EXPECT_NEAR(eigenvalues[1], 2.0, 1.0e-10);
    EXPECT_LE(solver.max_residual(), 1.0e-10);
}

TEST(DiagoPPCGPeriodicRitzTest, ClusteredGeneralizedSpectrumWithLockedPrefix)
{
    const int dimension = 16;
    const int bands = 4;
    std::vector<Real> spectrum(dimension);
    spectrum[0] = -1.0;
    spectrum[1] = 0.0;
    spectrum[2] = 1.0e-6;
    spectrum[3] = 1.0e-6;
    spectrum[4] = 2.0e-6;
    std::vector<Real> overlap(dimension);
    std::vector<Real> prec(dimension, 1.0);
    for (int i = 0; i < dimension; ++i)
    {
        overlap[i] = std::pow(10.0, 3.0 * i / (dimension - 1));
        if (i >= 5)
        {
            spectrum[i] = 1.0 + i;
        }
        prec[i] = overlap[i];
    }
    std::mt19937 random(392);
    std::uniform_real_distribution<double> uniform(-1.0, 1.0);
    std::vector<T> initial(dimension * bands, T(0));
    initial[0] = T(1.0 / std::sqrt(overlap[0]));
    for (int j = 1; j < bands; ++j)
    {
        for (int i = 1; i < dimension; ++i)
        {
            const double real = uniform(random);
            const double imaginary = uniform(random);
            initial[i + j * dimension] = T(real, imaginary) / std::sqrt(overlap[i]);
        }
    }
    const auto h_op = [&spectrum, &overlap, dimension](T* in, T* out, int ld, int ncol) {
        for (int j = 0; j < ncol; ++j)
        {
            for (int i = 0; i < dimension; ++i)
            {
                out[i + j * ld] = spectrum[i] * overlap[i] * in[i + j * ld];
            }
        }
    };
    const auto s_op = [&overlap, dimension](T* in, T* out, int ld, int ncol) {
        for (int j = 0; j < ncol; ++j)
        {
            for (int i = 0; i < dimension; ++i)
            {
                out[i + j * ld] = overlap[i] * in[i + j * ld];
            }
        }
    };
    const std::vector<double> thresholds(bands, 1.0e-10);
    const int intervals[] = {1, 16};
    for (const int interval : intervals)
    {
        SCOPED_TRACE(interval);
        std::vector<T> psi = initial;
        std::vector<Real> eigenvalues(bands, 0.0);
        hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
            1.0e-12, 600, 2, interval);
        solver.diag(h_op, s_op, dimension, bands, dimension,
                    psi.data(), eigenvalues.data(), thresholds, prec.data());
        EXPECT_TRUE(solver.converged());
        EXPECT_LE(solver.max_residual(), 1.0e-10);
        for (int j = 0; j < bands; ++j)
        {
            EXPECT_NEAR(eigenvalues[j], spectrum[j], 1.0e-10);
            for (int k = 0; k < bands; ++k)
            {
                T product = T(0);
                for (int i = 0; i < dimension; ++i)
                {
                    product += std::conj(psi[i + k * dimension]) * overlap[i] * psi[i + j * dimension];
                }
                const T expected = j == k ? T(1) : T(0);
                const double error = std::abs(product - expected);
                EXPECT_LT(error, 1.0e-10);
            }
        }
        const double prefix_weight = std::norm(psi[0]) * overlap[0];
        EXPECT_NEAR(prefix_weight, 1.0, 1.0e-12);
    }
}

TEST(DiagoPPCGStatusTest, ReportsActiveTailStart)
{
    const int n_dim = 4;
    const int nband = 3;
    const int ld = n_dim;
    std::vector<T> h_mat(n_dim * n_dim, T(0));
    for (int i = 0; i < n_dim; ++i)
    {
        h_mat[i + i * ld] = T(Real(i + 1), 0);
    }
    const Real inverse_sqrt_two = Real(1) / std::sqrt(Real(2));
    std::vector<T> psi(ld * nband, T(0));
    psi[0] = T(1, 0);
    psi[1 + ld] = T(1, 0);
    psi[2 + 2 * ld] = T(inverse_sqrt_two, 0);
    psi[3 + 2 * ld] = T(inverse_sqrt_two, 0);
    std::vector<Real> eval(nband, Real(0));
    const std::vector<double> ethr(nband, 1.0e-12);
    const std::vector<Real> prec = {1, 2, 3, 4};
    const auto h_op = [&h_mat, n_dim](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(h_mat.data(), n_dim, in, out, ld_in, ncol);
    };
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1.0e-12,
        /* max_iter = */ 0,
        /* sbsize   = */ 3,
        /* rr_step  = */ 1);

    solver.diag(h_op, nullptr, ld, nband, n_dim,
                psi.data(), eval.data(), ethr, prec.data());

    EXPECT_FALSE(solver.converged());
    EXPECT_EQ(solver.first_active_band(), 2);
    EXPECT_EQ(solver.active_band_count(), 1);
    EXPECT_DOUBLE_EQ(solver.first_active_threshold(), Real(1.0e-12));
    EXPECT_GT(solver.max_residual(), solver.first_active_threshold());
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
        /* rr_step  = */ 2);

    solver.diag(h_op, nullptr, ld, nband, n_dim,
                psi.data(), eval.data(), ethr, prec.data());

    EXPECT_FALSE(solver.converged());
    EXPECT_EQ(solver.active_band_count(), nband);
    EXPECT_FALSE(solver.iteration_limit_reached());
    EXPECT_EQ(solver.first_active_band(), 0);
    EXPECT_TRUE(std::isinf(solver.max_residual()));
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
        /* rr_step  = */ 2);

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
    const std::string trace_path = testing::TempDir() + "./ppcg_residual_trace_test.csv";
    std::remove(trace_path.c_str());

    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 100,
        /* sbsize   = */ 4,
        /* rr_step  = */ 4);

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    std::ostringstream running_log;
    auto& residual_timer = ModuleBase::timer::timer_pool["DiagoPPCG"]["max_generalized_residual"];
    const size_t before_disabled = residual_timer.calls;
    ASSERT_EQ(::unsetenv(env_name), 0);
    solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data(), running_log);
    const size_t disabled_calls = residual_timer.calls - before_disabled;
    EXPECT_TRUE(solver.converged());

    // An unset variable and a failed open must both skip trace-only work.
    // Opening the existing temporary directory fails without creating a file.
    const std::string directory_path = testing::TempDir();
    ASSERT_EQ(::setenv(env_name, directory_path.c_str(), 1), 0);
    psi_run = psi;
    const size_t before_failed_open = residual_timer.calls;
    solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data(), running_log);
    const size_t failed_open_calls = residual_timer.calls - before_failed_open;
    EXPECT_EQ(disabled_calls, failed_open_calls);
    EXPECT_TRUE(solver.converged());
    EXPECT_TRUE(running_log.str().empty());

    ASSERT_EQ(::setenv(env_name, trace_path.c_str(), 1), 0);
    psi_run = psi;
    const size_t before_enabled = residual_timer.calls;
    solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data(), running_log);
    const size_t enabled_calls = residual_timer.calls - before_enabled;
    EXPECT_GT(enabled_calls, disabled_calls);

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
    char* resolved_path = ::realpath(trace_path.c_str(), nullptr);
    ASSERT_NE(resolved_path, nullptr);
    const std::string full_path(resolved_path);
    std::free(resolved_path);
    const std::string expected_log = "PPCG residual trace written: " + full_path;
    EXPECT_NE(running_log.str().find(expected_log), std::string::npos);
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

// Two explicit startup subspaces exercise the same low-root recovery check.
class DiagoPPCGStartupCorrectionTest : public DiagoPPCGMatrixTest,
                                     public ::testing::WithParamInterface<bool>
{
};

TEST_P(DiagoPPCGStartupCorrectionTest, KeepsLowestRoots)
{
    n_dim = 8;
    nband = 4;
    ld = n_dim;
    const Real inv_sqrt_two = Real(1) / std::sqrt(Real(2));
    const int matrix_size = n_dim * n_dim;
    H_mat.assign(matrix_size, T(0));
    prec.resize(n_dim);
    for (int i = 0; i < n_dim; ++i)
    {
        const Real diagonal = Real(i + 1);
        H_mat[i + i * n_dim] = T(diagonal, 0);
        prec[i] = diagonal;
    }

    const int vector_size = ld * nband;
    psi.assign(vector_size, T(0));
    psi[0] = T(inv_sqrt_two, 0);
    psi[2] = T(inv_sqrt_two, 0);
    psi[1 + ld] = T(inv_sqrt_two, 0);
    psi[3 + ld] = T(inv_sqrt_two, 0);
    const bool all_vectors_active = GetParam();
    if (all_vectors_active)
    {
        // Every initial vector has a nonzero residual. Independent two-column
        // blocks would nevertheless retain roots 1, 2, 5, and 6, discarding the
        // third and fourth roots before a global comparison is possible.
        psi[4 + 2 * ld] = T(inv_sqrt_two, 0);
        psi[6 + 2 * ld] = T(inv_sqrt_two, 0);
        psi[5 + 3 * ld] = T(inv_sqrt_two, 0);
        psi[7 + 3 * ld] = T(inv_sqrt_two, 0);
    }
    else
    {
        // The first two columns contain the missing third and fourth roots. The
        // last two are exact high-energy eigenvectors and would otherwise lock
        // before their lower replacements can enter the active subspace.
        psi[6 + 2 * ld] = T(1, 0);
        psi[7 + 3 * ld] = T(1, 0);
    }

    ethr.assign(nband, 1e-12);
    hsolver::DiagoPPCG<T, base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 50,
        /* sbsize   = */ 2,
        /* rr_step  = */ 4);
    const double avg_iter = solve(solver);
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], Real(i + 1), 1e-10) << "root[" << i << "] mismatch";
    }
    EXPECT_TRUE(solver.converged());
    EXPECT_LE(avg_iter, 50.0);
}

INSTANTIATE_TEST_SUITE_P(
    StartupVectors,
    DiagoPPCGStartupCorrectionTest,
    ::testing::Bool(),
    [](const ::testing::TestParamInfo<bool>& info) {
        return info.param ? "AllVectorsActive" : "LockedHighRoots";
    });

TEST(DiagoPPCGRootSelectionTest, StartupCorrectionKeepsConvergedPrefixInactive)
{
    const int n_dim = 10;
    const int nband = 6;
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

    // The first two roots are already converged.  The active tail still needs
    // one global comparison: its first two vectors contain roots 3--6, while
    // its last two vectors are exact high-energy states.
    std::vector<T> psi(ld * nband, T(0));
    psi[0] = T(1, 0);
    psi[1 + ld] = T(1, 0);
    psi[2 + 2 * ld] = T(inv_sqrt_two, 0);
    psi[4 + 2 * ld] = T(inv_sqrt_two, 0);
    psi[3 + 3 * ld] = T(inv_sqrt_two, 0);
    psi[5 + 3 * ld] = T(inv_sqrt_two, 0);
    psi[8 + 4 * ld] = T(1, 0);
    psi[9 + 5 * ld] = T(1, 0);

    std::vector<Real> eval(nband, Real(0));
    const std::vector<double> ethr(nband, 1e-12);
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 1,
        /* sbsize   = */ 2,
        /* rr_step  = */ 4);
    const auto h_op = [&h_mat, n_dim](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(h_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    const double avg_iter = solver.diag(
        h_op, nullptr, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());

    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], Real(i + 1), 1e-10) << "root[" << i << "] mismatch";
    }
    const double expected_avg_iter = 1.0 + 4.0 / 6.0;
    EXPECT_TRUE(solver.converged());
    EXPECT_NEAR(avg_iter, expected_avg_iter, 1e-12);
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
        /* rr_step  = */ 4);
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
        /* rr_step  = */ 4);

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
            /* rr_step  = */ 1);
        const double avg_iter = solver.diag(
            h_op, nullptr, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());
        EXPECT_NEAR(eval[0], 1.0, 1e-8);
        EXPECT_NEAR(eval[1], 2.0, 1e-8);
        EXPECT_TRUE(solver.converged());
        EXPECT_LE(avg_iter, 80.0);
    }
}

TEST(DiagoPPCGRootSelectionTest, LaterCorrectionKeepsLowestRootsAcrossBlocks)
{
    const int n_dim = 4;
    const int nband = 2;
    const int ld = n_dim;
    const std::vector<Real> diagonal = {1.0, 2.0, 3.0, 100.0};
    const std::vector<Real> prec(n_dim, 1.0);
    const std::vector<double> ethr(nband, 1.0e-10);
    const Real coefficient = 1.0 / std::sqrt(3.0);
    const int block_sizes[] = {1, 2};
    const auto h_op = [&diagonal, n_dim](T* in, T* out, int ld_in, int ncol) {
        for (int column = 0; column < ncol; ++column)
        {
            for (int row = 0; row < n_dim; ++row)
            {
                out[row + column * ld_in] = diagonal[row] * in[row + column * ld_in];
            }
        }
    };

    for (const int block_size : block_sizes)
    {
        SCOPED_TRACE("block_size = " + std::to_string(block_size));
        // The first expanded subspace keeps a low mixture and the exact
        // third root. The missing second root needs a later correction;
        // convergence of the retained vectors alone does not certify it.
        std::vector<T> psi(ld * nband, T(0));
        psi[0] = T(coefficient);
        psi[1] = T(coefficient);
        psi[3] = T(coefficient);
        psi[2 + ld] = T(1);
        std::vector<Real> eval(nband, 0.0);
        hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
            1.0e-10, 100, block_size, 1);
        solver.diag(h_op, nullptr, ld, nband, n_dim,
                    psi.data(), eval.data(), ethr, prec.data());

        double max_residual = 0.0;
        for (int column = 0; column < nband; ++column)
        {
            double residual_squared = 0.0;
            for (int row = 0; row < n_dim; ++row)
            {
                const T residual = (diagonal[row] - eval[column]) * psi[row + column * ld];
                residual_squared += std::norm(residual);
            }
            const double residual_norm = std::sqrt(residual_squared);
            max_residual = std::max(max_residual, residual_norm);
            for (int other = 0; other < nband; ++other)
            {
                T overlap(0);
                for (int row = 0; row < n_dim; ++row)
                {
                    overlap += std::conj(psi[row + column * ld]) * psi[row + other * ld];
                }
                const T expected_overlap = column == other ? T(1) : T(0);
                const double overlap_error = std::abs(overlap - expected_overlap);
                EXPECT_LE(overlap_error, 1.0e-12);
            }
        }
        EXPECT_TRUE(solver.converged());
        EXPECT_LE(max_residual, 1.0e-10);
        EXPECT_NEAR(solver.max_residual(), max_residual, 1.0e-13);
        EXPECT_NEAR(eval[0], 1.0, 1.0e-9);
        EXPECT_NEAR(eval[1], 2.0, 1.0e-9)
            << "returned roots = " << eval[0] << ", " << eval[1]
            << "; fresh max residual = " << max_residual;
        if (block_size == 1)
        {
            EXPECT_GT(solver.root_recovery_count(), 0);
        }
        else
        {
            EXPECT_EQ(solver.root_recovery_count(), 0);
        }
    }
}

TEST(DiagoPPCGRootSelectionTest, LaterCorrectionsRecoverMultipleComplexGeneralizedRoots)
{
    const int dimension = 6;
    const int bands = 3;
    const Real spectrum[] = {1.0, 2.0, 3.0, 4.0, 5.0, 100.0};
    const Real overlap_spectrum[] = {1.0, 2.0, 4.0, 7.0, 3.0, 5.0};
    std::vector<T> unitary(dimension * dimension);
    const Real normalization = 1.0 / std::sqrt(Real(dimension));
    for (int column = 0; column < dimension; ++column)
    {
        for (int row = 0; row < dimension; ++row)
        {
            const Real angle = 2.0 * PPCG_PI * row * column / dimension;
            const T phase(std::cos(angle), std::sin(angle));
            unitary[row + column * dimension] = normalization * phase;
        }
    }
    std::vector<T> hamiltonian(dimension * dimension, T(0));
    std::vector<T> overlap(dimension * dimension, T(0));
    for (int column = 0; column < dimension; ++column)
    {
        for (int row = 0; row < dimension; ++row)
        {
            for (int state = 0; state < dimension; ++state)
            {
                const T projector = unitary[row + state * dimension]
                                    * std::conj(unitary[column + state * dimension]);
                hamiltonian[row + column * dimension] += projector * overlap_spectrum[state] * spectrum[state];
                overlap[row + column * dimension] += projector * overlap_spectrum[state];
            }
        }
    }
    const auto h_op = [&hamiltonian, dimension](T* in, T* out, int ld, int count) {
        dense_h_multiply(hamiltonian.data(), dimension, in, out, ld, count);
    };
    const auto s_op = [&overlap, dimension](T* in, T* out, int ld, int count) {
        dense_h_multiply(overlap.data(), dimension, in, out, ld, count);
    };
    const std::vector<Real> preconditioner(dimension, 1.0);
    const std::vector<double> thresholds(bands, 1.0e-9);
    const int block_sizes[] = {1, 2, 3};
    const int ritz_intervals[] = {1, 4};
    for (const int block_size : block_sizes)
    {
        for (const int ritz_interval : ritz_intervals)
        {
            SCOPED_TRACE("block_size = " + std::to_string(block_size)
                         + ", ritz_interval = " + std::to_string(ritz_interval));
            std::vector<T> psi(dimension * bands, T(0));
            // Four spectral components in the first vector require more
            // than the first [X,W] expansion to resolve all three low roots.
            const int mixed_states[] = {0, 1, 2, 5};
            for (int row = 0; row < dimension; ++row)
            {
                for (const int state : mixed_states)
                {
                    psi[row] += 0.5 * unitary[row + state * dimension]
                                / std::sqrt(overlap_spectrum[state]);
                }
                psi[row + dimension] = unitary[row + 3 * dimension] / std::sqrt(overlap_spectrum[3]);
                psi[row + 2 * dimension] = unitary[row + 4 * dimension] / std::sqrt(overlap_spectrum[4]);
            }
            std::vector<Real> eigenvalues(bands, 0.0);
            hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
                1.0e-9, 120, block_size, ritz_interval);
            solver.diag(h_op, s_op, dimension, bands, dimension,
                        psi.data(), eigenvalues.data(), thresholds, preconditioner.data());
            std::vector<T> hpsi(psi.size(), T(0));
            std::vector<T> spsi(psi.size(), T(0));
            h_op(psi.data(), hpsi.data(), dimension, bands);
            s_op(psi.data(), spsi.data(), dimension, bands);
            EXPECT_TRUE(solver.converged());
            for (int column = 0; column < bands; ++column)
            {
                EXPECT_NEAR(eigenvalues[column], spectrum[column], 1.0e-8);
                double residual_squared = 0.0;
                for (int row = 0; row < dimension; ++row)
                {
                    const int offset = row + column * dimension;
                    const T residual = hpsi[offset] - eigenvalues[column] * spsi[offset];
                    residual_squared += std::norm(residual);
                }
                const double residual_norm = std::sqrt(residual_squared);
                EXPECT_LE(residual_norm, 1.0e-9);
                for (int other = 0; other < bands; ++other)
                {
                    T inner_product(0);
                    for (int row = 0; row < dimension; ++row)
                    {
                        inner_product += std::conj(psi[row + column * dimension]) * spsi[row + other * dimension];
                    }
                    const T expected = column == other ? T(1) : T(0);
                    const double orthogonality_error = std::abs(inner_product - expected);
                    EXPECT_LE(orthogonality_error, 1.0e-11);
                }
            }
        }
    }
}

// Test fixture: 2×2 matrix — smallest non-trivial case
// H = [[2, 1], [1, 2]], eigenvalues: 1, 3
TEST_F(DiagoPPCGMatrixTest, TwoByTwo)
{
    initialize(2, 2, 123, 1e-10);
    // H = [[2, 1], [1, 2]]
    H_mat.assign(n_dim * n_dim, T(0));
    H_mat[0 + 0 * n_dim] = T(2.0, 0);
    H_mat[1 + 1 * n_dim] = T(2.0, 0);
    H_mat[0 + 1 * n_dim] = T(1.0, 0);
    H_mat[1 + 0 * n_dim] = T(1.0, 0);

    prec.assign(n_dim, 2.0);

    // λ₁ = 1, λ₂ = 3
    exact = {1.0, 3.0};
    check_eigenvalues(1e-12, 50, 2, 2, 1e-8);
}

TEST_F(DiagoPPCGMatrixTest, DegenerateEigenvalues)
{
    initialize(4, 4, 456, 1e-10);
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
    check_eigenvalues(1e-12, 100, 4, 4, 1e-8);
}

// Test fixture: dense matrix with known eigenvalues
// H = Q^T * D * Q where Q is a known orthogonal matrix (a Givens rotation
// repeated on different index pairs) and D is diagonal.
// For an 8×8 case: D = diag(1, 2, 3, 4, 5, 6, 7, 8), then apply several
// Givens rotations to mix all rows/cols. The exact eigenvalues remain 1..8.
TEST_F(DiagoPPCGMatrixTest, DenseKnownSpectrum)
{
    initialize(8, 4, 222, 1e-10);
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
        const Real angle = angle_dist(rng_dense);
        apply_givens(p, q, angle);
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
    check_eigenvalues(1e-12, 200, 4, 4, 1e-8);
}

// Euclidean column norm for independent residual checks
static Real column_norm(const T* x, int n_dim)
{
    Real nrm = 0;
    for (int i = 0; i < n_dim; ++i)
    {
        nrm += std::norm(x[i]);
    }
    return std::sqrt(nrm);
}

struct PPCGOverlapCase
{
    const char* name;
    unsigned seed;
    int bands;
    int leading_dimension;
    int max_iterations;
    bool tridiagonal;
    double residual_threshold;
};

static void PrintTo(const PPCGOverlapCase& config, std::ostream* out)
{
    *out << config.name;
}

// H is the six-dimensional Laplacian. Diagonal S also exercises padded vectors;
// tridiagonal S shares H's sine eigenvectors and has analytical eigenvalues.
class DiagoPPCGGeneralizedTest : public DiagoPPCGMatrixTest,
                                 public ::testing::WithParamInterface<PPCGOverlapCase>
{
  protected:
    void SetUp() override
    {
        const PPCGOverlapCase& config = GetParam();
        n_dim = 6;
        nband = config.bands;
        ld = config.leading_dimension;
        ethr.assign(nband, 1e-8);
        set_laplacian(1.0);
        const int matrix_size = n_dim * n_dim;
        S_mat.assign(matrix_size, T(0));
        const Real diagonal_overlap[] = {1.1, 1.0, 0.9, 1.0, 1.1, 1.0};
        for (int i = 0; i < n_dim; ++i)
        {
            const Real diagonal = config.tridiagonal ? 1.0 : diagonal_overlap[i];
            S_mat[i + i * n_dim] = diagonal;
            if (config.tridiagonal && i > 0)
            {
                S_mat[i + (i - 1) * n_dim] = 0.2;
                S_mat[i - 1 + i * n_dim] = 0.2;
            }
        }
        if (config.tridiagonal)
        {
            for (int k = 0; k < nband; ++k)
            {
                const Real theta = Real(k + 1) * PPCG_PI / Real(n_dim + 1);
                const Real cosine = std::cos(theta);
                exact[k] = (2.0 - 2.0 * cosine) / (1.0 + 0.4 * cosine);
            }
        }

        std::mt19937 rng(config.seed);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);
        const int vector_size = ld * nband;
        psi.assign(vector_size, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                const Real value = dist(rng);
                psi[i + j * ld] = value;
            }
        }
        std::vector<T> s_column(n_dim);
        for (int j = 0; j < nband; ++j)
        {
            T* column = psi.data() + j * ld;
            for (int k = 0; k < j; ++k)
            {
                const T* previous = psi.data() + k * ld;
                dense_h_multiply(S_mat.data(), n_dim, previous, s_column.data(), n_dim, 1);
                T overlap = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    overlap += std::conj(s_column[i]) * column[i];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    column[i] -= overlap * previous[i];
                }
            }
            dense_h_multiply(S_mat.data(), n_dim, column, s_column.data(), n_dim, 1);
            Real norm_squared = 0;
            for (int i = 0; i < n_dim; ++i)
            {
                const T norm_term = std::conj(column[i]) * s_column[i];
                norm_squared += std::real(norm_term);
            }
            ASSERT_GT(norm_squared, 0.0);
            const Real norm = std::sqrt(norm_squared);
            for (int i = 0; i < n_dim; ++i)
            {
                column[i] /= norm;
            }
        }
    }

    void apply_overlap(T* in, T* out, int leading_dimension, int columns) const
    {
        const bool tridiagonal = GetParam().tridiagonal;
        for (int j = 0; j < columns; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                const int offset = i + j * leading_dimension;
                T value = S_mat[i + i * n_dim] * in[offset];
                if (tridiagonal && i > 0)
                {
                    value += 0.2 * in[offset - 1];
                }
                if (tridiagonal && i + 1 < n_dim)
                {
                    value += 0.2 * in[offset + 1];
                }
                out[offset] = value;
            }
        }
    }

    std::vector<T> S_mat;
};

TEST_P(DiagoPPCGGeneralizedTest, EigenpairsAndSOrthogonality)
{
    const PPCGOverlapCase& config = GetParam();
    std::vector<Real> eval(nband, 0.0);
    hsolver::DiagoPPCG<T, base_device::DEVICE_CPU> solver(1e-10, config.max_iterations, 3, 3);
    const auto h_op = [this](T* in, T* out, int stride, int columns) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, stride, columns);
    };
    const auto s_op = [this](T* in, T* out, int stride, int columns) {
        apply_overlap(in, out, stride, columns);
    };
    const double iterations = solver.diag(
        h_op, s_op, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());
    if (config.tridiagonal)
    {
        EXPECT_TRUE(solver.converged());
        EXPECT_EQ(solver.active_band_count(), 0);
        EXPECT_FALSE(solver.iteration_limit_reached());
    }

    // Use fresh dense products, independently of the structured S callback.
    const int image_size = ld * nband;
    std::vector<T> hpsi(image_size);
    std::vector<T> spsi(image_size);
    std::vector<T> residual(n_dim);
    dense_h_multiply(H_mat.data(), n_dim, psi.data(), hpsi.data(), ld, nband);
    dense_h_multiply(S_mat.data(), n_dim, psi.data(), spsi.data(), ld, nband);
    for (int i = 0; i < nband; ++i)
    {
        if (config.tridiagonal)
        {
            EXPECT_NEAR(eval[i], exact[i], ethr[i]) << "band " << i;
        }
        else
        {
            EXPECT_GT(eval[i], 0.0) << "band " << i;
            EXPECT_LT(eval[i], 10.0) << "band " << i;
        }
        for (int j = 0; j < n_dim; ++j)
        {
            const int offset = j + i * ld;
            residual[j] = hpsi[offset] - eval[i] * spsi[offset];
        }
        const Real residual_norm = column_norm(residual.data(), n_dim);
        EXPECT_LE(residual_norm, config.residual_threshold) << "band " << i;
        for (int k = 0; k < nband; ++k)
        {
            T overlap = 0;
            for (int j = 0; j < n_dim; ++j)
            {
                overlap += std::conj(psi[j + k * ld]) * spsi[j + i * ld];
            }
            const T expected_overlap = (i == k) ? T(1) : T(0);
            const T difference = overlap - expected_overlap;
            const Real overlap_error = std::abs(difference);
            EXPECT_LE(overlap_error, 1e-10) << "S-overlap at (" << k << ", " << i << ")";
        }
    }
    const double iteration_budget = config.max_iterations;
    EXPECT_LE(iterations, iteration_budget);
}

static const PPCGOverlapCase generalized_overlap_cases[] = {
    // name, seed, bands, leading dimension, iteration budget, tridiagonal, residual tolerance
    {"Diagonal", 333, 3, 8, 100, false, 1e-6},
    {"Tridiagonal", 444, 2, 6, 150, true, 1e-8}};

INSTANTIATE_TEST_SUITE_P(
    OverlapOperators,
    DiagoPPCGGeneralizedTest,
    ::testing::ValuesIn(generalized_overlap_cases),
    [](const ::testing::TestParamInfo<PPCGOverlapCase>& info) { return info.param.name; });

// Test fixture: eigenvector quality — verify Hψ ≈ εψ and ψ^H ψ = I
// Uses the 10×10 tridiagonal Laplacian. After convergence, check:
//   1. ||Hψ_i - ε_i ψ_i|| < tol   for each band
//   2. |ψ_i^H ψ_j - δ_ij| < tol    for all i,j
TEST_F(DiagoPPCGMatrixTest, EigenvectorQuality)
{
    initialize(10, 3, 888, 1e-8);
    set_laplacian(1.0);
    check_eigenvalues(1e-12, 100, 3, 3, 1e-8);
    check_eigenvectors(1e-4, 1e-8);
}

// Test fixture: pentadiagonal Toeplitz (discrete biharmonic operator)
// H[i,i]=6, H[i,i±1]=-4, H[i,i±2]=1.  Eigenvalues:
// λ_k = 16·sin⁴(k·π / (2·(n+1))),  k = 1,...,n
// Wider bandwidth (5 vs 3) tests the solver with more off-diagonal coupling.
TEST_F(DiagoPPCGMatrixTest, SquaredLaplacian)
{
    initialize(8, 4, 505, 1e-10);
    // H = T² where T is the tridiagonal Laplacian (2 on diag, -1 on off-diag).
    // The corners of T² have diag=5 (not 6) since (T²)[0,0] = 2² + (-1)² = 5.
    // Interior: (T²)[i,i] = (-1)² + 2² + (-1)² = 6.
    H_mat.assign(n_dim * n_dim, T(0));
    for (int i = 0; i < n_dim; ++i)
    {
        const Real diagonal_value = (i == 0 || i == n_dim - 1) ? 5.0 : 6.0;
        H_mat[i + i * n_dim] = T(diagonal_value, 0);
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
        Real theta = Real(k + 1) * PPCG_PI / Real(2 * (n_dim + 1));
        Real s = std::sin(theta);
        exact[k] = Real(16) * s * s * s * s;
    }
    check_eigenvalues(1e-12, 150, 4, 4, 1e-8);
}

// Test fixture: gapped spectrum
// H = diag(0.1, 0.5, 5.0, 6.0, 10.0), nband=3.
// Large gaps between eigenvalue groups test the solver's band separation.
TEST_F(DiagoPPCGMatrixTest, GappedSpectrum)
{
    initialize(5, 3, 606, 1e-10);
    H_mat.assign(n_dim * n_dim, T(0));
    H_mat[0 + 0 * n_dim] = T(0.1, 0);
    H_mat[1 + 1 * n_dim] = T(0.5, 0);
    H_mat[2 + 2 * n_dim] = T(5.0, 0);
    H_mat[3 + 3 * n_dim] = T(6.0, 0);
    H_mat[4 + 4 * n_dim] = T(10.0, 0);

    prec.assign(n_dim, 1.0);

    exact = {0.1, 0.5, 5.0};
    check_eigenvalues(1e-12, 100, 3, 3, 1e-8);
}

// Test fixture: preconditioner stress test
// Uses the 10×10 tridiagonal Laplacian with a fixed, nonuniform positive
// preconditioner that differs from the exact diagonal. The test checks that
// the solver still converges and returns a small independently computed residual.
TEST_F(DiagoPPCGMatrixTest, NonuniformPreconditioner)
{
    initialize(10, 3, 707, 1e-10);
    set_laplacian(1.0);
    for (int i = 0; i < n_dim; ++i)
    {
        prec[i] = 0.6 + 0.1 * Real(i);
    }
    hsolver::DiagoPPCG<T, base_device::DEVICE_CPU> solver(1e-12, 200, 4, 4);
    const double iterations = solve(solver);
    EXPECT_TRUE(solver.converged());
    EXPECT_LE(iterations, 200.0);
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-8) << "band " << i;
    }
    check_eigenvectors(1e-8, 1e-8);
}

// Test fixture: n_dim = 1, nband = 1 — absolute minimum
// H is a 1×1 matrix [5.0], eigenvalue = 5.0
TEST_F(DiagoPPCGMatrixTest, ScalarMatrix)
{
    initialize(1, 1, 0, 1e-10);
    H_mat = {T(5.0, 0)};
    prec = {5.0};
    exact = {5.0};
    psi = {T(1.0, 0)}; // already normalized
    check_eigenvalues(1e-12, 10, 1, 1, 1e-8);
}

// Test fixture: Neumann boundary tridiagonal Laplacian
// H[0,0]=H[n-1,n-1]=1, interior diag=2, off-diag=-1.
// Eigenvalues: λ_k = 2 - 2·cos(k·π/n), k=0,1,...,n-1
TEST_F(DiagoPPCGMatrixTest, NeumannZeroEigenvalue)
{
    initialize(8, 4, 222, 1e-10);
    H_mat.assign(n_dim * n_dim, T(0));
    for (int i = 0; i < n_dim; ++i)
    {
        const Real diagonal_value = (i == 0 || i == n_dim - 1) ? 1.0 : 2.0;
        H_mat[i + i * n_dim] = T(diagonal_value, 0);
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
        exact[k] = 2.0 - 2.0 * std::cos(Real(k) * PPCG_PI / Real(n_dim));
    }
    check_eigenvalues(1e-12, 100, 4, 4, 1e-8);
}

// A deterministic sparse matrix exercises convergence beyond analytical fixtures.
TEST(DiagoPPCGRandomMatrixTest, ConvergesWithDiagonalPreconditioner)
{
    const int n_dim = 80;
    const int nband = 8;
    const int sparsity_pct = 60;
    const int matrix_size = n_dim * n_dim;
    const unsigned seed = n_dim * 100 + sparsity_pct;
    std::vector<T> h(matrix_size, T(0));
    std::mt19937 rng(seed);
    std::uniform_real_distribution<Real> dist(-1.0, 1.0);
    for (int i = 0; i < n_dim; ++i)
    {
        for (int j = i; j < n_dim; ++j)
        {
            if (i != j && (rng() % 100) < sparsity_pct)
            {
                continue;
            }
            const Real value = (i == j) ? std::abs(dist(rng)) * n_dim + 1.0 : dist(rng) * 0.5;
            h[i + j * n_dim] = T(value, 0);
            h[j + i * n_dim] = T(value, 0);
        }
    }
    std::vector<Real> prec(n_dim);
    for (int i = 0; i < n_dim; ++i)
    {
        const Real diagonal = std::real(h[i + i * n_dim]);
        prec[i] = std::max(diagonal, 1e-6);
    }
    std::vector<T> psi;
    initialize_vectors(psi, n_dim, nband, n_dim, 42);
    std::vector<Real> eval(nband, 0.0);
    const std::vector<double> ethr(nband, 1e-4);
    const auto h_op = [&h, n_dim](T* in, T* out, int ld, int ncol) {
        dense_h_multiply(h.data(), n_dim, in, out, ld, ncol);
    };
    hsolver::DiagoPPCG<T, base_device::DEVICE_CPU> solver(1e-8, 500, nband, 4);
    const double iterations = solver.diag(h_op, nullptr, n_dim, nband, n_dim,
                                          psi.data(), eval.data(), ethr, prec.data());
    EXPECT_TRUE(solver.converged());
    EXPECT_TRUE(std::isfinite(solver.max_residual()));
    EXPECT_LE(solver.max_residual(), 1e-4);
    EXPECT_LE(iterations, 500.0);
}

int main(int argc, char** argv)
{
#ifdef __MPI
    int nproc = 1;
    int myrank = 0;
    setupmpi(argc, argv, nproc, myrank);
    int nproc_in_pool = 1;
    int kpar = nproc;
    int mypool = 0;
    int rank_in_pool = 0;
    divide_pools(nproc, myrank, nproc_in_pool, kpar, mypool, rank_in_pool);
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    finishmpi();
#endif
    return result;
}
