// Iteration budgets, operator image reuse, convergence status and trace handling.
#include "diago_ppcg_test_support.h"
#include "source_base/timer.h"

#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <string>

using ppcg_test::DiagoPPCGMatrixTest;
using ppcg_test::Real;
using ppcg_test::T;
using ppcg_test::dense_h_multiply;

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
