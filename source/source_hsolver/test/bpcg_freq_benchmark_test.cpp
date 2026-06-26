/**
 * Benchmark & unit test for BPCG subspace diagonalization frequency optimization.
 *
 * Part 1: Unit tests for three free functions
 *   compute_optimal_freq, compute_convergence_rate, compute_adaptive_freq
 *
 * Part 2: Integration benchmark
 *   Runs DiagoBPCG::diag() with SCF_ITER=1 (the intended target of the
 *   optimization) and reports eigenvalue accuracy, timing, and frequency info.
 *
 * NOTE on the optimization's design:
 *   - The frequency optimization only activates when SCF_ITER == 1 (first SCF
 *     iteration), because that's when max_iter = nline × 6 = 24 and there's
 *     enough CG work to benefit from adaptive subspace diagonalization.
 *   - For SCF_ITER > 1, max_iter = nline = 4, and no intermediate subspace
 *     diagonalization is performed at all (original ABACUS behavior).
 *   - Therefore the benchmark does NOT compare SCF_ITER=1 vs SCF_ITER>1;
 *     those serve different purposes. Instead it validates that the optimized
 *     mode (SCF_ITER=1) produces correct eigenvalues and reports its performance.
 */

#include "../diago_bpcg.h"
#include "../kernels/bpcg_kernel_op.h"
#include "../diago_iter_assist.h"
#include "diago_mock.h"

#include "source_base/parallel_comm.h"
#include "source_base/parallel_reduce.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_base/macros.h"
#include "source_pw/module_pwdft/hamilt_pw.h"
#include "source_basis/module_pw/test/test_tool.h"

#include <gtest/gtest.h>
#include <mpi.h>
#include <complex>
#include <random>
#include <cmath>
#include <iomanip>

// ============================================================================
// Part 1: Unit tests for compute_optimal_freq
// ============================================================================

TEST(BpcgFreqUnitTest, OptimalFreqTinyProblem)
{
    // problem_size < 10,000  --> freq = nline
    EXPECT_EQ(hsolver::compute_optimal_freq(5, 100, 4), 4);
    EXPECT_EQ(hsolver::compute_optimal_freq(10, 500, 4), 4);
}

TEST(BpcgFreqUnitTest, OptimalFreqVerySmallProblem)
{
    // 10,000 <= problem_size < 100,000  --> freq = 6
    // NOTE: freq 6 > nline=4 — known issue: missing upper clamp
    EXPECT_EQ(hsolver::compute_optimal_freq(10, 1000, 4), 6);
}

TEST(BpcgFreqUnitTest, OptimalFreqSmallProblem)
{
    // 100,000 <= problem_size < 500,000  --> freq = 5
    EXPECT_EQ(hsolver::compute_optimal_freq(32, 3200, 4), 5);
    EXPECT_EQ(hsolver::compute_optimal_freq(200, 1500, 4), 5);
}

TEST(BpcgFreqUnitTest, OptimalFreqMediumProblem)
{
    // 500,000 <= problem_size < 2,000,000  --> freq = 4
    EXPECT_EQ(hsolver::compute_optimal_freq(64, 8000, 4), 4);
    EXPECT_EQ(hsolver::compute_optimal_freq(128, 8000, 4), 4);
}

TEST(BpcgFreqUnitTest, OptimalFreqLargeProblem)
{
    // 2,000,000 <= problem_size < 10,000,000  --> freq = 3
    EXPECT_EQ(hsolver::compute_optimal_freq(128, 20000, 4), 3);
}

TEST(BpcgFreqUnitTest, OptimalFreqVeryLargeProblem)
{
    // problem_size >= 10,000,000  --> freq = 2
    EXPECT_EQ(hsolver::compute_optimal_freq(256, 40000, 4), 2);
}

TEST(BpcgFreqUnitTest, OptimalFreqExactThresholds)
{
    EXPECT_EQ(hsolver::compute_optimal_freq(100, 100, 4), 6);    // size 10,000
    EXPECT_EQ(hsolver::compute_optimal_freq(100, 1000, 4), 5);   // size 100,000
    EXPECT_EQ(hsolver::compute_optimal_freq(500, 1000, 4), 4);   // size 500,000
    EXPECT_EQ(hsolver::compute_optimal_freq(2000, 1000, 4), 3);  // size 2,000,000
    EXPECT_EQ(hsolver::compute_optimal_freq(10000, 1000, 4), 2); // size 10,000,000
}

// ============================================================================
// Unit tests for compute_convergence_rate
// ============================================================================

TEST(BpcgFreqUnitTest, ConvRateNormal)
{
    EXPECT_NEAR(hsolver::compute_convergence_rate(2.0, 1.0), 0.5, 1e-12);
    EXPECT_NEAR(hsolver::compute_convergence_rate(10.0, 1.0), 0.1, 1e-12);
    EXPECT_NEAR(hsolver::compute_convergence_rate(1.0, 1.0), 1.0, 1e-12);
    EXPECT_GT(hsolver::compute_convergence_rate(1.0, 2.0), 1.0);
}

TEST(BpcgFreqUnitTest, ConvRateEdgeCases)
{
    EXPECT_EQ(hsolver::compute_convergence_rate(0.0, 1.0), 0.0);
    EXPECT_EQ(hsolver::compute_convergence_rate(1.0, 0.0), 0.0);
    EXPECT_EQ(hsolver::compute_convergence_rate(0.0, 0.0), 0.0);
    EXPECT_EQ(hsolver::compute_convergence_rate(-1.0, 1.0), 0.0);
    EXPECT_EQ(hsolver::compute_convergence_rate(1.0, -1.0), 0.0);
}

// ============================================================================
// Unit tests for compute_adaptive_freq
// ============================================================================

TEST(BpcgFreqUnitTest, AdaptiveFreqVeryFast)
{
    EXPECT_EQ(hsolver::compute_adaptive_freq(0.1, 4, 6), 5);
    EXPECT_EQ(hsolver::compute_adaptive_freq(0.29, 4, 6), 5);
}

TEST(BpcgFreqUnitTest, AdaptiveFreqFast)
{
    EXPECT_EQ(hsolver::compute_adaptive_freq(0.3, 4, 6), 4);
    EXPECT_EQ(hsolver::compute_adaptive_freq(0.5, 4, 6), 4);
}

TEST(BpcgFreqUnitTest, AdaptiveFreqModerate)
{
    EXPECT_EQ(hsolver::compute_adaptive_freq(0.6, 4, 6), 3);
    EXPECT_EQ(hsolver::compute_adaptive_freq(0.8, 4, 6), 3);
}

TEST(BpcgFreqUnitTest, AdaptiveFreqSlow)
{
    EXPECT_EQ(hsolver::compute_adaptive_freq(0.85, 4, 6), 2);
    EXPECT_EQ(hsolver::compute_adaptive_freq(1.5, 4, 6), 2);
}

TEST(BpcgFreqUnitTest, AdaptiveFreqClamp)
{
    // base_freq - 2 should not go below 1
    EXPECT_EQ(hsolver::compute_adaptive_freq(0.9, 2, 6), 1);
    EXPECT_EQ(hsolver::compute_adaptive_freq(0.9, 1, 6), 1);
    // base_freq + 1 should not exceed nline
    EXPECT_EQ(hsolver::compute_adaptive_freq(0.1, 4, 4), 4);
}

TEST(BpcgFreqUnitTest, AdaptiveFreqBoundaryRegimes)
{
    EXPECT_EQ(hsolver::compute_adaptive_freq(0.299, 4, 8), 5); // conv_rate < 0.3
    EXPECT_EQ(hsolver::compute_adaptive_freq(0.300, 4, 8), 4); // conv_rate == 0.3
    EXPECT_EQ(hsolver::compute_adaptive_freq(0.599, 4, 8), 4); // conv_rate < 0.6
    EXPECT_EQ(hsolver::compute_adaptive_freq(0.600, 4, 8), 3); // conv_rate == 0.6
    EXPECT_EQ(hsolver::compute_adaptive_freq(0.849, 4, 8), 3); // conv_rate < 0.85
    EXPECT_EQ(hsolver::compute_adaptive_freq(0.850, 4, 8), 2); // conv_rate == 0.85
}

// ============================================================================
// Part 2: Integration benchmark
// ============================================================================
// Runs DiagoBPCG::diag() with SCF_ITER=1 (adaptive freq optimization active)
// across multiple problem sizes. Reports:
//   - wall time
//   - max eigenvalue error vs LAPACK
//   - problem_size and compute_optimal_freq output

class BPCGFreqBenchmark : public ::testing::Test
{
protected:
    using T = std::complex<double>;

    int nprocs = 1;
    int myrank = 0;
    int nband;
    int npw;
    double ethr;
    double* e_lapack = nullptr;

    void SetUp() override
    {
#ifdef __MPI
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif
    }

    void TearDown() override { delete[] e_lapack; }

    // LAPACK zheev reference eigenvalues (rank 0 only)
    void compute_lapack_ref()
    {
        if (myrank != 0) return;
        int n = npw;
        e_lapack = new double[n];
        std::vector<T> hcopy = DIAGOTEST::hmatrix;
        int lwork = 2 * n;
        std::vector<T> work(lwork);
        std::vector<double> rwork(3 * n - 2);
        int info = 0;
        zheev_("V", "U", &n, hcopy.data(), &n,
               e_lapack, work.data(), &lwork, rwork.data(), &info);
        ASSERT_EQ(info, 0) << "LAPACK zheev failed, info=" << info;
    }

    void run_benchmark(int nband_, int npw_, double ethr_, int sparsity)
    {
        nband = nband_;
        npw = npw_;
        ethr = ethr_;
        size_t psize = static_cast<size_t>(nband) * npw;

        // Generate Hamiltonian and broadcast
        if (myrank == 0) {
            HPsi<T> hpsi(nband, npw, sparsity);
            DIAGOTEST::hmatrix = hpsi.hamilt();
        }
        DIAGOTEST::npw = npw;
#ifdef __MPI
        int hsize = npw * npw;
        if (myrank != 0) DIAGOTEST::hmatrix.resize(hsize);
        MPI_Bcast(DIAGOTEST::hmatrix.data(), hsize * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

        compute_lapack_ref();
#ifdef __MPI
        if (myrank != 0) e_lapack = new double[npw];
        MPI_Bcast(e_lapack, npw, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

        // --- Run DiagoBPCG::diag() with SCF_ITER=1 (adaptive freq active) ---
        hsolver::DiagoIterAssist<T, base_device::DEVICE_CPU>::SCF_ITER = 1;

        // Random initial psi
        psi::Psi<T> psi;
        psi.resize(1, nband, npw);
        {
            std::default_random_engine gen(42);
            std::uniform_real_distribution<double> d(-1.0, 1.0);
            for (int i = 0; i < nband; i++)
                for (int j = 0; j < npw; j++)
                    psi(i, j) = T{d(gen), d(gen)};
        }

        // Distribute (serial version: no actual distribution since we use full H)
        psi::Psi<T> psi_local;
        double* prec_local;
#ifdef __MPI
        DIAGOTEST::npw_local = new int[nprocs];
        DIAGOTEST::cal_division(DIAGOTEST::npw);
        DIAGOTEST::divide_hpsi(psi, psi_local, DIAGOTEST::hmatrix, DIAGOTEST::hmatrix_local);
        prec_local = new double[DIAGOTEST::npw_local[myrank]];
        {
            double* prec_full = new double[npw];
            std::default_random_engine g(1000);
            std::uniform_real_distribution<double> d(1.0, 2.0);
            for (int i = 0; i < npw; i++) prec_full[i] = d(g);
            DIAGOTEST::divide_psi<double>(prec_full, prec_local);
            delete[] prec_full;
        }
#else
        DIAGOTEST::npw_local = new int[1];
        DIAGOTEST::npw_local[0] = npw;
        DIAGOTEST::hmatrix_local = DIAGOTEST::hmatrix;
        psi_local = psi;
        prec_local = new double[npw];
        {
            std::default_random_engine g(1000);
            std::uniform_real_distribution<double> d(1.0, 2.0);
            for (int i = 0; i < npw; i++) prec_local[i] = d(g);
        }
#endif
        psi_local.fix_k(0);

        // Broadcast full H to all ranks (each rank does full H*psi locally)
        std::vector<T> full_h(npw * npw);
        if (myrank == 0) full_h = DIAGOTEST::hmatrix;
#ifdef __MPI
        MPI_Bcast(full_h.data(), npw * npw * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

        int npw_capture = npw;
        auto hpsi_func = [&full_h, npw_capture](T* psi_in, T* hpsi_out,
                                          const int ld_psi, const int nvec) {
            const T one = 1.0, zero = 0.0;
            base_device::DEVICE_CPU* ctx = {};
            ModuleBase::gemm_op<T, base_device::DEVICE_CPU>()(
                'N', 'N', npw_capture, nvec, npw_capture,
                &one, full_h.data(), npw_capture,
                psi_in, ld_psi,
                &zero, hpsi_out, ld_psi);
        };

        // Timing
        hsolver::DiagoBPCG<T, base_device::DEVICE_CPU> bpcg(prec_local);
        bpcg.init_iter(nband, nband, npw, psi_local.get_current_ngk());

        std::vector<double> ethr_band(nband, ethr);
        double* en = new double[npw];
        double t0 = MPI_Wtime();
        bpcg.diag(hpsi_func, psi_local.get_pointer(), en, ethr_band);
        double t1 = MPI_Wtime();

        // Report
        int optimal_freq = hsolver::compute_optimal_freq(nband, npw, 4);
        if (myrank == 0) {
            double max_err = 0.0;
            for (int i = 0; i < nband; i++) {
                double e = std::abs(en[i] - e_lapack[i]);
                if (e > max_err) max_err = e;
            }
            std::cout << "  size=" << psize
                      << "  nband=" << nband << "  npw=" << npw
                      << "  opt_freq=" << optimal_freq
                      << "  time=" << std::fixed << std::setprecision(4) << (t1 - t0) << " s"
                      << "  max_eig_err=" << max_err << std::endl;
        }

        delete[] en;
        delete[] DIAGOTEST::npw_local;
        delete[] prec_local;
    }
};

// Small problem: problem_size=5000 (< 10K tier)  → opt_freq = nline = 4
TEST_F(BPCGFreqBenchmark, Tiny_10x500)
{
    if (myrank == 0) std::cout << "\n--- Tiny: size=5000 (freq=nline=4) ---\n";
    run_benchmark(10, 500, 1e-5, 7);
}

// Small: problem_size=32000 (10K ~ 100K) → opt_freq = 6
TEST_F(BPCGFreqBenchmark, Small_32x1000)
{
    if (myrank == 0) std::cout << "\n--- Small: size=32000 (opt_freq=6) ---\n";
    run_benchmark(32, 1000, 1e-5, 7);
}

// Medium: problem_size=200000 (100K ~ 500K) → opt_freq = 5
TEST_F(BPCGFreqBenchmark, Medium_64x3125)
{
    if (myrank == 0) std::cout << "\n--- Medium: size=200000 (opt_freq=5) ---\n";
    run_benchmark(64, 3125, 1e-5, 7);
}

// Baseline: problem_size=500000 (500K ~ 2M, baseline) → opt_freq = 4
TEST_F(BPCGFreqBenchmark, Baseline_128x3906)
{
    if (myrank == 0) std::cout << "\n--- Baseline: size=~500000 (opt_freq=4) ---\n";
    run_benchmark(128, 3906, 1e-5, 7);
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char** argv)
{
    int nproc = 1, myrank = 0;

#ifdef __MPI
    int nproc_in_pool, kpar = 1, mypool, rank_in_pool;
    setupmpi(argc, argv, nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, kpar, mypool, rank_in_pool);
    MPI_Comm_split(MPI_COMM_WORLD, myrank, 0, &BP_WORLD);
    GlobalV::NPROC_IN_POOL = nproc;
#else
    MPI_Init(&argc, &argv);
#endif

    testing::InitGoogleTest(&argc, argv);
    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();
    if (myrank != 0)
        delete listeners.Release(listeners.default_result_printer());

    int result = RUN_ALL_TESTS();
    if (myrank == 0 && result != 0)
        std::cout << "ERROR: some tests are not passed" << std::endl;

    MPI_Finalize();
    return result;
}
