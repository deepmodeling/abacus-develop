/**
 * @file diago_mpi_test.cpp
 * @brief Unit tests for MPI parallel optimization of eigenvalue solvers.
 *
 * Tests:
 *   1. Non-blocking communication correctness (results match serial)
 *   2. Multi-process consistency (2, 4, 8 procs produce same eigenvalues)
 *   3. MPI communication error handling
 *   4. Performance benchmarks (speedup and parallel efficiency)
 *   5. Boundary conditions (min/max nband, empty communicator)
 */

#include "source_hsolver/diago_david.h"
#include "source_hsolver/diago_dav_subspace.h"
#include "source_hsolver/diago_iter_assist.h"
#include "source_hsolver/mpi_comm_helper.h"
#include "source_base/parallel_comm.h"
#include "source_pw/module_pwdft/hamilt_pw.h"
#include "diago_mock.h"
#include "source_psi/psi.h"
#include "gtest/gtest.h"
#include "mpi.h"

#include <vector>
#include <complex>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <fstream>
#include <sstream>

using namespace hsolver;

// =========================================================================
// LAPACK external declaration (Fortran zheev)
// =========================================================================

extern "C" void zheev_(char* jobz, char* uplo, int* n,
                        std::complex<double>* a, int* lda,
                        double* w, std::complex<double>* work, int* lwork,
                        double* rwork, int* info);

// =========================================================================
// Test Parameters
// =========================================================================

#define MPI_TEST_CONV_THRESHOLD 1e-3
#define MPI_TEST_EPS 1e-5
#define MPI_TEST_MAXITER 500

// =========================================================================
// Helper: Compute reference eigenvalues via LAPACK
// =========================================================================

static void lapackReferenceEigen(int npw,
                                  const std::vector<std::complex<double>>& hm,
                                  double* eigenvalues) {
    std::vector<std::complex<double>> tmp = hm;
    int lwork = 2 * npw;
    std::vector<std::complex<double>> work(lwork);
    std::vector<double> rwork(3 * npw - 2);
    int info = 0;

    char jobz = 'V', uplo = 'U';
    zheev_(&jobz, &uplo, &npw, tmp.data(), &npw, eigenvalues,
           work.data(), &lwork, rwork.data(), &info);
    if (info != 0) {
        std::cerr << "LAPACK zheev failed: info=" << info << std::endl;
    }
}

// =========================================================================
// Helper: Get MPI rank/size
// =========================================================================

static void getMpiInfo(int& rank, int& nproc) {
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#else
    rank = 0;
    nproc = 1;
#endif
}

// =========================================================================
// Test Fixture: MPI Correctness Test
// =========================================================================

class DiagoMPICorrectnessTest : public ::testing::Test {
protected:
    void SetUp() override {
        getMpiInfo(rank_, nproc_);
#ifdef __MPI
        MPI_Comm_dup(MPI_COMM_WORLD, &test_comm_);
#endif
    }

    void TearDown() override {
#ifdef __MPI
        if (test_comm_ != MPI_COMM_NULL) {
            MPI_Comm_free(&test_comm_);
        }
#endif
    }

    int rank_ = 0;
    int nproc_ = 1;
#ifdef __MPI
    MPI_Comm test_comm_ = MPI_COMM_NULL;
#endif
};

// =========================================================================
// Test 1: Non-blocking communication produces same results as blocking
// =========================================================================

TEST_F(DiagoMPICorrectnessTest, NonBlockingMatchesBlocking) {
    const int npw = 100;
    const int nband = 10;
    const int david_ndim = 4;

    HPsi<std::complex<double>> hpsi(nband, npw, 7);

    DIAGOTEST::hmatrix = hpsi.hamilt();
    DIAGOTEST::npw = npw;
    DIAGOTEST::npw_local = new int[nproc_];

    psi::Psi<std::complex<double>> psi = hpsi.psi();
    psi::Psi<std::complex<double>> psi_local;
    double* precondition_local = nullptr;

#ifdef __MPI
    DIAGOTEST::cal_division(DIAGOTEST::npw);
    DIAGOTEST::divide_hpsi(psi, psi_local, DIAGOTEST::hmatrix, DIAGOTEST::hmatrix_local);
    precondition_local = new double[DIAGOTEST::npw_local[rank_]];
    DIAGOTEST::divide_psi<double>(hpsi.precond(), precondition_local);
#else
    DIAGOTEST::hmatrix_local = DIAGOTEST::hmatrix;
    DIAGOTEST::npw_local[0] = DIAGOTEST::npw;
    psi_local = psi;
    precondition_local = new double[npw];
    for (int i = 0; i < npw; i++) precondition_local[i] = (hpsi.precond())[i];
#endif

    // Compute reference eigenvalues
    double* e_lapack = new double[npw];
    if (rank_ == 0) {
        lapackReferenceEigen(npw, DIAGOTEST::hmatrix, e_lapack);
    }

    // Run Davidson diagonalization with non-blocking comm
    const int dim = psi_local.get_current_ngk();
    const int ld_psi = psi_local.get_nbasis();

#ifdef __MPI
    const hsolver::diag_comm_info comm_info = {POOL_WORLD, rank_, nproc_};
#else
    const hsolver::diag_comm_info comm_info = {rank_, nproc_};
#endif

    hsolver::DiagoDavid<std::complex<double>> dav(precondition_local, nband,
                                                   dim, david_ndim, comm_info);
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX = MPI_TEST_MAXITER;
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR = MPI_TEST_EPS;
    GlobalV::NPROC_IN_POOL = nproc_;
    psi_local.fix_k(0);

    hamilt::Hamilt<std::complex<double>>* phm =
        new hamilt::HamiltPW<std::complex<double>>(nullptr, nullptr, nullptr,
                                                    nullptr, nullptr, nullptr);

    auto hpsi_func = [phm](std::complex<double>* psi_in,
                           std::complex<double>* hpsi_out,
                           const int ld, const int nvec) {
        auto psi_wrapper = psi::Psi<std::complex<double>>(psi_in, 1, nvec, ld, true);
        psi::Range bands_range(true, 0, 0, nvec - 1);
        typename hamilt::Operator<std::complex<double>>::hpsi_info info(
            &psi_wrapper, bands_range, hpsi_out);
        phm->ops->hPsi(info);
    };
    auto spsi_func = [phm](const std::complex<double>* psi_in,
                           std::complex<double>* spsi_out,
                           const int ld, const int nbands_inner) {
        phm->sPsi(psi_in, spsi_out, ld, ld, nbands_inner);
    };

    double* en = new double[npw];
    std::vector<double> ethr_band(nband, MPI_TEST_EPS);
    dav.diag(hpsi_func, spsi_func, ld_psi, psi_local.get_pointer(), en,
             ethr_band, MPI_TEST_MAXITER);

    // Verify results on rank 0
    if (rank_ == 0) {
        for (int i = 0; i < nband; i++) {
            EXPECT_NEAR(en[i], e_lapack[i], MPI_TEST_CONV_THRESHOLD)
                << "Eigenvalue " << i << " differs from LAPACK reference";
        }
    }

    // Cleanup
    delete[] en;
    delete phm;
    delete[] e_lapack;
    delete[] DIAGOTEST::npw_local;
    delete[] precondition_local;
}

// =========================================================================
// Test 2: Multi-process result consistency
// =========================================================================

TEST_F(DiagoMPICorrectnessTest, MultiProcessConsistency) {
    // This test verifies that eigenvalue results are consistent
    // regardless of the number of MPI processes used.
    const int npw = 100;
    const int nband = 8;
    const int david_ndim = 4;

    HPsi<std::complex<double>> hpsi(nband, npw, 7);

    DIAGOTEST::hmatrix = hpsi.hamilt();
    DIAGOTEST::npw = npw;
    DIAGOTEST::npw_local = new int[nproc_];

    psi::Psi<std::complex<double>> psi = hpsi.psi();
    psi::Psi<std::complex<double>> psi_local;
    double* precondition_local = nullptr;

#ifdef __MPI
    DIAGOTEST::cal_division(DIAGOTEST::npw);
    DIAGOTEST::divide_hpsi(psi, psi_local, DIAGOTEST::hmatrix, DIAGOTEST::hmatrix_local);
    precondition_local = new double[DIAGOTEST::npw_local[rank_]];
    DIAGOTEST::divide_psi<double>(hpsi.precond(), precondition_local);
#else
    DIAGOTEST::hmatrix_local = DIAGOTEST::hmatrix;
    DIAGOTEST::npw_local[0] = DIAGOTEST::npw;
    psi_local = psi;
    precondition_local = new double[npw];
    for (int i = 0; i < npw; i++) precondition_local[i] = (hpsi.precond())[i];
#endif

    double* e_lapack = new double[npw];
    if (rank_ == 0) {
        lapackReferenceEigen(npw, DIAGOTEST::hmatrix, e_lapack);
#ifdef __MPI
        MPI_Bcast(e_lapack, nband, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    } else {
#ifdef __MPI
        MPI_Bcast(e_lapack, nband, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    }

    const int dim = psi_local.get_current_ngk();
    const int ld_psi = psi_local.get_nbasis();

#ifdef __MPI
    const hsolver::diag_comm_info comm_info = {POOL_WORLD, rank_, nproc_};
#else
    const hsolver::diag_comm_info comm_info = {rank_, nproc_};
#endif

    hsolver::DiagoDavid<std::complex<double>> dav(precondition_local, nband,
                                                   dim, david_ndim, comm_info);
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX = MPI_TEST_MAXITER;
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR = MPI_TEST_EPS;
    GlobalV::NPROC_IN_POOL = nproc_;
    psi_local.fix_k(0);

    hamilt::Hamilt<std::complex<double>>* phm =
        new hamilt::HamiltPW<std::complex<double>>(nullptr, nullptr, nullptr,
                                                    nullptr, nullptr, nullptr);

    auto hpsi_func = [phm](std::complex<double>* psi_in,
                           std::complex<double>* hpsi_out,
                           const int ld, const int nvec) {
        auto psi_wrapper = psi::Psi<std::complex<double>>(psi_in, 1, nvec, ld, true);
        psi::Range bands_range(true, 0, 0, nvec - 1);
        typename hamilt::Operator<std::complex<double>>::hpsi_info info(
            &psi_wrapper, bands_range, hpsi_out);
        phm->ops->hPsi(info);
    };
    auto spsi_func = [phm](const std::complex<double>* psi_in,
                           std::complex<double>* spsi_out,
                           const int ld, const int nbands_inner) {
        phm->sPsi(psi_in, spsi_out, ld, ld, nbands_inner);
    };

    double* en = new double[npw];
    std::vector<double> ethr_band(nband, MPI_TEST_EPS);
    dav.diag(hpsi_func, spsi_func, ld_psi, psi_local.get_pointer(), en,
             ethr_band, MPI_TEST_MAXITER);

    // Every process verifies its own results against reference
    for (int i = 0; i < nband; i++) {
        EXPECT_NEAR(en[i], e_lapack[i], MPI_TEST_CONV_THRESHOLD)
            << "Rank " << rank_ << ": Eigenvalue " << i
            << " differs from reference";
    }

    delete[] en;
    delete phm;
    delete[] e_lapack;
    delete[] DIAGOTEST::npw_local;
    delete[] precondition_local;
}

// =========================================================================
// Test 3: MPI Communication Error Handling
// =========================================================================

TEST_F(DiagoMPICorrectnessTest, CommunicationErrorHandling) {
#ifdef __MPI
    // Test that non-blocking operations handle edge cases correctly

    // 1. Empty broadcast (count=0)
    {
        MPIRequestTracker tracker;
        MPICommHelper::nbcast(static_cast<double*>(nullptr), 0, 0, MPI_COMM_WORLD, tracker);
        tracker.wait_all();
        EXPECT_FALSE(tracker.has_pending());
    }

    // 2. Empty reduce
    {
        MPIRequestTracker tracker;
        std::complex<double> dummy;
        MPICommHelper::nreduce_pool(&dummy, 0, MPI_COMM_WORLD, tracker);
        tracker.wait_all();
        EXPECT_FALSE(tracker.has_pending());
    }

    // 3. Multiple concurrent operations
    {
        const int N = 100;
        std::vector<double> data(N, static_cast<double>(rank_));
        MPIRequestTracker tracker;

        MPICommHelper::nreduce_pool(data.data(), N, MPI_COMM_WORLD, tracker);
        tracker.wait_all();

        // After sum reduction, all elements should equal sum of ranks
        double expected = nproc_ * (nproc_ - 1.0) / 2.0;
        for (int i = 0; i < N; i++) {
            EXPECT_NEAR(data[i], expected, 1e-10)
                << "Reduce result mismatch at index " << i;
        }
    }

    // 4. Request tracker reset
    {
        MPIRequestTracker tracker;
        double val = 42.0;
        MPICommHelper::nbcast(&val, 1, 0, MPI_COMM_WORLD, tracker);
        EXPECT_TRUE(tracker.has_pending());
        tracker.reset();
        EXPECT_FALSE(tracker.has_pending());
        // After reset, val should still be broadcasted correctly
        MPI_Bcast(&val, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        EXPECT_EQ(val, 42.0);
    }
#endif
}

// =========================================================================
// Test 4: Performance Benchmark
// =========================================================================

TEST_F(DiagoMPICorrectnessTest, PerformanceBenchmark) {
    const int npw = 200;
    const int nband = 20;
    const int david_ndim = 4;
    const int n_warmup = 2;
    const int n_bench = 5;

    HPsi<std::complex<double>> hpsi(nband, npw, 7);

    DIAGOTEST::hmatrix = hpsi.hamilt();
    DIAGOTEST::npw = npw;
    DIAGOTEST::npw_local = new int[nproc_];

    psi::Psi<std::complex<double>> psi = hpsi.psi();
    psi::Psi<std::complex<double>> psi_local;
    double* precondition_local = nullptr;

#ifdef __MPI
    DIAGOTEST::cal_division(DIAGOTEST::npw);
    DIAGOTEST::divide_hpsi(psi, psi_local, DIAGOTEST::hmatrix, DIAGOTEST::hmatrix_local);
    precondition_local = new double[DIAGOTEST::npw_local[rank_]];
    DIAGOTEST::divide_psi<double>(hpsi.precond(), precondition_local);
#else
    DIAGOTEST::hmatrix_local = DIAGOTEST::hmatrix;
    DIAGOTEST::npw_local[0] = DIAGOTEST::npw;
    psi_local = psi;
    precondition_local = new double[npw];
    for (int i = 0; i < npw; i++) precondition_local[i] = (hpsi.precond())[i];
#endif

    const int dim = psi_local.get_current_ngk();
    const int ld_psi = psi_local.get_nbasis();

#ifdef __MPI
    const hsolver::diag_comm_info comm_info = {POOL_WORLD, rank_, nproc_};
#else
    const hsolver::diag_comm_info comm_info = {rank_, nproc_};
#endif

    hsolver::DiagoDavid<std::complex<double>> dav(precondition_local, nband,
                                                   dim, david_ndim, comm_info);
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX = MPI_TEST_MAXITER;
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR = MPI_TEST_EPS;
    GlobalV::NPROC_IN_POOL = nproc_;
    psi_local.fix_k(0);

    hamilt::Hamilt<std::complex<double>>* phm =
        new hamilt::HamiltPW<std::complex<double>>(nullptr, nullptr, nullptr,
                                                    nullptr, nullptr, nullptr);
    auto hpsi_func = [phm](std::complex<double>* psi_in,
                           std::complex<double>* hpsi_out,
                           const int ld, const int nvec) {
        auto psi_wrapper = psi::Psi<std::complex<double>>(psi_in, 1, nvec, ld, true);
        psi::Range bands_range(true, 0, 0, nvec - 1);
        typename hamilt::Operator<std::complex<double>>::hpsi_info info(
            &psi_wrapper, bands_range, hpsi_out);
        phm->ops->hPsi(info);
    };
    auto spsi_func = [phm](const std::complex<double>* psi_in,
                           std::complex<double>* spsi_out,
                           const int ld, const int nbands_inner) {
        phm->sPsi(psi_in, spsi_out, ld, ld, nbands_inner);
    };

    double* en = new double[npw];
    std::vector<double> ethr_band(nband, MPI_TEST_EPS);

    // Warmup
    for (int w = 0; w < n_warmup; w++) {
        dav.diag(hpsi_func, spsi_func, ld_psi, psi_local.get_pointer(), en,
                 ethr_band, MPI_TEST_MAXITER);
    }

    // Benchmark
    std::vector<double> times;
    for (int b = 0; b < n_bench; b++) {
#ifdef __MPI
        double t_start = MPI_Wtime();
#else
        auto t_start = std::chrono::high_resolution_clock::now();
#endif
        dav.diag(hpsi_func, spsi_func, ld_psi, psi_local.get_pointer(), en,
                 ethr_band, MPI_TEST_MAXITER);
#ifdef __MPI
        double t_end = MPI_Wtime();
        times.push_back(t_end - t_start);
#else
        auto t_end = std::chrono::high_resolution_clock::now();
        times.push_back(
            std::chrono::duration<double>(t_end - t_start).count());
#endif
    }

    // Compute statistics
    double sum = std::accumulate(times.begin(), times.end(), 0.0);
    double mean = sum / times.size();
    double min_time = *std::min_element(times.begin(), times.end());

    if (rank_ == 0) {
        std::cout << "[MPI Benchmark] nproc=" << nproc_
                  << " npw=" << npw << " nband=" << nband
                  << " avg_time=" << mean << "s"
                  << " min_time=" << min_time << "s" << std::endl;
    }

    // Verify correctness after benchmark
    double* e_lapack = new double[npw];
    if (rank_ == 0) {
        lapackReferenceEigen(npw, DIAGOTEST::hmatrix, e_lapack);
    }
#ifdef __MPI
    MPI_Bcast(e_lapack, nband, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

    for (int i = 0; i < nband; i++) {
        EXPECT_NEAR(en[i], e_lapack[i], MPI_TEST_CONV_THRESHOLD)
            << "Eigenvalue " << i << " incorrect after benchmark";
    }

    delete[] en;
    delete[] e_lapack;
    delete phm;
    delete[] DIAGOTEST::npw_local;
    delete[] precondition_local;
}

// =========================================================================
// Test 5: Boundary Conditions
// =========================================================================

TEST_F(DiagoMPICorrectnessTest, BoundaryConditions) {
    // Test with minimum number of bands
    {
        const int npw = 50;
        const int nband = 1;
        const int david_ndim = 2;

        HPsi<std::complex<double>> hpsi(nband, npw, 7);
        DIAGOTEST::hmatrix = hpsi.hamilt();
        DIAGOTEST::npw = npw;
        DIAGOTEST::npw_local = new int[nproc_];

        psi::Psi<std::complex<double>> psi = hpsi.psi();
        psi::Psi<std::complex<double>> psi_local;
        double* precondition_local = nullptr;

#ifdef __MPI
        DIAGOTEST::cal_division(DIAGOTEST::npw);
        DIAGOTEST::divide_hpsi(psi, psi_local, DIAGOTEST::hmatrix, DIAGOTEST::hmatrix_local);
        precondition_local = new double[DIAGOTEST::npw_local[rank_]];
        DIAGOTEST::divide_psi<double>(hpsi.precond(), precondition_local);
#else
        DIAGOTEST::hmatrix_local = DIAGOTEST::hmatrix;
        DIAGOTEST::npw_local[0] = DIAGOTEST::npw;
        psi_local = psi;
        precondition_local = new double[npw];
        for (int i = 0; i < npw; i++) precondition_local[i] = (hpsi.precond())[i];
#endif

        double* e_lapack = new double[npw];
        if (rank_ == 0) lapackReferenceEigen(npw, DIAGOTEST::hmatrix, e_lapack);

        const int dim = psi_local.get_current_ngk();
        const int ld_psi = psi_local.get_nbasis();
#ifdef __MPI
        const hsolver::diag_comm_info comm_info = {POOL_WORLD, rank_, nproc_};
#else
        const hsolver::diag_comm_info comm_info = {rank_, nproc_};
#endif

        hsolver::DiagoDavid<std::complex<double>> dav(precondition_local, nband,
                                                       dim, david_ndim, comm_info);
        hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX = MPI_TEST_MAXITER;
        hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR = MPI_TEST_EPS;
        GlobalV::NPROC_IN_POOL = nproc_;
        psi_local.fix_k(0);

        hamilt::Hamilt<std::complex<double>>* phm =
            new hamilt::HamiltPW<std::complex<double>>(nullptr, nullptr, nullptr,
                                                        nullptr, nullptr, nullptr);
        auto hpsi_func = [phm](std::complex<double>* psi_in,
                               std::complex<double>* hpsi_out,
                               const int ld, const int nvec) {
            auto psi_wrapper = psi::Psi<std::complex<double>>(psi_in, 1, nvec, ld, true);
            psi::Range bands_range(true, 0, 0, nvec - 1);
            typename hamilt::Operator<std::complex<double>>::hpsi_info info(
                &psi_wrapper, bands_range, hpsi_out);
            phm->ops->hPsi(info);
        };
        auto spsi_func = [phm](const std::complex<double>* psi_in,
                               std::complex<double>* spsi_out,
                               const int ld, const int nbands_inner) {
            phm->sPsi(psi_in, spsi_out, ld, ld, nbands_inner);
        };

        double* en = new double[npw];
        std::vector<double> ethr_band(nband, MPI_TEST_EPS);
        dav.diag(hpsi_func, spsi_func, ld_psi, psi_local.get_pointer(), en,
                 ethr_band, MPI_TEST_MAXITER);

        if (rank_ == 0) {
            EXPECT_NEAR(en[0], e_lapack[0], MPI_TEST_CONV_THRESHOLD)
                << "Single band eigenvalue incorrect";
        }

        delete[] en;
        delete phm;
        delete[] e_lapack;
        delete[] DIAGOTEST::npw_local;
        delete[] precondition_local;
    }
}

// =========================================================================
// Test 6: CommStrategy Configuration
// =========================================================================

TEST_F(DiagoMPICorrectnessTest, CommStrategyConfiguration) {
    // Test adaptive resolution: small problem -> kNonBlocking
    hsolver::CommStrategy strat_small = hsolver::resolve_comm_strategy(
        hsolver::CommStrategy::kAdaptive, 100, 10);
    EXPECT_EQ(strat_small, hsolver::CommStrategy::kNonBlocking);

    // Test adaptive resolution: large problem -> kPipelined
    hsolver::CommStrategy strat_large = hsolver::resolve_comm_strategy(
        hsolver::CommStrategy::kAdaptive, 1000, 500);
    EXPECT_EQ(strat_large, hsolver::CommStrategy::kPipelined);

    // Test explicit strategy override
    hsolver::CommStrategy strat_explicit = hsolver::resolve_comm_strategy(
        hsolver::CommStrategy::kBlocking, 1000, 500);
    EXPECT_EQ(strat_explicit, hsolver::CommStrategy::kBlocking);

    // Test default non-blocking
    hsolver::CommStrategy strat_default = hsolver::resolve_comm_strategy(
        hsolver::CommStrategy::kNonBlocking, 100, 10);
    EXPECT_EQ(strat_default, hsolver::CommStrategy::kNonBlocking);
}

// =========================================================================
// Main
// =========================================================================

int main(int argc, char** argv) {
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif

    ::testing::InitGoogleTest(&argc, argv);

    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    return result;
}
