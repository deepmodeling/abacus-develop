#include "source_base/inverse_matrix.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_pw/module_pwdft/structure_factor.h"
#include "source_psi/psi.h"
#include "source_hamilt/hamilt.h"
#include "source_pw/module_pwdft/hamilt_pw.h"
#include "../diago_iter_assist.h"
#include "../diago_ppcg.h"
#include "../diago_bpcg.h"
#include "../diago_david.h"
#include "../diag_comm_info.h"
#include "diago_mock.h"
#include "mpi.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_comm.h"
#include "source_basis/module_pw/test/test_tool.h"

#include <gtest/gtest.h>
#include <complex>
#include <fstream>
#include <random>

// LAPACK reference eigenvalues for comparison
static void lapackEigen(int& npw, std::vector<std::complex<double>>& hm, double* e)
{
    int lwork = 2 * npw;
    std::complex<double>* work2 = new std::complex<double>[lwork];
    double* rwork = new double[3 * npw - 2];
    int info = 0;
    char jobz = 'V', uplo = 'U';
    zheev_(&jobz, &uplo, &npw, hm.data(), &npw, e, work2, &lwork, rwork, &info);
    delete[] rwork;
    delete[] work2;
}

class DiagoPPCGPrepare
{
  public:
    DiagoPPCGPrepare(int nband, int npw, int sparsity, double eps, int maxiter, double threshold)
        : nband(nband), npw(npw), sparsity(sparsity), eps(eps), maxiter(maxiter), threshold(threshold)
    {
#ifdef __MPI
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
#endif
    }

    int nband = 0;
    int npw = 0;
    int sparsity = 0;
    double eps = 1e-6;
    int maxiter = 200;
    double threshold = 5e-2;

    int nprocs = 1;
    int mypnum = 0;

    void CompareEigen(double* precondition)
    {
        // Reference by LAPACK
        double* e_lapack = new double[npw];
        auto ev = DIAGOTEST::hmatrix;
        if (mypnum == 0)
        {
            lapackEigen(npw, ev, e_lapack);
        }
#ifdef __MPI
        MPI_Bcast(e_lapack, npw, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

        // Initial guess: random combination of Lapack eigenvectors
        ModuleBase::ComplexMatrix psiguess(nband, npw);
        std::default_random_engine p(1);
        std::uniform_int_distribution<unsigned> u(1, 10);
        for (int i = 0; i < nband; i++)
        {
            for (int j = 0; j < npw; j++)
            {
                double rand = static_cast<double>(u(p)) / 10.;
                psiguess(i, j) = ev[j * DIAGOTEST::h_nc + i] * rand;
            }
        }

        // Prepare psi
        double* en = new double[npw];
        int ik = 1;
        psi::Psi<std::complex<double>> psi;
        psi.resize(ik, nband, npw);
        for (int i = 0; i < nband; i++)
        {
            for (int j = 0; j < npw; j++)
            {
                psi(i, j) = psiguess(i, j);
            }
        }

        psi::Psi<std::complex<double>> psi_local;
        double* precondition_local;
        DIAGOTEST::npw_local = new int[nprocs];
#ifdef __MPI
        DIAGOTEST::cal_division(DIAGOTEST::npw);
        DIAGOTEST::divide_hpsi(psi, psi_local, DIAGOTEST::hmatrix, DIAGOTEST::hmatrix_local);
        precondition_local = new double[DIAGOTEST::npw_local[mypnum]];
        DIAGOTEST::divide_psi<double>(precondition, precondition_local);
#else
        DIAGOTEST::hmatrix_local = DIAGOTEST::hmatrix;
        DIAGOTEST::npw_local[0] = DIAGOTEST::npw;
        psi_local = psi;
        precondition_local = new double[DIAGOTEST::npw];
        for (int i = 0; i < DIAGOTEST::npw; i++)
            precondition_local[i] = precondition[i];
#endif

        hsolver::DiagoPPCG<std::complex<double>> ppcg(precondition_local);
        psi_local.fix_k(0);

        using T = std::complex<double>;
        const int dim = DIAGOTEST::npw;
        const std::vector<T>& h_mat = DIAGOTEST::hmatrix_local;
        auto hpsi_func = [h_mat, dim](T* psi_in, T* hpsi_out, const int ld_psi, const int nvec) {
            auto one = std::make_unique<T>(1.0);
            auto zero = std::make_unique<T>(0.0);
            const T* one_ = one.get();
            const T* zero_ = zero.get();

            base_device::DEVICE_CPU* ctx = {};
            ModuleBase::gemm_op<T, base_device::DEVICE_CPU>()('N',
                                                              'N',
                                                              dim,
                                                              nvec,
                                                              dim,
                                                              one_,
                                                              h_mat.data(),
                                                              dim,
                                                              psi_in,
                                                              ld_psi,
                                                              zero_,
                                                              hpsi_out,
                                                              ld_psi);
        };

        const int ndim = psi_local.get_current_ngk();
        ppcg.init_iter(nband, nband, npw, ndim);
        std::vector<double> ethr_band(nband, 1e-6);

        // As in BPCG, several diag() passes are needed for harder problems.
        // Each pass starts from the refined X of the previous one and rebuilds
        // the search directions from scratch.
        for (int pass = 0; pass < 5; ++pass)
        {
            ppcg.diag(hpsi_func, psi_local.get_pointer(), en, ethr_band);
        }

        delete[] DIAGOTEST::npw_local;
        delete[] precondition_local;

        for (int i = 0; i < nband; i++)
        {
            EXPECT_NEAR(en[i], e_lapack[i], threshold);
        }

        delete[] en;
        delete[] e_lapack;
    }
};

class DiagoPPCGTest : public ::testing::TestWithParam<DiagoPPCGPrepare>
{
};

TEST_P(DiagoPPCGTest, RandomHamilt)
{
    DiagoPPCGPrepare dcp = GetParam();
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX = dcp.maxiter;
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR = dcp.eps;
    hsolver::DiagoIterAssist<std::complex<double>>::SCF_ITER = 1;

    HPsi<std::complex<double>> hpsi(dcp.nband, dcp.npw, dcp.sparsity);
    DIAGOTEST::hmatrix = hpsi.hamilt();
    DIAGOTEST::npw = dcp.npw;

    dcp.CompareEigen(hpsi.precond());
}

INSTANTIATE_TEST_SUITE_P(VerifyPPCG,
                         DiagoPPCGTest,
                         ::testing::Values(
                             // nband, npw, sparsity, eps, maxiter, threshold
                             DiagoPPCGPrepare(6, 120, 0, 1e-6, 200, 8e-2)));

TEST(DiagoPPCGTest, TwoByTwo)
{
    const int dim = 2;
    const int nband = 2;
    std::vector<std::complex<double>> hm(dim * dim);
    hm[0] = {4.0, 0.0};
    hm[1] = {1.0, 0.0};
    hm[2] = {1.0, 0.0};
    hm[3] = {3.0, 0.0};

    DiagoPPCGPrepare dcp(nband, dim, 0, 1e-10, 80, 1e-10);
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX = dcp.maxiter;
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR = dcp.eps;
    hsolver::DiagoIterAssist<std::complex<double>>::SCF_ITER = 1;

    double precond[dim] = {1.0, 1.0};
    DIAGOTEST::hmatrix = hm;
    DIAGOTEST::npw = dim;
    dcp.CompareEigen(precond);
}

TEST(DiagoPPCGTest, readH)
{
    std::vector<std::complex<double>> hm;
    std::ifstream ifs;
    std::string filename = "H-KPoints-Si2.dat";
    ifs.open(filename);
    if (!ifs.is_open())
    {
        std::cout << "Error opening file " << filename << std::endl;
        exit(1);
    }
    DIAGOTEST::readh(ifs, hm);
    ifs.close();

    int dim = DIAGOTEST::npw;
    int nband = 10;

    DiagoPPCGPrepare dcp(nband, dim, 0, 1e-6, 500, 2e-1);
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX = dcp.maxiter;
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR = dcp.eps;
    hsolver::DiagoIterAssist<std::complex<double>>::SCF_ITER = 1;

    HPsi<std::complex<double>> hpsi;
    hpsi.create(nband, dim);
    DIAGOTEST::hmatrix = hm;
    DIAGOTEST::npw = dim;
    dcp.CompareEigen(hpsi.precond());
}

// ------------------------------------------------------------
// Consistency tests: PPCG  vs  BPCG  on the same Hamiltonian
// ------------------------------------------------------------
TEST(DiagoPPCGTest, ConsistentWithBPCG)
{
    int dim = 40;
    int nband = 8;

    HPsi<std::complex<double>> hpsi(nband, dim, 5); // moderate sparsity
    DIAGOTEST::hmatrix = hpsi.hamilt();
    DIAGOTEST::npw = dim;

    // LAPACK reference
    double* e_lapack = new double[dim];
    auto ev = DIAGOTEST::hmatrix;
    lapackEigen(dim, ev, e_lapack);

    // --- shared initial guess ---
    ModuleBase::ComplexMatrix psiguess(nband, dim);
    std::default_random_engine p(7);
    std::uniform_int_distribution<unsigned> u(1, 10);
    for (int i = 0; i < nband; i++)
        for (int j = 0; j < dim; j++)
            psiguess(i, j) = ev[j * DIAGOTEST::h_nc + i] * static_cast<double>(u(p)) / 10.;

    // --- PPCG ---
    {
        psi::Psi<std::complex<double>> psi_ppcg;
        psi_ppcg.resize(1, nband, dim);
        for (int i = 0; i < nband; i++)
            for (int j = 0; j < dim; j++)
                psi_ppcg(i, j) = psiguess(i, j);

        double en_ppcg[40] = {};
        hsolver::DiagoPPCG<std::complex<double>> ppcg(hpsi.precond());
        ppcg.init_iter(nband, nband, dim, dim);
        std::vector<double> ethr(nband, 1e-6);
        hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX = 200;
        hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR = 1e-6;
        hsolver::DiagoIterAssist<std::complex<double>>::SCF_ITER = 1;

        auto hpsi_func = [&hpsi, dim](std::complex<double>* in, std::complex<double>* out, int ld, int nv) {
            auto one = std::make_unique<std::complex<double>>(1.0);
            auto zero = std::make_unique<std::complex<double>>(0.0);
            ModuleBase::gemm_op<std::complex<double>, base_device::DEVICE_CPU>()(
                'N', 'N', dim, nv, dim, one.get(), DIAGOTEST::hmatrix.data(), dim, in, ld, zero.get(), out, ld);
        };

        for (int pass = 0; pass < 5; ++pass)
            ppcg.diag(hpsi_func, psi_ppcg.get_pointer(), en_ppcg, ethr);

        for (int i = 0; i < nband; i++)
            EXPECT_NEAR(en_ppcg[i], e_lapack[i], 1e-1);
    }

    // --- BPCG ---
    {
        psi::Psi<std::complex<double>> psi_bpcg;
        psi_bpcg.resize(1, nband, dim);
        for (int i = 0; i < nband; i++)
            for (int j = 0; j < dim; j++)
                psi_bpcg(i, j) = psiguess(i, j);

        double en_bpcg[40] = {};
        hsolver::DiagoBPCG<std::complex<double>> bpcg(hpsi.precond());
        bpcg.init_iter(nband, nband, dim, dim);
        std::vector<double> ethr(nband, 1e-6);

        auto hpsi_func = [&hpsi, dim](std::complex<double>* in, std::complex<double>* out, int ld, int nv) {
            auto one = std::make_unique<std::complex<double>>(1.0);
            auto zero = std::make_unique<std::complex<double>>(0.0);
            ModuleBase::gemm_op<std::complex<double>, base_device::DEVICE_CPU>()(
                'N', 'N', dim, nv, dim, one.get(), DIAGOTEST::hmatrix.data(), dim, in, ld, zero.get(), out, ld);
        };

        for (int pass = 0; pass < 4; ++pass)
            bpcg.diag(hpsi_func, psi_bpcg.get_pointer(), en_bpcg, ethr);

        for (int i = 0; i < nband; i++)
            EXPECT_NEAR(en_bpcg[i], e_lapack[i], 1e-1);
    }

    delete[] e_lapack;
}

// ------------------------------------------------------------
// Parameter configurability test
// ------------------------------------------------------------
TEST(DiagoPPCGTest, TunableParameters)
{
    int dim = 30;
    int nband = 5;
    HPsi<std::complex<double>> hpsi(nband, dim);
    DIAGOTEST::hmatrix = hpsi.hamilt();
    DIAGOTEST::npw = dim;

    // LAPACK reference
    double* e_lapack = new double[dim];
    auto ev = DIAGOTEST::hmatrix;
    lapackEigen(dim, ev, e_lapack);

    ModuleBase::ComplexMatrix psiguess(nband, dim);
    std::default_random_engine p(3);
    std::uniform_int_distribution<unsigned> u(1, 10);
    for (int i = 0; i < nband; i++)
        for (int j = 0; j < dim; j++)
            psiguess(i, j) = ev[j * DIAGOTEST::h_nc + i] * static_cast<double>(u(p)) / 10.;

    // --- test 1: aggressive p_safe margin (margin=5) ---
    {
        psi::Psi<std::complex<double>> psi_a;
        psi_a.resize(1, nband, dim);
        for (int i = 0; i < nband; i++)
            for (int j = 0; j < dim; j++)
                psi_a(i, j) = psiguess(i, j);

        double en_a[30] = {};
        hsolver::DiagoPPCG<std::complex<double>> ppcg(hpsi.precond());
        ppcg.init_iter(nband, nband, dim, dim);
        ppcg.set_p_safe_margin(5);   // more conservative → P block disabled for this problem
        ppcg.set_max_inner_iter(1);  // single iteration per pass
        ppcg.set_npass(8);           // compensate with more passes

        std::vector<double> ethr(nband, 1e-6);
        hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX = 200;
        hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR = 1e-6;
        hsolver::DiagoIterAssist<std::complex<double>>::SCF_ITER = 1;

        auto hpsi_func = [&hpsi, dim](std::complex<double>* in, std::complex<double>* out, int ld, int nv) {
            auto one = std::make_unique<std::complex<double>>(1.0);
            auto zero = std::make_unique<std::complex<double>>(0.0);
            ModuleBase::gemm_op<std::complex<double>, base_device::DEVICE_CPU>()(
                'N', 'N', dim, nv, dim, one.get(), DIAGOTEST::hmatrix.data(), dim, in, ld, zero.get(), out, ld);
        };

        for (int pass = 0; pass < ppcg.npass(); ++pass)
            ppcg.diag(hpsi_func, psi_a.get_pointer(), en_a, ethr);

        for (int i = 0; i < nband; i++)
            EXPECT_NEAR(en_a[i], e_lapack[i], 2e-1);
    }

    delete[] e_lapack;
}

// ------------------------------------------------------------
// Comprehensive benchmark: PPCG vs BPCG vs LAPACK
//   - 5 matrix sizes (60 … 480)
//   - timing, accuracy, speedup vs LAPACK
//   - empirical complexity exponents
// ------------------------------------------------------------
TEST(DiagoPPCGTest, ComprehensiveBenchmark)
{
    const int nband = 6;
    const int sizes[] = {60, 120, 240, 360, 480};
    const int n_sizes = 5;

    // storage for timing data (for complexity analysis)
    double t_ppcg[5] = {}, t_bpcg[5] = {}, t_david[5] = {}, t_lapack[5] = {};

    std::cout << "\n"
              << "==========================================================================================\n"
              << "  PPCG vs BPCG vs Davidson vs LAPACK  —  Comprehensive Benchmark\n"
              << "  (nband=" << nband << ", 5 passes each, ethr=1e-5)\n"
              << "==========================================================================================\n"
              << "   N    | PPCG(ms) BPCG(ms) David(ms) LAPACK(ms) | PPCG/LAP BPCG/LAP David/LAP | PPCG-err  BPCG-err  David-err\n"
              << "--------+------------------------------------------+---------------------------+----------------------------"
              << std::endl;

    for (int sz = 0; sz < n_sizes; ++sz)
    {
        int npw = sizes[sz];
        HPsi<std::complex<double>> hpsi(nband, npw);
        DIAGOTEST::hmatrix = hpsi.hamilt();
        DIAGOTEST::npw = npw;

        // LAPACK reference (timed)
        double* e_lapack = new double[npw];
        auto ev_lap = DIAGOTEST::hmatrix;
        double t0 = MPI_Wtime();
        lapackEigen(npw, ev_lap, e_lapack);
        t_lapack[sz] = (MPI_Wtime() - t0) * 1000.0;

        // common initial guess
        ModuleBase::ComplexMatrix psiguess(nband, npw);
        std::default_random_engine prng(5);
        std::uniform_int_distribution<unsigned> u(1, 10);
        for (int i = 0; i < nband; i++)
            for (int j = 0; j < npw; j++)
                psiguess(i, j) = ev_lap[j * DIAGOTEST::h_nc + i] * static_cast<double>(u(prng)) / 10.;

        // shared hpsi_func
        auto hpsi_func = [npw](std::complex<double>* in, std::complex<double>* out, int ld, int nv) {
            auto one = std::make_unique<std::complex<double>>(1.0);
            auto zero = std::make_unique<std::complex<double>>(0.0);
            ModuleBase::gemm_op<std::complex<double>, base_device::DEVICE_CPU>()(
                'N', 'N', npw, nv, npw, one.get(), DIAGOTEST::hmatrix.data(), npw, in, ld, zero.get(), out, ld);
        };

        std::vector<double> ethr(nband, 1e-5);
        hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX = 200;
        hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR = 1e-5;
        hsolver::DiagoIterAssist<std::complex<double>>::SCF_ITER = 1;

        // ---- PPCG ----
        double en_ppcg[500] = {};
        {
            psi::Psi<std::complex<double>> psi_ppcg;
            psi_ppcg.resize(1, nband, npw);
            for (int i = 0; i < nband; i++)
                for (int j = 0; j < npw; j++)
                    psi_ppcg(i, j) = psiguess(i, j);

            hsolver::DiagoPPCG<std::complex<double>> ppcg(hpsi.precond());
            ppcg.init_iter(nband, nband, npw, npw);
            double start = MPI_Wtime();
            for (int pass = 0; pass < 5; ++pass)
                ppcg.diag(hpsi_func, psi_ppcg.get_pointer(), en_ppcg, ethr);
            t_ppcg[sz] = (MPI_Wtime() - start) * 1000.0;
        }

        // ---- BPCG ----
        double en_bpcg[500] = {};
        {
            psi::Psi<std::complex<double>> psi_bpcg;
            psi_bpcg.resize(1, nband, npw);
            for (int i = 0; i < nband; i++)
                for (int j = 0; j < npw; j++)
                    psi_bpcg(i, j) = psiguess(i, j);

            hsolver::DiagoBPCG<std::complex<double>> bpcg(hpsi.precond());
            bpcg.init_iter(nband, nband, npw, npw);
            double start = MPI_Wtime();
            for (int pass = 0; pass < 4; ++pass)
                bpcg.diag(hpsi_func, psi_bpcg.get_pointer(), en_bpcg, ethr);
            t_bpcg[sz] = (MPI_Wtime() - start) * 1000.0;
        }

        // ---- Davidson ----
        double en_david[500] = {};
        {
            psi::Psi<std::complex<double>> psi_dav;
            psi_dav.resize(1, nband, npw);
            for (int i = 0; i < nband; i++)
                for (int j = 0; j < npw; j++)
                    psi_dav(i, j) = psiguess(i, j);

            hsolver::diag_comm_info comm_info(
#ifdef __MPI
                MPI_COMM_SELF,
#endif
                0, 1);
            hsolver::DiagoDavid<std::complex<double>> david(hpsi.precond(), nband, npw, 4, comm_info);
            auto spsi_func = [npw](const std::complex<double>* in, std::complex<double>* out, int ld, int nv) {
                for (int ib = 0; ib < nv; ib++)
                    for (int i = 0; i < npw; i++)
                        out[ib * ld + i] = in[ib * ld + i];
            };
            double start = MPI_Wtime();
            david.diag(hpsi_func, spsi_func, npw, psi_dav.get_pointer(), en_david, ethr, 200);
            t_david[sz] = (MPI_Wtime() - start) * 1000.0;
        }

        // errors
        double err_ppcg = std::abs(en_ppcg[0] - e_lapack[0]);
        double err_bpcg = std::abs(en_bpcg[0] - e_lapack[0]);
        double err_david = std::abs(en_david[0] - e_lapack[0]);

        double s_ppcg = t_lapack[sz] / std::max(t_ppcg[sz], 1e-6);
        double s_bpcg = t_lapack[sz] / std::max(t_bpcg[sz], 1e-6);
        double s_david = t_lapack[sz] / std::max(t_david[sz], 1e-6);

        char buf[256];
        snprintf(buf, sizeof(buf),
                 " %5d  | %7.1f  %7.1f  %8.1f  %8.1f | %7.1fx  %7.1fx  %7.1fx  | %8.1e %8.1e %8.1e",
                 npw, t_ppcg[sz], t_bpcg[sz], t_david[sz], t_lapack[sz],
                 s_ppcg, s_bpcg, s_david,
                 err_ppcg, err_bpcg, err_david);
        std::cout << buf << std::endl;

        EXPECT_NEAR(en_ppcg[0], e_lapack[0], std::abs(e_lapack[0]) * 0.1 + 0.5);
        EXPECT_NEAR(en_bpcg[0], e_lapack[0], std::abs(e_lapack[0]) * 0.1 + 0.5);
        EXPECT_NEAR(en_david[0], e_lapack[0], std::abs(e_lapack[0]) * 0.1 + 0.5);

        delete[] e_lapack;
    }

    // ---- empirical complexity analysis ----
    // Fit  t = C * N^k  →  log(t) = log(C) + k * log(N)
    // Use adjacent pairs to estimate local exponent: k ≈ log(t2/t1) / log(N2/N1)
    std::cout << "\n--- Empirical complexity exponents (k in  t ∝ N^k) ---\n";
    for (int sz = 1; sz < n_sizes; ++sz)
    {
        double ratio_N = std::log((double)sizes[sz] / sizes[sz - 1]);
        double k_ppcg  = std::log(std::max(t_ppcg[sz], 1e-9) / std::max(t_ppcg[sz - 1], 1e-9)) / ratio_N;
        double k_bpcg  = std::log(std::max(t_bpcg[sz], 1e-9) / std::max(t_bpcg[sz - 1], 1e-9)) / ratio_N;
        double k_david = std::log(std::max(t_david[sz], 1e-9) / std::max(t_david[sz - 1], 1e-9)) / ratio_N;
        double k_lap   = std::log(std::max(t_lapack[sz], 1e-9) / std::max(t_lapack[sz - 1], 1e-9)) / ratio_N;
        printf("  %4d→%4d :  PPCG k=%.2f   BPCG k=%.2f   David k=%.2f   LAPACK k=%.2f\n",
               sizes[sz - 1], sizes[sz], k_ppcg, k_bpcg, k_david, k_lap);
    }

    // average speedup
    double avg_ppcg_vs_lap = 0, avg_bpcg_vs_lap = 0, avg_david_vs_lap = 0, avg_ppcg_vs_bpcg = 0, avg_ppcg_vs_david = 0;
    for (int sz = 0; sz < n_sizes; ++sz)
    {
        avg_ppcg_vs_lap   += t_lapack[sz] / std::max(t_ppcg[sz], 1e-6);
        avg_bpcg_vs_lap   += t_lapack[sz] / std::max(t_bpcg[sz], 1e-6);
        avg_david_vs_lap  += t_lapack[sz] / std::max(t_david[sz], 1e-6);
        avg_ppcg_vs_bpcg  += t_bpcg[sz] / std::max(t_ppcg[sz], 1e-6);
        avg_ppcg_vs_david += t_david[sz] / std::max(t_ppcg[sz], 1e-6);
    }
    avg_ppcg_vs_lap /= n_sizes;
    avg_bpcg_vs_lap /= n_sizes;
    avg_david_vs_lap /= n_sizes;
    avg_ppcg_vs_bpcg /= n_sizes;
    avg_ppcg_vs_david /= n_sizes;

    std::cout << "\n--- Average speedup (geometric mean over 5 sizes) ---\n"
              << "  PPCG  vs LAPACK  : " << avg_ppcg_vs_lap << "x\n"
              << "  BPCG  vs LAPACK  : " << avg_bpcg_vs_lap << "x\n"
              << "  David vs LAPACK  : " << avg_david_vs_lap << "x\n"
              << "  PPCG  vs BPCG    : " << avg_ppcg_vs_bpcg << "x\n"
              << "  PPCG  vs David   : " << avg_ppcg_vs_david << "x\n"
              << std::endl;
}

int main(int argc, char** argv)
{
    int nproc = 1, myrank = 0;

#ifdef __MPI
    int nproc_in_pool, kpar = 1, mypool, rank_in_pool;
    setupmpi(argc, argv, nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, kpar, mypool, rank_in_pool);
    MPI_Comm_dup(POOL_WORLD, &BP_WORLD);
    GlobalV::NPROC_IN_POOL = nproc;
#else
    MPI_Init(&argc, &argv);
#endif

    testing::InitGoogleTest(&argc, argv);
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
    if (myrank != 0)
    {
        delete listeners.Release(listeners.default_result_printer());
    }

    int result = RUN_ALL_TESTS();
    if (myrank == 0 && result != 0)
    {
        std::cout << "ERROR:some tests are not passed" << std::endl;
        return result;
    }

    MPI_Finalize();
    return 0;
}
