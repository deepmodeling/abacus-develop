#include "source_base/inverse_matrix.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_pw/module_pwdft/structure_factor.h"
#include "source_psi/psi.h"
#include "source_hamilt/hamilt.h"
#include "source_pw/module_pwdft/hamilt_pw.h"
#include "../diago_iter_assist.h"
#include "../diago_ppcg.h"
#include "diago_mock.h"
#include "mpi.h"
#include "source_basis/module_pw/test/test_tool.h"

#include <gtest/gtest.h>
#include <complex>
#include <random>

/************************************************
 *  unit test of functions in DiagoPPCG
 ***********************************************/

/**
 * Class DiagoPPCG is an approach for eigenvalue problems
 * using the Polak-Ribiere conjugate gradient with soft-locking.
 * This unittest tests DiagoPPCG::diag() for complex<double>, Device=cpu
 * with random Hermitian matrices of varying size and sparsity.
 *
 * The test is passed when the eigenvalues are close to those
 * calculated by LAPACK (zheev).
 */

void lapackEigen(int &npw, std::vector<std::complex<double>> &hm, double *e, bool outtime = false)
{
    clock_t start, end;
    start = clock();
    int lwork = 2 * npw;
    std::complex<double> *work2 = new std::complex<double>[lwork];
    double *rwork = new double[3 * npw - 2];
    int info = 0;
    char tmp_c1 = 'V', tmp_c2 = 'U';
    zheev_(&tmp_c1, &tmp_c2, &npw, hm.data(), &npw, e, work2, &lwork, rwork, &info);
    end = clock();
    if (outtime) {
        std::cout << "Lapack Run time: " << (double)(end - start) / CLOCKS_PER_SEC << " S" << std::endl;
    }
    delete[] rwork;
    delete[] work2;
}

class DiagoPPCGPrepare
{
  public:
    DiagoPPCGPrepare(int nband, int npw, int sparsity, bool reorder, double eps, int maxiter, double threshold)
        : nband(nband), npw(npw), sparsity(sparsity), reorder(reorder), eps(eps), maxiter(maxiter),
          threshold(threshold)
    {
#ifdef __MPI
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
#endif
    }

    int nband, npw, sparsity, maxiter, notconv;
    double eps, avg_iter;
    bool reorder;
    double threshold;
    int nprocs=1, mypnum=0;

    void CompareEigen(double *precondition)
    {
        // calculate eigenvalues by LAPACK
        double *e_lapack = new double[npw];
        auto ev = DIAGOTEST::hmatrix;
        if (mypnum == 0) {
            lapackEigen(npw, ev, e_lapack, false);
        }
        // initial guess of psi by perturbing lapack psi
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
        // run ppcg
        double *en = new double[npw];
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
        for (int i = 0; i < DIAGOTEST::npw; i++) precondition_local[i] = precondition[i];
#endif
        hsolver::DiagoPPCG<std::complex<double>> ppcg(precondition_local);
        psi_local.fix_k(0);
        double start, end;
        start = MPI_Wtime();
        using T = std::complex<double>;
        const int dim = DIAGOTEST::npw;
        const std::vector<T> &h_mat = DIAGOTEST::hmatrix_local;
        auto hpsi_func = [h_mat, dim](T *psi_in, T *hpsi_out,
                                const int ld_psi, const int nvec) {
            auto one = std::make_unique<T>(1.0);
            auto zero = std::make_unique<T>(0.0);
            const T *one_ = one.get();
            const T *zero_ = zero.get();

            base_device::DEVICE_CPU *ctx = {};
            ModuleBase::gemm_op<T, base_device::DEVICE_CPU>()(
                'N', 'N',
                dim, nvec, dim,
                one_,
                h_mat.data(), dim,
                psi_in, ld_psi,
                zero_,
                hpsi_out, ld_psi);
        };
        const int ndim = psi_local.get_current_ngk();
        ppcg.init_iter(nband, nband, npw, ndim);
        std::vector<double> ethr_band(nband, 1e-5);
        ppcg.diag(hpsi_func, psi_local.get_pointer(), en, ethr_band);
        ppcg.diag(hpsi_func, psi_local.get_pointer(), en, ethr_band);
        ppcg.diag(hpsi_func, psi_local.get_pointer(), en, ethr_band);
        ppcg.diag(hpsi_func, psi_local.get_pointer(), en, ethr_band);
        end = MPI_Wtime();

        delete [] DIAGOTEST::npw_local;
        delete [] precondition_local;

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
    HPsi<std::complex<double>> hpsi(dcp.nband, dcp.npw, dcp.sparsity);
    DIAGOTEST::hmatrix = hpsi.hamilt();
    DIAGOTEST::npw = dcp.npw;
    dcp.CompareEigen(hpsi.precond());
}

INSTANTIATE_TEST_SUITE_P(VerifyPPCG,
                         DiagoPPCGTest,
                         ::testing::Values(
                             DiagoPPCGPrepare(10, 500, 0, true, 1e-5, 300, 5e-2),
                             DiagoPPCGPrepare(20, 500, 6, true, 1e-5, 300, 5e-2),
                             DiagoPPCGPrepare(20, 1000, 8, true, 1e-5, 300, 5e-2),
                             DiagoPPCGPrepare(40, 1000, 8, true, 1e-6, 300, 5e-2)));

TEST(DiagoPPCGTest, SoftLock)
{
    // Verify that PPCG soft-locking converges on a small problem
    int dim = 100;
    int nbnd = 10;
    HPsi<std::complex<double>> hpsi(nbnd, dim);
    DIAGOTEST::hmatrix = hpsi.hamilt();
    DIAGOTEST::npw = dim;

    double *e_lapack = new double[dim];
    auto ev = DIAGOTEST::hmatrix;
    lapackEigen(dim, ev, e_lapack, false);

    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX = 300;
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR = 1e-5;

    ModuleBase::ComplexMatrix psiguess(nbnd, dim);
    std::default_random_engine p(1);
    std::uniform_int_distribution<unsigned> u(1, 10);
    for (int i = 0; i < nbnd; i++) {
        for (int j = 0; j < dim; j++) {
            double rand = static_cast<double>(u(p)) / 10.;
            psiguess(i, j) = ev[j * DIAGOTEST::h_nc + i] * rand;
        }
    }

    psi::Psi<std::complex<double>> psi;
    psi.resize(1, nbnd, dim);
    for (int i = 0; i < nbnd; i++) {
        for (int j = 0; j < dim; j++) {
            psi(i, j) = psiguess(i, j);
        }
    }

    double* precondition_local = new double[dim];
    for (int i = 0; i < dim; i++) precondition_local[i] = hpsi.precond()[i];
    DIAGOTEST::hmatrix_local = DIAGOTEST::hmatrix;
    DIAGOTEST::npw_local = new int[1];
    DIAGOTEST::npw_local[0] = dim;
    psi.fix_k(0);

    double *en = new double[dim];
    using T = std::complex<double>;
    auto hpsi_func = [dim](T *psi_in, T *hpsi_out,
                           const int ld_psi, const int nvec) {
        auto one = std::make_unique<T>(1.0);
        auto zero = std::make_unique<T>(0.0);
        const T *one_ = one.get();
        const T *zero_ = zero.get();
        base_device::DEVICE_CPU *ctx = {};
        ModuleBase::gemm_op<T, base_device::DEVICE_CPU>()(
            'N', 'N', dim, nvec, dim,
            one_, DIAGOTEST::hmatrix.data(), dim,
            psi_in, ld_psi, zero_, hpsi_out, ld_psi);
    };

    hsolver::DiagoPPCG<std::complex<double>> ppcg(precondition_local);
    ppcg.init_iter(nbnd, nbnd, dim, dim);
    std::vector<double> ethr_band(nbnd, 1e-5);
    // Run multiple diag calls to exercise soft-locking across iterations
    for (int iter = 0; iter < 4; iter++) {
        ppcg.diag(hpsi_func, psi.get_pointer(), en, ethr_band);
    }

    for (int i = 0; i < nbnd; i++) {
        EXPECT_NEAR(en[i], e_lapack[i], 5e-2);
    }

    delete[] DIAGOTEST::npw_local;
    delete[] precondition_local;
    delete[] en;
    delete[] e_lapack;
}

int main(int argc, char **argv)
{
    int nproc = 1, myrank = 0;

#ifdef __MPI
    int nproc_in_pool, kpar=1, mypool, rank_in_pool;
    setupmpi(argc, argv, nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, kpar, mypool, rank_in_pool);
    MPI_Comm_split(MPI_COMM_WORLD, myrank, 0, &BP_WORLD);
    GlobalV::NPROC_IN_POOL = nproc;
#else
    MPI_Init(&argc, &argv);
#endif

    testing::InitGoogleTest(&argc, argv);
    ::testing::TestEventListeners &listeners = ::testing::UnitTest::GetInstance()->listeners();
    if (myrank != 0) {
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
