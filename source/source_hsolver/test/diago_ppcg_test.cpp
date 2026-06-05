#include "source_base/inverse_matrix.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_hamilt/hamilt.h"
#include "source_psi/psi.h"
#include "../diago_iter_assist.h"
#include "../diago_ppcg.h"
#include "diago_mock.h"
#include "source_basis/module_pw/test/test_tool.h"

#include <gtest/gtest.h>
#include <complex>
#include <random>

#include "mpi.h"

/************************************************
 *  unit test of functions in DiagoPPCG
 ***********************************************/

/**
 * Test the PPCG (Projected Preconditioned Conjugate Gradient) method
 * for eigenvalue problems.
 *
 * The test generates a random Hermitian matrix, computes eigenvalues
 * using LAPACK as reference, and compares with PPCG results.
 */

// call lapack to get reference eigenvalues
void lapackEigen(int &npw, std::vector<std::complex<double>> &hm, double *e)
{
    int lwork = 2 * npw;
    std::complex<double> *work2 = new std::complex<double>[lwork];
    double *rwork = new double[3 * npw - 2];
    int info = 0;
    char tmp_c1 = 'V', tmp_c2 = 'U';
    zheev_(&tmp_c1, &tmp_c2, &npw, hm.data(), &npw, e, work2, &lwork, rwork, &info);
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

    int nband, npw, sparsity, maxiter;
    double eps, avg_iter;
    double threshold;
    int nprocs = 1, mypnum = 0;

    void CompareEigen(double *precondition)
    {
        // Calculate eigenvalues by LAPACK
        double *e_lapack = new double[npw];
        auto ev = DIAGOTEST::hmatrix;
        if (mypnum == 0)
        {
            lapackEigen(npw, ev, e_lapack);
        }

        // Initial guess of psi by perturbing lapack psi
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

        double *en = new double[npw];
        hamilt::Hamilt<std::complex<double>> *ha = new hamilt::HamiltPW<std::complex<double>>(
            nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);

        psi::Psi<std::complex<double>> psi;
        psi.resize(1, nband, npw);
        for (int i = 0; i < nband; i++)
        {
            for (int j = 0; j < npw; j++)
            {
                psi(i, j) = psiguess(i, j);
            }
        }

        psi::Psi<std::complex<double>> psi_local;
        double *precondition_local;
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

        // Run PPCG
        psi_local.fix_k(0);
        using T = std::complex<double>;
        const int dim = DIAGOTEST::npw;
        const std::vector<T> &h_mat = DIAGOTEST::hmatrix_local;

        auto hpsi_func = [h_mat, dim](T *psi_in, T *hpsi_out,
                                      const int ld_psi, const int nvec)
        {
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

        // For PPCG (non-generalized), S = I
        auto spsi_func = [](T *psi_in, T *spsi_out,
                            const int ld_psi, const int nvec)
        {
            for (int i = 0; i < ld_psi * nvec; i++)
            {
                spsi_out[i] = psi_in[i];
            }
        };

        hsolver::DiagoPPCG<T> ppcg("pw", "scf");
        std::vector<double> ethr_band(nband, 1e-5);

        double avg_iter = ppcg.diag(hpsi_func,
                                    spsi_func,
                                    npw,
                                    nband,
                                    dim,
                                    psi_local.get_pointer(),
                                    en,
                                    ethr_band,
                                    precondition_local);

        // Run multiple times to ensure convergence
        for (int r = 0; r < 3; r++)
        {
            avg_iter = ppcg.diag(hpsi_func,
                                 spsi_func,
                                 npw,
                                 nband,
                                 dim,
                                 psi_local.get_pointer(),
                                 en,
                                 ethr_band,
                                 precondition_local);
        }

        // Compare with LAPACK reference
        for (int i = 0; i < nband; i++)
        {
            EXPECT_NEAR(en[i], e_lapack[i], threshold);
        }

        delete[] DIAGOTEST::npw_local;
        delete[] precondition_local;
        delete[] en;
        delete[] e_lapack;
        delete ha;
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
                             // nband, npw, sparsity, eps, maxiter, threshold
                             DiagoPPCGPrepare(4, 100, 0, 1e-5, 300, 5e-2)
                             ));

int main(int argc, char **argv)
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
    ::testing::TestEventListeners &listeners = ::testing::UnitTest::GetInstance()->listeners();
    if (myrank != 0)
    {
        delete listeners.Release(listeners.default_result_printer());
    }

    int result = RUN_ALL_TESTS();
    if (myrank == 0 && result != 0)
    {
        std::cout << "Some tests failed. Check the output above for details." << std::endl;
    }

    MPI_Finalize();
    return result;
}
