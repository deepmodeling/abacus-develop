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

        // A few passes for robustness on random problems
        for (int pass = 0; pass < 2; ++pass)
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
