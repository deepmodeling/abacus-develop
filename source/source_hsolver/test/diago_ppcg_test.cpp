/**
 * PPCG (Projected Preconditioned Conjugate Gradient) solver tests.
 *
 * Test cases:
 *   Fixed4x4Matrix     — fixed 4x4 Hermitian matrix with known eigenvalues
 *   SmallDense         — random 40x40 dense, 4 bands
 *   MediumDense        — random 100x100 dense, 10 bands
 *   MediumSparse       — random 100x100 sparse (60%), 10 bands
 *   LargeSparse        — random 200x200 sparse (80%), 20 bands
 *
 * Each test generates a random Hermitian matrix via HPsi, computes reference
 * eigenvalues with LAPACK zheev, runs PPCG with a perturbed initial guess,
 * and asserts the results match within tolerance.
 */
#include "gtest/gtest.h"

#include "../diago_iter_assist.h"
#include "../diago_ppcg.h"
#include "diago_mock.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_basis/module_pw/test/test_tool.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_hamilt/hamilt.h"
#include "source_pw/module_pwdft/hamilt_pw.h"
#include "source_psi/psi.h"

#include <complex>
#include <random>
#include <vector>

namespace
{

/// Compute all eigenvalues of a Hermitian matrix using LAPACK zheev.
void lapackEigen(const int npw, std::vector<std::complex<double>>& hm, double* e)
{
    int lwork = 2 * npw;
    std::vector<std::complex<double>> work(lwork);
    std::vector<double> rwork(3 * npw - 2);
    int info = 0;
    char jobz = 'V';
    char uplo = 'U';
    zheev_(&jobz, &uplo, &npw, hm.data(), &npw, e, work.data(), &lwork, rwork.data(), &info);
    ASSERT_EQ(info, 0);
}

/// Common PPCG test runner: generate random H, compare PPCG eigenvalues with LAPACK.
void runPPCGTest(const int nband, const int npw, const int sparsity, const double tolerance)
{
    int nprocs = 1;
    int mypnum = 0;
#ifdef __MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
#endif

    // Generate random Hermitian matrix + precondition via HPsi
    HPsi<std::complex<double>> hpsi_mock(nband, npw, sparsity);
    DIAGOTEST::hmatrix = hpsi_mock.hamilt();
    DIAGOTEST::npw = npw;

    // Reference eigenvalues from LAPACK
    std::vector<double> e_lapack(npw, 0.0);
    auto h_lapack = DIAGOTEST::hmatrix;
    lapackEigen(npw, h_lapack, e_lapack.data());
#ifdef __MPI
    MPI_Bcast(e_lapack.data(), npw, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

    // Initial psi: perturb LAPACK eigenvectors to simulate a poor initial guess
    psi::Psi<std::complex<double>> psi;
    psi.resize(1, nband, npw);
    std::default_random_engine engine(7);
    std::uniform_real_distribution<double> dist(0.2, 1.0);
    for (int ib = 0; ib < nband; ++ib)
    {
        for (int ig = 0; ig < npw; ++ig)
        {
            psi(ib, ig) = h_lapack[ig + ib * npw] * dist(engine);
        }
    }

    // Distribute data across MPI processes
    psi::Psi<std::complex<double>> psi_local;
    DIAGOTEST::npw_local = new int[nprocs];
    double* precondition_local = nullptr;
#ifdef __MPI
    DIAGOTEST::cal_division(DIAGOTEST::npw);
    DIAGOTEST::divide_hpsi(psi, psi_local, DIAGOTEST::hmatrix, DIAGOTEST::hmatrix_local);
    precondition_local = new double[DIAGOTEST::npw_local[mypnum]];
    DIAGOTEST::divide_psi<double>(hpsi_mock.precond(), precondition_local);
#else
    DIAGOTEST::hmatrix_local = DIAGOTEST::hmatrix;
    DIAGOTEST::npw_local[0] = DIAGOTEST::npw;
    psi_local = psi;
    precondition_local = new double[DIAGOTEST::npw];
    for (int ig = 0; ig < DIAGOTEST::npw; ++ig)
    {
        precondition_local[ig] = hpsi_mock.precond()[ig];
    }
#endif

    psi_local.fix_k(0);
    using T = std::complex<double>;
    const int dim = DIAGOTEST::npw;
    const std::vector<T>& h_mat = DIAGOTEST::hmatrix_local;
    auto hpsi_func = [h_mat, dim](T* psi_in, T* hpsi_out, const int ld_psi, const int nvec) {
        const T one(1.0);
        const T zero(0.0);
        ModuleBase::gemm_op<T, base_device::DEVICE_CPU>()(
            'N', 'N',
            dim, nvec, dim,
            &one,
            h_mat.data(), dim,
            psi_in, ld_psi,
            &zero,
            hpsi_out, ld_psi);
    };

    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX = 80;
    hsolver::DiagoPPCG<std::complex<double>> ppcg(precondition_local);
    ppcg.init_iter(nband, nband, npw, psi_local.get_current_ngk());

    std::vector<double> eigen(nband, 0.0);
    std::vector<double> ethr_band(nband, 1e-7);
    ppcg.diag(hpsi_func, psi_local.get_pointer(), eigen.data(), ethr_band);

    for (int ib = 0; ib < nband; ++ib)
    {
        EXPECT_NEAR(eigen[ib], e_lapack[ib], tolerance);
    }

    delete[] DIAGOTEST::npw_local;
    delete[] precondition_local;
}

} // namespace

// ====== Fixed matrix tests ======

TEST(DiagoPPCGTest, Fixed4x4Matrix)
{
    const int nband = 2;
    const int npw = 4;
    const int sparsity = 0;

    int nprocs = 1;
    int mypnum = 0;
#ifdef __MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
#endif

    // clang-format off
    std::vector<std::complex<double>> h_fixed(16);
    h_fixed[0]  = {4.0, 0.0}; h_fixed[1]  = {1.0, 0.0}; h_fixed[2]  = {1.0, 0.0}; h_fixed[3]  = {0.0, 0.0};
    h_fixed[4]  = {1.0, 0.0}; h_fixed[5]  = {3.0, 0.0}; h_fixed[6]  = {0.0, 0.0}; h_fixed[7]  = {1.0, 0.0};
    h_fixed[8]  = {1.0, 0.0}; h_fixed[9]  = {0.0, 0.0}; h_fixed[10] = {2.0, 0.0}; h_fixed[11] = {1.0, 0.0};
    h_fixed[12] = {0.0, 0.0}; h_fixed[13] = {1.0, 0.0}; h_fixed[14] = {1.0, 0.0}; h_fixed[15] = {5.0, 0.0};
    // clang-format on

    HPsi<std::complex<double>> hpsi_mock(nband, npw, sparsity);
    DIAGOTEST::hmatrix = h_fixed;
    DIAGOTEST::npw = npw;

    std::vector<double> e_lapack(npw, 0.0);
    auto h_lapack = h_fixed;
    lapackEigen(npw, h_lapack, e_lapack.data());
#ifdef __MPI
    MPI_Bcast(e_lapack.data(), npw, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

    psi::Psi<std::complex<double>> psi;
    psi.resize(1, nband, npw);
    std::default_random_engine engine(42);
    std::uniform_real_distribution<double> dist(0.1, 1.0);
    for (int ib = 0; ib < nband; ++ib)
    {
        for (int ig = 0; ig < npw; ++ig)
        {
            psi(ib, ig) = h_lapack[ig + ib * npw] * dist(engine);
        }
    }

    psi::Psi<std::complex<double>> psi_local;
    DIAGOTEST::npw_local = new int[nprocs];
    double* precondition_local = nullptr;
#ifdef __MPI
    DIAGOTEST::cal_division(DIAGOTEST::npw);
    DIAGOTEST::divide_hpsi(psi, psi_local, DIAGOTEST::hmatrix, DIAGOTEST::hmatrix_local);
    precondition_local = new double[DIAGOTEST::npw_local[mypnum]];
    DIAGOTEST::divide_psi<double>(hpsi_mock.precond(), precondition_local);
#else
    DIAGOTEST::hmatrix_local = DIAGOTEST::hmatrix;
    DIAGOTEST::npw_local[0] = DIAGOTEST::npw;
    psi_local = psi;
    precondition_local = new double[DIAGOTEST::npw];
    for (int ig = 0; ig < DIAGOTEST::npw; ++ig)
    {
        precondition_local[ig] = hpsi_mock.precond()[ig];
    }
#endif

    psi_local.fix_k(0);
    using T = std::complex<double>;
    const int dim = DIAGOTEST::npw;
    const std::vector<T>& h_mat = DIAGOTEST::hmatrix_local;
    auto hpsi_func = [h_mat, dim](T* psi_in, T* hpsi_out, const int ld_psi, const int nvec) {
        const T one(1.0);
        const T zero(0.0);
        ModuleBase::gemm_op<T, base_device::DEVICE_CPU>()(
            'N', 'N',
            dim, nvec, dim,
            &one,
            h_mat.data(), dim,
            psi_in, ld_psi,
            &zero,
            hpsi_out, ld_psi);
    };

    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX = 50;
    hsolver::DiagoPPCG<std::complex<double>> ppcg(precondition_local);
    ppcg.init_iter(nband, nband, npw, psi_local.get_current_ngk());

    std::vector<double> eigen(nband, 0.0);
    std::vector<double> ethr_band(nband, 1e-7);
    ppcg.diag(hpsi_func, psi_local.get_pointer(), eigen.data(), ethr_band);

    for (int ib = 0; ib < nband; ++ib)
    {
        EXPECT_NEAR(eigen[ib], e_lapack[ib], 1e-2);
    }

    delete[] DIAGOTEST::npw_local;
    delete[] precondition_local;
}

// ====== Random Hermitian matrix tests ======

TEST(DiagoPPCGTest, SmallDense)
{
    runPPCGTest(4, 40, 0, 1e-2);
}

TEST(DiagoPPCGTest, MediumDense)
{
    runPPCGTest(10, 100, 0, 5e-2);
}

TEST(DiagoPPCGTest, MediumSparse)
{
    runPPCGTest(10, 100, 6, 5e-2);
}

TEST(DiagoPPCGTest, LargeSparse)
{
    runPPCGTest(20, 200, 8, 5e-2);
}


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
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
    if (myrank != 0)
    {
        delete listeners.Release(listeners.default_result_printer());
    }

    int result = RUN_ALL_TESTS();
    if (myrank == 0 && result != 0)
    {
        std::cout << "ERROR: some tests are not passed" << std::endl;
        return result;
    }

    MPI_Finalize();
    return 0;
}
