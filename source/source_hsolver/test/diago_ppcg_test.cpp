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

} // namespace

TEST(DiagoPPCGTest, RandomHermitianEigenvalues)
{
    const int nband = 4;
    const int npw = 60;
    const int sparsity = 0;

    int nprocs = 1;
    int mypnum = 0;
#ifdef __MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
#endif

    HPsi<std::complex<double>> hpsi_mock(nband, npw, sparsity);
    DIAGOTEST::hmatrix = hpsi_mock.hamilt();
    DIAGOTEST::npw = npw;

    std::vector<double> e_lapack(npw, 0.0);
    auto h_lapack = DIAGOTEST::hmatrix;
    lapackEigen(npw, h_lapack, e_lapack.data());
#ifdef __MPI
    MPI_Bcast(e_lapack.data(), npw, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

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
        ModuleBase::gemm_op<T, base_device::DEVICE_CPU>()('N',
                                                          'N',
                                                          dim,
                                                          nvec,
                                                          dim,
                                                          &one,
                                                          h_mat.data(),
                                                          dim,
                                                          psi_in,
                                                          ld_psi,
                                                          &zero,
                                                          hpsi_out,
                                                          ld_psi);
    };

    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX = 80;
    hsolver::DiagoPPCG<std::complex<double>> ppcg(precondition_local);
    ppcg.init_iter(nband, nband, npw, psi_local.get_current_ngk());

    std::vector<double> eigen(nband, 0.0);
    std::vector<double> ethr_band(nband, 1e-7);
    ppcg.diag(hpsi_func, psi_local.get_pointer(), eigen.data(), ethr_band);
    ppcg.diag(hpsi_func, psi_local.get_pointer(), eigen.data(), ethr_band);

    for (int ib = 0; ib < nband; ++ib)
    {
        EXPECT_NEAR(eigen[ib], e_lapack[ib], 5e-2);
    }

    delete[] DIAGOTEST::npw_local;
    delete[] precondition_local;
}
