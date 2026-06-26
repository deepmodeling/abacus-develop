/**
 * OpenMP consistency test for eigenvalue solvers.
 * Verifies that BPCG and Davidson produce identical results
 * across different OMP_NUM_THREADS values.
 */
#include <cstdlib>
#include "source_base/module_external/lapack_connector.h"
#include "source_pw/module_pwdft/hamilt_pw.h"
#include "source_psi/psi.h"
#include "source_hamilt/hamilt.h"
#include "../diago_iter_assist.h"
#include "../diago_bpcg.h"
#include "../diago_david.h"
#include "diago_mock.h"
#include "source_basis/module_pw/test/test_tool.h"

#include <gtest/gtest.h>
#include <complex>
#include <random>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{

void lapackEigen(int& npw, std::vector<std::complex<double>>& hm, double* e)
{
    int lwork = 2 * npw;
    std::complex<double>* work2 = new std::complex<double>[lwork];
    double* rwork = new double[3 * npw - 2];
    int info = 0;
    char tmp_c1 = 'V', tmp_c2 = 'U';
    zheev_(&tmp_c1, &tmp_c2, &npw, hm.data(), &npw, e, work2, &lwork, rwork, &info);
    delete[] rwork;
    delete[] work2;
}

// Run BPCG with given matrix and precondition, return first nband eigenvalues
std::vector<double> run_bpcg(int nband, int npw,
                             const std::vector<std::complex<double>>& hmatrix,
                             const std::vector<double>& precondition)
{
    DIAGOTEST::hmatrix = hmatrix;
    DIAGOTEST::npw = npw;
    DIAGOTEST::npw_local = new int[1];
    DIAGOTEST::npw_local[0] = npw;
    DIAGOTEST::hmatrix_local = hmatrix;

    psi::Psi<std::complex<double>> psi;
    psi.resize(1, nband, npw);
    std::default_random_engine p(1);
    std::uniform_int_distribution<unsigned> u(1, 10);
    for (int i = 0; i < nband; i++)
    {
        for (int j = 0; j < npw; j++)
        {
            psi(i, j) = static_cast<double>(u(p)) / 10.0;
        }
    }

    double* precondition_local = new double[npw];
    for (int i = 0; i < npw; i++) precondition_local[i] = precondition[i];

    hsolver::DiagoBPCG<std::complex<double>> bpcg(precondition_local);
    psi.fix_k(0);
    const int dim = npw;
    using T = std::complex<double>;
    auto hpsi_func = [hmatrix, dim](T* psi_in, T* hpsi_out,
                                    const int ld_psi, const int nvec) {
        T one(1.0), zero(0.0);
        ModuleBase::gemm_op<T, base_device::DEVICE_CPU>()(
            'N', 'N', dim, nvec, dim, &one,
            hmatrix.data(), dim, psi_in, ld_psi,
            &zero, hpsi_out, ld_psi);
    };

    bpcg.init_iter(nband, nband, npw, npw);
    std::vector<double> eigen(nband, 0.0);
    std::vector<double> ethr_band(nband, 1e-5);
    bpcg.diag(hpsi_func, psi.get_pointer(), eigen.data(), ethr_band);

    delete[] precondition_local;
    delete[] DIAGOTEST::npw_local;
    return eigen;
}

// Run Davidson with given matrix and precondition, return first nband eigenvalues
std::vector<double> run_davidson(int nband, int npw,
                                 const std::vector<std::complex<double>>& hmatrix,
                                 const std::vector<double>& precondition)
{
    DIAGOTEST::hmatrix = hmatrix;
    DIAGOTEST::npw = npw;
    DIAGOTEST::npw_local = new int[1];
    DIAGOTEST::npw_local[0] = npw;
    DIAGOTEST::hmatrix_local = hmatrix;

    psi::Psi<std::complex<double>> psi;
    psi.resize(1, nband, npw);
    std::default_random_engine p(1);
    std::uniform_int_distribution<unsigned> u(1, 10);
    for (int i = 0; i < nband; i++)
    {
        for (int j = 0; j < npw; j++)
        {
            psi(i, j) = static_cast<double>(u(p)) / 10.0;
        }
    }

    double* precondition_local = new double[npw];
    for (int i = 0; i < npw; i++) precondition_local[i] = precondition[i];

#ifdef __MPI
    hsolver::diag_comm_info comm_info(MPI_COMM_WORLD, 0, 1);
#else
    hsolver::diag_comm_info comm_info(0, 1);
#endif
    hsolver::DiagoDavid<std::complex<double>> dav(precondition_local, nband, npw, 4, comm_info);
    psi.fix_k(0);
    const int dim = npw;
    using T = std::complex<double>;
    auto hpsi_func = [hmatrix, dim](T* psi_in, T* hpsi_out,
                                    const int ld_psi, const int nvec) {
        T one(1.0), zero(0.0);
        ModuleBase::gemm_op<T, base_device::DEVICE_CPU>()(
            'N', 'N', dim, nvec, dim, &one,
            hmatrix.data(), dim, psi_in, ld_psi,
            &zero, hpsi_out, ld_psi);
    };
    auto spsi_func = [](T* psi_in, T* spsi_out,
                        const int ld_psi, const int nvec) {
        std::copy(psi_in, psi_in + static_cast<size_t>(ld_psi) * nvec, spsi_out);
    };

    std::vector<double> eigen(nband, 0.0);
    std::vector<double> ethr_band(nband, 1e-12);
    dav.diag(hpsi_func, spsi_func, npw, psi.get_pointer(), eigen.data(), ethr_band, 500);

    delete[] precondition_local;
    delete[] DIAGOTEST::npw_local;
    return eigen;
}

} // namespace

class OpenMPConsistencyTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Ensure consistent random state
        std::srand(42);
    }
};

TEST_F(OpenMPConsistencyTest, BPCG_ThreadConsistency)
{
    const int npw = 200;
    const int nband = 20;
    const int sparsity = 5;

    HPsi<std::complex<double>> hpsi(nband, npw, sparsity);
    auto hmatrix = hpsi.hamilt();
    std::vector<double> precondition(npw);
    for (int i = 0; i < npw; i++) precondition[i] = hpsi.precond()[i];

    // Reference eigenvalues with 1 thread
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    auto ref_eigen = run_bpcg(nband, npw, hmatrix, precondition);

    // Test with 2 and 4 threads
    for (int nthreads : {2, 4})
    {
#ifdef _OPENMP
        omp_set_num_threads(nthreads);
#endif
        auto test_eigen = run_bpcg(nband, npw, hmatrix, precondition);

        for (int i = 0; i < nband; i++)
        {
            EXPECT_NEAR(test_eigen[i], ref_eigen[i], 1e-5)
                << "BPCG eigenvalue mismatch at band " << i
                << " with threads=" << nthreads;
        }
    }
}

TEST_F(OpenMPConsistencyTest, Davidson_ThreadConsistency)
{
    const int npw = 200;
    const int nband = 20;
    const int sparsity = 5;

    HPsi<std::complex<double>> hpsi(nband, npw, sparsity);
    auto hmatrix = hpsi.hamilt();
    std::vector<double> precondition(npw);
    for (int i = 0; i < npw; i++) precondition[i] = hpsi.precond()[i];

    // Reference eigenvalues with 1 thread
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    auto ref_eigen = run_davidson(nband, npw, hmatrix, precondition);

    // Test with 2 and 4 threads
    for (int nthreads : {2, 4})
    {
#ifdef _OPENMP
        omp_set_num_threads(nthreads);
#endif
        auto test_eigen = run_davidson(nband, npw, hmatrix, precondition);

        for (int i = 0; i < nband; i++)
        {
            EXPECT_NEAR(test_eigen[i], ref_eigen[i], 1e-5)
                << "Davidson eigenvalue mismatch at band " << i
                << " with threads=" << nthreads;
        }
    }
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
    int result = RUN_ALL_TESTS();
#ifdef __MPI
    finishmpi();
#else
    MPI_Finalize();
#endif
    return result;
}
