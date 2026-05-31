#include "../diago_iter_assist.h"
#include "../diago_lobpcg.h"
#include "source_base/global_variable.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_base/parallel_comm.h"
#include "source_basis/module_pw/test/test_tool.h"

#ifdef __MPI
#include "mpi.h"
#endif

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <random>
#include <vector>

/************************************************
 *  unit test of DiagoLobpcg (NC, CPU-only)
 *
 *  Validates eigenvalues, orthonormality, and
 *  residual against LAPACK zheev for random
 *  well-conditioned Hermitian matrices.
 ***********************************************/

using TestT      = std::complex<double>;
using TestDevice = base_device::DEVICE_CPU;
using TestReal   = double;

/// Reference eigenvalues via LAPACK zheev (eigenvalues only).
static int lapackEigen(int npw, std::vector<TestT>& hm, TestReal* e)
{
    int lwork = std::max(1, 2 * npw - 1);
    std::vector<TestT> work(lwork);
    std::vector<TestReal> rwork(std::max(1, 3 * npw - 2));
    int info = 0;
    char jobz = 'N', uplo = 'U';
    zheev_(&jobz, &uplo, &npw, hm.data(), &npw, e,
           work.data(), &lwork, rwork.data(), &info);
    return info;
}

static void run_generalized_lobpcg()
{
    const int npw = 8, nband = 3;
    std::vector<TestT> psi(nband * npw, {0.0, 0.0});
    for (int ib = 0; ib < nband; ++ib)
        psi[ib * npw + ib] = {1.0, 0.0};
    std::vector<TestReal> eigens(nband, 0.0);
    std::vector<TestReal> prec(npw, 1.0);
    std::vector<double> ethr(nband, 1e-6);

    auto hpsi_func = [](TestT* psi_in, TestT* hpsi_out,
                        int ld_psi, int nvec) {
        std::copy(psi_in, psi_in + ld_psi * nvec, hpsi_out);
    };
    auto spsi_func = [](const TestT* psi_in, TestT* spsi_out,
                        int ld_psi, int nvec) {
        std::copy(psi_in, psi_in + ld_psi * nvec, spsi_out);
        spsi_out[0] *= 2.0;
    };

    hsolver::DiagoLobpcg<TestT, TestDevice> lobpcg(prec.data());
    lobpcg.init_iter(nband, nband, npw, npw);
    lobpcg.diag(hpsi_func, spsi_func, psi.data(), eigens.data(), ethr);
}

class DiagoLobpcgTest : public ::testing::Test
{
  protected:
    static int idx(int row, int col, int ld) { return col * ld + row; }

    void run_and_validate(int npw, int nband,
                          const std::vector<TestT>& hmat,
                          const std::vector<TestReal>& prec,
                          const std::vector<TestReal>& e_ref,
                          double eig_tol, double orth_tol, double res_tol)
    {
        const int ld = npw;

        // ---- random orthonormal initial guess ----
        std::vector<TestT> psi(nband * npw, {0.0, 0.0});
        {
            std::mt19937 gen(123);
            std::uniform_real_distribution<TestReal> dist(-1.0, 1.0);
            for (int ib = 0; ib < nband; ib++)
                for (int ig = 0; ig < npw; ig++)
                    psi[ib * npw + ig] = {dist(gen), dist(gen)};

            for (int ib = 0; ib < nband; ib++)
            {
                for (int jb = 0; jb < ib; jb++)
                {
                    TestT dot = {0.0, 0.0};
                    for (int ig = 0; ig < npw; ig++)
                        dot += std::conj(psi[jb * npw + ig]) * psi[ib * npw + ig];
                    for (int ig = 0; ig < npw; ig++)
                        psi[ib * npw + ig] -= dot * psi[jb * npw + ig];
                }
                TestReal norm = 0.0;
                for (int ig = 0; ig < npw; ig++)
                    norm += std::norm(psi[ib * npw + ig]);
                norm = std::sqrt(norm);
                for (int ig = 0; ig < npw; ig++)
                    psi[ib * npw + ig] /= norm;
            }
        }

        auto hpsi_func = [&](TestT* psi_in, TestT* hpsi_out,
                              int ld_psi, int nvec) {
            for (int iv = 0; iv < nvec; iv++)
                for (int i = 0; i < npw; i++)
                {
                    TestT sum = {0.0, 0.0};
                    for (int j = 0; j < npw; j++)
                        sum += hmat[idx(i, j, ld)] * psi_in[iv * ld_psi + j];
                    hpsi_out[iv * ld_psi + i] = sum;
                }
        };

        // ---- run LOBPCG ----
        std::vector<TestReal> eigens(nband, 0.0);
        std::vector<double> ethr(nband, 1e-6);

        // SCF_ITER = 1 triggers periodic R-R restart that clears P, which
        // naturally limits subspace ill-conditioning.  Use moderate nline.
        const int old_scf = hsolver::DiagoIterAssist<TestT, TestDevice>::SCF_ITER;
        hsolver::DiagoIterAssist<TestT, TestDevice>::SCF_ITER = 1;

        hsolver::DiagoLobpcg<TestT, TestDevice> lobpcg(prec.data());
        lobpcg.init_iter(nband, nband, npw, npw);
        lobpcg.set_nline(4);
        auto spsi_func = [](const TestT* psi_in, TestT* spsi_out,
                            int ld_psi, int nvec) {
            std::copy(psi_in, psi_in + ld_psi * nvec, spsi_out);
        };
        lobpcg.diag(hpsi_func, spsi_func, psi.data(), eigens.data(), ethr);

        hsolver::DiagoIterAssist<TestT, TestDevice>::SCF_ITER = old_scf;

        // ---- validate eigenvalues ----
        for (int ib = 0; ib < nband; ib++)
            ASSERT_NEAR(eigens[ib], e_ref[ib], eig_tol)
                << "eigenvalue[" << ib << "] mismatch: "
                << eigens[ib] << " vs ref " << e_ref[ib];

        // ---- validate orthonormality ----
        for (int i = 0; i < nband; i++)
            for (int j = 0; j < nband; j++)
            {
                TestT dot = {0.0, 0.0};
                for (int ig = 0; ig < npw; ig++)
                    dot += std::conj(psi[i * npw + ig]) * psi[j * npw + ig];
                if (i == j)
                {
                    EXPECT_NEAR(std::real(dot), 1.0, orth_tol)
                        << "psi^H psi diag[" << i << "] = " << std::real(dot);
                    EXPECT_NEAR(std::imag(dot), 0.0, orth_tol)
                        << "psi^H psi diag[" << i << "] imag = " << std::imag(dot);
                }
                else
                    EXPECT_NEAR(std::abs(dot), 0.0, orth_tol)
                        << "psi^H psi (" << i << "," << j << ") = " << std::abs(dot);
            }

        // ---- validate residual: ||H*psi_i - eig_i * psi_i|| ----
        for (int ib = 0; ib < nband; ib++)
        {
            TestReal res2 = 0.0;
            for (int i = 0; i < npw; i++)
            {
                TestT hxi = {0.0, 0.0};
                for (int j = 0; j < npw; j++)
                    hxi += hmat[idx(i, j, ld)] * psi[ib * npw + j];
                const auto r = hxi - eigens[ib] * psi[ib * npw + i];
                res2 += std::norm(r);
            }
            EXPECT_LT(std::sqrt(res2), res_tol)
                << "residual[" << ib << "] = " << std::sqrt(res2);
        }
    }

    /// Build a strongly diagonally-dominant random Hermitian matrix.
    static void build_well_conditioned(int npw, std::vector<TestT>& hmat,
                                       std::vector<TestReal>& prec,
                                       std::vector<TestReal>& e_ref)
    {
        const int ld = npw;
        hmat.assign(npw * npw, {0.0, 0.0});

        std::mt19937 gen(42);
        std::uniform_real_distribution<TestReal> dist(-1.0, 1.0);

        for (int i = 0; i < npw; i++)
        {
            for (int j = i; j < npw; j++)
            {
                TestReal re = dist(gen) * 0.5;
                TestReal im = (i != j) ? dist(gen) * 0.5 : 0.0;
                hmat[idx(i, j, ld)] = {re, im};
                hmat[idx(j, i, ld)] = {re, -im};
            }
            hmat[idx(i, i, ld)] += TestT(
                static_cast<TestReal>(2.0 * (i + 1) * (i + 1)), 0.0);
        }

        // Reference
        auto hcopy = hmat;
        e_ref.resize(npw);
        ASSERT_EQ(lapackEigen(npw, hcopy, e_ref.data()), 0);

        // Preconditioner
        prec.resize(npw);
        for (int i = 0; i < npw; i++)
            prec[i] = std::max(static_cast<TestReal>(1.0),
                               std::real(hmat[idx(i, i, ld)]));
    }
};

// ============================================================================
// Test cases: various matrix sizes and band counts
// ============================================================================

TEST_F(DiagoLobpcgTest, SmallMatrix)
{
    const int npw = 50, nband = 10;
    std::vector<TestT> hmat;
    std::vector<TestReal> prec, e_ref;
    build_well_conditioned(npw, hmat, prec, e_ref);
    run_and_validate(npw, nband, hmat, prec, e_ref, 1e-5, 1e-8, 1e-4);
}

TEST_F(DiagoLobpcgTest, MediumMatrix)
{
    const int npw = 200, nband = 20;
    std::vector<TestT> hmat;
    std::vector<TestReal> prec, e_ref;
    build_well_conditioned(npw, hmat, prec, e_ref);
    run_and_validate(npw, nband, hmat, prec, e_ref, 1e-4, 1e-6, 2e-4);
}

TEST_F(DiagoLobpcgTest, LargerMatrixFewBands)
{
    const int npw = 400, nband = 12;
    std::vector<TestT> hmat;
    std::vector<TestReal> prec, e_ref;
    build_well_conditioned(npw, hmat, prec, e_ref);
    run_and_validate(npw, nband, hmat, prec, e_ref, 1e-4, 1e-6, 3e-4);
}

TEST_F(DiagoLobpcgTest, NonIdentityOverlapNotImplemented)
{
    EXPECT_EXIT(
        run_generalized_lobpcg(),
        ::testing::ExitedWithCode(1),
        ".*"
    );
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
