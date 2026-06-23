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
#include <cstdlib>
#include <limits>
#include <random>
#include <string>
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

/// Reference generalized eigenvalues via LAPACK zhegvd (eigenvalues only).
static int lapackGeneralizedEigen(int npw,
                                  std::vector<TestT>& hm,
                                  std::vector<TestT>& sm,
                                  TestReal* e)
{
    int info = 0;
    int itype = 1;
    char jobz = 'N', uplo = 'U';
    int lwork = -1, lrwork = -1, liwork = -1;
    TestT work_query = {0.0, 0.0};
    TestReal rwork_query = 0.0;
    int iwork_query = 0;
    zhegvd_(&itype, &jobz, &uplo, &npw, hm.data(), &npw, sm.data(), &npw, e,
            &work_query, &lwork, &rwork_query, &lrwork, &iwork_query, &liwork, &info);
    if (info != 0)
        return info;

    lwork = std::max(1, static_cast<int>(std::real(work_query)));
    lrwork = std::max(1, static_cast<int>(rwork_query));
    liwork = std::max(1, iwork_query);
    std::vector<TestT> work(lwork);
    std::vector<TestReal> rwork(lrwork);
    std::vector<int> iwork(liwork);
    zhegvd_(&itype, &jobz, &uplo, &npw, hm.data(), &npw, sm.data(), &npw, e,
            work.data(), &lwork, rwork.data(), &lrwork, iwork.data(), &liwork, &info);
    return info;
}

class DiagoLobpcgTest : public ::testing::Test
{
  protected:
    struct HegvdMetrics
    {
        int info = 0;
        TestReal max_rel_residual = 0.0;
        TestReal max_s_orth_error = 0.0;
        TestReal max_abs_coeff = 0.0;
    };

    static int idx(int row, int col, int ld) { return col * ld + row; }

    static HegvdMetrics solve_generalized_eigenvectors(
        int dim,
        std::vector<TestT> hmat,
        std::vector<TestT> smat)
    {
        const auto h_orig = hmat;
        const auto s_orig = smat;
        std::vector<TestReal> eval(dim, 0.0);

        HegvdMetrics metrics;
        int itype = 1;
        char jobz = 'V', uplo = 'U';
        int lwork = -1, lrwork = -1, liwork = -1;
        TestT work_query = {0.0, 0.0};
        TestReal rwork_query = 0.0;
        int iwork_query = 0;
        zhegvd_(&itype, &jobz, &uplo, &dim, hmat.data(), &dim, smat.data(), &dim,
                eval.data(), &work_query, &lwork, &rwork_query, &lrwork,
                &iwork_query, &liwork, &metrics.info);
        if (metrics.info != 0)
            return metrics;

        lwork = std::max(1, static_cast<int>(std::real(work_query)));
        lrwork = std::max(1, static_cast<int>(rwork_query));
        liwork = std::max(1, iwork_query);
        std::vector<TestT> work(lwork);
        std::vector<TestReal> rwork(lrwork);
        std::vector<int> iwork(liwork);
        zhegvd_(&itype, &jobz, &uplo, &dim, hmat.data(), &dim, smat.data(), &dim,
                eval.data(), work.data(), &lwork, rwork.data(), &lrwork,
                iwork.data(), &liwork, &metrics.info);
        if (metrics.info != 0)
            return metrics;

        std::vector<TestT> av(dim), bv(dim);
        for (int iv = 0; iv < dim; ++iv)
        {
            const TestT* vec = hmat.data() + iv * dim;
            for (int i = 0; i < dim; ++i)
            {
                metrics.max_abs_coeff = std::max(
                    metrics.max_abs_coeff,
                    static_cast<TestReal>(std::abs(vec[i])));

                TestT ah = {0.0, 0.0};
                TestT bs = {0.0, 0.0};
                for (int j = 0; j < dim; ++j)
                {
                    ah += h_orig[idx(i, j, dim)] * vec[j];
                    bs += s_orig[idx(i, j, dim)] * vec[j];
                }
                av[i] = ah;
                bv[i] = bs;
            }

            TestReal res2 = 0.0;
            TestReal av2 = 0.0;
            TestReal bv2 = 0.0;
            for (int i = 0; i < dim; ++i)
            {
                res2 += std::norm(av[i] - eval[iv] * bv[i]);
                av2 += std::norm(av[i]);
                bv2 += std::norm(bv[i]);
            }
            const TestReal denom = std::max(
                static_cast<TestReal>(1.0),
                std::sqrt(av2) + std::abs(eval[iv]) * std::sqrt(bv2));
            metrics.max_rel_residual = std::max(metrics.max_rel_residual,
                                                std::sqrt(res2) / denom);

            for (int jv = 0; jv < dim; ++jv)
            {
                const TestT* vec_j = hmat.data() + jv * dim;
                TestT dot = {0.0, 0.0};
                for (int i = 0; i < dim; ++i)
                {
                    TestT bvi = {0.0, 0.0};
                    for (int j = 0; j < dim; ++j)
                        bvi += s_orig[idx(i, j, dim)] * vec_j[j];
                    dot += std::conj(vec[i]) * bvi;
                }
                const TestT target = (iv == jv) ? TestT(1.0, 0.0) : TestT(0.0, 0.0);
                metrics.max_s_orth_error = std::max(
                    metrics.max_s_orth_error,
                    static_cast<TestReal>(std::abs(dot - target)));
            }
        }
        return metrics;
    }

    static void build_nearly_dependent_overlap_problem(
        int dim,
        TestReal delta,
        std::vector<TestT>& hmat,
        std::vector<TestT>& smat,
        std::vector<TestReal>& prec,
        std::vector<TestReal>& e_ref)
    {
        build_well_conditioned(dim, hmat, prec, e_ref);
        smat.assign(dim * dim, {0.0, 0.0});
        for (int i = 0; i < dim; ++i)
            smat[idx(i, i, dim)] = {1.0, 0.0};
        smat[idx(0, 1, dim)] = {1.0 - delta, 0.0};
        smat[idx(1, 0, dim)] = {1.0 - delta, 0.0};

        auto hcopy = hmat;
        auto scopy = smat;
        ASSERT_EQ(lapackGeneralizedEigen(dim, hcopy, scopy, e_ref.data()), 0);
    }

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
        auto spsi_func = [&](const TestT* psi_in, TestT* spsi_out,
                              int ld_psi, int nvec) {
            for (int iv = 0; iv < nvec; iv++)
            {
                for (int i = 0; i < npw; i++)
                    spsi_out[iv * ld_psi + i] = psi_in[iv * ld_psi + i];
                for (int i = npw; i < ld_psi; i++)
                    spsi_out[iv * ld_psi + i] = {0.0, 0.0};
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

    static void matvec(const std::vector<TestT>& mat,
                       const TestT* psi_in,
                       TestT* out,
                       int npw,
                       int ld_psi,
                       int nvec)
    {
        const int ld = npw;
        for (int iv = 0; iv < nvec; iv++)
        {
            for (int i = 0; i < npw; i++)
            {
                TestT sum = {0.0, 0.0};
                for (int j = 0; j < npw; j++)
                    sum += mat[idx(i, j, ld)] * psi_in[iv * ld_psi + j];
                out[iv * ld_psi + i] = sum;
            }
            for (int i = npw; i < ld_psi; i++)
                out[iv * ld_psi + i] = {123.0, -456.0};
        }
    }

    static void build_generalized_problem(int npw,
                                          TestReal overlap_diag,
                                          TestReal overlap_scale,
                                          int seed,
                                          std::vector<TestT>& hmat,
                                          std::vector<TestT>& smat,
                                          std::vector<TestReal>& prec,
                                          std::vector<TestReal>& e_ref)
    {
        build_well_conditioned(npw, hmat, prec, e_ref);
        smat.assign(npw * npw, {0.0, 0.0});

        std::mt19937 gen(seed);
        std::uniform_real_distribution<TestReal> dist(-1.0, 1.0);
        for (int i = 0; i < npw; i++)
        {
            for (int j = i; j < npw; j++)
            {
                TestReal re = dist(gen) * overlap_scale;
                TestReal im = (i != j) ? dist(gen) * overlap_scale : 0.0;
                smat[idx(i, j, npw)] = {re, im};
                smat[idx(j, i, npw)] = {re, -im};
            }
            smat[idx(i, i, npw)] += TestT(overlap_diag, 0.0);
        }

        auto hcopy = hmat;
        auto scopy = smat;
        e_ref.resize(npw);
        ASSERT_EQ(lapackGeneralizedEigen(npw, hcopy, scopy, e_ref.data()), 0);
    }

    void run_generalized_and_validate(int npw,
                                      int nband,
                                      int ld_psi,
                                      TestReal overlap_diag,
                                      TestReal overlap_scale,
                                      int seed,
                                      TestReal min_s_minus_identity,
                                      double eig_tol,
                                      double orth_tol,
                                      double res_tol)
    {
        std::vector<TestT> hmat, smat;
        std::vector<TestReal> prec, e_ref;
        build_generalized_problem(npw, overlap_diag, overlap_scale, seed,
                                  hmat, smat, prec, e_ref);

        TestReal max_s_minus_identity = 0.0;
        for (int j = 0; j < npw; j++)
            for (int i = 0; i < npw; i++)
            {
                const TestT identity = (i == j) ? TestT(1.0, 0.0) : TestT(0.0, 0.0);
                max_s_minus_identity = std::max(
                    max_s_minus_identity,
                    static_cast<TestReal>(std::abs(smat[idx(i, j, npw)] - identity)));
            }
        ASSERT_GT(max_s_minus_identity, min_s_minus_identity);

        std::vector<TestT> psi(nband * ld_psi, {0.0, 0.0});
        for (int ib = 0; ib < nband; ib++)
        {
            psi[ib * ld_psi + ib] = {1.0, 0.0};
            for (int ig = npw; ig < ld_psi; ig++)
                psi[ib * ld_psi + ig] = {99.0, -77.0};
        }

        auto hpsi_func = [&](TestT* psi_in, TestT* hpsi_out,
                              int ld_in, int nvec) {
            matvec(hmat, psi_in, hpsi_out, npw, ld_in, nvec);
        };
        auto spsi_func = [&](const TestT* psi_in, TestT* spsi_out,
                              int ld_in, int nvec) {
            matvec(smat, psi_in, spsi_out, npw, ld_in, nvec);
        };

        std::vector<TestReal> eigens(nband, 0.0);
        std::vector<double> ethr(nband, res_tol);
        const int old_scf = hsolver::DiagoIterAssist<TestT, TestDevice>::SCF_ITER;
        hsolver::DiagoIterAssist<TestT, TestDevice>::SCF_ITER = 1;

        hsolver::DiagoLobpcg<TestT, TestDevice> lobpcg(prec.data());
        lobpcg.init_iter(nband, nband, ld_psi, npw);
        lobpcg.set_nline(10);
        lobpcg.diag(hpsi_func, spsi_func, psi.data(), eigens.data(), ethr);

        hsolver::DiagoIterAssist<TestT, TestDevice>::SCF_ITER = old_scf;

        for (int ib = 0; ib < nband; ib++)
            ASSERT_NEAR(eigens[ib], e_ref[ib], eig_tol);

        std::vector<TestT> hpsi(nband * ld_psi), spsi(nband * ld_psi);
        matvec(hmat, psi.data(), hpsi.data(), npw, ld_psi, nband);
        matvec(smat, psi.data(), spsi.data(), npw, ld_psi, nband);
        for (int i = 0; i < nband; i++)
        {
            for (int ig = npw; ig < ld_psi; ig++)
                EXPECT_EQ(psi[i * ld_psi + ig], TestT(0.0, 0.0));

            for (int j = 0; j < nband; j++)
            {
                TestT dot = {0.0, 0.0};
                for (int ig = 0; ig < npw; ig++)
                    dot += std::conj(psi[i * ld_psi + ig]) * spsi[j * ld_psi + ig];
                EXPECT_NEAR(
                    std::abs(dot - (i == j ? TestT(1.0, 0.0) : TestT(0.0, 0.0))),
                    0.0, orth_tol);
            }

            TestReal res2 = 0.0;
            for (int ig = 0; ig < npw; ig++)
                res2 += std::norm(hpsi[i * ld_psi + ig] - eigens[i] * spsi[i * ld_psi + ig]);
            EXPECT_LT(std::sqrt(res2), res_tol);
        }
    }

#ifdef __MPI
    void run_generalized_band_parallel_and_validate()
    {
        int nproc = 1;
        int rank = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (nproc < 2)
            GTEST_SKIP() << "band-parallel LOBPCG test requires at least 2 MPI ranks";

        const int npw = 18;
        const int nband = 11;
        const int ld_psi = npw + 3;
        const int nband_l = nband / nproc + (rank < nband % nproc ? 1 : 0);
        const int band_start = nband / nproc * rank + std::min(rank, nband % nproc);

        std::vector<TestT> hmat, smat;
        std::vector<TestReal> prec, e_ref;
        build_generalized_problem(npw, 1.7, 0.03, 101, hmat, smat, prec, e_ref);

        std::vector<TestT> psi(nband_l * ld_psi, {0.0, 0.0});
        for (int ib = 0; ib < nband_l; ++ib)
        {
            const int global_band = band_start + ib;
            psi[ib * ld_psi + global_band] = {1.0, 0.0};
            for (int ig = npw; ig < ld_psi; ++ig)
                psi[ib * ld_psi + ig] = {99.0, -77.0};
        }

        auto hpsi_func = [&](TestT* psi_in, TestT* hpsi_out,
                              int ld_in, int nvec) {
            matvec(hmat, psi_in, hpsi_out, npw, ld_in, nvec);
        };
        auto spsi_func = [&](const TestT* psi_in, TestT* spsi_out,
                              int ld_in, int nvec) {
            matvec(smat, psi_in, spsi_out, npw, ld_in, nvec);
        };

        std::vector<TestReal> eigens(nband_l, 0.0);
        std::vector<double> ethr(nband_l, 1e-8);
        const int old_scf = hsolver::DiagoIterAssist<TestT, TestDevice>::SCF_ITER;
        hsolver::DiagoIterAssist<TestT, TestDevice>::SCF_ITER = 1;

        hsolver::DiagoLobpcg<TestT, TestDevice> lobpcg(prec.data());
        lobpcg.init_iter(nband, nband_l, ld_psi, npw);
        lobpcg.set_nline(10);
        lobpcg.set_max_iter(80);
        lobpcg.set_diag_context("parallel-unit-rank-compression");
        lobpcg.diag(hpsi_func, spsi_func, psi.data(), eigens.data(), ethr);

        hsolver::DiagoIterAssist<TestT, TestDevice>::SCF_ITER = old_scf;

        for (int ib = 0; ib < nband_l; ++ib)
        {
            const int global_band = band_start + ib;
            ASSERT_NEAR(eigens[ib], e_ref[global_band], 2e-5)
                << "global_band=" << global_band;
        }

        std::vector<TestT> hpsi(nband_l * ld_psi), spsi(nband_l * ld_psi);
        matvec(hmat, psi.data(), hpsi.data(), npw, ld_psi, nband_l);
        matvec(smat, psi.data(), spsi.data(), npw, ld_psi, nband_l);

        for (int ib = 0; ib < nband_l; ++ib)
        {
            const int global_band = band_start + ib;
            for (int ig = npw; ig < ld_psi; ++ig)
                EXPECT_EQ(psi[ib * ld_psi + ig], TestT(0.0, 0.0));

            TestReal res2 = 0.0;
            for (int ig = 0; ig < npw; ++ig)
                res2 += std::norm(hpsi[ib * ld_psi + ig]
                                - eigens[ib] * spsi[ib * ld_psi + ig]);
            EXPECT_LT(std::sqrt(res2), 2e-4)
                << "global_band=" << global_band;
        }

        std::vector<TestT> global_psi(nband * ld_psi, {0.0, 0.0});
        std::vector<int> counts(nproc, 0), displs(nproc, 0);
        for (int ip = 0; ip < nproc; ++ip)
        {
            const int nlocal = nband / nproc + (ip < nband % nproc ? 1 : 0);
            counts[ip] = nlocal * ld_psi;
            displs[ip] = (nband / nproc * ip + std::min(ip, nband % nproc)) * ld_psi;
        }
        MPI_Allgatherv(psi.data(), nband_l * ld_psi, MPI_DOUBLE_COMPLEX,
                       global_psi.data(), counts.data(), displs.data(),
                       MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);

        std::vector<TestT> global_spsi(nband * ld_psi, {0.0, 0.0});
        matvec(smat, global_psi.data(), global_spsi.data(), npw, ld_psi, nband);
        for (int i = 0; i < nband; ++i)
        {
            for (int j = 0; j < nband; ++j)
            {
                TestT dot = {0.0, 0.0};
                for (int ig = 0; ig < npw; ++ig)
                    dot += std::conj(global_psi[i * ld_psi + ig])
                         * global_spsi[j * ld_psi + ig];
                const TestT target = (i == j) ? TestT(1.0, 0.0) : TestT(0.0, 0.0);
                EXPECT_NEAR(std::abs(dot - target), 0.0, 2e-6)
                    << "S-orth(" << i << "," << j << ")";
            }
        }
    }

    void run_generalized_band_parallel_operator_count()
    {
        int nproc = 1;
        int rank = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (nproc < 2)
            GTEST_SKIP() << "band-parallel LOBPCG test requires at least 2 MPI ranks";

        const int npw = 18;
        const int nband = 10;
        const int ld_psi = npw + 3;
        const int nband_l = nband / nproc + (rank < nband % nproc ? 1 : 0);
        const int band_start = nband / nproc * rank + std::min(rank, nband % nproc);

        std::vector<TestT> hmat, smat;
        std::vector<TestReal> prec, e_ref;
        build_generalized_problem(npw, 1.7, 0.03, 103, hmat, smat, prec, e_ref);

        std::vector<TestT> psi(nband_l * ld_psi, {0.0, 0.0});
        for (int ib = 0; ib < nband_l; ++ib)
        {
            const int global_band = band_start + ib;
            psi[ib * ld_psi + global_band] = {1.0, 0.0};
        }

        int hpsi_calls = 0;
        int spsi_calls = 0;
        auto hpsi_func = [&](TestT* psi_in, TestT* hpsi_out,
                              int ld_in, int nvec) {
            ++hpsi_calls;
            matvec(hmat, psi_in, hpsi_out, npw, ld_in, nvec);
        };
        auto spsi_func = [&](const TestT* psi_in, TestT* spsi_out,
                              int ld_in, int nvec) {
            ++spsi_calls;
            matvec(smat, psi_in, spsi_out, npw, ld_in, nvec);
        };

        std::vector<TestReal> eigens(nband_l, 0.0);
        std::vector<double> ethr(nband_l, 0.0);
        const int old_scf = hsolver::DiagoIterAssist<TestT, TestDevice>::SCF_ITER;
        hsolver::DiagoIterAssist<TestT, TestDevice>::SCF_ITER = 1;

        hsolver::DiagoLobpcg<TestT, TestDevice> lobpcg(prec.data());
        lobpcg.init_iter(nband, nband_l, ld_psi, npw);
        lobpcg.set_max_iter(2);
        lobpcg.set_diag_context("parallel-unit-operator-count");
        lobpcg.diag(hpsi_func, spsi_func, psi.data(), eigens.data(), ethr);

        hsolver::DiagoIterAssist<TestT, TestDevice>::SCF_ITER = old_scf;

        EXPECT_EQ(hpsi_calls, 2);
        EXPECT_EQ(spsi_calls, 2);
    }

#endif
};

// ============================================================================
// Test cases: various matrix sizes and band counts
// ============================================================================

TEST(DiagoLobpcgDetailTest, GeneralizedResidualGuardAllowsObservedBp2Growth)
{
    using hsolver::lobpcg_detail::compressed_guard_is_acceptable;
    using hsolver::lobpcg_detail::generalized_residual_growth_limit;
    using hsolver::lobpcg_detail::residual_guard_limit;
    using hsolver::lobpcg_detail::should_reject_residual_update;

    const double growth_limit = generalized_residual_growth_limit<double>();
    EXPECT_DOUBLE_EQ(growth_limit, 10.0);

    EXPECT_FALSE(should_reject_residual_update(5, 0.023084, 5, 0.025190, growth_limit));
    EXPECT_TRUE(compressed_guard_is_acceptable(5, 0.023084, 5, 0.025190, growth_limit));
    EXPECT_GT(residual_guard_limit(0.023084, growth_limit), 0.025190);

    EXPECT_TRUE(should_reject_residual_update(5, 0.023084, 5, 0.250000, growth_limit));
    EXPECT_FALSE(compressed_guard_is_acceptable(5, 0.023084, 5, 0.250000, growth_limit));
}

TEST_F(DiagoLobpcgTest, RejectsNonLocalEthrBandSize)
{
    const int npw = 20, nband = 6;
    std::vector<TestT> hmat;
    std::vector<TestReal> prec, e_ref;
    build_well_conditioned(npw, hmat, prec, e_ref);

    std::vector<TestT> psi(nband * npw, {0.0, 0.0});
    for (int ib = 0; ib < nband; ++ib)
        psi[ib * npw + ib] = {1.0, 0.0};

    auto hpsi_func = [&](TestT* psi_in, TestT* hpsi_out,
                          int ld_psi, int nvec) {
        matvec(hmat, psi_in, hpsi_out, npw, ld_psi, nvec);
    };
    auto spsi_func = [](const TestT* psi_in, TestT* spsi_out,
                         int ld_psi, int nvec) {
        std::copy(psi_in, psi_in + nvec * ld_psi, spsi_out);
    };

    std::vector<TestReal> eigens(nband, 0.0);
    std::vector<double> ethr(nband + 1, 1e-6);
    hsolver::DiagoLobpcg<TestT, TestDevice> lobpcg(prec.data());
    lobpcg.init_iter(nband, nband, npw, npw);
    lobpcg.set_diag_context("k=1/1, npw=20, nbands=6");

    try {
        lobpcg.diag(hpsi_func, spsi_func, psi.data(), eigens.data(), ethr);
        FAIL() << "Expected invalid_argument for non-local ethr_band size";
    } catch (const std::invalid_argument& e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find("local ethr_band size mismatch"), std::string::npos);
        EXPECT_NE(msg.find("size=7"), std::string::npos);
        EXPECT_NE(msg.find("required local bands=6"), std::string::npos);
        EXPECT_NE(msg.find("global bands=6"), std::string::npos);
        EXPECT_NE(msg.find("context={k=1/1, npw=20, nbands=6}"), std::string::npos);
    }
}

TEST_F(DiagoLobpcgTest, HpsiNonFiniteErrorIncludesContext)
{
    const int npw = 8, nband = 2;
    std::vector<TestReal> prec(npw, 1.0);
    std::vector<TestT> psi(nband * npw, {0.0, 0.0});
    for (int ib = 0; ib < nband; ++ib)
        psi[ib * npw + ib] = {1.0, 0.0};

    auto hpsi_func = [](TestT* psi_in, TestT* hpsi_out,
                         int ld_psi, int nvec) {
        std::copy(psi_in, psi_in + nvec * ld_psi, hpsi_out);
        hpsi_out[0] = {std::numeric_limits<double>::quiet_NaN(), 0.0};
    };
    auto spsi_func = [](const TestT* psi_in, TestT* spsi_out,
                         int ld_psi, int nvec) {
        std::copy(psi_in, psi_in + nvec * ld_psi, spsi_out);
    };

    std::vector<TestReal> eigens(nband, 0.0);
    std::vector<double> ethr(nband, 1e-6);
    hsolver::DiagoLobpcg<TestT, TestDevice> lobpcg(prec.data());
    lobpcg.init_iter(nband, nband, npw, npw);
    lobpcg.set_diag_context("k=4/4, npw=8, nbands=2, use_uspp=0");

    try {
        lobpcg.diag(hpsi_func, spsi_func, psi.data(), eigens.data(), ethr);
        FAIL() << "Expected runtime_error for non-finite hPsi";
    } catch (const std::runtime_error& e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find("LOBPCG hPsi produced non-finite values"), std::string::npos);
        EXPECT_NE(msg.find("context={k=4/4, npw=8, nbands=2, use_uspp=0}"), std::string::npos);
    }
}

TEST_F(DiagoLobpcgTest, GeneralizedNonPositiveOverlapFails)
{
    const int npw = 20, nband = 6;
    std::vector<TestT> hmat;
    std::vector<TestReal> prec, e_ref;
    build_well_conditioned(npw, hmat, prec, e_ref);

    std::vector<TestT> smat(npw * npw, {0.0, 0.0});
    for (int i = 0; i < npw; ++i)
        smat[idx(i, i, npw)] = {-1.0, 0.0};

    std::vector<TestT> psi(nband * npw, {0.0, 0.0});
    for (int ib = 0; ib < nband; ++ib)
        psi[ib * npw + ib] = {1.0, 0.0};

    auto hpsi_func = [&](TestT* psi_in, TestT* hpsi_out,
                          int ld_psi, int nvec) {
        matvec(hmat, psi_in, hpsi_out, npw, ld_psi, nvec);
    };
    auto spsi_func = [&](const TestT* psi_in, TestT* spsi_out,
                          int ld_psi, int nvec) {
        matvec(smat, psi_in, spsi_out, npw, ld_psi, nvec);
    };

    std::vector<TestReal> eigens(nband, 0.0);
    std::vector<double> ethr(nband, 1e-6);
    const int old_scf = hsolver::DiagoIterAssist<TestT, TestDevice>::SCF_ITER;
    hsolver::DiagoIterAssist<TestT, TestDevice>::SCF_ITER = 1;

    hsolver::DiagoLobpcg<TestT, TestDevice> lobpcg(prec.data());
    lobpcg.init_iter(nband, nband, npw, npw);
    EXPECT_THROW(lobpcg.diag(hpsi_func, spsi_func, psi.data(), eigens.data(), ethr),
                 std::runtime_error);

    hsolver::DiagoIterAssist<TestT, TestDevice>::SCF_ITER = old_scf;
}

TEST_F(DiagoLobpcgTest, GeneralizedUsppLikeOverlap)
{
    const int npw = 40, nband = 8;
    run_generalized_and_validate(npw, nband, npw + 7,
                                 1.5, 0.02, 7, 0.1,
                                 1e-5, 1e-7, 1e-4);
}

TEST_F(DiagoLobpcgTest, GeneralizedNearlyIdentityOverlap)
{
    const int npw = 60, nband = 10;
    run_generalized_and_validate(npw, nband, npw,
                                 1.05, 0.005, 17, 0.01,
                                 1e-5, 1e-7, 1e-4);
}

TEST_F(DiagoLobpcgTest, GeneralizedModerateCouplingOverlap)
{
    const int npw = 80, nband = 12;
    run_generalized_and_validate(npw, nband, npw + 5,
                                 2.0, 0.05, 29, 0.5,
                                 2e-5, 2e-7, 2e-4);
}

TEST_F(DiagoLobpcgTest, HegvdIsAccurateForWellConditionedOverlap)
{
    const int npw = 40;
    std::vector<TestT> hmat, smat;
    std::vector<TestReal> prec, e_ref;
    build_generalized_problem(npw, 1.5, 0.02, 7, hmat, smat, prec, e_ref);

    const HegvdMetrics metrics = solve_generalized_eigenvectors(npw, hmat, smat);
    EXPECT_EQ(metrics.info, 0);
    EXPECT_LT(metrics.max_rel_residual, 1e-12);
    EXPECT_LT(metrics.max_s_orth_error, 1e-12);
    EXPECT_LT(metrics.max_abs_coeff, 2.0);
}

TEST_F(DiagoLobpcgTest, HegvdAccuracyDegradesForNearlySingularOverlap)
{
    const int npw = 20;
    std::vector<TestT> hmat, smat;
    std::vector<TestReal> prec, e_ref;
    build_nearly_dependent_overlap_problem(npw, 1e-10, hmat, smat, prec, e_ref);

    const HegvdMetrics metrics = solve_generalized_eigenvectors(npw, hmat, smat);
    EXPECT_EQ(metrics.info, 0);
    EXPECT_GT(metrics.max_rel_residual, 1e-10);
    EXPECT_LT(metrics.max_rel_residual, 1e-5);
    EXPECT_GT(metrics.max_s_orth_error, 1e-8);
    EXPECT_GT(metrics.max_abs_coeff, 1e4);
}

#ifdef __MPI
TEST_F(DiagoLobpcgTest, GeneralizedBandParallelRankCompressedSubspace)
{
    run_generalized_band_parallel_and_validate();
}

TEST_F(DiagoLobpcgTest, BandParallelReusesProjectedSearchDirectionProducts)
{
    run_generalized_band_parallel_operator_count();
}

#endif

int main(int argc, char** argv)
{
    int nproc = 1, myrank = 0;

#ifdef __MPI
    int nproc_in_pool, kpar = 1, mypool, rank_in_pool;
    setupmpi(argc, argv, nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, kpar, mypool, rank_in_pool);
    const bool use_band_parallel_world = std::getenv("ABACUS_LOBPCG_TEST_BNDPAR") != nullptr;
    MPI_Comm_split(MPI_COMM_WORLD, use_band_parallel_world ? 0 : myrank, 0, &BP_WORLD);
    if (use_band_parallel_world)
    {
        GlobalV::MY_BNDGROUP = myrank;
        GlobalV::NPROC_IN_BNDGROUP = nproc;
        MPI_Comm_free(&POOL_WORLD);
        MPI_Comm_split(MPI_COMM_WORLD, myrank, 0, &POOL_WORLD);
    }
    GlobalV::NPROC_IN_POOL = nproc;
#endif

    testing::InitGoogleTest(&argc, argv);
    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();
    if (myrank != 0)
        delete listeners.Release(listeners.default_result_printer());

    int result = RUN_ALL_TESTS();
    if (myrank == 0 && result != 0)
        std::cout << "ERROR: some tests are not passed" << std::endl;

#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
