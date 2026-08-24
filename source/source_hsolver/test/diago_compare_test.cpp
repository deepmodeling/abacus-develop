/**
 * diago_compare_test.cpp — head-to-head comparison of the iterative
 * diagonalization solvers available in source_hsolver, on identical
 * random Hermitian matrices:
 *   - PPCG   (DiagoPPCG,  BLOCK_SUBSPACE)
 *   - CG     (DiagoCG,     band-by-band Polak-Ribiere)
 *   - BPCG   (DiagoBPCG,   block PCG)
 *   - Davidson (DiagoDavid)
 *
 * Every solver is fed the SAME Hamiltonian, the SAME initial guess and the
 * SAME per-band convergence threshold, so wall-clock time and the eigenvalue
 * error vs. a LAPACK reference are directly comparable.
 *
 * This is a benchmark/audit aid, not a correctness unit test: it is DISABLED
 * by default and must be run explicitly.
 */

#include "../diago_ppcg.h"
#include "../diago_cg.h"
#include "../diago_bpcg.h"
#include "../diago_david.h"

#include "source_base/module_external/lapack_connector.h"
#include "source_base/parallel_comm.h"
#include "source_base/global_variable.h"
#include "source_basis/module_pw/test/test_tool.h"

#include "mpi.h"

#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <malloc.h>
#include <random>
#include <string>
#include <vector>

using T = std::complex<double>;
using Real = double;

// Total heap memory currently allocated (bytes).  Used to compare the peak
// working memory of the solvers: PPCG keeps a bounded subspace, while
// Davidson grows its basis with the number of iterations.
static long heap_bytes()
{
    struct mallinfo2 mi = mallinfo2();
    return static_cast<long>(mi.uordblks) + static_cast<long>(mi.hblkhd);
}

extern "C" void zgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k, const T* alpha,
                       const T* a, const int* lda, const T* b, const int* ldb, const T* beta, T* c, const int* ldc);

static void dense_h_multiply(const T* H, int n, const T* in, T* out, int ld, int ncol)
{
    const T one(1.0, 0.0);
    const T zero(0.0, 0.0);
    zgemm_("N", "N", &n, &ncol, &n, &one, H, &n, in, &ld, &zero, out, &ld);
}

static void identity_s(const T* in, T* out, int ld, int ncol)
{
    for (int j = 0; j < ncol; ++j)
    {
        for (int i = 0; i < ld; ++i)
        {
            out[i + j * ld] = in[i + j * ld];
        }
    }
}

// Reference eigenvalues via LAPACK zheev (H is Hermitian, S = I).
static void ref_eigen(const T* H, int n, Real* e)
{
    std::vector<T> a(H, H + n * n);
    int lwork = 2 * n;
    std::vector<T> work(lwork);
    std::vector<Real> rwork(3 * n - 2);
    int info = 0;
    char jobz = 'N', uplo = 'U';
    zheev_(&jobz, &uplo, &n, a.data(), &n, e, work.data(), &lwork, rwork.data(), &info);
}

// Diagonal-dominant random Hermitian matrix (same recipe as the PPCG benchmark).
static void make_H(int n, int sparsity_pct, std::vector<T>& H, std::vector<Real>& prec)
{
    H.assign(n * n, T(0));
    std::mt19937 rng(unsigned(n * 100 + sparsity_pct));
    std::uniform_real_distribution<Real> dist(-1.0, 1.0);
    for (int i = 0; i < n; ++i)
    {
        for (int j = i; j < n; ++j)
        {
            if (i != j && (rng() % 100) < sparsity_pct)
            {
                continue;
            }
            Real val = (i == j) ? std::abs(dist(rng)) * n + 1.0 : dist(rng) * 0.5;
            H[i + j * n] = T(val, 0);
            if (i != j)
            {
                H[j + i * n] = T(val, 0);
            }
        }
    }
    prec.resize(n);
    for (int i = 0; i < n; ++i)
    {
        prec[i] = std::max(std::real(H[i + i * n]), 1e-6);
    }
}

// Random orthonormalized initial guess (identical for every solver).
static void make_psi(int n, int nband, std::vector<T>& psi)
{
    int ld = n;
    psi.assign(ld * nband, T(0));
    std::mt19937 rng(42);
    std::uniform_real_distribution<Real> dist(-1.0, 1.0);
    for (int j = 0; j < nband; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            psi[i + j * ld] = T(dist(rng), 0.0);
        }
    }
    for (int j = 0; j < nband; ++j)
    {
        for (int k = 0; k < j; ++k)
        {
            T d = 0;
            for (int i = 0; i < n; ++i)
            {
                d += std::conj(psi[i + k * ld]) * psi[i + j * ld];
            }
            for (int i = 0; i < n; ++i)
            {
                psi[i + j * ld] -= d * psi[i + k * ld];
            }
        }
        Real nr = 0;
        for (int i = 0; i < n; ++i)
        {
            nr += std::norm(psi[i + j * ld]);
        }
        nr = std::sqrt(nr);
        for (int i = 0; i < n; ++i)
        {
            psi[i + j * ld] /= nr;
        }
    }
}

// Rayleigh-Ritz subspace diagonalization used as CG's subspace_func.
static void rr_subspace(const T* H, int n, T* psi_in, T* psi_out, int ld, int nband)
{
    std::vector<T> hpsi(size_t(n) * nband, T(0));
    dense_h_multiply(H, n, psi_in, hpsi.data(), n, nband);

    // S_sub = Psi^H Psi (S = I), H_sub = Psi^H H Psi
    std::vector<T> s_sub(nband * nband, T(0)), h_sub(nband * nband, T(0));
    for (int i = 0; i < nband; ++i)
    {
        for (int j = 0; j < nband; ++j)
        {
            T s = 0, h = 0;
            for (int k = 0; k < n; ++k)
            {
                T pk = psi_in[k + i * ld];
                s += std::conj(pk) * psi_in[k + j * ld];
                h += std::conj(pk) * hpsi[k + j * n];
            }
            s_sub[i + j * nband] = s;
            h_sub[i + j * nband] = h;
        }
    }

    // Generalized Hermitian eigenproblem: H_sub C = S_sub C Lambda
    int lwork = 2 * nband;
    std::vector<T> work(lwork);
    std::vector<Real> rwork(3 * nband - 2);
    std::vector<Real> w(nband);
    int info = 0, itype = 1, nn = nband;
    char jobz = 'V', uplo = 'U';
    zhegv_(&itype, &jobz, &uplo, &nn, h_sub.data(), &nn, s_sub.data(), &nn, w.data(), work.data(), &lwork, rwork.data(),
           &info);

    // psi_out = psi_in * C  (C now holds the eigenvectors in h_sub)
    for (int j = 0; j < nband; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            T acc = 0;
            for (int c = 0; c < nband; ++c)
            {
                acc += psi_in[i + c * ld] * h_sub[c + j * nband];
            }
            psi_out[i + j * ld] = acc;
        }
    }
}

struct Result
{
    double wall_s = 0.0;
    double avg_iter = -1.0; // -1 when the solver does not report it
    double max_err = 0.0;   // max |eval_i - ref_i| over the requested bands
    long mem_bytes = 0;     // peak heap memory allocated by the solver
    bool ok = false;
};

static Result run_ppcg(const std::vector<T>& H, int n, int nband, const std::vector<Real>& prec,
                       const std::vector<T>& psi0, const std::vector<double>& ethr, const Real* ref)
{
    Result r;
    std::vector<T> psi = psi0;
    std::vector<Real> eval(nband, 0.0);
    long mem0 = heap_bytes();
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(1e-8, 500, nband, std::min(nband, 4), false,
                                                                   hsolver::PpcgStrategy::BLOCK_SUBSPACE);
    auto h_op = [&H, n](T* in, T* out, int ld, int nc) { dense_h_multiply(H.data(), n, in, out, ld, nc); };
    auto t0 = std::chrono::high_resolution_clock::now();
    double avg = solver.diag(h_op, nullptr, n, nband, n, psi.data(), eval.data(), ethr, prec.data());
    auto t1 = std::chrono::high_resolution_clock::now();
    r.wall_s = std::chrono::duration<double>(t1 - t0).count();
    r.avg_iter = avg;
    r.mem_bytes = heap_bytes() - mem0;
    for (int i = 0; i < nband; ++i)
    {
        r.max_err = std::max(r.max_err, std::abs(eval[i] - ref[i]));
    }
    r.ok = true;
    return r;
}

static Result run_cg(const std::vector<T>& H, int n, int nband, const std::vector<Real>& prec,
                     const std::vector<T>& psi0, const std::vector<double>& ethr, const Real* ref)
{
    Result r;
    std::vector<T> psi = psi0;
    std::vector<Real> eval(nband, 0.0);
    auto subspace_func = [&H, n](T* psi_in, T* psi_out, int ld, int nband, bool) {
        rr_subspace(H.data(), n, psi_in, psi_out, ld, nband);
    };
    long mem0 = heap_bytes();
    hsolver::DiagoCG<T, hsolver::base_device::DEVICE_CPU> cg("pw", "scf", true, subspace_func, 1e-8, 500, 1);
    auto h_op = [&H, n](T* in, T* out, int ld, int nc) { dense_h_multiply(H.data(), n, in, out, ld, nc); };
    auto s_op = [](T* in, T* out, int ld, int nc) { identity_s(in, out, ld, nc); };
    auto t0 = std::chrono::high_resolution_clock::now();
    double avg = cg.diag(h_op, s_op, n, nband, n, psi.data(), eval.data(), ethr, prec.data());
    auto t1 = std::chrono::high_resolution_clock::now();
    r.wall_s = std::chrono::duration<double>(t1 - t0).count();
    r.avg_iter = avg;
    r.mem_bytes = heap_bytes() - mem0;
    for (int i = 0; i < nband; ++i)
    {
        r.max_err = std::max(r.max_err, std::abs(eval[i] - ref[i]));
    }
    r.ok = true;
    return r;
}

static Result run_bpcg(const std::vector<T>& H, int n, int nband, const std::vector<Real>& prec,
                       const std::vector<T>& psi0, const std::vector<double>& ethr, const Real* ref)
{
    Result r;
    std::vector<T> psi = psi0;
    std::vector<Real> eval(nband, 0.0);
    long mem0 = heap_bytes();
    hsolver::DiagoBPCG<T, hsolver::base_device::DEVICE_CPU> bpcg(prec.data());
    bpcg.init_iter(nband, nband, n, n);
    auto h_op = [&H, n](T* in, T* out, int ld, int nc) { dense_h_multiply(H.data(), n, in, out, ld, nc); };
    // BPCG::diag() is a single block-CG sweep; iterate until convergence.
    int it = 0;
    auto t0 = std::chrono::high_resolution_clock::now();
    for (; it < 200; ++it)
    {
        bpcg.diag(h_op, psi.data(), eval.data(), ethr);
        double err = 0.0;
        for (int i = 0; i < nband; ++i)
        {
            err = std::max(err, std::abs(eval[i] - ref[i]));
        }
        if (err < ethr[0])
        {
            break;
        }
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    r.wall_s = std::chrono::duration<double>(t1 - t0).count();
    r.avg_iter = it;
    r.mem_bytes = heap_bytes() - mem0;
    for (int i = 0; i < nband; ++i)
    {
        r.max_err = std::max(r.max_err, std::abs(eval[i] - ref[i]));
    }
    r.ok = true;
    return r;
}

static Result run_dav(const std::vector<T>& H, int n, int nband, const std::vector<Real>& prec,
                      const std::vector<T>& psi0, const std::vector<double>& ethr, const Real* ref)
{
    Result r;
    std::vector<T> psi = psi0;
    std::vector<Real> eval(nband, 0.0);
    hsolver::diag_comm_info comm(MPI_COMM_WORLD, 0, 1);
    long mem0 = heap_bytes();
    hsolver::DiagoDavid<T, hsolver::base_device::DEVICE_CPU> dav(prec.data(), nband, n, 4, comm);
    auto h_op = [&H, n](T* in, T* out, int ld, int nc) { dense_h_multiply(H.data(), n, in, out, ld, nc); };
    auto s_op = [](T* in, T* out, int ld, int nc) { identity_s(in, out, ld, nc); };
    auto t0 = std::chrono::high_resolution_clock::now();
    int it = dav.diag(h_op, s_op, n, psi.data(), eval.data(), ethr, 500);
    auto t1 = std::chrono::high_resolution_clock::now();
    r.wall_s = std::chrono::duration<double>(t1 - t0).count();
    r.avg_iter = it;
    r.mem_bytes = heap_bytes() - mem0;
    for (int i = 0; i < nband; ++i)
    {
        r.max_err = std::max(r.max_err, std::abs(eval[i] - ref[i]));
    }
    r.ok = true;
    return r;
}

int main(int argc, char** argv)
{
    int nproc = 1, myrank = 0;
    int nproc_in_pool, kpar = 1, mypool, rank_in_pool;
    setupmpi(argc, argv, nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, kpar, mypool, rank_in_pool);
    MPI_Comm_split(MPI_COMM_WORLD, myrank, 0, &BP_WORLD);

    struct Case
    {
        int n;
        int nband;
        int sparsity;
    };
    // Without arguments a small default grid is used.  To benchmark a single
    // (possibly large) problem, pass:  <n> <nband> <sparsity_pct>
    std::vector<Case> cases;
    if (argc >= 4)
    {
        cases.push_back({std::atoi(argv[1]), std::atoi(argv[2]), std::atoi(argv[3])});
    }
    else
    {
        cases = {
            {50, 10, 0}, {50, 10, 60}, {100, 10, 60}, {200, 10, 80}, {500, 10, 80},
        };
    }

    std::printf("\n=== Solver comparison (identical H, psi0, ethr) ===\n");
    std::printf("%-5s %-5s %-6s %-10s %-14s %-12s %-10s %-12s\n", "n", "nband", "spars", "solver", "wall_time(s)", "avg_iter",
                "max_err", "mem(MB)");
    std::printf("---------------------------------------------------------------------------\n");

    for (const auto& c : cases)
    {
        std::vector<T> H;
        std::vector<Real> prec;
        make_H(c.n, c.sparsity, H, prec);
        std::vector<Real> ref(c.n, 0.0);
        ref_eigen(H.data(), c.n, ref.data());
        std::vector<T> psi0;
        make_psi(c.n, c.nband, psi0);
        std::vector<double> ethr(c.nband, 1e-6);

        Result r_ppcg = run_ppcg(H, c.n, c.nband, prec, psi0, ethr, ref.data());
        Result r_cg = run_cg(H, c.n, c.nband, prec, psi0, ethr, ref.data());
        Result r_bpcg = run_bpcg(H, c.n, c.nband, prec, psi0, ethr, ref.data());
        Result r_dav = run_dav(H, c.n, c.nband, prec, psi0, ethr, ref.data());

        std::printf("%-5d %-5d %-6d %-10s %-14.5f %-12.1f %-10.2e %-12.2f\n", c.n, c.nband, c.sparsity, "PPCG", r_ppcg.wall_s,
                    r_ppcg.avg_iter, r_ppcg.max_err, r_ppcg.mem_bytes / 1048576.0);
        std::printf("%-5s %-5s %-6s %-10s %-14.5f %-12.1f %-10.2e %-12.2f\n", "", "", "", "CG", r_cg.wall_s, r_cg.avg_iter,
                    r_cg.max_err, r_cg.mem_bytes / 1048576.0);
        std::printf("%-5s %-5s %-6s %-10s %-14.5f %-12.1f %-10.2e %-12.2f\n", "", "", "", "BPCG", r_bpcg.wall_s,
                    r_bpcg.avg_iter, r_bpcg.max_err, r_bpcg.mem_bytes / 1048576.0);
        std::printf("%-5s %-5s %-6s %-10s %-14.5f %-12.1f %-10.2e %-12.2f\n", "", "", "", "Davidson", r_dav.wall_s,
                    r_dav.avg_iter, r_dav.max_err, r_dav.mem_bytes / 1048576.0);
        std::printf("---------------------------------------------------------------------------\n");
    }

    MPI_Finalize();
    return 0;
}
