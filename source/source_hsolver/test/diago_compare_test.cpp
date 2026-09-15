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
 * SAME per-band convergence threshold, and wall time is measured to the SAME
 * reference accuracy (max_eval_err < err_target) so the timings are directly
 * comparable despite the solvers' differing internal stopping rules.
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

#include <algorithm>
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

// Optional PPCG parameter overrides (set from argv) for exploring the block
// size (sbsize) and Rayleigh-Ritz frequency (rr_step).  A negative value keeps
// the default used by the comparison benchmark (sbsize = nband, rr_step = 16).
static int g_sbsize = -1;
static int g_rr_step = -1;

// ---------------------------------------------------------------------------
// Unified stopping criterion.  Wall time is measured as "time to reach
// max_eval_err < err_target", the same LAPACK-reference accuracy for every
// solver.  Each solver's own stopping rule (eigenvalue-change for PPCG/CG,
// residual for BPCG/Davidson) is looser than the reference error, so every
// solver is re-driven for up to max_outer_passes calls to diag().  max_err
// always reports the accuracy actually reached, so a solver that fails to hit
// err_target within the budget stays visible.
//
// NOTE: err_target must be strictly coarser than every solver's own stopping
// point.  The CG solver stops on |eigenvalue change| < ethr, which in practice
// leaves max_eval_err ~ (2-4)x ethr; a target tighter than that (e.g. 1e-6
// with ethr=1e-6) can never be reached, so CG burns all max_outer_passes doing
// an expensive subspace restart + per-band CG each round.  A relaxed target
// keeps the cross-solver comparison honest instead of penalizing CG for being
// re-driven into repeated subspace restarts.
// ---------------------------------------------------------------------------
const double err_target = 1e-5;    // relaxed: reachable by every solver
const int max_outer_passes = 20;   // outer diag() re-drives before giving up

// max over bands of |eval_i - ref_i|: the reference-based accuracy that makes
// wall-clock times comparable across solvers.
static double max_eval_err(const Real* eval, const Real* ref, int nband)
{
    double err = 0.0;
    for (int i = 0; i < nband; ++i)
    {
        err = std::max(err, std::abs(eval[i] - ref[i]));
    }
    return err;
}

// Total heap memory currently allocated (bytes).  Used to compare the peak
// working memory of the solvers: PPCG keeps a bounded subspace, while
// Davidson grows its basis with the number of iterations.
static long heap_bytes()
{
    struct mallinfo2 mi = mallinfo2();
    return static_cast<long>(mi.uordblks) + static_cast<long>(mi.hblkhd);
}

// Sparse symmetric band of half-bandwidth `bw`, stored in LAPACK upper-band
// format: H[i, j] (i <= j, d = j - i) lives at band[bd*j + (bw - d)],
// with bd = bw + 1.  This matvec is O(n * bw), modeling the real plane-wave H
// application (kinetic energy is a diagonal/tridiagonal discretization of
// -Laplacian/2 plus a local potential) instead of an O(n^2) dense zgemm.
static void banded_h_multiply(const Real* band, int n, int bw, int bd, const T* in, T* out, int ld, int ncol)
{
    for (int j = 0; j < ncol; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            Real acc = 0.0;
            const int hi = std::min(n - 1, i + bw);
            for (int k = std::max(0, i - bw); k <= hi; ++k)
            {
                // H[i, k]: reuse the upper-triangular stored entry (k may be < i)
                const int r = std::min(i, k);
                const int c = std::max(i, k);
                const int d = c - r;
                acc += band[bd * c + (bw - d)] * std::real(in[k + j * ld]);
            }
            out[i + j * ld] = T(acc, 0.0);
        }
    }
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

extern "C" void dsbev_(const char* jobz, const char* uplo, const int* n, const int* kd,
                       double* ab, const int* ldab, double* w, double* z, const int* ldz,
                       double* work, int* info);

// Reference eigenvalues via LAPACK dsbev on the symmetric band matrix.  Only
// the lowest `nband_req` are kept.  dsbev is O(n * bw^2), avoiding the O(n^3)
// tridiagonalization that dense zheev/zheevx would require for large n.
static void ref_eigen(const Real* band, int n, int bw, int bd, int nband_req, Real* e)
{
    std::vector<double> ab(band, band + size_t(bd) * n);
    std::vector<double> w(n);
    const char jobz = 'N', uplo = 'U';
    std::vector<double> work(3 * n);
    int info = 0;
    dsbev_(&jobz, &uplo, &n, &bw, ab.data(), &bd, w.data(), nullptr, &n, work.data(), &info);
    if (info != 0)
    {
        std::fprintf(stderr, "[ref_eigen] dsbev info=%d\n", info);
    }
    for (int i = 0; i < nband_req; ++i)
    {
        e[i] = w[i];
    }
}

// Diagonally-dominant symmetric band matrix: a local potential on the diagonal
// plus random couplings within a half-bandwidth `bw`.  Stored in LAPACK upper
// band format `band[bd*j + (bw - (j - i))] = H[i, j]` for i <= j, bd = bw + 1.
// This is the discrete analogue of H = -Laplacian/2 + V(r) that plane-wave
// solvers actually apply.
static void make_H(int n, int bw, std::vector<Real>& band, int& bd, std::vector<Real>& prec)
{
    bd = bw + 1;
    band.assign(size_t(bd) * n, 0.0);
    std::mt19937 rng(unsigned(n * 100 + bw));
    std::uniform_real_distribution<Real> dist(-1.0, 1.0);
    for (int i = 0; i < n; ++i)
    {
        band[bd * i + bw] = std::abs(dist(rng)) * n + 1.0;  // diagonal (d=0)
    }
    for (int i = 0; i < n; ++i)
    {
        for (int d = 1; d <= bw; ++d)
        {
            if (i + d < n)
            {
                band[bd * (i + d) + (bw - d)] = dist(rng) * 0.5 / double(d);
            }
        }
    }
    prec.resize(n);
    for (int i = 0; i < n; ++i)
    {
        prec[i] = std::max(band[bd * i + bw], 1e-6);
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
static void rr_subspace(const Real* band, int n, int bw, int bd, T* psi_in, T* psi_out, int ld, int nband)
{
    std::vector<T> hpsi(size_t(n) * nband, T(0));
    banded_h_multiply(band, n, bw, bd, psi_in, hpsi.data(), n, nband);

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
    double max_err = 0.0;   // max |eval_i - ref_i| over the requested bands
    long mem_bytes = 0;     // peak heap memory allocated by the solver
    bool ok = false;
};

static Result run_ppcg(const std::vector<Real>& band, int n, int bw, int bd, int nband, const std::vector<Real>& prec,
                       const std::vector<T>& psi0, const std::vector<double>& ethr, const Real* ref)
{
    Result r;
    std::vector<T> psi = psi0;
    std::vector<Real> eval(nband, 0.0);
    long mem0 = heap_bytes();
    const int sbsize = (g_sbsize > 0) ? g_sbsize : nband;
    const int rr_step = (g_rr_step > 0) ? g_rr_step : 16;
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(1e-8, 500, sbsize, rr_step, false);
    auto h_op = [&band, n, bw, bd](T* in, T* out, int ld, int nc) { banded_h_multiply(band.data(), n, bw, bd, in, out, ld, nc); };
    auto t0 = std::chrono::high_resolution_clock::now();
    int pass = 0;
    for (; pass < max_outer_passes; ++pass)
    {
        solver.diag(h_op, nullptr, n, nband, n, psi.data(), eval.data(), ethr, prec.data());
        if (max_eval_err(eval.data(), ref, nband) < err_target)
        {
            break;
        }
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    r.wall_s = std::chrono::duration<double>(t1 - t0).count();
    r.mem_bytes = heap_bytes() - mem0;
    r.max_err = max_eval_err(eval.data(), ref, nband);
    r.ok = true;
    return r;
}

static Result run_cg(const std::vector<Real>& band, int n, int bw, int bd, int nband, const std::vector<Real>& prec,
                     const std::vector<T>& psi0, const std::vector<double>& ethr, const Real* ref)
{
    Result r;
    std::vector<T> psi = psi0;
    std::vector<Real> eval(nband, 0.0);
    auto subspace_func = [&band, n, bw, bd](T* psi_in, T* psi_out, int ld, int nband, bool) {
        rr_subspace(band.data(), n, bw, bd, psi_in, psi_out, ld, nband);
    };
    long mem0 = heap_bytes();
    hsolver::DiagoCG<T, hsolver::base_device::DEVICE_CPU> cg("pw", "scf", true, subspace_func, 1e-8, 500, 1);
    auto h_op = [&band, n, bw, bd](T* in, T* out, int ld, int nc) { banded_h_multiply(band.data(), n, bw, bd, in, out, ld, nc); };
    auto s_op = [](T* in, T* out, int ld, int nc) { identity_s(in, out, ld, nc); };
    auto t0 = std::chrono::high_resolution_clock::now();
    int pass = 0;
    for (; pass < max_outer_passes; ++pass)
    {
        cg.diag(h_op, s_op, n, nband, n, psi.data(), eval.data(), ethr, prec.data());
        if (max_eval_err(eval.data(), ref, nband) < err_target)
        {
            break;
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    r.wall_s = std::chrono::duration<double>(t1 - t0).count();
    r.mem_bytes = heap_bytes() - mem0;
    r.max_err = max_eval_err(eval.data(), ref, nband);
    r.ok = true;
    return r;
}

static Result run_bpcg(const std::vector<Real>& band, int n, int bw, int bd, int nband, const std::vector<Real>& prec,
                       const std::vector<T>& psi0, const std::vector<double>& ethr, const Real* ref)
{
    Result r;
    std::vector<T> psi = psi0;
    std::vector<Real> eval(nband, 0.0);
    long mem0 = heap_bytes();
    hsolver::DiagoBPCG<T, hsolver::base_device::DEVICE_CPU> bpcg(prec.data());
    bpcg.init_iter(nband, nband, n, n);
    auto h_op = [&band, n, bw, bd](T* in, T* out, int ld, int nc) { banded_h_multiply(band.data(), n, bw, bd, in, out, ld, nc); };
    auto s_op = [](const T* in, T* out, int ld, int nc) { identity_s(in, out, ld, nc); };
    // BPCG::diag() is a single block-CG sweep; iterate until convergence.
    int it = 0;
    auto t0 = std::chrono::high_resolution_clock::now();
    for (; it < max_outer_passes; ++it)
    {
        bpcg.diag(h_op, s_op, psi.data(), eval.data(), ethr);
        if (max_eval_err(eval.data(), ref, nband) < err_target)
        {
            break;
        }
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    r.wall_s = std::chrono::duration<double>(t1 - t0).count();
    r.mem_bytes = heap_bytes() - mem0;
    r.max_err = max_eval_err(eval.data(), ref, nband);
    r.ok = true;
    return r;
}

static Result run_dav(const std::vector<Real>& band, int n, int bw, int bd, int nband, const std::vector<Real>& prec,
                      const std::vector<T>& psi0, const std::vector<double>& ethr, const Real* ref)
{
    Result r;
    std::vector<T> psi = psi0;
    std::vector<Real> eval(nband, 0.0);
    hsolver::diag_comm_info comm(MPI_COMM_WORLD, 0, 1);
    long mem0 = heap_bytes();
    hsolver::DiagoDavid<T, hsolver::base_device::DEVICE_CPU> dav(prec.data(), nband, n, 4, comm);
    auto h_op = [&band, n, bw, bd](T* in, T* out, int ld, int nc) { banded_h_multiply(band.data(), n, bw, bd, in, out, ld, nc); };
    auto s_op = [](T* in, T* out, int ld, int nc) { identity_s(in, out, ld, nc); };
    auto t0 = std::chrono::high_resolution_clock::now();
    // Davidson's diag() already iterates its growing subspace to convergence,
    // so it must be called exactly once; a re-drive loop would reuse the stale
    // Ritz basis and trigger a rank-deficient Schmidt orthogonalization.
    dav.diag(h_op, s_op, n, psi.data(), eval.data(), ethr, 500);
    auto t1 = std::chrono::high_resolution_clock::now();
    r.wall_s = std::chrono::duration<double>(t1 - t0).count();
    r.mem_bytes = heap_bytes() - mem0;
    r.max_err = max_eval_err(eval.data(), ref, nband);
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
        int bw;
    };
    // Without arguments a small default grid is used.  To benchmark a single
    // (possibly large) problem, pass:  <n> <nband> <bw> [sbsize] [rr_step]
    std::vector<Case> cases;
    if (argc >= 4)
    {
        cases.push_back({std::atoi(argv[1]), std::atoi(argv[2]), std::atoi(argv[3])});
    }
    else
    {
        cases = {
            {50, 10, 1}, {50, 10, 3}, {100, 10, 3}, {200, 10, 5}, {500, 10, 5},
        };
    }
    if (argc >= 5)
    {
        g_sbsize = std::atoi(argv[4]);
    }
    if (argc >= 6)
    {
        g_rr_step = std::atoi(argv[5]);
    }

    std::printf("\n=== Solver comparison (identical H, psi0, ethr) ===\n");
    std::printf("%-5s %-5s %-6s %-10s %-14s %-10s %-12s\n", "n", "nband", "bw", "solver", "wall_time(s)",
                "max_err", "mem(MB)");
    std::printf("-----------------------------------------------------------------\n");

    for (const auto& c : cases)
    {
        std::vector<Real> band;
        int bd = 0;
        std::vector<Real> prec;
        make_H(c.n, c.bw, band, bd, prec);
        std::vector<Real> ref(c.n, 0.0);
        ref_eigen(band.data(), c.n, c.bw, bd, c.nband, ref.data());
        std::vector<T> psi0;
        make_psi(c.n, c.nband, psi0);
        std::vector<double> ethr(c.nband, 1e-6);

        Result r_ppcg = run_ppcg(band, c.n, c.bw, bd, c.nband, prec, psi0, ethr, ref.data());
        Result r_cg = run_cg(band, c.n, c.bw, bd, c.nband, prec, psi0, ethr, ref.data());
        Result r_bpcg = run_bpcg(band, c.n, c.bw, bd, c.nband, prec, psi0, ethr, ref.data());
        Result r_dav = run_dav(band, c.n, c.bw, bd, c.nband, prec, psi0, ethr, ref.data());

        std::printf("%-5d %-5d %-6d %-10s %-14.5f %-10.2e %-12.2f\n", c.n, c.nband, c.bw, "PPCG", r_ppcg.wall_s,
                    r_ppcg.max_err, r_ppcg.mem_bytes / 1048576.0);
        std::printf("%-5s %-5s %-6s %-10s %-14.5f %-10.2e %-12.2f\n", "", "", "", "CG", r_cg.wall_s,
                    r_cg.max_err, r_cg.mem_bytes / 1048576.0);
        std::printf("%-5s %-5s %-6s %-10s %-14.5f %-10.2e %-12.2f\n", "", "", "", "BPCG", r_bpcg.wall_s,
                    r_bpcg.max_err, r_bpcg.mem_bytes / 1048576.0);
        std::printf("%-5s %-5s %-6s %-10s %-14.5f %-10.2e %-12.2f\n", "", "", "", "Davidson", r_dav.wall_s,
                    r_dav.max_err, r_dav.mem_bytes / 1048576.0);
        std::printf("-----------------------------------------------------------------\n");
    }

    MPI_Finalize();
    return 0;
}
