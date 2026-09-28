/**
 * diago_ppcg_parallel_test.cpp — MPI parallel test for DiagoPPCG.
 *
 * Distributes the rows of a diagonal matrix across MPI processes: each process
 * owns a slice of the diagonal and the corresponding rows of psi, computes the
 * partial Gram matrix / residual locally, and relies on the solver's pooled
 * MPI reductions (reduce_pool) to sum the partial results.  The eigenvalues of
 * the global diagonal matrix must be recovered identically on every process.
 *
 * Run with: mpirun -np <N> ./MODULE_HSOLVER_ppcg_parallel
 */

#include "../diago_ppcg.h"

#include "source_base/parallel_comm.h"
#include "source_base/parallel_global.h"
#include "source_basis/module_pw/test/test_tool.h"

#include "mpi.h"

#include <cmath>
#include <complex>
#include <cstdio>
#include <random>
#include <vector>

int main(int argc, char** argv)
{
    int nproc = 1;
    int myrank = 0;
    setupmpi(argc, argv, nproc, myrank);
    int nproc_in_pool = 0;
    int kpar = 1;
    int mypool = 0;
    int rank_in_pool = 0;
    divide_pools(nproc, myrank, nproc_in_pool, kpar, mypool, rank_in_pool);
    MPI_Comm_split(MPI_COMM_WORLD, myrank, 0, &BP_WORLD);

    using T    = std::complex<double>;
    using Real = double;

    const int nband = 3;
    const int n_local = 2 * nband;  // each process owns this many rows
    const int n_dim_total = nproc * n_local;

    // Diagonal H with entries 1, 2, ..., n_dim_total; process owns a slice.
    std::vector<Real> diag_local(n_local);
    std::vector<Real> prec(n_local);
    for (int i = 0; i < n_local; ++i)
    {
        diag_local[i] = Real(myrank * n_local + i + 1);
        prec[i] = diag_local[i];
    }

    std::vector<double> ethr(nband, 1e-8);

    // Random initial guess (fixed seed for reproducibility).
    std::mt19937 rng(42);
    std::uniform_real_distribution<Real> dist(-1.0, 1.0);
    std::vector<T> psi(n_local * nband, T(0.0, 0.0));
    for (int j = 0; j < nband; ++j)
    {
        for (int i = 0; i < n_local; ++i)
        {
            psi[i + j * n_local] = T(dist(rng), 0.0);
        }
    }

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 100,
        /* sbsize   = */ nband,
        /* rr_step  = */ nband,
        /* gamma_g0 = */ false);

    auto h_op = [&](T* in, T* out, int ld, int ncol) {
        for (int j = 0; j < ncol; ++j)
        {
            for (int i = 0; i < n_local; ++i)
            {
                out[i + j * ld] = diag_local[i] * in[i + j * ld];
            }
        }
    };

    std::vector<Real> eval(nband, 0.0);
    solver.diag(h_op, nullptr, n_local, nband, n_local,
                psi.data(), eval.data(), ethr, prec.data());

    // The lowest nband eigenvalues of the global diagonal matrix are 1..nband.
    int ok = 1;
    for (int i = 0; i < nband; ++i)
    {
        if (std::abs(eval[i] - Real(i + 1)) > 1e-6)
        {
            std::printf("rank %d: eval[%d] = %.12f != %d\n", myrank, i, eval[i], i + 1);
            ok = 0;
        }
    }

    int global_ok = 0;
    MPI_Allreduce(&ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    MPI_Finalize();
    if (myrank == 0)
    {
        if (global_ok == 1)
        {
            std::printf("PPCG MPI parallel test PASSED\n");
            return 0;
        }
        std::printf("PPCG MPI parallel test FAILED\n");
        return 1;
    }
    return global_ok == 1 ? 0 : 1;
}
