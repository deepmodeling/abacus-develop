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
#include "source_base/parallel_reduce.h"
#include "source_basis/module_pw/test/test_tool.h"

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// A local block can discover two low roots while another block already holds
// an exact higher root. The discarded low direction must survive global root
// selection, including when its rows and projection norms are distributed.
static bool check_competing_roots(const int rank_in_pool, const int nproc_in_pool)
{
    using Scalar = std::complex<double>;
    const int bands = 2;
    const int local_rows = nproc_in_pool == 1 ? 4 : 2;
    const int size = local_rows * bands;
    std::vector<double> diagonal(local_rows);
    std::vector<double> preconditioner(local_rows, 1.0);
    std::vector<double> thresholds(bands, 1.0e-10);
    std::vector<Scalar> vectors(size, Scalar(0));
    const double scale = 1.0 / std::sqrt(3.0);
    for (int row = 0; row < local_rows; ++row)
    {
        const int global_row = rank_in_pool * local_rows + row;
        diagonal[row] = global_row < 3 ? global_row + 1.0 : 97.0 + global_row;
        if (global_row == 0 || global_row == 1 || global_row == 3)
        {
            vectors[row] = Scalar(scale);
        }
        if (global_row == 2)
        {
            vectors[row + local_rows] = Scalar(1);
        }
    }
    bool operator_capacity_ok = true;
    const auto h_op = [&](Scalar* in, Scalar* out, int ld, int ncol) {
        operator_capacity_ok = operator_capacity_ok && ncol <= bands;
        for (int column = 0; column < ncol; ++column)
        {
            for (int row = 0; row < local_rows; ++row)
            {
                out[row + column * ld] = diagonal[row] * in[row + column * ld];
            }
        }
    };
    hsolver::DiagoPPCG<Scalar, base_device::DEVICE_CPU> solver(1.0e-10, 80, 1, 4);
    std::vector<double> eigenvalues(bands);
    solver.diag(h_op, nullptr, local_rows, bands, local_rows,
                vectors.data(), eigenvalues.data(), thresholds, preconditioner.data());

    std::vector<double> residuals(bands, 0.0);
    std::vector<Scalar> overlap(bands * bands, Scalar(0));
    for (int column = 0; column < bands; ++column)
    {
        for (int row = 0; row < local_rows; ++row)
        {
            const Scalar value = vectors[row + column * local_rows];
            const Scalar residual = (diagonal[row] - eigenvalues[column]) * value;
            residuals[column] += std::norm(residual);
            for (int other = 0; other < bands; ++other)
            {
                overlap[other + column * bands] += std::conj(vectors[row + other * local_rows]) * value;
            }
        }
    }
    Parallel_Reduce::reduce_pool(residuals.data(), bands);
    const int overlap_size = bands * bands;
    Parallel_Reduce::reduce_pool(overlap.data(), overlap_size);
    bool correct = solver.converged() && operator_capacity_ok && solver.root_recovery_count() > 0;
    for (int column = 0; column < bands; ++column)
    {
        const double expected = column + 1.0;
        correct = correct && std::isfinite(eigenvalues[column]);
        correct = correct && std::abs(eigenvalues[column] - expected) < 1.0e-8;
        correct = correct && std::isfinite(residuals[column]) && residuals[column] < 1.0e-16;
        for (int other = 0; other < bands; ++other)
        {
            const Scalar expected_overlap = other == column ? Scalar(1) : Scalar(0);
            const double error = std::abs(overlap[other + column * bands] - expected_overlap);
            correct = correct && std::isfinite(error) && error < 1.0e-8;
        }
    }
    if (!correct)
    {
        std::fprintf(stderr, "pool rank %d: competing-root regression returned %.12g, %.12g\n",
                     rank_in_pool, eigenvalues[0], eigenvalues[1]);
    }
    return correct;
}

int main(int argc, char** argv)
{
    int nproc = 1;
    int myrank = 0;
    setupmpi(argc, argv, nproc, myrank);

    int nproc_in_pool = 0;
    int kpar = argc > 1 ? std::atoi(argv[1]) : 1;
    int mypool = 0;
    int rank_in_pool = 0;
    divide_pools(nproc, myrank, nproc_in_pool, kpar, mypool, rank_in_pool);
    const bool world_root_is_correct = Parallel_Reduce::is_world_root() == (myrank == 0);
    BP_WORLD = MPI_COMM_WORLD;

    using T    = std::complex<double>;
    using Real = double;

    const int nband = 3;
    const int n_local = 2 * nband;  // each process owns this many rows
    const int n_dim_total = nproc_in_pool * n_local;

    // With one pool, only rank zero requests a trace to exercise inconsistent
    // process environments. With multiple pools, every rank requests the same
    // path to verify that only the MPI world root writes it.
    const char* trace_environment = "ABACUS_PPCG_RESIDUAL_TRACE";
    const std::string trace_path = "ppcg_residual_trace_parallel_test.csv";
    int trace_environment_ok = 1;
    if (myrank == 0)
    {
        std::remove(trace_path.c_str());
        trace_environment_ok = ::setenv(trace_environment, trace_path.c_str(), 1) == 0 ? 1 : 0;
    }
    else if (kpar > 1)
    {
        trace_environment_ok = ::setenv(trace_environment, trace_path.c_str(), 1) == 0 ? 1 : 0;
    }
    else
    {
        trace_environment_ok = ::unsetenv(trace_environment) == 0 ? 1 : 0;
    }

    // Diagonal H with entries 1, 2, ..., n_dim_total; process owns a slice.
    std::vector<Real> diag_local(n_local);
    std::vector<Real> prec(n_local);
    for (int i = 0; i < n_local; ++i)
    {
        diag_local[i] = Real(rank_in_pool * n_local + i + 1);
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
        /* sbsize   = */ 2,
        /* rr_step  = */ nband);

    int h_calls = 0;
    auto h_op = [&](T* in, T* out, int ld, int ncol) {
        ++h_calls;
        for (int j = 0; j < ncol; ++j)
        {
            for (int i = 0; i < n_local; ++i)
            {
                out[i + j * ld] = diag_local[i] * in[i + j * ld];
            }
        }
    };

    std::vector<Real> eval(nband, 0.0);
    // One invalid rank must stop every rank in its pool before H/S collectives.
    // Otherwise the valid ranks would hang while the invalid rank unwinds.
    std::vector<Real> invalid_prec = prec;
    if (rank_in_pool == 0)
    {
        invalid_prec[0] = std::numeric_limits<Real>::quiet_NaN();
    }
    bool collective_validation_failed = false;
    try
    {
        solver.diag(h_op, nullptr, n_local, nband, n_local,
                    psi.data(), eval.data(), ethr, invalid_prec.data());
    }
    catch (const std::invalid_argument&)
    {
        collective_validation_failed = true;
    }
    const bool validation_before_operator = collective_validation_failed && h_calls == 0;

    // Only one rank advertises an unrepresentable local [X,W,P] workspace.
    // The wavefunction product still fits, and the tiny real buffers must not
    // be read on that rank or passed to H on any of its peers.
    int checked_ld = n_local;
    if (rank_in_pool == 0)
    {
        checked_ld = std::numeric_limits<int>::max() / nband;
    }
    bool collective_size_failed = false;
    try
    {
        solver.diag(h_op, nullptr, checked_ld, nband, n_local,
                    psi.data(), eval.data(), ethr, prec.data());
    }
    catch (const std::invalid_argument&)
    {
        collective_size_failed = true;
    }
    const bool size_validation_before_operator = collective_size_failed && h_calls == 0;

    solver.diag(h_op, nullptr, n_local, nband, n_local,
                psi.data(), eval.data(), ethr, prec.data());

    // The lowest nband eigenvalues of the global diagonal matrix are 1..nband.
    int ok = 1;
    if (!world_root_is_correct
        || !validation_before_operator
        || !size_validation_before_operator
        || trace_environment_ok != 1)
    {
        ok = 0;
    }
    if (!world_root_is_correct)
    {
        std::fprintf(stderr, "rank %d: MPI world-root detection is inconsistent\n", myrank);
    }
    if (!validation_before_operator)
    {
        std::fprintf(stderr, "rank %d: invalid input was not rejected before operator calls\n", myrank);
    }
    if (!size_validation_before_operator)
    {
        std::fprintf(stderr, "rank %d: workspace overflow was not rejected before operator calls\n", myrank);
    }
    if (!solver.converged()
        || solver.active_band_count() != 0
        || solver.iteration_limit_reached()
        || !std::isfinite(solver.max_residual())
        || solver.max_residual() > 1.0e-8)
    {
        const int converged = solver.converged() ? 1 : 0;
        const int active_band_count = solver.active_band_count();
        const Real max_residual = solver.max_residual();
        std::fprintf(stderr,
                     "rank %d: solver status is inconsistent (converged=%d, active=%d, max_residual=%.12e)\n",
                     myrank,
                     converged,
                     active_band_count,
                     max_residual);
        ok = 0;
    }
    for (int i = 0; i < nband; ++i)
    {
        if (std::abs(eval[i] - Real(i + 1)) > 1e-6)
        {
            std::printf("rank %d: eval[%d] = %.12f != %d\n", myrank, i, eval[i], i + 1);
            ok = 0;
        }
    }

    if (myrank == 0)
    {
        std::ifstream trace(trace_path.c_str());
        std::string line;
        int header_count = 0;
        int record_count = 0;
        bool has_final_record = false;
        while (std::getline(trace, line))
        {
            if (line.find("iteration,stage,max_residual") == 0)
            {
                ++header_count;
            }
            else if (!line.empty())
            {
                ++record_count;
                if (line.find(",final,") != std::string::npos)
                {
                    has_final_record = true;
                }
            }
        }
        if (!trace.eof() || header_count != 1 || record_count == 0 || !has_final_record)
        {
            const int final_record = has_final_record ? 1 : 0;
            std::fprintf(stderr,
                         "rank 0: residual trace is incomplete (headers=%d, records=%d, final=%d)\n",
                         header_count,
                         record_count,
                         final_record);
            ok = 0;
        }
        std::remove(trace_path.c_str());
    }

    // Complete the single-trace assertions before this additional solve writes
    // its own trace. Each pool must independently recover both lowest roots.
    const int trace_unset_status = ::unsetenv(trace_environment);
    const bool competing_roots_ok = check_competing_roots(rank_in_pool, nproc_in_pool);
    if (trace_unset_status != 0 || !competing_roots_ok)
    {
        ok = 0;
    }

    int global_ok = ok;
    Parallel_Reduce::reduce_min(global_ok);

    if (myrank == 0)
    {
        if (global_ok == 1)
        {
            std::printf("PPCG MPI parallel test PASSED\n");
        }
        else
        {
            std::printf("PPCG MPI parallel test FAILED\n");
        }
    }
    const int exit_code = global_ok == 1 ? 0 : 1;
    finishmpi();
    return exit_code;
}
