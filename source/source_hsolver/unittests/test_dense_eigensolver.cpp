#include "source_hsolver/dense_eigensolver.h"
#include "source_base/parallel_cell.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

namespace
{
using Z = std::complex<double>;

// Column-major H/S with the exact generalized spectrum 1,...,n.
struct Pencil
{
    int n;
    std::vector<Z> h;
    std::vector<Z> s;
    std::vector<double> exact;
    std::vector<double> values;
    std::vector<Z> vectors;

    explicit Pencil(const int dimension) : n(dimension)
    {
        h.resize(n * n);
        s.resize(n * n);
        exact.resize(n);
        values.resize(n);
        vectors.resize(n * n);
        for (int i = 0; i < n; ++i)
        {
            s[i + i * n] = 2.0;
            h[i + i * n] = 2.0 * (i + 1);
            exact[i] = i + 1;
        }
    }
};

// Check against the preserved input matrices using explicit host sums.
void check_pairs(const Pencil& input, const Pencil& output, const int wanted)
{
    EXPECT_EQ(output.h, input.h);
    EXPECT_EQ(output.s, input.s);
    const int n = input.n;
    std::vector<Z> sx(n * wanted);
    for (int j = 0; j < wanted; ++j)
    {
        EXPECT_NEAR(output.values[j], input.exact[j], 2e-12);
        for (int i = 0; i < n; ++i)
        {
            Z hx(0);
            for (int k = 0; k < n; ++k)
            {
                hx += input.h[i + k * n] * output.vectors[k + j * n];
                sx[i + j * n] += input.s[i + k * n] * output.vectors[k + j * n];
            }
            EXPECT_LT(std::abs(hx - output.values[j] * sx[i + j * n]), 2e-10);
        }
        for (int i = 0; i < wanted; ++i)
        {
            Z overlap(0);
            for (int k = 0; k < n; ++k)
            {
                overlap += std::conj(output.vectors[k + i * n]) * sx[k + j * n];
            }
            const double identity = i == j ? 1 : 0;
            EXPECT_LT(std::abs(overlap - identity), 2e-10);
        }
    }
}

// H=F diag(2*lambda) F^H, S=2I, with unitary Fourier F.
Pencil clustered_pencil()
{
    Pencil p(16);
    std::fill(p.h.begin(), p.h.end(), Z(0));
    for (int k = 0; k < p.n; ++k)
    {
        p.exact[k] = k / 4;
        for (int j = 0; j < p.n; ++j)
        {
            for (int i = 0; i < p.n; ++i)
            {
                const double phase = 2 * std::acos(-1.0) * (i - j) * k / p.n;
                p.h[i + j * p.n] += 2 * p.exact[k] * Z(std::cos(phase), std::sin(phase)) / double(p.n);
            }
        }
    }
    return p;
}

#if !defined(__CUDA) && !defined(__ROCM)
void solve_serial(hsolver::dense::SerialGvWorkspace& solver, Pencil& p, const int wanted)
{
    const Pencil input(p);
    ASSERT_EQ(solver.solve(p.h.data(), p.s.data(), p.n, wanted, p.values.data(), p.vectors.data()), 0);
    check_pairs(input, p, wanted);
}

TEST(SerialDenseEigensolver, PartialAndFullSpectrum)
{
    hsolver::dense::SerialGvWorkspace solver;
    Pencil p = clustered_pencil();
    solve_serial(solver, p, 6);
    solve_serial(solver, p, p.n);
}

TEST(SerialDenseEigensolver, ReusesWorkspace)
{
    hsolver::dense::SerialGvWorkspace solver;
    Pencil p(16);
    solve_serial(solver, p, 6);
    const std::size_t capacity = solver.storage_bytes();
    solve_serial(solver, p, 6);
    EXPECT_EQ(solver.storage_bytes(), capacity);
    Pencil smaller(4);
    solve_serial(solver, smaller, 2);
    EXPECT_EQ(solver.storage_bytes(), capacity);
}

TEST(SerialDenseEigensolver, InvalidMetricAndRetry)
{
    hsolver::dense::SerialGvWorkspace solver;
    Pencil p(4);
    p.s[0] = -1;
    EXPECT_GT(solver.solve(p.h.data(), p.s.data(), p.n, 2, p.values.data(), p.vectors.data()), p.n);
    EXPECT_EQ(p.s[0], Z(-1));
    p.s[0] = 2;
    solve_serial(solver, p, 2);
}

#ifdef __MPI
void solve_parallel(hsolver::dense::ParallelDenseEigensolver& solver, Pencil& p, const int wanted,
                    const MPI_Comm comm, const int block)
{
    const Pencil input(p);
    solver.solve_generalized(p.h.data(), p.s.data(), p.n, wanted,
                             p.values.data(), p.vectors.data(), comm, block);
    ModuleBase::CommunicationDomain domain;
    domain.initialize(comm);
    if (domain.rank() == 0) { check_pairs(input, p, wanted); }
}

TEST(ParallelDenseEigensolver, PartialAndFullClusteredSpectrum)
{
    hsolver::dense::ParallelDenseEigensolver solver;
    Pencil p = clustered_pencil();
    solve_parallel(solver, p, 6, MPI_COMM_WORLD, 4); // Cuts a fourfold cluster.
    solve_parallel(solver, p, p.n, MPI_COMM_WORLD, 4);
}

TEST(ParallelDenseEigensolver, FullyDegenerateGeneralizedSpectrum)
{
    Pencil p(16);
    // A positive, nondiagonal S; H=2.5S gives a fully degenerate spectrum.
    for (int i = 1; i < p.n; ++i)
    {
        p.s[i - 1 + i * p.n] = Z(.2, .1);
        p.s[i + (i - 1) * p.n] = Z(.2, -.1);
    }
    for (int i = 0; i < p.n * p.n; ++i) { p.h[i] = 2.5 * p.s[i]; }
    p.exact.assign(p.n, 2.5);
    hsolver::dense::ParallelDenseEigensolver solver;
    solve_parallel(solver, p, 6, MPI_COMM_WORLD, 4);
    solve_parallel(solver, p, p.n, MPI_COMM_WORLD, 4);
}

TEST(ParallelDenseEigensolver, WorkspaceGrowthAndReuse)
{
    hsolver::dense::ParallelDenseEigensolver solver;
    Pencil p(16);
    solve_parallel(solver, p, 6, MPI_COMM_WORLD, 4);
    const std::size_t capacity = solver.storage_bytes();
    solve_parallel(solver, p, 6, MPI_COMM_WORLD, 4);
    EXPECT_EQ(solver.storage_bytes(), capacity);
    Pencil larger(32);
    solve_parallel(solver, larger, 12, MPI_COMM_WORLD, 4);
    EXPECT_GE(solver.storage_bytes(), capacity);
    Pencil tiny(1);
    solve_parallel(solver, tiny, 1, MPI_COMM_WORLD, 4);
}

TEST(ParallelDenseEigensolver, RebindsCommunicator)
{
    hsolver::dense::ParallelDenseEigensolver solver;
    Pencil p(16);
    solve_parallel(solver, p, 6, MPI_COMM_WORLD, 4);
    solve_parallel(solver, p, 6, MPI_COMM_SELF, 4);
    solve_parallel(solver, p, 6, MPI_COMM_WORLD, 4);
}

TEST(ParallelDenseEigensolver, InvalidMetricAndRetry)
{
    hsolver::dense::ParallelDenseEigensolver solver;
    Pencil p(16);
    p.s[0] = -1;
    EXPECT_THROW(solver.solve_generalized(p.h.data(), p.s.data(), p.n, 6,
                                         p.values.data(), p.vectors.data(), MPI_COMM_WORLD, 4), std::runtime_error);
    EXPECT_EQ(p.s[0], Z(-1));
    p.s[0] = 2;
    solve_parallel(solver, p, 6, MPI_COMM_WORLD, 4);
}
#endif
#endif

#ifdef __CUDA
void solve_gpu(hsolver::dense::GpuDenseEigensolver& solver, Pencil& p, const int wanted)
{
    const Pencil input(p);
    ASSERT_TRUE(solver.solve_generalized(p.h.data(), p.s.data(), p.n, wanted,
                                         p.values.data(), p.vectors.data()));
    check_pairs(input, p, wanted);
}

TEST(GpuDenseEigensolver, PartialAndFullClusteredSpectrum)
{
    hsolver::dense::GpuDenseEigensolver solver;
    Pencil p = clustered_pencil();
    solve_gpu(solver, p, 6);
    solve_gpu(solver, p, p.n);
}

TEST(GpuDenseEigensolver, WorkspaceGrowthAndReuse)
{
    hsolver::dense::GpuDenseEigensolver solver;
    Pencil p(16);
    solve_gpu(solver, p, 6);
    const std::size_t capacity = solver.storage_bytes();
    solve_gpu(solver, p, 6);
    EXPECT_EQ(solver.storage_bytes(), capacity);
    solve_gpu(solver, p, 8);
    Pencil larger(32);
    solve_gpu(solver, larger, 12);
    EXPECT_GE(solver.storage_bytes(), capacity);
    const std::size_t grown = solver.storage_bytes();
    solve_gpu(solver, larger, 12);
    EXPECT_EQ(solver.storage_bytes(), grown);
}

TEST(GpuDenseEigensolver, SingleEigenpair)
{
    hsolver::dense::GpuDenseEigensolver solver;
    Pencil p(1);
    solve_gpu(solver, p, 1);
}

TEST(GpuDenseEigensolver, InvalidMetricAndRetry)
{
    hsolver::dense::GpuDenseEigensolver solver;
    Pencil p(16);
    solve_gpu(solver, p, 6);
    p.s[0] = -1;
    EXPECT_FALSE(solver.solve_generalized(p.h.data(), p.s.data(), p.n, 6,
                                         p.values.data(), p.vectors.data()));
    EXPECT_EQ(p.s[0], Z(-1));
    p.s[0] = 2;
    solve_gpu(solver, p, 6);
}
#endif
} // namespace
