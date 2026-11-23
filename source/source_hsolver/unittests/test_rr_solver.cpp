#include "source_hsolver/rr_solver.h"
#ifdef __MPI
#include "source_base/parallel_blkcyc.h"
#endif
#include <gtest/gtest.h>
#include <cmath>
#include <complex>
#include <random>
#include <stdexcept>
#include <vector>

namespace
{
using Z = std::complex<double>;

struct Pencil
{
    int n;
    std::vector<Z> h;
    std::vector<Z> s;
    std::vector<double> exact;

    explicit Pencil(const int dimension) : n(dimension)
    {
        h.resize(n * n);
        s.resize(n * n);
        for (int i = 0; i < n; ++i)
        {
            s[i + i * n] = 2.0;
        }
    }
};

// Verify eigenpairs against the input H/S matrices using explicit host sums.
template <typename Solve>
void check_pairs(Pencil& p, const int wanted, const double eigen_error,
                 const double vector_error, const bool owns_output, Solve solve)
{
    SCOPED_TRACE(::testing::Message() << "n=" << p.n << " wanted=" << wanted);
    const auto h = p.h;
    const auto s = p.s;
    const Z guard(713, 29);
    std::vector<Z> vectors(p.n * wanted + 2, guard);
    std::vector<double> values(p.n + 2, 713);
    Z* x = vectors.data() + 1;
    double* w = values.data() + 1;
    solve(w, x);
    EXPECT_EQ(p.h, h);
    EXPECT_EQ(p.s, s);
    EXPECT_EQ(vectors.front(), guard);
    EXPECT_EQ(vectors.back(), guard);
    EXPECT_EQ(values.front(), 713);
    EXPECT_EQ(values.back(), 713);
    if (!owns_output) { return; }
    std::vector<Z> sx(p.n * wanted);
    for (int j = 0; j < wanted; ++j)
    {
        if (!p.exact.empty()) { EXPECT_NEAR(w[j], p.exact[j], eigen_error); }
        for (int i = 0; i < p.n; ++i)
        {
            Z hx(0);
            for (int k = 0; k < p.n; ++k)
            {
                hx += h[i + k * p.n] * x[k + j * p.n];
                sx[i + j * p.n] += s[i + k * p.n] * x[k + j * p.n];
            }
            EXPECT_LT(std::abs(hx - w[j] * sx[i + j * p.n]), vector_error);
        }
        for (int i = 0; i < wanted; ++i)
        {
            Z overlap(0);
            for (int k = 0; k < p.n; ++k)
            {
                overlap += std::conj(x[k + i * p.n]) * sx[k + j * p.n];
            }
            const double identity = i == j ? 1 : 0;
            EXPECT_LT(std::abs(overlap - identity), vector_error);
        }
    }
}

Pencil diagonal_pencil(const int n)
{
    Pencil p(n);
    for (int i = 0; i < n; ++i)
    {
        p.h[i + i * n] = 2.0 * (i + 1);
        p.exact.push_back(i + 1);
    }
    return p;
}

#if defined(__MPI) || defined(__CUDA)
Pencil fourier_pencil(const std::vector<double>& spectrum)
{
    const int n = spectrum.size();
    Pencil p(n);
    const double pi = std::acos(-1.0);
    for (int k = 0; k < n; ++k)
    {
        p.exact.push_back(spectrum[k] / 2);
        for (int j = 0; j < n; ++j)
        {
            for (int i = 0; i < n; ++i)
            {
                const double phase = 2 * pi * (i - j) * k / n;
                p.h[i + j * n] += spectrum[k] * Z(std::cos(phase), std::sin(phase)) / double(n);
            }
        }
    }
    return p;
}
#endif

#ifdef __CUDA
void check_gpu(hsolver::rr::GpuRREigensolver& solver, const int n, const int wanted)
{
    std::vector<double> spectrum(n);
    for (int k = 0; k < n; ++k) { spectrum[k] = k < 3 ? -2.0 : 0.2 * k; }
    auto p = fourier_pencil(spectrum);
    check_pairs(p, wanted, 1e-10, 1e-10, true, [&](double* w, Z* x) {
        ASSERT_TRUE(solver.solve_generalized(p.h.data(), p.s.data(), n, wanted, w, x));
        EXPECT_EQ(w[wanted], 713); // GPU copies only the requested eigenvalues.
    });
}
#endif

#if defined(__MPI) && !defined(__CUDA) && !defined(__ROCM)
Pencil random_pencil(const int n, const int seed)
{
    Pencil p(n);
    auto& h = p.h;
    auto& s = p.s;
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> uniform(-1.0, 1.0);
    std::vector<Z> factor(n * n);
    h.assign(n * n, Z(0));
    s.assign(n * n, Z(0));
    for (int i = 0; i < n; ++i)
    {
        for (int j = i; j < n; ++j)
        {
            const double real = uniform(generator);
            const double imag = uniform(generator);
            h[i * n + j] = Z(real, imag);
            h[j * n + i] = i == j ? Z(real, 0) : Z(real, -imag);
            const double sr = uniform(generator);
            const double si = uniform(generator);
            factor[i * n + j] = Z(sr, si);
        }
    }
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            for (int k = 0; k < n; ++k)
            {
                s[i * n + j] += factor[i * n + k] * std::conj(factor[j * n + k]);
            }
            if (i == j)
            {
                s[i * n + j] += 2.0;
            }
        }
    }
    return p;
}

Pencil tridiagonal_pencil(const int n, const double diagonal, const double width)
{
    Pencil p(n);
    for (int j = 0; j < n; ++j)
    {
        p.h[j + j * n] = diagonal;
        if (j > 0)
        {
            p.h[j - 1 + j * n] = width * Z(0.3, 0.4);
            p.h[j + (j - 1) * n] = width * Z(0.3, -0.4);
        }
        const double angle = std::acos(-1.0) * (j + 1) / (n + 1);
        p.exact.push_back((diagonal - width * std::cos(angle)) / 2);
    }
    return p;
}

void check_parallel(hsolver::rr::ParallelRREigensolver& solver, Pencil& p, const int wanted,
                    const MPI_Comm comm, const int block, const double eigen_error, const double vector_error)
{
    Parallel_Blkcyclic ownership;
    ownership.init(p.n, p.n, p.n, comm, false);
    const bool root = ownership.blacs_in_this_processor(0, 0);
    check_pairs(p, wanted, eigen_error, vector_error, root, [&](double* w, Z* x) {
        solver.solve_generalized(p.h.data(), p.s.data(), p.n, wanted, w, x, comm, block);
    });
}
#endif
} // namespace

#ifdef __CUDA
TEST(GpuRREigensolver, ClusterPartialFullAndShrinkingWorkspace)
{
    ASSERT_TRUE(hsolver::rr::GpuRREigensolver::available());
    hsolver::rr::GpuRREigensolver solver;
    EXPECT_EQ(solver.storage_bytes(), 0);
    check_gpu(solver, 48, 7);
    const std::size_t partial = solver.storage_bytes();
    check_gpu(solver, 48, 7);
    EXPECT_EQ(solver.storage_bytes(), partial);
    check_gpu(solver, 48, 48);
    const std::size_t warm = solver.storage_bytes();
    check_gpu(solver, 12, 12);
    EXPECT_EQ(solver.storage_bytes(), warm);
    check_gpu(solver, 1, 1);
    check_gpu(solver, 48, 2); // The requested boundary cuts a degenerate cluster.
    EXPECT_EQ(solver.storage_bytes(), warm);
    check_gpu(solver, 96, 9);
    const std::size_t grown = solver.storage_bytes();
    EXPECT_GT(grown, warm);
    check_gpu(solver, 48, 7);
    EXPECT_EQ(solver.storage_bytes(), grown);
}

TEST(GpuRREigensolver, QueryInvalidationAndFailureReuse)
{
    hsolver::rr::GpuRREigensolver solver;
    check_gpu(solver, 16, 5);
    check_gpu(solver, 16, 7);
    auto p = diagonal_pencil(4);
    std::vector<Z> vectors(16);
    std::vector<double> values(4);
    p.s[0] = -1.0;
    EXPECT_FALSE(solver.solve_generalized(p.h.data(), p.s.data(), 4, 4, values.data(), vectors.data()));
    const std::size_t capacity = solver.storage_bytes();
    check_gpu(solver, 16, 5);
    EXPECT_EQ(solver.storage_bytes(), capacity);
}
#endif

#if !defined(__CUDA) && !defined(__ROCM)
TEST(SerialRR, PartialSpectrumInputPreservationAndFailureReuse)
{
    hsolver::rr::SerialGvWorkspace solver;
    for (const int n : {12, 12, 6, 24})
    {
        auto p = diagonal_pencil(n);
        const int wanted = n / 2;
        check_pairs(p, wanted, 1e-12, 1e-12, true, [&](double* w, Z* x) {
            ASSERT_EQ(solver.solve(p.h.data(), p.s.data(), n, wanted, w, x), 0);
        });
    }
    auto p = diagonal_pencil(12);
    check_pairs(p, 6, 1e-12, 1e-12, true, [&](double* w, Z* x) {
        p.s[0] = -1.0;
        EXPECT_GT(solver.solve(p.h.data(), p.s.data(), 12, 6, w, x), 12);
        p.s[0] = 2.0;
        ASSERT_EQ(solver.solve(p.h.data(), p.s.data(), 12, 6, w, x), 0);
    });
}

#ifdef __MPI
TEST(ParallelRR, ReusesWorkspaceRebindsAndRestoresInputs)
{
    hsolver::rr::ParallelRREigensolver solver;
    std::size_t previous_capacity = 0;
    for (int step = 0; step < 4; ++step)
    {
        const int n = 32;
        const int seed = 100 + step;
        auto p = random_pencil(n, seed);
        // Rebind at the same dimension, then restore the original communicator.
        const MPI_Comm comm = step == 2 ? MPI_COMM_SELF : MPI_COMM_WORLD;
        check_parallel(solver, p, n, comm, 4, 1e-10, 1e-10);
        if (step == 1) { EXPECT_EQ(solver.storage_bytes(), previous_capacity); }
        previous_capacity = solver.storage_bytes();
    }
}

TEST(ParallelRR, GeneralizedSpectrumAndChangingLayouts)
{
    hsolver::rr::ParallelRREigensolver solver;
    const int dimensions[] = {1, 17, 32, 32, 33};
    const int blocks[] = {4, 16, 16, 16, 32};
    std::size_t previous_capacity = 0;
    for (int step = 0; step < 5; ++step)
    {
        const int n = dimensions[step];
        const double diagonal = 2.0 + step;
        const double width = step == 4 ? 1e-8 : 1.0;
        auto p = tridiagonal_pencil(n, diagonal, width);
        check_parallel(solver, p, n, MPI_COMM_WORLD, blocks[step], 1e-12, 1e-11);
        if (step == 3) { EXPECT_EQ(solver.storage_bytes(), previous_capacity); }
        previous_capacity = solver.storage_bytes();
    }
}

TEST(ParallelRR, RequestedPairsAndWorkspaceChanges)
{
    hsolver::rr::ParallelRREigensolver solver;
    const int dimensions[] = {17, 33, 64, 64, 64, 8, 64};
    const int requested[] = {6, 11, 21, 21, 30, 3, 21};
    for (int step = 0; step < 7; ++step)
    {
        auto p = tridiagonal_pencil(dimensions[step], 2.0, 1.0);
        const std::size_t capacity = solver.storage_bytes();
        check_parallel(solver, p, requested[step], MPI_COMM_WORLD, 16, 1e-12, 1e-11);
        if (step == 3) { EXPECT_EQ(solver.storage_bytes(), capacity); }
    }
}

TEST(ParallelRR, RequestedBoundarySplitsDegenerateCluster)
{
    std::vector<double> spectrum(33);
    for (int k = 0; k < 33; ++k)
    {
        // wanted=8 cuts the fivefold cluster k=5,...,9.
        spectrum[k] = k / 5;
        if (k >= 10) { spectrum[k] += 1e-12 * k; }
    }
    auto p = fourier_pencil(spectrum);
    hsolver::rr::ParallelRREigensolver solver;
    check_parallel(solver, p, 8, MPI_COMM_WORLD, 16, 2e-12, 2e-11);
}

TEST(ParallelRR, FullDegenerateClusterWithPartialAndFullSpectrum)
{
    auto p = random_pencil(33, 42);
    for (int i = 0; i < p.n * p.n; ++i) { p.h[i] = 2.5 * p.s[i]; }
    p.exact.assign(p.n, 2.5); // All N eigenvalues belong to the same cluster.
    hsolver::rr::ParallelRREigensolver solver;
    for (const int wanted : {8, p.n})
    {
        check_parallel(solver, p, wanted, MPI_COMM_WORLD, 16, 2e-12, 2e-10);
    }
}

TEST(ParallelRR, NonPositiveMetricFailurePreservesInputsAndAllowsReuse)
{
    auto p = diagonal_pencil(16);
    for (int i = 0; i < p.n; ++i)
    {
        p.h[i + i * p.n] /= 2;
        p.s[i + i * p.n] = 1;
    }
    p.s[0] = -1;
    const auto h = p.h;
    const auto s = p.s;
    std::vector<Z> x(p.n * 5);
    std::vector<double> w(p.n);
    hsolver::rr::ParallelRREigensolver solver;
    EXPECT_THROW(solver.solve_generalized(p.h.data(), p.s.data(), p.n, 5, w.data(), x.data(),
                                         MPI_COMM_WORLD, 4), std::runtime_error);
    EXPECT_EQ(p.h, h);
    EXPECT_EQ(p.s, s);
    p.s[0] = 1;
    check_parallel(solver, p, 5, MPI_COMM_WORLD, 4, 1e-12, 1e-12);
}
#endif
#endif
