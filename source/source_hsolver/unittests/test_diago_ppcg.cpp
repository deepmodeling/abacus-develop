// Analytical and generalized matrix checks for DiagoPPCG, plus the shared test main.
// Numerical regressions are grouped under ppcg/{roots,iteration,subspace}; CMake
// links all four sources into MODULE_HSOLVER_ppcg so CTest runs the complete set.

#include "diago_ppcg_test_support.h"
#ifdef __MPI
#include "source_basis/module_pw/test/test_tool.h"
#endif

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <ostream>
#include <random>
#include <vector>
#include <complex>

using ppcg_test::DiagoPPCGMatrixTest;
using ppcg_test::PPCG_PI;
using ppcg_test::Real;
using ppcg_test::T;
using ppcg_test::dense_h_multiply;
using ppcg_test::initialize_vectors;

// Like the CG test's preparation table, each row describes one complete solve.
struct PPCGLaplacianCase
{
    const char* name;
    int dimension;
    int bands;
    unsigned seed;
    double band_threshold;
    double matrix_scale;
    double diag_threshold;
    int max_iterations;
    int block_size;
    int ritz_interval;
    double eigenvalue_tolerance;
};

static void PrintTo(const PPCGLaplacianCase& config, std::ostream* out)
{
    *out << config.name;
}

class DiagoPPCGLaplacianTest : public DiagoPPCGMatrixTest,
                               public ::testing::WithParamInterface<PPCGLaplacianCase>
{
};

TEST_P(DiagoPPCGLaplacianTest, MatchesAnalyticalEigenvalues)
{
    const PPCGLaplacianCase& config = GetParam();
    initialize(config.dimension, config.bands, config.seed, config.band_threshold);
    set_laplacian(config.matrix_scale);
    check_eigenvalues(config.diag_threshold, config.max_iterations, config.block_size,
                      config.ritz_interval, config.eigenvalue_tolerance);
}

static const PPCGLaplacianCase laplacian_cases[] = {
    // name, dimension, bands, seed, band threshold, scale, diag threshold, budget, block, Ritz interval, tolerance
    {"Tridiagonal10", 10, 3, 42, 1e-10, 1.0, 1e-12, 100, 4, 4, 1e-8},
    {"SingleBand", 5, 1, 42, 1e-10, 1.0, 1e-12, 50, 1, 1, 1e-8},
    {"AllBands", 3, 3, 101, 1e-10, 1.0, 1e-12, 100, 3, 3, 1e-8},
    {"Laplacian15", 15, 4, 202, 1e-10, 1.0, 1e-12, 120, 4, 4, 1e-8},
    {"Laplacian20", 20, 5, 789, 1e-10, 1.0, 1e-12, 150, 5, 5, 1e-8},
    {"ScaledLaplacian", 8, 3, 808, 1e-8, 100.0, 1e-10, 120, 4, 4, 1e-6},
    {"ManyBands", 12, 4, 909, 1e-10, 1.0, 1e-12, 150, 4, 4, 1e-8},
    {"RitzEveryIteration", 8, 3, 111, 1e-10, 1.0, 1e-12, 100, 3, 1, 1e-8},
    {"TightThreshold", 6, 2, 333, 1e-14, 1.0, 1e-14, 200, 3, 3, 1e-8}};

INSTANTIATE_TEST_SUITE_P(
    Laplacians,
    DiagoPPCGLaplacianTest,
    ::testing::ValuesIn(laplacian_cases),
    [](const ::testing::TestParamInfo<PPCGLaplacianCase>& info) { return info.param.name; });

// Test fixture: 2×2 matrix — smallest non-trivial case
// H = [[2, 1], [1, 2]], eigenvalues: 1, 3
TEST_F(DiagoPPCGMatrixTest, TwoByTwo)
{
    initialize(2, 2, 123, 1e-10);
    // H = [[2, 1], [1, 2]]
    H_mat.assign(n_dim * n_dim, T(0));
    H_mat[0 + 0 * n_dim] = T(2.0, 0);
    H_mat[1 + 1 * n_dim] = T(2.0, 0);
    H_mat[0 + 1 * n_dim] = T(1.0, 0);
    H_mat[1 + 0 * n_dim] = T(1.0, 0);

    prec.assign(n_dim, 2.0);

    // λ₁ = 1, λ₂ = 3
    exact = {1.0, 3.0};
    check_eigenvalues(1e-12, 50, 2, 2, 1e-8);
}

TEST_F(DiagoPPCGMatrixTest, DegenerateEigenvalues)
{
    initialize(4, 4, 456, 1e-10);
    // H = I + J where J is the all-ones matrix
    H_mat.assign(n_dim * n_dim, T(0));
    for (int i = 0; i < n_dim; ++i)
    {
        for (int j = 0; j < n_dim; ++j)
        {
            H_mat[i + j * n_dim] = T(1.0, 0); // all-ones J
        }
    }
    for (int i = 0; i < n_dim; ++i)
    {
        H_mat[i + i * n_dim] += T(1.0, 0); // J → I+J
    }

    // Preconditioner: diagonal = 2
    prec.assign(n_dim, 2.0);

    // λ = {1, 1, 1, 5}
    exact = {1.0, 1.0, 1.0, 5.0};
    check_eigenvalues(1e-12, 100, 4, 4, 1e-8);
}

// Test fixture: dense matrix with known eigenvalues
// H = Q^T * D * Q where Q is a known orthogonal matrix (a Givens rotation
// repeated on different index pairs) and D is diagonal.
// For an 8×8 case: D = diag(1, 2, 3, 4, 5, 6, 7, 8), then apply several
// Givens rotations to mix all rows/cols. The exact eigenvalues remain 1..8.
TEST_F(DiagoPPCGMatrixTest, DenseKnownSpectrum)
{
    initialize(8, 4, 222, 1e-10);
    // Start with diagonal matrix
    std::vector<Real> dense(n_dim * n_dim, Real(0));
    for (int i = 0; i < n_dim; ++i)
    {
        dense[i + i * n_dim] = Real(i + 1);
    }

    // Apply several Givens rotations to make it dense while preserving
    // eigenvalues. Each rotation: A' = G(i,j,θ)^T * A * G(i,j,θ)
    auto apply_givens = [&](int p, int q, Real theta) {
        Real c = std::cos(theta);
        Real s = std::sin(theta);
        // Apply to columns
        for (int i = 0; i < n_dim; ++i)
        {
            Real aip = dense[i + p * n_dim];
            Real aiq = dense[i + q * n_dim];
            dense[i + p * n_dim] = c * aip + s * aiq;
            dense[i + q * n_dim] = -s * aip + c * aiq;
        }
        // Apply to rows
        for (int j = 0; j < n_dim; ++j)
        {
            Real apj = dense[p + j * n_dim];
            Real aqj = dense[q + j * n_dim];
            dense[p + j * n_dim] = c * apj + s * aqj;
            dense[q + j * n_dim] = -s * apj + c * aqj;
        }
    };

    // Several rotations with different angles to create a genuinely
    // dense matrix (all off-diagonals become non-zero)
    std::mt19937 rng_dense(111);
    std::uniform_real_distribution<Real> angle_dist(Real(0.2), Real(1.3));
    for (int k = 0; k < 20; ++k)
    {
        int p = k % (n_dim - 1);
        int q = p + 1 + (k / (n_dim - 1)) % (n_dim - 1 - p);
        if (q >= n_dim)
        {
            q = n_dim - 1;
        }
        if (p == q)
        {
            continue;
        }
        const Real angle = angle_dist(rng_dense);
        apply_givens(p, q, angle);
    }

    // Copy to complex H_mat
    H_mat.assign(n_dim * n_dim, T(0));
    for (int i = 0; i < n_dim * n_dim; ++i)
    {
        H_mat[i] = T(dense[i], 0);
    }

    // Preconditioner: use diagonal of the rotated H
    prec.resize(n_dim);
    for (int i = 0; i < n_dim; ++i)
    {
        prec[i] = std::real(H_mat[i + i * n_dim]);
    }

    // Lowest 4 eigenvalues: 1, 2, 3, 4
    exact = {1.0, 2.0, 3.0, 4.0};
    check_eigenvalues(1e-12, 200, 4, 4, 1e-8);
}

// Euclidean column norm for independent residual checks
static Real column_norm(const T* x, int n_dim)
{
    Real nrm = 0;
    for (int i = 0; i < n_dim; ++i)
    {
        nrm += std::norm(x[i]);
    }
    return std::sqrt(nrm);
}

struct PPCGOverlapCase
{
    const char* name;
    unsigned seed;
    int bands;
    int leading_dimension;
    int max_iterations;
    bool tridiagonal;
    double residual_threshold;
};

static void PrintTo(const PPCGOverlapCase& config, std::ostream* out)
{
    *out << config.name;
}

// H is the six-dimensional Laplacian. Diagonal S also exercises padded vectors;
// tridiagonal S shares H's sine eigenvectors and has analytical eigenvalues.
class DiagoPPCGGeneralizedTest : public DiagoPPCGMatrixTest,
                                 public ::testing::WithParamInterface<PPCGOverlapCase>
{
  protected:
    void SetUp() override
    {
        const PPCGOverlapCase& config = GetParam();
        n_dim = 6;
        nband = config.bands;
        ld = config.leading_dimension;
        ethr.assign(nband, 1e-8);
        set_laplacian(1.0);
        const int matrix_size = n_dim * n_dim;
        S_mat.assign(matrix_size, T(0));
        const Real diagonal_overlap[] = {1.1, 1.0, 0.9, 1.0, 1.1, 1.0};
        for (int i = 0; i < n_dim; ++i)
        {
            const Real diagonal = config.tridiagonal ? 1.0 : diagonal_overlap[i];
            S_mat[i + i * n_dim] = diagonal;
            if (config.tridiagonal && i > 0)
            {
                S_mat[i + (i - 1) * n_dim] = 0.2;
                S_mat[i - 1 + i * n_dim] = 0.2;
            }
        }
        if (config.tridiagonal)
        {
            for (int k = 0; k < nband; ++k)
            {
                const Real theta = Real(k + 1) * PPCG_PI / Real(n_dim + 1);
                const Real cosine = std::cos(theta);
                exact[k] = (2.0 - 2.0 * cosine) / (1.0 + 0.4 * cosine);
            }
        }

        std::mt19937 rng(config.seed);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);
        const int vector_size = ld * nband;
        psi.assign(vector_size, T(0));
        for (int j = 0; j < nband; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                const Real value = dist(rng);
                psi[i + j * ld] = value;
            }
        }
        std::vector<T> s_column(n_dim);
        for (int j = 0; j < nband; ++j)
        {
            T* column = psi.data() + j * ld;
            for (int k = 0; k < j; ++k)
            {
                const T* previous = psi.data() + k * ld;
                dense_h_multiply(S_mat.data(), n_dim, previous, s_column.data(), n_dim, 1);
                T overlap = 0;
                for (int i = 0; i < n_dim; ++i)
                {
                    overlap += std::conj(s_column[i]) * column[i];
                }
                for (int i = 0; i < n_dim; ++i)
                {
                    column[i] -= overlap * previous[i];
                }
            }
            dense_h_multiply(S_mat.data(), n_dim, column, s_column.data(), n_dim, 1);
            Real norm_squared = 0;
            for (int i = 0; i < n_dim; ++i)
            {
                const T norm_term = std::conj(column[i]) * s_column[i];
                norm_squared += std::real(norm_term);
            }
            ASSERT_GT(norm_squared, 0.0);
            const Real norm = std::sqrt(norm_squared);
            for (int i = 0; i < n_dim; ++i)
            {
                column[i] /= norm;
            }
        }
    }

    void apply_overlap(T* in, T* out, int leading_dimension, int columns) const
    {
        const bool tridiagonal = GetParam().tridiagonal;
        for (int j = 0; j < columns; ++j)
        {
            for (int i = 0; i < n_dim; ++i)
            {
                const int offset = i + j * leading_dimension;
                T value = S_mat[i + i * n_dim] * in[offset];
                if (tridiagonal && i > 0)
                {
                    value += 0.2 * in[offset - 1];
                }
                if (tridiagonal && i + 1 < n_dim)
                {
                    value += 0.2 * in[offset + 1];
                }
                out[offset] = value;
            }
        }
    }

    std::vector<T> S_mat;
};

TEST_P(DiagoPPCGGeneralizedTest, EigenpairsAndSOrthogonality)
{
    const PPCGOverlapCase& config = GetParam();
    std::vector<Real> eval(nband, 0.0);
    hsolver::DiagoPPCG<T, base_device::DEVICE_CPU> solver(1e-10, config.max_iterations, 3, 3);
    const auto h_op = [this](T* in, T* out, int stride, int columns) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, stride, columns);
    };
    const auto s_op = [this](T* in, T* out, int stride, int columns) {
        apply_overlap(in, out, stride, columns);
    };
    const double iterations = solver.diag(
        h_op, s_op, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());
    if (config.tridiagonal)
    {
        EXPECT_TRUE(solver.converged());
        EXPECT_EQ(solver.active_band_count(), 0);
        EXPECT_FALSE(solver.iteration_limit_reached());
    }

    // Use fresh dense products, independently of the structured S callback.
    const int image_size = ld * nband;
    std::vector<T> hpsi(image_size);
    std::vector<T> spsi(image_size);
    std::vector<T> residual(n_dim);
    dense_h_multiply(H_mat.data(), n_dim, psi.data(), hpsi.data(), ld, nband);
    dense_h_multiply(S_mat.data(), n_dim, psi.data(), spsi.data(), ld, nband);
    for (int i = 0; i < nband; ++i)
    {
        if (config.tridiagonal)
        {
            EXPECT_NEAR(eval[i], exact[i], ethr[i]) << "band " << i;
        }
        else
        {
            EXPECT_GT(eval[i], 0.0) << "band " << i;
            EXPECT_LT(eval[i], 10.0) << "band " << i;
        }
        for (int j = 0; j < n_dim; ++j)
        {
            const int offset = j + i * ld;
            residual[j] = hpsi[offset] - eval[i] * spsi[offset];
        }
        const Real residual_norm = column_norm(residual.data(), n_dim);
        EXPECT_LE(residual_norm, config.residual_threshold) << "band " << i;
        for (int k = 0; k < nband; ++k)
        {
            T overlap = 0;
            for (int j = 0; j < n_dim; ++j)
            {
                overlap += std::conj(psi[j + k * ld]) * spsi[j + i * ld];
            }
            const T expected_overlap = (i == k) ? T(1) : T(0);
            const T difference = overlap - expected_overlap;
            const Real overlap_error = std::abs(difference);
            EXPECT_LE(overlap_error, 1e-10) << "S-overlap at (" << k << ", " << i << ")";
        }
    }
    const double iteration_budget = config.max_iterations;
    EXPECT_LE(iterations, iteration_budget);
}

static const PPCGOverlapCase generalized_overlap_cases[] = {
    // name, seed, bands, leading dimension, iteration budget, tridiagonal, residual tolerance
    {"Diagonal", 333, 3, 8, 100, false, 1e-6},
    {"Tridiagonal", 444, 2, 6, 150, true, 1e-8}};

INSTANTIATE_TEST_SUITE_P(
    OverlapOperators,
    DiagoPPCGGeneralizedTest,
    ::testing::ValuesIn(generalized_overlap_cases),
    [](const ::testing::TestParamInfo<PPCGOverlapCase>& info) { return info.param.name; });

// Test fixture: eigenvector quality — verify Hψ ≈ εψ and ψ^H ψ = I
// Uses the 10×10 tridiagonal Laplacian. After convergence, check:
//   1. ||Hψ_i - ε_i ψ_i|| < tol   for each band
//   2. |ψ_i^H ψ_j - δ_ij| < tol    for all i,j
TEST_F(DiagoPPCGMatrixTest, EigenvectorQuality)
{
    initialize(10, 3, 888, 1e-8);
    set_laplacian(1.0);
    check_eigenvalues(1e-12, 100, 3, 3, 1e-8);
    check_eigenvectors(1e-4, 1e-8);
}

// Test fixture: pentadiagonal Toeplitz (discrete biharmonic operator)
// H[i,i]=6, H[i,i±1]=-4, H[i,i±2]=1.  Eigenvalues:
// λ_k = 16·sin⁴(k·π / (2·(n+1))),  k = 1,...,n
// Wider bandwidth (5 vs 3) tests the solver with more off-diagonal coupling.
TEST_F(DiagoPPCGMatrixTest, SquaredLaplacian)
{
    initialize(8, 4, 505, 1e-10);
    // H = T² where T is the tridiagonal Laplacian (2 on diag, -1 on off-diag).
    // The corners of T² have diag=5 (not 6) since (T²)[0,0] = 2² + (-1)² = 5.
    // Interior: (T²)[i,i] = (-1)² + 2² + (-1)² = 6.
    H_mat.assign(n_dim * n_dim, T(0));
    for (int i = 0; i < n_dim; ++i)
    {
        const Real diagonal_value = (i == 0 || i == n_dim - 1) ? 5.0 : 6.0;
        H_mat[i + i * n_dim] = T(diagonal_value, 0);
        if (i >= 1)
        {
            H_mat[i + (i - 1) * n_dim] = T(-4.0, 0);
        }
        if (i < n_dim - 1)
        {
            H_mat[i + (i + 1) * n_dim] = T(-4.0, 0);
        }
        if (i >= 2)
        {
            H_mat[i + (i - 2) * n_dim] = T(1.0, 0);
        }
        if (i < n_dim - 2)
        {
            H_mat[i + (i + 2) * n_dim] = T(1.0, 0);
        }
    }

    prec.assign(n_dim, 6.0);
    prec[0] = 5.0;
    prec[n_dim - 1] = 5.0;

    exact.resize(nband);
    for (int k = 0; k < nband; ++k)
    {
        Real theta = Real(k + 1) * PPCG_PI / Real(2 * (n_dim + 1));
        Real s = std::sin(theta);
        exact[k] = Real(16) * s * s * s * s;
    }
    check_eigenvalues(1e-12, 150, 4, 4, 1e-8);
}

// Test fixture: gapped spectrum
// H = diag(0.1, 0.5, 5.0, 6.0, 10.0), nband=3.
// Large gaps between eigenvalue groups test the solver's band separation.
TEST_F(DiagoPPCGMatrixTest, GappedSpectrum)
{
    initialize(5, 3, 606, 1e-10);
    H_mat.assign(n_dim * n_dim, T(0));
    H_mat[0 + 0 * n_dim] = T(0.1, 0);
    H_mat[1 + 1 * n_dim] = T(0.5, 0);
    H_mat[2 + 2 * n_dim] = T(5.0, 0);
    H_mat[3 + 3 * n_dim] = T(6.0, 0);
    H_mat[4 + 4 * n_dim] = T(10.0, 0);

    prec.assign(n_dim, 1.0);

    exact = {0.1, 0.5, 5.0};
    check_eigenvalues(1e-12, 100, 3, 3, 1e-8);
}

// Test fixture: preconditioner stress test
// Uses the 10×10 tridiagonal Laplacian with a fixed, nonuniform positive
// preconditioner that differs from the exact diagonal. The test checks that
// the solver still converges and returns a small independently computed residual.
TEST_F(DiagoPPCGMatrixTest, NonuniformPreconditioner)
{
    initialize(10, 3, 707, 1e-10);
    set_laplacian(1.0);
    for (int i = 0; i < n_dim; ++i)
    {
        prec[i] = 0.6 + 0.1 * Real(i);
    }
    hsolver::DiagoPPCG<T, base_device::DEVICE_CPU> solver(1e-12, 200, 4, 4);
    const double iterations = solve(solver);
    EXPECT_TRUE(solver.converged());
    EXPECT_LE(iterations, 200.0);
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], exact[i], 1e-8) << "band " << i;
    }
    check_eigenvectors(1e-8, 1e-8);
}

// Test fixture: n_dim = 1, nband = 1 — absolute minimum
// H is a 1×1 matrix [5.0], eigenvalue = 5.0
TEST_F(DiagoPPCGMatrixTest, ScalarMatrix)
{
    initialize(1, 1, 0, 1e-10);
    H_mat = {T(5.0, 0)};
    prec = {5.0};
    exact = {5.0};
    psi = {T(1.0, 0)}; // already normalized
    check_eigenvalues(1e-12, 10, 1, 1, 1e-8);
}

// Test fixture: Neumann boundary tridiagonal Laplacian
// H[0,0]=H[n-1,n-1]=1, interior diag=2, off-diag=-1.
// Eigenvalues: λ_k = 2 - 2·cos(k·π/n), k=0,1,...,n-1
TEST_F(DiagoPPCGMatrixTest, NeumannZeroEigenvalue)
{
    initialize(8, 4, 222, 1e-10);
    H_mat.assign(n_dim * n_dim, T(0));
    for (int i = 0; i < n_dim; ++i)
    {
        const Real diagonal_value = (i == 0 || i == n_dim - 1) ? 1.0 : 2.0;
        H_mat[i + i * n_dim] = T(diagonal_value, 0);
        if (i > 0)
        {
            H_mat[i + (i - 1) * n_dim] = T(-1.0, 0);
        }
        if (i < n_dim - 1)
        {
            H_mat[i + (i + 1) * n_dim] = T(-1.0, 0);
        }
    }
    prec.assign(n_dim, 2.0);
    prec[0] = 1.0;
    prec[n_dim - 1] = 1.0;
    exact.resize(nband);
    for (int k = 0; k < nband; ++k)
    {
        exact[k] = 2.0 - 2.0 * std::cos(Real(k) * PPCG_PI / Real(n_dim));
    }
    check_eigenvalues(1e-12, 100, 4, 4, 1e-8);
}

// A deterministic sparse matrix exercises convergence beyond analytical fixtures.
TEST(DiagoPPCGRandomMatrixTest, ConvergesWithDiagonalPreconditioner)
{
    const int n_dim = 80;
    const int nband = 8;
    const int sparsity_pct = 60;
    const int matrix_size = n_dim * n_dim;
    const unsigned seed = n_dim * 100 + sparsity_pct;
    std::vector<T> h(matrix_size, T(0));
    std::mt19937 rng(seed);
    std::uniform_real_distribution<Real> dist(-1.0, 1.0);
    for (int i = 0; i < n_dim; ++i)
    {
        for (int j = i; j < n_dim; ++j)
        {
            if (i != j && (rng() % 100) < sparsity_pct)
            {
                continue;
            }
            const Real value = (i == j) ? std::abs(dist(rng)) * n_dim + 1.0 : dist(rng) * 0.5;
            h[i + j * n_dim] = T(value, 0);
            h[j + i * n_dim] = T(value, 0);
        }
    }
    std::vector<Real> prec(n_dim);
    for (int i = 0; i < n_dim; ++i)
    {
        const Real diagonal = std::real(h[i + i * n_dim]);
        prec[i] = std::max(diagonal, 1e-6);
    }
    std::vector<T> psi;
    initialize_vectors(psi, n_dim, nband, n_dim, 42);
    std::vector<Real> eval(nband, 0.0);
    const std::vector<double> ethr(nband, 1e-4);
    const auto h_op = [&h, n_dim](T* in, T* out, int ld, int ncol) {
        dense_h_multiply(h.data(), n_dim, in, out, ld, ncol);
    };
    hsolver::DiagoPPCG<T, base_device::DEVICE_CPU> solver(1e-8, 500, nband, 4);
    const double iterations = solver.diag(h_op, nullptr, n_dim, nband, n_dim,
                                          psi.data(), eval.data(), ethr, prec.data());
    EXPECT_TRUE(solver.converged());
    EXPECT_TRUE(std::isfinite(solver.max_residual()));
    EXPECT_LE(solver.max_residual(), 1e-4);
    EXPECT_LE(iterations, 500.0);
}

int main(int argc, char** argv)
{
#ifdef __MPI
    int nproc = 1;
    int myrank = 0;
    setupmpi(argc, argv, nproc, myrank);
    int nproc_in_pool = 1;
    int kpar = nproc;
    int mypool = 0;
    int rank_in_pool = 0;
    divide_pools(nproc, myrank, nproc_in_pool, kpar, mypool, rank_in_pool);
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    finishmpi();
#endif
    return result;
}
