// Local subspace transformations, rank truncation and input/storage contracts.
#include "../../diago_ppcg_test_support.h"

#include <algorithm>
#include <limits>
#include <stdexcept>

using ppcg_test::DiagoPPCGMatrixTest;
using ppcg_test::PPCG_PI;
using ppcg_test::Real;
using ppcg_test::T;
using ppcg_test::dense_h_multiply;
using ppcg_test::initialize_vectors;

namespace hsolver
{

// Narrow friend access exercises cancellation in the projection helper without
// relying on a particular sequence of LAPACK rotations in an end-to-end solve.
class DiagoPPCGProjectionTest : public ::testing::Test
{
protected:
    void check_direction_rotation(const bool ritz)
    {
        const int dimension = 3;
        const int leading_dimension = 4;
        const int bands = 2;
        const int size = leading_dimension * bands;
        const std::vector<T> h = {T(1), T(0, -0.4), T(0),
                                  T(0, 0.4), T(4), T(0), T(0), T(0), T(10)};
        const std::vector<T> s = {T(1), T(0), T(0), T(0), T(2), T(0), T(0), T(0), T(3)};
        std::vector<T> x(size, T(0));
        x[0] = T(1);
        x[leading_dimension] = T(0.3);
        x[leading_dimension + 1] = T(1);
        DiagoPPCG<T, base_device::DEVICE_CPU> solver(1e-12, 20, 1, 4);
        solver.ld_psi_ = leading_dimension;
        solver.n_dim_ = dimension;
        solver.n_band_ = bands;
        solver.hpsi_.assign(size, T(0));
        solver.spsi_.assign(size, T(0));
        solver.w_.assign(size, T(0));
        solver.sw_.assign(size, T(0));
        solver.hw_.assign(size, T(0));
        solver.p_.assign(size, T(0));
        solver.p_[2] = T(0.4, 0.3);
        solver.p_[leading_dimension + 2] = T(-0.2, 0.6);
        solver.hp_.assign(size, T(0));
        solver.sp_.assign(size, T(0));
        solver.rr_eval_.resize(bands);
        dense_h_multiply(h.data(), dimension, x.data(), solver.hpsi_.data(), leading_dimension, bands);
        dense_h_multiply(s.data(), dimension, x.data(), solver.spsi_.data(), leading_dimension, bands);
        dense_h_multiply(h.data(), dimension, solver.p_.data(), solver.hp_.data(), leading_dimension, bands);
        dense_h_multiply(s.data(), dimension, solver.p_.data(), solver.sp_.data(), leading_dimension, bands);
        const std::vector<T> original_p = solver.p_;
        if (ritz)
        {
            std::vector<Real> eigenvalues(bands);
            const std::vector<double> thresholds(bands, 1e-12);
            std::vector<int> active;
            ASSERT_TRUE(solver.rayleigh_ritz(x.data(), eigenvalues.data(), active, thresholds));
        }
        else
        {
            ASSERT_TRUE(solver.orthonormalize(x.data()));
        }
        double change = 0.0;
        for (int column = 0; column < bands; ++column)
        {
            // Recover C from X_new = X_old C using the known triangular
            // leading 2-by-2 block, independently of the solver's workspace.
            const int first = column * leading_dimension;
            const T c0 = x[first] - 0.3 * x[first + 1];
            const T c1 = x[first + 1];
            for (int row = 0; row < dimension; ++row)
            {
                const int index = first + row;
                const T expected = original_p[row] * c0 + original_p[row + leading_dimension] * c1;
                const T direction_difference = solver.p_[index] - expected;
                const double error = std::abs(direction_difference);
                EXPECT_LT(error, 1e-13);
                const T direction_change = solver.p_[index] - original_p[index];
                change += std::norm(direction_change);
            }
            EXPECT_EQ(solver.p_[first + dimension], T(0));
        }
        EXPECT_GT(change, 1e-3);
        std::vector<T> hp(size, T(0));
        std::vector<T> sp(size, T(0));
        dense_h_multiply(h.data(), dimension, solver.p_.data(), hp.data(), leading_dimension, bands);
        dense_h_multiply(s.data(), dimension, solver.p_.data(), sp.data(), leading_dimension, bands);
        for (int index = 0; index < size; ++index)
        {
            const T h_difference = hp[index] - solver.hp_[index];
            const T s_difference = sp[index] - solver.sp_[index];
            const double h_error = std::abs(h_difference);
            const double s_error = std::abs(s_difference);
            EXPECT_LT(h_error, 1e-13);
            EXPECT_LT(s_error, 1e-13);
        }
    }

    void check_rank_truncated_outputs(const int kept_count)
    {
        const int dimension = 8;
        const int rank = 7;
        const int matrix_size = dimension * dimension;
        DiagoPPCG<T, base_device::DEVICE_CPU> solver(1.0e-12, 20, kept_count, 1);
        DiagoPPCG<T, base_device::DEVICE_CPU>::SmallSubspace subspace;
        subspace.k.assign(matrix_size, T(0));
        subspace.m.assign(matrix_size, T(0));
        subspace.eval.assign(dimension, 0.0);
        for (int column = 0; column < rank; ++column)
        {
            const double weight = column + 1.0;
            subspace.m[column + column * dimension] = T(weight);
            subspace.k[column + column * dimension] = T(weight * weight);
            if (column + 1 < rank)
            {
                const double scale = std::sqrt(weight * (weight + 1.0));
                const T coupling = T(-0.3, 0.2) * scale;
                subspace.k[column + (column + 1) * dimension] = coupling;
                subspace.k[column + 1 + column * dimension] = std::conj(coupling);
            }
        }
        const std::vector<T> hamiltonian = subspace.k;
        const std::vector<T> overlap = subspace.m;
        const int valid_count = solver.solve_small_generalized(dimension, kept_count, subspace);
        ASSERT_EQ(valid_count, rank);
        // Check every exported vector, including those not used for the local
        // update: the recovery path must receive valid, normalized Ritz pairs.
        for (int column = 0; column < valid_count; ++column)
        {
            double residual_squared = 0.0;
            for (int row = 0; row < dimension; ++row)
            {
                T residual(0);
                for (int inner = 0; inner < dimension; ++inner)
                {
                    const int offset = row + inner * dimension;
                    residual += (hamiltonian[offset] - subspace.eval[column] * overlap[offset])
                                * subspace.k[inner + column * dimension];
                }
                residual_squared += std::norm(residual);
            }
            const double residual_norm = std::sqrt(residual_squared);
            EXPECT_LT(residual_norm, 1.0e-12);
            for (int other = 0; other < valid_count; ++other)
            {
                T inner_product(0);
                for (int row = 0; row < rank; ++row)
                {
                    inner_product += std::conj(subspace.k[row + column * dimension])
                                     * (row + 1.0) * subspace.k[row + other * dimension];
                }
                const T expected = column == other ? T(1) : T(0);
                const double error = std::abs(inner_product - expected);
                EXPECT_LT(error, 1.0e-12);
            }
        }
    }

    void check_projection(const bool nearly_parallel, const bool zero_direction)
    {
        const int dimension = 2;
        const int bands = 1;
        const std::vector<Real> h_diagonal = {1.0, 1000.0};
        const std::vector<Real> s_diagonal = {1.0, 3.0};
        const std::vector<T> x = {T(0.5), T(0.5)};
        const double scale = 1.0 / std::sqrt(12.0);
        const std::vector<T> z = {T(3.0 * scale), T(-scale)};
        DiagoPPCG<T, base_device::DEVICE_CPU> solver(1.0e-12, 20, bands, 1);
        solver.ld_psi_ = dimension;
        solver.n_dim_ = dimension;
        solver.n_band_ = bands;
        solver.hpsi_.resize(dimension);
        solver.spsi_.resize(dimension);
        solver.p_.resize(dimension);
        solver.hp_.resize(dimension);
        solver.sp_.resize(dimension);
        for (int i = 0; i < dimension; ++i)
        {
            solver.hpsi_[i] = h_diagonal[i] * x[i];
            solver.spsi_[i] = s_diagonal[i] * x[i];
            solver.p_[i] = z[i];
            if (nearly_parallel)
            {
                solver.p_[i] = x[i] + 1.0e-10 * z[i];
            }
            solver.hp_[i] = h_diagonal[i] * solver.p_[i];
            solver.sp_[i] = s_diagonal[i] * solver.p_[i];
            if (zero_direction)
            {
                solver.p_[i] = T(0);
                solver.hp_[i] = T(1.0e-15);
                solver.sp_[i] = T(1.0e-15);
            }
        }
        int h_calls = 0;
        int s_calls = 0;
        const auto h_op = [&h_calls, &h_diagonal, dimension](T* in, T* out, int ld, int ncol) {
            ++h_calls;
            for (int j = 0; j < ncol; ++j)
            {
                for (int i = 0; i < dimension; ++i)
                {
                    out[i + j * ld] = h_diagonal[i] * in[i + j * ld];
                }
            }
        };
        solver.spsi_func_ = [&s_calls, &s_diagonal, dimension](T* in, T* out, int ld, int ncol) {
            ++s_calls;
            for (int j = 0; j < ncol; ++j)
            {
                for (int i = 0; i < dimension; ++i)
                {
                    out[i + j * ld] = s_diagonal[i] * in[i + j * ld];
                }
            }
        };
        const std::vector<int> columns = {0};
        solver.project_search_directions(h_op, x.data(), columns, columns);
        double direction_norm_squared = 0.0;
        double h_error_squared = 0.0;
        double s_error_squared = 0.0;
        for (int i = 0; i < dimension; ++i)
        {
            direction_norm_squared += std::norm(solver.p_[i]);
            h_error_squared += std::norm(solver.hp_[i] - h_diagonal[i] * solver.p_[i]);
            s_error_squared += std::norm(solver.sp_[i] - s_diagonal[i] * solver.p_[i]);
        }
        if (zero_direction)
        {
            EXPECT_DOUBLE_EQ(direction_norm_squared, 0.0);
            EXPECT_DOUBLE_EQ(h_error_squared, 0.0);
            EXPECT_DOUBLE_EQ(s_error_squared, 0.0);
            EXPECT_EQ(h_calls, 0);
            EXPECT_EQ(s_calls, 0);
            return;
        }
        ASSERT_GT(direction_norm_squared, 0.0);
        const double h_error = std::sqrt(h_error_squared / direction_norm_squared);
        const double s_error = std::sqrt(s_error_squared / direction_norm_squared);
        EXPECT_LT(h_error, 1.0e-12);
        EXPECT_LT(s_error, 1.0e-12);
        const int expected_calls = nearly_parallel ? 1 : 0;
        EXPECT_EQ(h_calls, expected_calls);
        EXPECT_EQ(s_calls, expected_calls);
    }
};

TEST_F(DiagoPPCGProjectionTest, RefreshesImagesAfterSevereProjectionCancellation)
{
    check_projection(true, false);
}

TEST_F(DiagoPPCGProjectionTest, RitzRotationPreservesDirectionColumnAssociation)
{
    check_direction_rotation(true);
}

TEST_F(DiagoPPCGProjectionTest, OrthonormalizationRotatesDirectionImagesTogether)
{
    check_direction_rotation(false);
}

TEST_F(DiagoPPCGProjectionTest, ReusesImagesForWellSeparatedDirection)
{
    check_projection(false, false);
}

TEST_F(DiagoPPCGProjectionTest, ClearsStaleImagesForZeroDirection)
{
    check_projection(false, true);
}

TEST_F(DiagoPPCGProjectionTest, RankTruncationExportsKeptAndDiscardedRitzVectors)
{
    const int kept_counts[] = {2, 3, 7};
    for (const int kept_count : kept_counts)
    {
        SCOPED_TRACE(kept_count);
        check_rank_truncated_outputs(kept_count);
    }
}

} // namespace hsolver

// Test fixture: diagonal matrix (simplest possible Hamiltonian)
class DiagoPPCGDiagonalTest : public DiagoPPCGMatrixTest
{
  protected:
    void SetUp() override
    {
        initialize(5, 3, 42, 1e-10);
        const int matrix_size = n_dim * n_dim;
        H_mat.assign(matrix_size, T(0));
        prec.resize(n_dim);
        for (int i = 0; i < n_dim; ++i)
        {
            prec[i] = Real(i + 1);
            H_mat[i + i * n_dim] = prec[i];
        }
        exact = {1.0, 2.0, 3.0};
    }
};

TEST_F(DiagoPPCGDiagonalTest, BlockSubspace)
{
    check_eigenvalues(1e-12, 50, 3, 3, 1e-8);
}

TEST(DiagoPPCGRankTruncationTest, ComplexHermitianStartupSubspace)
{
    const int n_dim = 6;
    const int nband = 3;
    const int ld = n_dim;
    const Real inv_sqrt_two = Real(1) / std::sqrt(Real(2));

    // Build H = U diag(1,...,6) U^H.  The first two columns of U form a
    // complex rotation, so the reduced Hamiltonian needs a conjugate
    // transpose when it is projected into the retained overlap subspace.
    std::vector<T> unitary(n_dim * n_dim, T(0));
    unitary[0] = T(inv_sqrt_two, 0);
    unitary[1] = T(0, inv_sqrt_two);
    unitary[n_dim] = T(0, inv_sqrt_two);
    unitary[1 + n_dim] = T(inv_sqrt_two, 0);
    for (int i = 2; i < n_dim; ++i)
    {
        unitary[i + i * n_dim] = T(1, 0);
    }

    std::vector<T> h_mat(n_dim * n_dim, T(0));
    std::vector<Real> eigenvalues(n_dim);
    for (int state = 0; state < n_dim; ++state)
    {
        eigenvalues[state] = Real(state + 1);
        for (int col = 0; col < n_dim; ++col)
        {
            for (int row = 0; row < n_dim; ++row)
            {
                h_mat[row + col * n_dim] += unitary[row + state * n_dim]
                                            * eigenvalues[state]
                                            * std::conj(unitary[col + state * n_dim]);
            }
        }
    }
    EXPECT_GT(std::abs(std::imag(h_mat[n_dim])), Real(0.1));

    std::vector<Real> prec(n_dim);
    for (int i = 0; i < n_dim; ++i)
    {
        const Real diagonal_value = std::real(h_mat[i + i * n_dim]);
        prec[i] = std::max(Real(1), diagonal_value);
    }

    // The third column is an exact high-energy state.  Its zero residual
    // makes [psi, w] rank deficient during the first global correction,
    // while the first two columns and their residuals span the three lowest
    // roots.  This exercises Q^H K Q for a complex Hermitian matrix.
    std::vector<T> psi(ld * nband, T(0));
    for (int i = 0; i < n_dim; ++i)
    {
        psi[i] = (unitary[i] + unitary[i + 3 * n_dim]) * inv_sqrt_two;
        psi[i + ld] = (unitary[i + n_dim] + unitary[i + 2 * n_dim]) * inv_sqrt_two;
        psi[i + 2 * ld] = unitary[i + 5 * n_dim];
    }

    std::vector<Real> eval(nband, Real(0));
    std::vector<double> ethr(nband, 1e-12);
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 30,
        /* sbsize   = */ 2,
        /* rr_step  = */ 4);

    auto h_op = [&h_mat, n_dim](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(h_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    const double avg_iter = solver.diag(
        h_op, nullptr, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], Real(i + 1), 1e-10) << "root[" << i << "] mismatch";
    }
    EXPECT_TRUE(solver.converged());
    EXPECT_LE(avg_iter, 30.0);
}

TEST_F(DiagoPPCGDiagonalTest, EmptyHOperatorThrows)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 50,
        /* sbsize   = */ 3,
        /* rr_step  = */ 3);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU>::HPsiFunc h_op;
    EXPECT_THROW(solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data()),
                 std::invalid_argument);
}

TEST_F(DiagoPPCGDiagonalTest, NonFiniteInputThrows)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 50,
        /* sbsize   = */ 3,
        /* rr_step  = */ 3);

    auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    std::vector<double> bad_ethr = ethr;
    bad_ethr[0] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), bad_ethr, prec.data()),
                 std::invalid_argument);

    std::vector<Real> bad_prec = prec;
    bad_prec[0] = std::numeric_limits<Real>::infinity();
    EXPECT_THROW(solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, bad_prec.data()),
                 std::invalid_argument);
}

TEST(DiagoPPCGWorkspaceTest, RejectsOversizedWorkspacesBeforeUsingInputBuffers)
{
    struct Shape
    {
        const char* description;
        int leading_dimension;
        int bands;
    };
    const int int_max = std::numeric_limits<int>::max();
    const Shape shapes[] = {
        {"wavefunction storage", int_max / 2 + 1, 2},
        {"local basis storage", int_max / 3 + 1, 2},
        {"projected matrix storage", 1, 32768},
        {"LAPACK workspace", 1, 20000},
        {"root recovery workspace", 1, 12000},
    };

    for (const Shape& shape : shapes)
    {
        SCOPED_TRACE(shape.description);
        T psi_value(1.0, 0.0);
        Real eigenvalue = -1.0;
        const Real preconditioner = 1.0;
        // The size diagnostic must precede threshold/buffer validation. An
        // empty threshold vector also keeps this regression bounded if that
        // diagnostic is accidentally removed: no huge allocations are needed.
        const std::vector<double> thresholds;
        int h_calls = 0;
        int s_calls = 0;
        const auto h_op = [&h_calls](T*, T*, int, int) { ++h_calls; };
        const auto s_op = [&s_calls](T*, T*, int, int) { ++s_calls; };
        hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(1e-12, 1, 1, 1);

        try
        {
            solver.diag(h_op,
                        s_op,
                        shape.leading_dimension,
                        shape.bands,
                        1,
                        &psi_value,
                        &eigenvalue,
                        thresholds,
                        &preconditioner);
            FAIL() << "Oversized PPCG workspace was accepted";
        }
        catch (const std::invalid_argument& error)
        {
            EXPECT_STREQ(error.what(), "PPCG: workspace dimensions exceed integer indexing or LAPACK limits.");
        }
        EXPECT_EQ(h_calls, 0);
        EXPECT_EQ(s_calls, 0);
        EXPECT_EQ(psi_value, T(1.0, 0.0));
        EXPECT_DOUBLE_EQ(eigenvalue, -1.0);
    }
}

TEST_F(DiagoPPCGDiagonalTest, HandlesSubblockSizeAboveBandCount)
{
    std::vector<T> psi_run = psi;
    std::vector<Real> eval(nband, 0.0);
    const int subblock_size = std::numeric_limits<int>::max();
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(1e-12, 100, subblock_size, 3);
    const auto h_op = [this](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    solver.diag(h_op, nullptr, ld, nband, n_dim, psi_run.data(), eval.data(), ethr, prec.data());

    EXPECT_TRUE(solver.converged());
    for (int band = 0; band < nband; ++band)
    {
        EXPECT_NEAR(eval[band], exact[band], ethr[band]);
    }
}

TEST_F(DiagoPPCGDiagonalTest, BlockSubspaceWithPadding)
{
    ld = 8;
    initialize_vectors(psi, n_dim, nband, ld, 7);
    for (int band = 0; band < nband; ++band)
    {
        for (int row = n_dim; row < ld; ++row)
        {
            psi[row + band * ld] = T(17.0, -3.0);
        }
    }
    check_eigenvalues(1e-12, 80, 2, 3, 1e-8);
}

TEST(DiagoPPCGComplexHermitianTest, DefaultKeepsImaginaryProjection)
{
    const int n_dim = 2;
    const int nband = 2;
    const int ld = n_dim;

    // H = [[2, i], [-i, 3]].  Dropping Im(<x|Hy>) would incorrectly
    // diagonalize diag(2, 3); the Hermitian eigenvalues are 2.5 +/- sqrt(1.25).
    std::vector<T> H_mat(n_dim * n_dim, T(0));
    H_mat[0 + 0 * n_dim] = T(2.0, 0.0);
    H_mat[1 + 1 * n_dim] = T(3.0, 0.0);
    H_mat[0 + 1 * n_dim] = T(0.0, 1.0);
    H_mat[1 + 0 * n_dim] = T(0.0, -1.0);

    std::vector<T> psi(ld * nband, T(0));
    psi[0 + 0 * ld] = T(1.0, 0.0);
    psi[1 + 1 * ld] = T(1.0, 0.0);

    std::vector<Real> prec(n_dim, 2.0);
    std::vector<double> ethr(nband, 1e-12);
    std::vector<Real> eval(nband, 0.0);

    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 10,
        /* sbsize   = */ 2,
        /* rr_step  = */ 1);

    auto h_op = [&H_mat, n_dim](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(H_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    solver.diag(h_op, nullptr, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());

    EXPECT_TRUE(std::isfinite(eval[0]));
    EXPECT_TRUE(std::isfinite(eval[1]));
    const Real delta = std::sqrt(1.25);
    EXPECT_NEAR(eval[0], 2.5 - delta, 1e-10);
    EXPECT_NEAR(eval[1], 2.5 + delta, 1e-10);
}

TEST(DiagoPPCGComplexHermitianTest, ExactImaginaryLeadingCoefficientNeedsNoProjection)
{
    const int n_dim = 3;
    const int nband = 1;
    const int ld = n_dim;
    const int matrix_size = n_dim * n_dim;
    const Real diagonal[n_dim] = {Real(1), Real(4), Real(9)};
    std::vector<T> h_mat(matrix_size, T(0));
    for (int i = 0; i < n_dim; ++i)
    {
        h_mat[i + i * n_dim] = T(diagonal[i], 0);
    }

    // A valid full-complex eigenvector may have an imaginary first
    // coefficient. Clearing that coefficient would destroy its norm and
    // rank. A Ritz solve is free to choose another overall phase, so check
    // the eigenpair and normalization rather than a particular phase.
    std::vector<T> psi(ld, T(0));
    psi[0] = T(0, 1);
    std::vector<Real> eval(nband, Real(0));
    const std::vector<double> ethr(nband, 1e-12);
    const std::vector<Real> prec(n_dim, Real(1));
    const auto h_op = [&h_mat, n_dim](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(h_mat.data(), n_dim, in, out, ld_in, ncol);
    };
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 0,
        /* sbsize   = */ 1,
        /* rr_step  = */ 4);

    solver.diag(
        h_op, nullptr, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());

    Real norm_squared = Real(0);
    Real residual_squared = Real(0);
    for (int i = 0; i < n_dim; ++i)
    {
        const T residual = diagonal[i] * psi[i] - eval[0] * psi[i];
        norm_squared += std::norm(psi[i]);
        residual_squared += std::norm(residual);
    }
    const Real residual_norm = std::sqrt(residual_squared);
    EXPECT_TRUE(solver.converged());
    EXPECT_FALSE(solver.iteration_limit_reached());
    EXPECT_NEAR(eval[0], diagonal[0], 1e-12);
    EXPECT_NEAR(norm_squared, Real(1), 1e-12);
    EXPECT_LE(residual_norm, Real(1e-12));
}

// Test fixture: degenerate eigenvalues
// H = I + J  (identity plus all-ones), 4×4.
// J has eigenvector [1,1,1,1]^T with eigenvalue 4.
// All vectors orthogonal to [1,1,1,1]^T are eigenvectors with eigenvalue 0.
// So H = I + J has: λ₁ = 1 (multiplicity 3), λ₄ = 5.
