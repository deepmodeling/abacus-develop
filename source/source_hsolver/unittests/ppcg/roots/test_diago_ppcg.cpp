// Low-root recovery with startup correction, locking and changed Hamiltonians.
#include "../../diago_ppcg_test_support.h"

#include <algorithm>
#include <string>

using ppcg_test::DiagoPPCGMatrixTest;
using ppcg_test::Real;
using ppcg_test::T;
using ppcg_test::dense_h_multiply;

// Two explicit startup subspaces exercise the same low-root recovery check.
class DiagoPPCGStartupCorrectionTest : public DiagoPPCGMatrixTest,
                                     public ::testing::WithParamInterface<bool>
{
};

TEST_P(DiagoPPCGStartupCorrectionTest, KeepsLowestRoots)
{
    n_dim = 8;
    nband = 4;
    ld = n_dim;
    const Real inv_sqrt_two = Real(1) / std::sqrt(Real(2));
    const int matrix_size = n_dim * n_dim;
    H_mat.assign(matrix_size, T(0));
    prec.resize(n_dim);
    for (int i = 0; i < n_dim; ++i)
    {
        const Real diagonal = Real(i + 1);
        H_mat[i + i * n_dim] = T(diagonal, 0);
        prec[i] = diagonal;
    }

    const int vector_size = ld * nband;
    psi.assign(vector_size, T(0));
    psi[0] = T(inv_sqrt_two, 0);
    psi[2] = T(inv_sqrt_two, 0);
    psi[1 + ld] = T(inv_sqrt_two, 0);
    psi[3 + ld] = T(inv_sqrt_two, 0);
    const bool all_vectors_active = GetParam();
    if (all_vectors_active)
    {
        // Every initial vector has a nonzero residual. Independent two-column
        // blocks would nevertheless retain roots 1, 2, 5, and 6, discarding the
        // third and fourth roots before a global comparison is possible.
        psi[4 + 2 * ld] = T(inv_sqrt_two, 0);
        psi[6 + 2 * ld] = T(inv_sqrt_two, 0);
        psi[5 + 3 * ld] = T(inv_sqrt_two, 0);
        psi[7 + 3 * ld] = T(inv_sqrt_two, 0);
    }
    else
    {
        // The first two columns contain the missing third and fourth roots. The
        // last two are exact high-energy eigenvectors and would otherwise lock
        // before their lower replacements can enter the active subspace.
        psi[6 + 2 * ld] = T(1, 0);
        psi[7 + 3 * ld] = T(1, 0);
    }

    ethr.assign(nband, 1e-12);
    hsolver::DiagoPPCG<T, base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 50,
        /* sbsize   = */ 2,
        /* rr_step  = */ 4);
    const double avg_iter = solve(solver);
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], Real(i + 1), 1e-10) << "root[" << i << "] mismatch";
    }
    EXPECT_TRUE(solver.converged());
    EXPECT_LE(avg_iter, 50.0);
}

INSTANTIATE_TEST_SUITE_P(
    StartupVectors,
    DiagoPPCGStartupCorrectionTest,
    ::testing::Bool(),
    [](const ::testing::TestParamInfo<bool>& info) {
        return info.param ? "AllVectorsActive" : "LockedHighRoots";
    });

TEST(DiagoPPCGRootSelectionTest, StartupCorrectionKeepsConvergedPrefixInactive)
{
    const int n_dim = 10;
    const int nband = 6;
    const int ld = n_dim;
    const Real inv_sqrt_two = Real(1) / std::sqrt(Real(2));

    std::vector<T> h_mat(n_dim * n_dim, T(0));
    std::vector<Real> prec(n_dim);
    for (int i = 0; i < n_dim; ++i)
    {
        const Real diagonal = Real(i + 1);
        h_mat[i + i * n_dim] = T(diagonal, 0);
        prec[i] = diagonal;
    }

    // The first two roots are already converged.  The active tail still needs
    // one global comparison: its first two vectors contain roots 3--6, while
    // its last two vectors are exact high-energy states.
    std::vector<T> psi(ld * nband, T(0));
    psi[0] = T(1, 0);
    psi[1 + ld] = T(1, 0);
    psi[2 + 2 * ld] = T(inv_sqrt_two, 0);
    psi[4 + 2 * ld] = T(inv_sqrt_two, 0);
    psi[3 + 3 * ld] = T(inv_sqrt_two, 0);
    psi[5 + 3 * ld] = T(inv_sqrt_two, 0);
    psi[8 + 4 * ld] = T(1, 0);
    psi[9 + 5 * ld] = T(1, 0);

    std::vector<Real> eval(nband, Real(0));
    const std::vector<double> ethr(nband, 1e-12);
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 1,
        /* sbsize   = */ 2,
        /* rr_step  = */ 4);
    const auto h_op = [&h_mat, n_dim](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(h_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    const double avg_iter = solver.diag(
        h_op, nullptr, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());

    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], Real(i + 1), 1e-10) << "root[" << i << "] mismatch";
    }
    const double expected_avg_iter = 1.0 + 4.0 / 6.0;
    EXPECT_TRUE(solver.converged());
    EXPECT_NEAR(avg_iter, expected_avg_iter, 1e-12);
}

TEST(DiagoPPCGRootSelectionTest, WarmStartAfterHamiltonianChangeKeepsLowestRoots)
{
    const int n_dim = 8;
    const int nband = 4;
    const int ld = n_dim;
    const Real inv_sqrt_two = Real(1) / std::sqrt(Real(2));

    std::vector<T> h_mat(n_dim * n_dim, T(0));
    std::vector<Real> prec(n_dim, Real(1));
    std::vector<T> psi(ld * nband, T(0));
    const int first_indices[nband] = {0, 1, 4, 5};
    const int second_indices[nband] = {2, 3, 6, 7};
    for (int state = 0; state < nband; ++state)
    {
        const int first = first_indices[state];
        const int second = second_indices[state];
        const Real low_eigenvalue = Real(state + 1);
        const Real high_eigenvalue = Real(state + 5);
        const Real diagonal = Real(0.5) * (low_eigenvalue + high_eigenvalue);
        const Real off_diagonal = Real(0.5) * (low_eigenvalue - high_eigenvalue);
        h_mat[first + first * n_dim] = T(diagonal, 0);
        h_mat[second + second * n_dim] = T(diagonal, 0);
        h_mat[first + second * n_dim] = T(off_diagonal, 0);
        h_mat[second + first * n_dim] = T(off_diagonal, 0);
        psi[first + state * ld] = T(inv_sqrt_two, 0);
        psi[second + state * ld] = T(inv_sqrt_two, 0);
        prec[first] = diagonal;
        prec[second] = diagonal;
    }

    std::vector<Real> eval(nband, Real(0));
    std::vector<double> ethr(nband, 1e-12);
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 50,
        /* sbsize   = */ 2,
        /* rr_step  = */ 4);
    auto h_op = [&h_mat, n_dim](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(h_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    solver.diag(
        h_op, nullptr, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());
    ASSERT_TRUE(solver.converged());
    for (int i = 0; i < nband; ++i)
    {
        ASSERT_NEAR(eval[i], Real(i + 1), 1e-10);
    }

    // Model a later SCF call: the previous eigenvectors are a valid warm start,
    // but the changed Hamiltonian makes each one mix a low and a high root.
    std::fill(h_mat.begin(), h_mat.end(), T(0));
    for (int i = 0; i < n_dim; ++i)
    {
        const Real diagonal = Real(i + 1);
        h_mat[i + i * n_dim] = T(diagonal, 0);
        prec[i] = diagonal;
    }

    solver.diag(
        h_op, nullptr, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());
    EXPECT_TRUE(solver.converged());
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], Real(i + 1), 1e-10) << "root[" << i << "] mismatch";
    }
}

TEST(DiagoPPCGRootSelectionTest, StartupCorrectionKeepsLowestGeneralizedRoots)
{
    const int n_dim = 8;
    const int nband = 4;
    const int ld = n_dim;
    const Real inv_sqrt_two = Real(1) / std::sqrt(Real(2));

    std::vector<T> h_mat(n_dim * n_dim, T(0));
    std::vector<T> s_mat(n_dim * n_dim, T(0));
    std::vector<Real> s_diagonal(n_dim);
    std::vector<Real> prec(n_dim);
    for (int i = 0; i < n_dim; ++i)
    {
        const Real eigenvalue = Real(i + 1);
        const Real overlap_value = Real(1) + Real(0.1) * Real(i % 3);
        const Real hamiltonian_value = eigenvalue * overlap_value;
        h_mat[i + i * n_dim] = T(hamiltonian_value, 0);
        s_mat[i + i * n_dim] = T(overlap_value, 0);
        s_diagonal[i] = overlap_value;
        prec[i] = hamiltonian_value;
    }

    // H = S diag(1,...,8), so the generalized eigenvectors are the
    // coordinate vectors.  The initial columns create the same cross-block
    // root-selection hazard as the identity-overlap regression above.
    std::vector<T> psi(ld * nband, T(0));
    std::vector<Real> normalized_coefficients(n_dim);
    for (int i = 0; i < n_dim; ++i)
    {
        normalized_coefficients[i] = inv_sqrt_two / std::sqrt(s_diagonal[i]);
    }
    psi[0] = T(normalized_coefficients[0], 0);
    psi[2] = T(normalized_coefficients[2], 0);
    psi[1 + ld] = T(normalized_coefficients[1], 0);
    psi[3 + ld] = T(normalized_coefficients[3], 0);
    psi[4 + 2 * ld] = T(normalized_coefficients[4], 0);
    psi[6 + 2 * ld] = T(normalized_coefficients[6], 0);
    psi[5 + 3 * ld] = T(normalized_coefficients[5], 0);
    psi[7 + 3 * ld] = T(normalized_coefficients[7], 0);

    std::vector<Real> eval(nband, Real(0));
    std::vector<double> ethr(nband, 1e-12);
    hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
        /* diag_thr = */ 1e-12,
        /* max_iter = */ 50,
        /* sbsize   = */ 2,
        /* rr_step  = */ 4);

    auto h_op = [&h_mat, n_dim](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(h_mat.data(), n_dim, in, out, ld_in, ncol);
    };
    auto s_op = [&s_mat, n_dim](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(s_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    const double avg_iter = solver.diag(
        h_op, s_op, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eval[i], Real(i + 1), 1e-10) << "root[" << i << "] mismatch";
    }
    EXPECT_TRUE(solver.converged());
    EXPECT_LE(avg_iter, 50.0);
}

TEST(DiagoPPCGRootSelectionTest, LockingKeepsRequestedLowRoots)
{
    const int n_dim = 4;
    const int nband = 2;
    const int ld = n_dim;
    std::vector<T> h_mat(n_dim * n_dim, T(0));
    h_mat[0] = T(1, 0);
    h_mat[1 + ld] = T(2, 0);
    h_mat[2 + 2 * ld] = T(4, 0);
    h_mat[3 + 3 * ld] = T(3, 0);
    const std::vector<Real> prec = {1, 2, 4, 3};
    const std::vector<double> ethr(nband, 1e-12);
    const Real epsilons[] = {1e-4, 1e-3, 1e-2, 1e-1};

    const auto h_op = [&h_mat, n_dim](T* in, T* out, int ld_in, int ncol) {
        dense_h_multiply(h_mat.data(), n_dim, in, out, ld_in, ncol);
    };

    for (const Real epsilon : epsilons)
    {
        SCOPED_TRACE("epsilon = " + std::to_string(epsilon));
        std::vector<T> psi(ld * nband, T(0));
        const Real norm = std::sqrt(Real(2) + epsilon * epsilon);
        psi[0] = T(Real(1) / norm, 0);
        psi[1] = T(epsilon / norm, 0);
        psi[2] = T(Real(1) / norm, 0);
        psi[3 + ld] = T(1, 0);
        std::vector<Real> eval(nband, Real(0));

        hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
            /* diag_thr = */ 1e-12,
            /* max_iter = */ 80,
            /* sbsize   = */ 2,
            /* rr_step  = */ 1);
        const double avg_iter = solver.diag(
            h_op, nullptr, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());
        EXPECT_NEAR(eval[0], 1.0, 1e-8);
        EXPECT_NEAR(eval[1], 2.0, 1e-8);
        EXPECT_TRUE(solver.converged());
        EXPECT_LE(avg_iter, 80.0);
    }
}

TEST(DiagoPPCGRootSelectionTest, LaterCorrectionKeepsLowestRootsAcrossBlocks)
{
    const int n_dim = 4;
    const int nband = 2;
    const int ld = n_dim;
    const std::vector<Real> diagonal = {1.0, 2.0, 3.0, 100.0};
    const std::vector<Real> prec(n_dim, 1.0);
    const std::vector<double> ethr(nband, 1.0e-10);
    const Real coefficient = 1.0 / std::sqrt(3.0);
    const int block_sizes[] = {1, 2};
    const auto h_op = [&diagonal, n_dim](T* in, T* out, int ld_in, int ncol) {
        for (int column = 0; column < ncol; ++column)
        {
            for (int row = 0; row < n_dim; ++row)
            {
                out[row + column * ld_in] = diagonal[row] * in[row + column * ld_in];
            }
        }
    };

    for (const int block_size : block_sizes)
    {
        SCOPED_TRACE("block_size = " + std::to_string(block_size));
        // The first expanded subspace keeps a low mixture and the exact
        // third root. The missing second root needs a later correction;
        // convergence of the retained vectors alone does not certify it.
        std::vector<T> psi(ld * nband, T(0));
        psi[0] = T(coefficient);
        psi[1] = T(coefficient);
        psi[3] = T(coefficient);
        psi[2 + ld] = T(1);
        std::vector<Real> eval(nband, 0.0);
        hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
            1.0e-10, 100, block_size, 1);
        solver.diag(h_op, nullptr, ld, nband, n_dim,
                    psi.data(), eval.data(), ethr, prec.data());

        double max_residual = 0.0;
        for (int column = 0; column < nband; ++column)
        {
            double residual_squared = 0.0;
            for (int row = 0; row < n_dim; ++row)
            {
                const T residual = (diagonal[row] - eval[column]) * psi[row + column * ld];
                residual_squared += std::norm(residual);
            }
            const double residual_norm = std::sqrt(residual_squared);
            max_residual = std::max(max_residual, residual_norm);
            for (int other = 0; other < nband; ++other)
            {
                T overlap(0);
                for (int row = 0; row < n_dim; ++row)
                {
                    overlap += std::conj(psi[row + column * ld]) * psi[row + other * ld];
                }
                const T expected_overlap = column == other ? T(1) : T(0);
                const double overlap_error = std::abs(overlap - expected_overlap);
                EXPECT_LE(overlap_error, 1.0e-12);
            }
        }
        EXPECT_TRUE(solver.converged());
        EXPECT_LE(max_residual, 1.0e-10);
        EXPECT_NEAR(solver.max_residual(), max_residual, 1.0e-13);
        EXPECT_NEAR(eval[0], 1.0, 1.0e-9);
        EXPECT_NEAR(eval[1], 2.0, 1.0e-9)
            << "returned roots = " << eval[0] << ", " << eval[1]
            << "; fresh max residual = " << max_residual;
        if (block_size == 1)
        {
            EXPECT_GT(solver.root_recovery_count(), 0);
        }
        else
        {
            EXPECT_EQ(solver.root_recovery_count(), 0);
        }
    }
}

TEST(DiagoPPCGRootSelectionTest, LaterCorrectionsRecoverMultipleComplexGeneralizedRoots)
{
    const Real pi = 3.14159265358979323846;
    const int dimension = 6;
    const int bands = 3;
    const Real spectrum[] = {1.0, 2.0, 3.0, 4.0, 5.0, 100.0};
    const Real overlap_spectrum[] = {1.0, 2.0, 4.0, 7.0, 3.0, 5.0};
    std::vector<T> unitary(dimension * dimension);
    const Real normalization = 1.0 / std::sqrt(Real(dimension));
    for (int column = 0; column < dimension; ++column)
    {
        for (int row = 0; row < dimension; ++row)
        {
            const Real angle = 2.0 * pi * row * column / dimension;
            const T phase(std::cos(angle), std::sin(angle));
            unitary[row + column * dimension] = normalization * phase;
        }
    }
    std::vector<T> hamiltonian(dimension * dimension, T(0));
    std::vector<T> overlap(dimension * dimension, T(0));
    for (int column = 0; column < dimension; ++column)
    {
        for (int row = 0; row < dimension; ++row)
        {
            for (int state = 0; state < dimension; ++state)
            {
                const T projector = unitary[row + state * dimension]
                                    * std::conj(unitary[column + state * dimension]);
                hamiltonian[row + column * dimension] += projector * overlap_spectrum[state] * spectrum[state];
                overlap[row + column * dimension] += projector * overlap_spectrum[state];
            }
        }
    }
    const auto h_op = [&hamiltonian, dimension](T* in, T* out, int ld, int count) {
        dense_h_multiply(hamiltonian.data(), dimension, in, out, ld, count);
    };
    const auto s_op = [&overlap, dimension](T* in, T* out, int ld, int count) {
        dense_h_multiply(overlap.data(), dimension, in, out, ld, count);
    };
    const std::vector<Real> preconditioner(dimension, 1.0);
    const std::vector<double> thresholds(bands, 1.0e-9);
    const int block_sizes[] = {1, 2, 3};
    const int ritz_intervals[] = {1, 4};
    for (const int block_size : block_sizes)
    {
        for (const int ritz_interval : ritz_intervals)
        {
            SCOPED_TRACE("block_size = " + std::to_string(block_size)
                         + ", ritz_interval = " + std::to_string(ritz_interval));
            std::vector<T> psi(dimension * bands, T(0));
            // Four spectral components in the first vector require more
            // than the first [X,W] expansion to resolve all three low roots.
            const int mixed_states[] = {0, 1, 2, 5};
            for (int row = 0; row < dimension; ++row)
            {
                for (const int state : mixed_states)
                {
                    psi[row] += 0.5 * unitary[row + state * dimension]
                                / std::sqrt(overlap_spectrum[state]);
                }
                psi[row + dimension] = unitary[row + 3 * dimension] / std::sqrt(overlap_spectrum[3]);
                psi[row + 2 * dimension] = unitary[row + 4 * dimension] / std::sqrt(overlap_spectrum[4]);
            }
            std::vector<Real> eigenvalues(bands, 0.0);
            hsolver::DiagoPPCG<T, hsolver::base_device::DEVICE_CPU> solver(
                1.0e-9, 120, block_size, ritz_interval);
            solver.diag(h_op, s_op, dimension, bands, dimension,
                        psi.data(), eigenvalues.data(), thresholds, preconditioner.data());
            std::vector<T> hpsi(psi.size(), T(0));
            std::vector<T> spsi(psi.size(), T(0));
            h_op(psi.data(), hpsi.data(), dimension, bands);
            s_op(psi.data(), spsi.data(), dimension, bands);
            EXPECT_TRUE(solver.converged());
            for (int column = 0; column < bands; ++column)
            {
                EXPECT_NEAR(eigenvalues[column], spectrum[column], 1.0e-8);
                double residual_squared = 0.0;
                for (int row = 0; row < dimension; ++row)
                {
                    const int offset = row + column * dimension;
                    const T residual = hpsi[offset] - eigenvalues[column] * spsi[offset];
                    residual_squared += std::norm(residual);
                }
                const double residual_norm = std::sqrt(residual_squared);
                EXPECT_LE(residual_norm, 1.0e-9);
                for (int other = 0; other < bands; ++other)
                {
                    T inner_product(0);
                    for (int row = 0; row < dimension; ++row)
                    {
                        inner_product += std::conj(psi[row + column * dimension]) * spsi[row + other * dimension];
                    }
                    const T expected = column == other ? T(1) : T(0);
                    const double orthogonality_error = std::abs(inner_product - expected);
                    EXPECT_LE(orthogonality_error, 1.0e-11);
                }
            }
        }
    }
}
