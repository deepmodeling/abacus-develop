#ifndef ABACUS_HSOLVER_UNITTESTS_DIAGO_PPCG_TEST_SUPPORT_H
#define ABACUS_HSOLVER_UNITTESTS_DIAGO_PPCG_TEST_SUPPORT_H

#include "../diago_ppcg.h"

#include <gtest/gtest.h>
#include <cmath>
#include <complex>
#include <random>
#include <vector>

namespace ppcg_test
{
static constexpr double PPCG_PI = 3.14159265358979323846;

using T = std::complex<double>;
using Real = double;

// -----------------------------------------------------------------------------
// Helper: dense H-matrix times a set of column vectors
// H is stored column-major: H(row, col) = H_data[row + col * n_dim]
// -----------------------------------------------------------------------------
inline void dense_h_multiply(const T* H_data, int n_dim, const T* in, T* out, int ld, int ncol)
{
    for (int j = 0; j < ncol; ++j)
    {
        for (int i = 0; i < n_dim; ++i)
        {
            T sum = 0;
            for (int k = 0; k < n_dim; ++k)
            {
                sum += H_data[i + k * n_dim] * in[k + j * ld];
            }
            out[i + j * ld] = sum;
        }
    }
}

// Use the same seeded real initial vectors for the analytical matrix fixtures.
inline void initialize_vectors(std::vector<T>& psi, int n_dim, int nband, int ld, unsigned seed)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<Real> dist(-1.0, 1.0);
    const int size = ld * nband;
    psi.assign(size, T(0));
    for (int j = 0; j < nband; ++j)
    {
        for (int i = 0; i < n_dim; ++i)
        {
            psi[i + j * ld] = T(dist(rng), 0.0);
        }
    }
    for (int j = 0; j < nband; ++j)
    {
        for (int k = 0; k < j; ++k)
        {
            T dot = 0;
            for (int i = 0; i < n_dim; ++i)
            {
                dot += std::conj(psi[i + k * ld]) * psi[i + j * ld];
            }
            for (int i = 0; i < n_dim; ++i)
            {
                psi[i + j * ld] -= dot * psi[i + k * ld];
            }
        }
        Real norm_squared = 0;
        for (int i = 0; i < n_dim; ++i)
        {
            norm_squared += std::norm(psi[i + j * ld]);
        }
        const Real norm = std::sqrt(norm_squared);
        for (int i = 0; i < n_dim; ++i)
        {
            psi[i + j * ld] /= norm;
        }
    }
}

// Shared storage and independent checks for analytical matrix cases. Numerical
// regressions supply their own deliberately constructed vectors and callbacks.
class DiagoPPCGMatrixTest : public ::testing::Test
{
  protected:
    void initialize(int dimension, int bands, unsigned seed, double threshold)
    {
        n_dim = dimension;
        nband = bands;
        ld = dimension;
        ethr.assign(nband, threshold);
        initialize_vectors(psi, n_dim, nband, ld, seed);
    }

    void set_laplacian(double scale)
    {
        const int size = n_dim * n_dim;
        H_mat.assign(size, T(0));
        for (int i = 0; i < n_dim; ++i)
        {
            H_mat[i + i * n_dim] = 2.0 * scale;
            if (i > 0)
            {
                H_mat[i + (i - 1) * n_dim] = -scale;
                H_mat[i - 1 + i * n_dim] = -scale;
            }
        }
        const double diagonal = 2.0 * scale;
        prec.assign(n_dim, diagonal);
        exact.resize(nband);
        for (int k = 0; k < nband; ++k)
        {
            const double angle = Real(k + 1) * PPCG_PI / Real(n_dim + 1);
            exact[k] = scale * (2.0 - 2.0 * std::cos(angle));
        }
    }

    double solve(hsolver::DiagoPPCG<T, base_device::DEVICE_CPU>& solver)
    {
        eval.assign(nband, 0.0);
        const auto h_op = [this](T* in, T* out, int leading_dimension, int columns) {
            dense_h_multiply(H_mat.data(), n_dim, in, out, leading_dimension, columns);
        };
        return solver.diag(h_op, nullptr, ld, nband, n_dim, psi.data(), eval.data(), ethr, prec.data());
    }

    void check_eigenvalues(double diag_threshold, int max_iterations, int block_size,
                           int ritz_interval, double tolerance)
    {
        hsolver::DiagoPPCG<T, base_device::DEVICE_CPU> solver(
            diag_threshold, max_iterations, block_size, ritz_interval);
        const double iterations = solve(solver);
        for (int i = 0; i < nband; ++i)
        {
            EXPECT_NEAR(eval[i], exact[i], tolerance) << "eigenvalue " << i;
        }
        const double iteration_budget = max_iterations;
        EXPECT_LE(iterations, iteration_budget);
    }

    // Recompute operator images independently of the solver's cached vectors.
    void check_eigenvectors(double residual_tolerance, double orthogonality_tolerance)
    {
        const int image_size = ld * nband;
        std::vector<T> hpsi(image_size);
        dense_h_multiply(H_mat.data(), n_dim, psi.data(), hpsi.data(), ld, nband);
        for (int band = 0; band < nband; ++band)
        {
            Real residual_squared = 0;
            for (int row = 0; row < n_dim; ++row)
            {
                const int offset = row + band * ld;
                const T residual = hpsi[offset] - eval[band] * psi[offset];
                residual_squared += std::norm(residual);
            }
            const Real residual_norm = std::sqrt(residual_squared);
            EXPECT_TRUE(std::isfinite(residual_norm));
            EXPECT_LT(residual_norm, residual_tolerance) << "band " << band;
            for (int other = 0; other < nband; ++other)
            {
                T overlap = 0;
                for (int row = 0; row < n_dim; ++row)
                {
                    overlap += std::conj(psi[row + other * ld]) * psi[row + band * ld];
                }
                const T expected_overlap = (other == band) ? T(1) : T(0);
                const T difference = overlap - expected_overlap;
                const Real overlap_error = std::abs(difference);
                EXPECT_LT(overlap_error, orthogonality_tolerance) << "bands " << other << ", " << band;
            }
        }
    }

    int n_dim;
    int nband;
    int ld;
    std::vector<T> H_mat;
    std::vector<Real> prec;
    std::vector<Real> exact;
    std::vector<Real> eval;
    std::vector<double> ethr;
    std::vector<T> psi;
};

} // namespace ppcg_test

#endif
