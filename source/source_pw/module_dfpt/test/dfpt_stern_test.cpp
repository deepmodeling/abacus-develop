#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <vector>
#include "source_pw/module_dfpt/dfpt_stern.h"

/************************************************
 *  unit test of DFPT_Stern (C2)
 ***********************************************/

/**
 * - Tested Functions:
 *   - DFPT_Stern::solve - projected conjugate-gradient solution of the
 *     Sternheimer equation (H(k+q)-eps) P_c |dpsi> = -P_c |dV psi>.
 *   - apply_pv (implicitly) - the occupied-subspace projection.
 *
 * References are fully analytic:
 *   1. a diagonal (plane-wave-kinetic-like) operator with an exact
 *      closed-form complement solution;
 *   2. a dense Hermitian operator built from its known eigenbasis
 *      H = U D U^dagger (rotations x phases), the matrix analogue of a
 *      harmonic-oscillator Sternheimer problem in its eigenbasis, solved
 *      against the analytic spectral expansion;
 *   3. projection properties and degenerate right-hand sides.
 */

namespace {

unsigned g_seed = 20260814u;
double test_rand()
{
    g_seed = g_seed * 1664525u + 1013904223u;
    return ((g_seed >> 8) & 0xffffff) / 16777216.0 * 2.0 - 1.0;
}

std::complex<double> crand()
{
    return std::complex<double>(test_rand(), test_rand());
}

// Re <a|b>
double vdot(const std::vector<std::complex<double>>& a, const std::vector<std::complex<double>>& b)
{
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
    {
        s += a[i].real() * b[i].real() + a[i].imag() * b[i].imag();
    }
    return s;
}

// diagonal shifted operator: y_i = (d_i - eps) x_i
class DiagonalOperator : public ModuleDFPT::DFPT_Stern::LinearOperator
{
  public:
    DiagonalOperator(std::vector<double> d, double eps) : d_(std::move(d)), eps_(eps)
    {
    }
    int dimension() const override
    {
        return static_cast<int>(d_.size());
    }
    void apply(const std::complex<double>* x, std::complex<double>* y) const override
    {
        for (size_t i = 0; i < d_.size(); ++i)
        {
            y[i] = (d_[i] - eps_) * x[i];
        }
    }

  private:
    std::vector<double> d_;
    double eps_;
};

// dense Hermitian operator with a known eigenbasis: apply (H - eps I),
// H = U D U^dagger
class EigenbasisOperator : public ModuleDFPT::DFPT_Stern::LinearOperator
{
  public:
    EigenbasisOperator(const std::vector<std::vector<std::complex<double>>>& u,
                       const std::vector<double>& lambda,
                       double eps)
        : lambda_(lambda), eps_(eps)
    {
        const int n = static_cast<int>(lambda.size());
        h_.assign(n, std::vector<std::complex<double>>(n, std::complex<double>(0.0, 0.0)));
        for (int j = 0; j < n; ++j)
        {
            for (int i = 0; i < n; ++i)
            {
                for (int k = 0; k < n; ++k)
                {
                    // H(i,k) += lambda_j u(i,j) conj(u(k,j))
                    h_[i][k] += lambda_[j] * u[i][j] * std::conj(u[k][j]);
                }
            }
        }
    }
    int dimension() const override
    {
        return static_cast<int>(lambda_.size());
    }
    void apply(const std::complex<double>* x, std::complex<double>* y) const override
    {
        const int n = static_cast<int>(lambda_.size());
        for (int i = 0; i < n; ++i)
        {
            std::complex<double> s(0.0, 0.0);
            for (int k = 0; k < n; ++k)
            {
                s += h_[i][k] * x[k];
            }
            y[i] = s - eps_ * x[i];
        }
    }

  private:
    std::vector<std::vector<std::complex<double>>> h_;
    std::vector<double> lambda_;
    double eps_;
};

} // namespace

TEST(DFPTSternTest, SolvesDiagonalSystemExactlyOnTheComplement)
{
    const int n = 40;
    const int nocc = 5;
    std::vector<double> d(n);
    for (int i = 0; i < n; ++i)
    {
        d[i] = 1.0 + 0.37 * i;
    }
    const double eps = 2.0; // below every complement eigenvalue d_i - eps > 0
    DiagonalOperator aop(d, eps);

    std::vector<std::vector<std::complex<double>>> occ(nocc, std::vector<std::complex<double>>(n));
    for (int m = 0; m < nocc; ++m)
    {
        for (int i = 0; i < n; ++i)
        {
            occ[m][i] = (i == m) ? std::complex<double>(1.0, 0.0) : std::complex<double>(0.0, 0.0);
        }
    }
    std::vector<std::complex<double>> b(n);
    for (int i = 0; i < n; ++i)
    {
        b[i] = crand();
    }
    // analytic reference: x_i = (P_c b)_i / (d_i - eps), zero inside the occ space
    std::vector<std::complex<double>> ref = b;
    for (int m = 0; m < nocc; ++m)
    {
        for (int i = 0; i < n; ++i)
        {
            ref[i] -= occ[m][i] * std::conj(occ[m][i]) * b[i];
        }
    }
    for (int i = 0; i < nocc; ++i)
    {
        ref[i] = std::complex<double>(0.0, 0.0);
    }
    for (int i = nocc; i < n; ++i)
    {
        ref[i] /= (d[i] - eps);
    }

    ModuleDFPT::DFPT_Stern stern;
    std::vector<std::complex<double>> dpsi;
    double residual = 1.0;
    const int used = stern.solve(aop, occ, b, 200, 1.0e-10, dpsi, residual);
    EXPECT_GT(used, 0);
    EXPECT_LT(residual, 1.0e-10);
    ASSERT_EQ(dpsi.size(), static_cast<size_t>(n));
    double err2 = 0.0;
    double ref2 = 0.0;
    for (int i = 0; i < n; ++i)
    {
        err2 += std::norm(dpsi[i] - ref[i]);
        ref2 += std::norm(ref[i]);
    }
    EXPECT_LT(std::sqrt(err2 / ref2), 1.0e-8);
}

TEST(DFPTSternTest, SolvesDenseHermitianSystemAgainstSpectralReference)
{
    const int n = 24;
    const int nocc = 4;
    // unitary U = (Givens rotations) * diag(phases)
    std::vector<std::vector<std::complex<double>>> u(n, std::vector<std::complex<double>>(n));
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            u[i][j] = (i == j) ? std::complex<double>(1.0, 0.0) : std::complex<double>(0.0, 0.0);
        }
    }
    for (int j = 0; j < n; ++j)
    {
        const std::complex<double> ph = std::polar(1.0, 0.31 * j + 0.11);
        for (int i = 0; i < n; ++i)
        {
            u[i][j] *= ph;
        }
    }
    for (int j = 0; j + 1 < n; ++j)
    {
        // rotate the (j, j+1) plane of every column
        const double theta = 0.23 + 0.05 * j;
        const double c = std::cos(theta);
        const double s = std::sin(theta);
        for (int i = 0; i < n; ++i)
        {
            const std::complex<double> a = u[i][j];
            const std::complex<double> b = u[i][j + 1];
            u[i][j] = c * a - s * b;
            u[i][j + 1] = s * a + c * b;
        }
    }
    std::vector<double> lambda(n);
    for (int j = 0; j < n; ++j)
    {
        lambda[j] = 1.5 + 0.4 * j;
    }
    const double eps = 1.7; // inside the occ spectrum: indefinite/null directions
    EigenbasisOperator aop(u, lambda, eps);

    // occupied set = first nocc eigen-columns (machine-orthonormal)
    std::vector<std::vector<std::complex<double>>> occ(nocc, std::vector<std::complex<double>>(n));
    for (int m = 0; m < nocc; ++m)
    {
        for (int i = 0; i < n; ++i)
        {
            occ[m][i] = u[i][m];
        }
    }
    std::vector<std::complex<double>> b(n);
    for (int i = 0; i < n; ++i)
    {
        b[i] = crand();
    }
    // spectral reference: x = sum_{j >= nocc} u_j (u_j^dag b) / (lambda_j - eps)
    std::vector<std::complex<double>> ref(n, std::complex<double>(0.0, 0.0));
    for (int j = nocc; j < n; ++j)
    {
        std::complex<double> c(0.0, 0.0);
        for (int i = 0; i < n; ++i)
        {
            c += std::conj(u[i][j]) * b[i];
        }
        for (int i = 0; i < n; ++i)
        {
            ref[i] += u[i][j] * c / (lambda[j] - eps);
        }
    }

    ModuleDFPT::DFPT_Stern stern;
    std::vector<std::complex<double>> dpsi;
    double residual = 1.0;
    const int used = stern.solve(aop, occ, b, 500, 1.0e-11, dpsi, residual);
    EXPECT_GT(used, 0);
    EXPECT_LT(residual, 1.0e-10);
    double err2 = 0.0;
    double ref2 = 0.0;
    for (int i = 0; i < n; ++i)
    {
        err2 += std::norm(dpsi[i] - ref[i]);
        ref2 += std::norm(ref[i]);
    }
    EXPECT_LT(std::sqrt(err2 / ref2), 1.0e-7);
}

TEST(DFPTSternTest, SolutionStaysOrthogonalToOccupiedStates)
{
    const int n = 32;
    const int nocc = 6;
    std::vector<std::vector<std::complex<double>>> occ;
    for (int m = 0; m < nocc; ++m)
    {
        std::vector<std::complex<double>> v(n);
        for (int i = 0; i < n; ++i)
        {
            v[i] = crand();
        }
        // orthonormalize against the previous ones
        for (size_t k = 0; k < occ.size(); ++k)
        {
            std::complex<double> c(0.0, 0.0);
            for (int i = 0; i < n; ++i)
            {
                c += std::conj(occ[k][i]) * v[i];
            }
            for (int i = 0; i < n; ++i)
            {
                v[i] -= c * occ[k][i];
            }
        }
        const double nrm = std::sqrt(vdot(v, v));
        for (int i = 0; i < n; ++i)
        {
            v[i] /= nrm;
        }
        occ.push_back(v);
    }
    // random Hermitian positive definite operator H = A A^dag / n + I
    std::vector<std::vector<std::complex<double>>> a(n, std::vector<std::complex<double>>(n));
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            a[i][j] = crand() / std::sqrt(static_cast<double>(n));
        }
    }
    class PDOperator : public ModuleDFPT::DFPT_Stern::LinearOperator
    {
      public:
        explicit PDOperator(std::vector<std::vector<std::complex<double>>> a) : a_(std::move(a))
        {
            const int n = static_cast<int>(a_.size());
            h_.assign(n, std::vector<std::complex<double>>(n));
            for (int i = 0; i < n; ++i)
            {
                for (int k = 0; k < n; ++k)
                {
                    std::complex<double> s(0.0, 0.0);
                    for (int j = 0; j < n; ++j)
                    {
                        s += a_[i][j] * std::conj(a_[k][j]);
                    }
                    h_[i][k] = s + ((i == k) ? std::complex<double>(1.0, 0.0) : std::complex<double>(0.0, 0.0));
                }
            }
        }
        int dimension() const override
        {
            return static_cast<int>(a_.size());
        }
        void apply(const std::complex<double>* x, std::complex<double>* y) const override
        {
            const int n = static_cast<int>(a_.size());
            for (int i = 0; i < n; ++i)
            {
                std::complex<double> s(0.0, 0.0);
                for (int k = 0; k < n; ++k)
                {
                    s += h_[i][k] * x[k];
                }
                y[i] = s;
            }
        }

      private:
        std::vector<std::vector<std::complex<double>>> a_;
        std::vector<std::vector<std::complex<double>>> h_;
    };
    PDOperator aop(a);

    std::vector<std::complex<double>> b(n);
    for (int i = 0; i < n; ++i)
    {
        b[i] = crand();
    }
    ModuleDFPT::DFPT_Stern stern;
    std::vector<std::complex<double>> dpsi;
    double residual = 1.0;
    const int used = stern.solve(aop, occ, b, 500, 1.0e-10, dpsi, residual);
    EXPECT_GT(used, 0);
    EXPECT_LT(residual, 1.0e-9);
    ASSERT_EQ(dpsi.size(), static_cast<size_t>(n));
    for (int m = 0; m < nocc; ++m)
    {
        std::complex<double> c(0.0, 0.0);
        for (int i = 0; i < n; ++i)
        {
            c += std::conj(occ[m][i]) * dpsi[i];
        }
        EXPECT_LT(std::abs(c), 1.0e-9);
    }
}

TEST(DFPTSternTest, DegenerateRightHandSideInOccupiedSubspace)
{
    // b fully inside the occ subspace: dpsi must be exactly zero
    const int n = 16;
    const int nocc = 3;
    std::vector<std::vector<std::complex<double>>> occ(nocc, std::vector<std::complex<double>>(n));
    for (int m = 0; m < nocc; ++m)
    {
        for (int i = 0; i < n; ++i)
        {
            occ[m][i] = (i == m) ? std::complex<double>(1.0, 0.0) : std::complex<double>(0.0, 0.0);
        }
    }
    std::vector<double> d(n);
    for (int i = 0; i < n; ++i)
    {
        d[i] = 1.0 + i;
    }
    DiagonalOperator aop(d, 2.5);
    std::vector<std::complex<double>> b(n, std::complex<double>(0.0, 0.0));
    for (int m = 0; m < nocc; ++m)
    {
        for (int i = 0; i < n; ++i)
        {
            b[i] += (0.3 * m + 0.1) * occ[m][i];
        }
    }
    ModuleDFPT::DFPT_Stern stern;
    std::vector<std::complex<double>> dpsi;
    double residual = 1.0;
    const int used = stern.solve(aop, occ, b, 100, 1.0e-10, dpsi, residual);
    EXPECT_EQ(used, 0);
    EXPECT_EQ(residual, 0.0);
    ASSERT_EQ(dpsi.size(), static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
    {
        EXPECT_EQ(dpsi[i], std::complex<double>(0.0, 0.0));
    }
}

TEST(DFPTSternTest, ZeroRightHandSideReturnsImmediately)
{
    DiagonalOperator aop(std::vector<double>(8, 1.0), 0.5);
    std::vector<std::complex<double>> b(8, std::complex<double>(0.0, 0.0));
    ModuleDFPT::DFPT_Stern stern;
    std::vector<std::complex<double>> dpsi;
    double residual = 1.0;
    const int used = stern.solve(aop, std::vector<std::vector<std::complex<double>>>(), b, 50, 1.0e-8, dpsi, residual);
    EXPECT_EQ(used, 0);
    EXPECT_EQ(residual, 0.0);
    ASSERT_EQ(dpsi.size(), 8u);
    for (int i = 0; i < 8; ++i)
    {
        EXPECT_EQ(dpsi[i], std::complex<double>(0.0, 0.0));
    }
}
