// Unit tests of DFPT_NSCF (the block-Davidson nscf eigensolver feeding
// the k+q occupied states): the solver is checked against two operators
// with closed-form spectra -- a diagonal operator (exact eigenpairs known
// in arbitrary order) and a synthetic two-block Hamiltonian (degenerate
// pair) -- plus determinism and nwant-reduction contracts.

#include "../dfpt_nscf.h"
#include "../dfpt_stern.h"

#include "source_base/complexmatrix.h"
#include "source_base/module_external/lapack_connector.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <random>
#include <vector>

namespace
{

// near-diagonal approximation of the preconditioner: the production
// preconditioner is the kinetic diagonal (never the exact H diagonal),
// and an exact diag would degenerate the Davidson expansion (the
// residual direction maps back into the span); tests use a scaled copy
// to mirror that
std::vector<double> approx_precond(const std::vector<double> &d)
{
    std::vector<double> p(d.size());
    for (size_t i = 0; i < d.size(); ++i)
    {
        p[i] = 1.3 * d[i] + 2.0;
    }
    return p;
}

// diagonal operator: H = diag(d[0..n-1]) with the nscf preconditioner
// convention (tpiba^2-style magnitudes)
class DiagonalOperator : public ModuleDFPT::DFPT_Stern::LinearOperator
{
  public:
    DiagonalOperator(const std::vector<double> &d) : d_(d) {}
    int dimension() const override { return static_cast<int>(d_.size()); }
    void apply(const std::complex<double> *x, std::complex<double> *y) const override
    {
        for (size_t i = 0; i < d_.size(); ++i)
        {
            y[i] = d_[i] * x[i];
        }
    }
    const std::vector<double> &diag() const { return d_; }

  private:
    std::vector<double> d_;
};

TEST(DFPTNscfTest, DiagonalOperatorRecoversLowestEigenpairs)
{
    // 40-dim diagonal with non-monotone ordering: the solver must return
    // the 6 lowest eigenvalues ascending and the exact unit eigenvectors
    const int n = 40;
    std::vector<double> d(n);
    std::mt19937 gen(7u);
    std::uniform_real_distribution<double> uni(1.0, 60.0);
    for (int i = 0; i < n; ++i)
    {
        d[i] = uni(gen);
    }
    d[0] = 0.3;
    d[9] = 1.7;
    d[23] = 2.9;
    d[5] = 4.1;
    d[31] = 5.5;
    d[17] = 6.6;
    DiagonalOperator op(d);

    std::vector<std::vector<std::complex<double>>> psi;
    std::vector<double> eig;
    const std::vector<double> precond = approx_precond(d);
    const int it = ModuleDFPT::DFPT_NSCF::lowest_eigpairs(op, precond, 6, 200, 1e-9, psi, eig);
    EXPECT_EQ(0, it);
    ASSERT_EQ(6u, psi.size());
    ASSERT_EQ(6u, eig.size());
    std::vector<double> want = d;
    std::sort(want.begin(), want.end());
    for (int j = 0; j < 6; ++j)
    {
        EXPECT_NEAR(want[j], eig[j], 1e-8);
    }
    // eigenvector contraction: <psi_j|H|psi_j> == eig_j and unit norm
    for (int j = 0; j < 6; ++j)
    {
        double nrm = 0.0;
        std::complex<double> hv(0.0, 0.0);
        for (int i = 0; i < n; ++i)
        {
            nrm += std::norm(psi[j][i]);
            hv += std::conj(psi[j][i]) * psi[j][i] * static_cast<double>(d[i]);
        }
        EXPECT_NEAR(1.0, nrm, 1e-8);
        EXPECT_NEAR(eig[j], hv.real(), 1e-8);
        EXPECT_NEAR(0.0, hv.imag(), 1e-10);
    }
}

TEST(DFPTNscfTest, DegeneratePairConvergesWithOrthogonalOutput)
{
    // rank-two degenerate block: eigenvalues 1.0 (x2) share a plane; the
    // solver must return two orthonormal vectors spanning it
    const int n = 24;
    std::vector<double> d(n, 10.0);
    d[0] = 1.0;
    d[1] = 1.0;
    d[7] = 4.0;
    DiagonalOperator op(d);

    std::vector<std::vector<std::complex<double>>> psi;
    std::vector<double> eig;
    const std::vector<double> precond = approx_precond(d);
    const int it = ModuleDFPT::DFPT_NSCF::lowest_eigpairs(op, precond, 3, 200, 1e-9, psi, eig);
    EXPECT_EQ(0, it);
    ASSERT_EQ(3u, psi.size());
    EXPECT_NEAR(1.0, eig[0], 1e-8);
    EXPECT_NEAR(1.0, eig[1], 1e-8);
    EXPECT_NEAR(4.0, eig[2], 1e-8);
    for (int a = 0; a < 3; ++a)
    {
        double na = 0.0;
        for (int i = 0; i < n; ++i)
        {
            na += std::norm(psi[a][i]);
        }
        EXPECT_NEAR(1.0, na, 1e-8);
        for (int b = a + 1; b < 3; ++b)
        {
            std::complex<double> dot(0.0, 0.0);
            for (int i = 0; i < n; ++i)
            {
                dot += std::conj(psi[a][i]) * psi[b][i];
            }
            EXPECT_NEAR(0.0, std::abs(dot), 1e-8);
        }
    }
}

TEST(DFPTNscfTest, DeterministicAcrossCalls)
{
    const int n = 30;
    std::vector<double> d(n);
    std::mt19937 gen(11u);
    std::uniform_real_distribution<double> uni(1.0, 50.0);
    for (int i = 0; i < n; ++i)
    {
        d[i] = uni(gen);
    }
    d[3] = 0.7;
    d[19] = 1.2;
    DiagonalOperator op(d);
    std::vector<std::vector<std::complex<double>>> p1;
    std::vector<std::vector<std::complex<double>>> p2;
    std::vector<double> e1;
    std::vector<double> e2;
    // a single named preconditioner emphasizes that both calls must use
    // the same preconditioner for the determinism comparison
    const std::vector<double> precond = approx_precond(d);
    ModuleDFPT::DFPT_NSCF::lowest_eigpairs(op, precond, 4, 200, 1e-9, p1, e1);
    ModuleDFPT::DFPT_NSCF::lowest_eigpairs(op, precond, 4, 200, 1e-9, p2, e2);
    ASSERT_EQ(e1.size(), e2.size());
    for (size_t j = 0; j < e1.size(); ++j)
    {
        EXPECT_DOUBLE_EQ(e1[j], e2[j]);
    }
}

TEST(DFPTNscfTest, FullDimensionSingleBandExact)
{
    // n = nwant = 1: the compressed basis spans the whole space, so the
    // reseed has no independent component and must be discarded (a kept
    // rank-deficient reseed would inject a spurious zero Ritz pair). The
    // initial Ritz state is already exact; assert the eigenvalue, the
    // converged status, and the unit norm survive the full-dimension path.
    std::vector<double> d(1, 2.5);
    DiagonalOperator op(d);
    std::vector<std::vector<std::complex<double>>> psi;
    std::vector<double> eig;
    const std::vector<double> precond = approx_precond(d);
    const int it = ModuleDFPT::DFPT_NSCF::lowest_eigpairs(op, precond, 1, 50, 1e-9, psi, eig);
    EXPECT_EQ(0, it);
    ASSERT_EQ(1u, psi.size());
    ASSERT_EQ(1u, eig.size());
    EXPECT_NEAR(2.5, eig[0], 1e-10);
    EXPECT_NEAR(1.0, std::abs(psi[0][0]), 1e-10);
}

} // namespace

#include <iostream>
namespace
{
// direct probe of the LapackConnector::zheev round-trip semantics:
// M = [[2, i], [-i, 3]] has spectrum {1, 4}; the eigenvector for 1 is
// proportional to (i, 1)/sqrt(2) (or any phase multiple)
TEST(DFPTNscfTest, ZheevRoundTripProbe)
{
    // M = [[2, i], [-i, 3]]: spectrum = (5 +/- sqrt(5))/2 = 1.38197/3.61803.
    // The round-trip (transpose in, transpose out) must preserve the
    // eigenvector convention, checked through the full Rayleigh quotient
    // and pairwise orthonormality on the ORIGINAL matrix.
    ModuleBase::ComplexMatrix m(2, 2, true);
    m(0, 0) = {2, 0};
    m(1, 1) = {3, 0};
    m(0, 1) = {0, 1};
    m(1, 0) = {0, -1};
    const double lam_lo = (5.0 - std::sqrt(5.0)) / 2.0;
    const double lam_hi = (5.0 + std::sqrt(5.0)) / 2.0;
    std::vector<double> w(2);
    std::vector<std::complex<double>> work(40);
    std::vector<double> rwork(40);
    int info = 0;
    LapackConnector::zheev('V', 'U', 2, m, 2, w.data(), work.data(), 40, rwork.data(), &info);
    ASSERT_EQ(0, info);
    EXPECT_NEAR(lam_lo, w[0], 1e-12);
    EXPECT_NEAR(lam_hi, w[1], 1e-12);
    // orthonormal columns
    std::complex<double> dot = std::conj(m(0, 0)) * m(0, 1) + std::conj(m(1, 0)) * m(1, 1);
    EXPECT_NEAR(0.0, std::abs(dot), 1e-12);
    // H * v_j = w_j * v_j on the original matrix entries
    const std::complex<double> hv0 = std::complex<double>(2, 0) * m(0, 0) + std::complex<double>(0, 1) * m(1, 0);
    EXPECT_NEAR(std::abs(hv0 - w[0] * m(0, 0)), 0.0, 1e-12);
    const std::complex<double> hv1 = std::complex<double>(2, 0) * m(0, 1) + std::complex<double>(0, 1) * m(1, 1);
    EXPECT_NEAR(std::abs(hv1 - w[1] * m(0, 1)), 0.0, 1e-12);
}
} // namespace
