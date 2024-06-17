#include "module_base/math_sphbes.h"
#include "module_base/matrix3.h"
#include "module_base/vector3.h"
#include "module_io/numerical_basis_jyjy.h"

#include "gtest/gtest.h"

using namespace NumericalBasis;
using ModuleBase::Sphbes;

/***********************************************************
 *      Unit test of class "TwoCenterIntegrator"
 ***********************************************************/
/**
 * Tested functions:
 *
 * - indexgen
 *   generates a vector of composite index (itype, iatom, l, m)
 *   from given number of atoms and maximum angular momentum
 *
 * - cal_overlap_Sq (static function)
 *   computes <jy|jy> and <jy|-\nabla^2|jy> with two-center integration
 *
 */
class NumericalBasisTest : public ::testing::Test
{
  protected:
    void SetUp();
    void TearDown();
};

void NumericalBasisTest::SetUp()
{
}

void NumericalBasisTest::TearDown()
{
}

TEST_F(NumericalBasisTest, indexgen)
{
    const std::vector<int> natom = {1, 2, 3};
    const std::vector<int> lmax = {2, 1, 2};

    const auto index = indexgen(natom, lmax);
    size_t idx = 0;

    for (size_t it = 0; it < natom.size(); ++it)
    {
        for (int ia = 0; ia != natom[it]; ++ia)
        {
            for (int l = 0; l <= lmax[it]; ++l)
            {
                for (int M = 0; M < 2 * l + 1; ++M)
                {
                    // convert the "abacus M" to the conventional m
                    const int m = (M % 2 == 0) ? -M / 2 : (M + 1) / 2;
                    EXPECT_EQ(index[idx], std::make_tuple(int(it), ia, l, m));
                    ++idx;
                }
            }
        }
    }
}

TEST_F(NumericalBasisTest, cal_overlap_Sq)
{
    const int lmax = 3;
    const int nbes = 7;
    const double rcut = 7.0;

    const ModuleBase::Matrix3 latvec(30.0, 0.0, 0.0, 0.0, 30.0, 0.0, 0.0, 0.0, 30.0);

    // atom-0 and 1 have overlap with each other, but neither of them
    // has overlap with atom-2
    std::vector<std::vector<ModuleBase::Vector3<double>>> tau_cart;
    tau_cart.emplace_back();
    tau_cart[0].emplace_back(0.0, 0.0, 0.0);
    tau_cart[0].emplace_back(1.0, 1.0, 1.0);
    tau_cart[0].emplace_back(0.0, 0.0, 15.0);

    const std::vector<int> natom = {int(tau_cart[0].size())};
    const std::vector<int> lmax_v = {lmax};
    const auto index = indexgen(natom, lmax_v);
    const int nao = index.size();

    // zeros of sperical Bessel functions
    std::vector<double> zeros(nbes * (lmax + 1));
    ModuleBase::Sphbes::sphbes_zeros(lmax, nbes, zeros.data(), true);

    // <jy|jy> and <jy|-\nabla^2|jy>
    const double S_thr = 1e-5;
    const double T_thr = 1e-3;
    const auto S = cal_overlap_Sq('S', lmax, nbes, rcut, tau_cart, latvec, index);
    const auto T = cal_overlap_Sq('T', lmax, nbes, rcut, tau_cart, latvec, index);
    int t1, a1, l1, m1, t2, a2, l2, m2; // values will be updated in the loop
    for (int i = 0; i < nao; ++i)
    {
        std::tie(t1, a1, l1, m1) = index[i];
        for (int j = 0; j < nao; ++j)
        {
            std::tie(t2, a2, l2, m2) = index[j];
            for (int q1 = 0; q1 < nbes; ++q1)
            {
                for (int q2 = 0; q2 < nbes; ++q2)
                {
                    if (i == j)
                    {
                        if (q1 == q2)
                        {
                            // diagonal elements have analytical results
                            double S_ref
                                = std::pow(rcut, 3) * 0.5 * std::pow(Sphbes::sphbesj(l1 + 1, zeros[l1 * nbes + q1]), 2);
                            EXPECT_NEAR(S(i, i, q1, q1).real(), S_ref, S_thr);
                            EXPECT_EQ(S(i, i, q1, q1).imag(), 0);

                            double T_ref
                                = 0.5 * rcut
                                  * std::pow(zeros[l1 * nbes + q1] * Sphbes::sphbesj(l1 + 1, zeros[l1 * nbes + q1]), 2);
                            EXPECT_NEAR(T(i, i, q1, q1).real(), T_ref, T_thr);
                            EXPECT_EQ(T(i, i, q1, q1).imag(), 0);
                        }
                        else
                        {
                            // off-diagonal elements should be zero due to orthogonality
                            EXPECT_NEAR(S(i, i, q1, q2).real(), 0, S_thr);
                            EXPECT_EQ(S(i, i, q1, q2).imag(), 0);

                            EXPECT_NEAR(T(i, i, q1, q2).real(), 0, T_thr);
                            EXPECT_EQ(T(i, i, q1, q2).imag(), 0);
                        }
                    }

                    if ((a1 == 2) != (a2 == 2))
                    {
                        // atom-2 has no overlap with atom-0/1
                        EXPECT_EQ(S(i, j, q1, q2).real(), 0);
                        EXPECT_EQ(S(i, j, q1, q2).imag(), 0);

                        EXPECT_EQ(T(i, j, q1, q2).real(), 0);
                        EXPECT_EQ(T(i, j, q1, q2).imag(), 0);
                    }

                    if (a1 != a2 && a1 != 2 && a2 != 2)
                    {
                        // overlap between atom-0 and atom-1 orbitals should be non-zero
                        EXPECT_NE(S(i, j, q1, q2).real(), 0);
                        EXPECT_EQ(S(i, j, q1, q2).imag(), 0);

                        EXPECT_NE(T(i, j, q1, q2).real(), 0);
                        EXPECT_EQ(T(i, j, q1, q2).imag(), 0);
                    }
                }
            }
        }
    }

    // verifies that the calculated matrix is invariant w.r.t. translation
    // across the cell boundary
    tau_cart[0][0].x -= 0.5;
    tau_cart[0][1].x -= 0.5;
    tau_cart[0][2].x -= 0.5;
    const auto S_wrap = cal_overlap_Sq('S', lmax, nbes, rcut, tau_cart, latvec, index);

    for (int i = 0; i < nao; ++i)
    {
        std::tie(t1, a1, l1, m1) = index[i];
        for (int j = 0; j < nao; ++j)
        {
            std::tie(t2, a2, l2, m2) = index[j];
            for (int q1 = 0; q1 < nbes; ++q1)
            {
                for (int q2 = 0; q2 < nbes; ++q2)
                {
                    EXPECT_EQ(S(i, j, q1, q2), S_wrap(i, j, q1, q2));
                }
            }
        }
    }
}

TEST_F(NumericalBasisTest, neighbor_vec)
{
    const ModuleBase::Vector3<double> d0(2.5, 0.0, 0.0);
    const ModuleBase::Matrix3 latvec(3.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 5.0);

    const double r1 = 2.0;
    auto d1 = neighbor_vec(d0, latvec, r1);
    EXPECT_EQ(d1.size(), 1);
    EXPECT_EQ(d1[0], ModuleBase::Vector3<double>(-0.5, 0.0, 0.0));

    const double r2 = 3.6;
    auto d2 = neighbor_vec(d0, latvec, r2);
    EXPECT_EQ(d2.size(), 3);
    EXPECT_EQ(d2[0], ModuleBase::Vector3<double>(-3.5, 0.0, 0.0));
    EXPECT_EQ(d2[1], ModuleBase::Vector3<double>(-0.5, 0.0, 0.0));
    EXPECT_EQ(d2[2], ModuleBase::Vector3<double>(2.5, 0.0, 0.0));

    const double r3 = 4.1;
    auto d3 = neighbor_vec(d0, latvec, r3);
    EXPECT_EQ(d3.size(), 5);
    EXPECT_EQ(d3[0], ModuleBase::Vector3<double>(-3.5, 0.0, 0.0));
    EXPECT_EQ(d3[1], ModuleBase::Vector3<double>(-0.5, -4.0, 0.0));
    EXPECT_EQ(d3[2], ModuleBase::Vector3<double>(-0.5, 0.0, 0.0));
    EXPECT_EQ(d3[3], ModuleBase::Vector3<double>(-0.5, 4.0, 0.0));
    EXPECT_EQ(d3[4], ModuleBase::Vector3<double>(2.5, 0.0, 0.0));
}
