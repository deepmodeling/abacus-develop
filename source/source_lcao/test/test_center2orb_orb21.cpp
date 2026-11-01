#include "source_lcao/center2orb_orb11.h"
#include "source_lcao/center2orb_orb21.h"

#include "gtest/gtest.h"
#include <cmath>
#include <vector>

namespace
{
Numerical_Orbital_Lm make_orbital(
    const int l,
    const int index,
    const int nr,
    const int nk)
{
    const double dr = 0.02;
    const double cutoff = (nr - 1) * dr;
    std::vector<double> radial(nr);
    std::vector<double> rab(nr, dr);
    std::vector<double> values(nr);
    for (int ir = 0; ir < nr; ++ir)
    {
        const double r = ir * dr;
        radial[ir] = r;
        values[ir] = std::pow(r, l) * std::pow(cutoff - r, 2) * std::exp(-(0.7 + 0.03 * index) * r);
    }
    Numerical_Orbital_Lm orbital;
    const double *rab_ptr = rab.data();
    const double *radial_ptr = radial.data();
    const double *values_ptr = values.data();
    orbital.set_orbital_info(
        "X",
        0,
        l,
        index,
        nr,
        rab_ptr,
        radial_ptr,
        Numerical_Orbital_Lm::Psi_Type::Psi,
        values_ptr,
        nk,
        0.1,
        0.01,
        false,
        true,
        true);
    return orbital;
}

} // namespace

TEST(Center2Orb11, AllMagneticComponentsMatchScalarEvaluation)
{
    ORB_gaunt_table gaunt;
    gaunt.init_Gaunt_CH(6);
    gaunt.init_Gaunt(6);
    ModuleBase::Sph_Bessel_Recursive::D2 bessel;
    bessel.set_dx(0.001);
    bessel.cal_jlx(10, 405, 65);
    const ModuleBase::Vector3<double> origin(0.0, 0.0, 0.0);
    for (int la : {0,
             1,
             3,
             6})
    {
        for (int lb : {0, 2, 3})
        {
            const auto left = make_orbital(
                la,
                0,
                101,
                65);
            const auto right = make_orbital(
                lb,
                1,
                81,
                65);
            Center2_Orb::Orb11 orbital(
                left,
                right,
                &bessel,
                gaunt);
            orbital.init_radial_table();
            for (double distance : {0.0,
                     1e-12,
                     5e-11,
                     2e-10,
                     0.037,
                     0.73,
                     1.91,
                     3.57,
                     3.8})
            {
                const double x = 0.8 * distance;
                const double y = 0.36 * distance;
                const double z = 0.48 * distance;
                const ModuleBase::Vector3<double> position(x, y, z);
                const auto values = orbital.cal_overlap(origin, position);
                const auto gradients = orbital.cal_grad_overlap(origin, position);
                const int count = (2 * la + 1) * (2 * lb + 1);
                ASSERT_EQ(values.size(), count);
                ASSERT_EQ(gradients.size(), count);
                for (int ma = 0; ma <= 2 * la; ++ma)
                {
                    for (int mb = 0; mb <= 2 * lb; ++mb)
                    {
                        const int offset = ma * (2 * lb + 1) + mb;
                        const double expected = orbital.cal_overlap(
                            origin,
                            position,
                            ma,
                            mb);
                        const double tolerance = 2e-11 * (1.0 + std::abs(expected));
                        EXPECT_NEAR(values[offset], expected, tolerance);
                        const auto gradient = orbital.cal_grad_overlap(
                            origin,
                            position,
                            ma,
                            mb);
                        const double error = (gradients[offset] - gradient).norm();
                        const double gradient_tolerance = 2e-10 * (1.0 + gradient.norm());
                        EXPECT_LT(error, gradient_tolerance);
                    }
                }
            }
        }
    }
}

TEST(Center2Orb21, AllMagneticComponentsMatchScalarEvaluation)
{
    // The l=6 auxiliary channel and f orbital produce channels through l=9.
    const auto left = make_orbital(
        6,
        0,
        101,
        65);
    const auto middle = make_orbital(
        3,
        1,
        101,
        65);
    const auto right = make_orbital(
        3,
        2,
        81,
        65);
    ORB_gaunt_table gaunt;
    gaunt.init_Gaunt_CH(9);
    gaunt.init_Gaunt(9);
    ModuleBase::Sph_Bessel_Recursive::D2 bessel;
    bessel.set_dx(0.001);
    bessel.cal_jlx(13, 405, 65);
    Center2_Orb::Orb21 orbital(
        left,
        middle,
        right,
        &bessel,
        gaunt);
    orbital.init_radial_table();
    const ModuleBase::Vector3<double> origin(0.0, 0.0, 0.0);
    for (double distance : {0.0,
             1e-12,
             5e-11,
             2e-10,
             0.037,
             0.73,
             1.91,
             3.57,
             3.8})
    {
        const double x = 0.8 * distance;
        const double y = 0.36 * distance;
        const double z = 0.48 * distance;
        const ModuleBase::Vector3<double> position(x, y, z);
        const auto values = orbital.cal_overlap(origin, position);
        const auto gradients = orbital.cal_grad_overlap(origin, position);
        ASSERT_EQ(values.size(), 13 * 7 * 7);
        ASSERT_EQ(gradients.size(), values.size());
        for (int ma = 0; ma < 13; ++ma)
        {
            for (int mb = 0; mb < 7; ++mb)
            {
                for (int mc = 0; mc < 7; ++mc)
                {
                    const int offset = (ma * 7 + mb) * 7 + mc;
                    const double expected = orbital.cal_overlap(
                        origin,
                        position,
                        ma,
                        mb,
                        mc);
                    const double tolerance = 2e-11 * (1.0 + std::abs(expected));
                    EXPECT_NEAR(values[offset], expected, tolerance);
                    const auto gradient = orbital.cal_grad_overlap(
                        origin,
                        position,
                        ma,
                        mb,
                        mc);
                    const double error = (gradients[offset] - gradient).norm();
                    const double gradient_tolerance = 2e-10 * (1.0 + gradient.norm());
                    EXPECT_LT(error, gradient_tolerance);
                }
            }
        }
    }
}

TEST(Center2Orb21, ProductGradientMatchesFiniteDifference)
{
    const auto left = make_orbital(
        2,
        0,
        101,
        65);
    const auto middle = make_orbital(
        1,
        1,
        101,
        65);
    const auto right = make_orbital(
        3,
        2,
        101,
        65);
    ORB_gaunt_table gaunt;
    gaunt.init_Gaunt_CH(3);
    gaunt.init_Gaunt(3);
    ModuleBase::Sph_Bessel_Recursive::D2 bessel;
    bessel.set_dx(0.001);
    bessel.cal_jlx(7, 405, 65);
    Center2_Orb::Orb21 orbital(
        left,
        middle,
        right,
        &bessel,
        gaunt);
    orbital.init_radial_table();
    const ModuleBase::Vector3<double> origin(0.0, 0.0, 0.0);
    const ModuleBase::Vector3<double> position(0.813, 0.321, 0.473);
    const double step = 1e-5;
    const auto gradients = orbital.cal_grad_overlap(origin, position);
    for (int axis = 0; axis < 3; ++axis)
    {
        auto plus = position;
        auto minus = position;
        plus[axis] += step;
        minus[axis] -= step;
        const auto upper = orbital.cal_overlap(origin, plus);
        const auto lower = orbital.cal_overlap(origin, minus);
        ASSERT_EQ(upper.size(), gradients.size());
        ASSERT_EQ(lower.size(), gradients.size());
        for (std::size_t i = 0; i < gradients.size(); ++i)
        {
            const double finite_difference = (upper[i] - lower[i]) / (2 * step);
            const double tolerance = 2e-7 * (1.0 + std::abs(finite_difference));
            EXPECT_NEAR(gradients[i][axis], finite_difference, tolerance);
        }
    }
}
