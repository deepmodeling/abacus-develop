#include "gtest/gtest.h"
#include "source_basis/module_ao/orb_atomic_lm.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

TEST(NumericalOrbitalLmUniformGrid, DerivativeMatchesInterpolatedValues)
{
    const int nr = 101;
    const double dr = 0.02;
    const double fine_dr = 0.001;
    const double cutoff = 2.0;
    std::vector<double> radial(nr);
    std::vector<double> rab(nr, dr);
    std::vector<double> values(nr);
    for (int ir = 0; ir < nr; ++ir)
    {
        const double r = ir * dr;
        const double remaining = cutoff - r;
        const double envelope = std::exp(-r);
        radial[ir] = r;
        values[ir] = r * r * remaining * remaining * envelope;
    }
    const double* radial_data = radial.data();
    const double* rab_data = rab.data();
    const double* values_data = values.data();
    double max_error_all = 0.0;
    // Cover every explicit spline boundary-condition branch and the default.
    for (int angular_momentum = 0; angular_momentum <= 5; ++angular_momentum)
    {
        SCOPED_TRACE(angular_momentum);
        Numerical_Orbital_Lm orbital;
        orbital.set_orbital_info("X", 0, angular_momentum, 0, nr, rab_data, radial_data,
                                Numerical_Orbital_Lm::Psi_Type::Psi, values_data,
                                65, 0.1, fine_dr, false, true, false);
        const double* f = orbital.getPsiuniform();
        const double* df = orbital.getDpsiuniform();
        double max_error = 0.0;
        // All five samples lie inside one original spline interval. This
        // fourth-order centered difference differentiates a cubic exactly.
        for (int ir = 10; ir < 1990; ir += 20)
        {
            const double fd = (f[ir - 2] - 8.0 * f[ir - 1]
                               + 8.0 * f[ir + 1] - f[ir + 2]) / (12.0 * fine_dr);
            const double error = std::abs(fd - df[ir]);
            const bool finite_error = std::isfinite(error);
            ASSERT_TRUE(finite_error);
            max_error = std::max(max_error, error);
        }
        EXPECT_LT(max_error, 5e-10);
        max_error_all = std::max(max_error_all, max_error);
        int padding_points = 0;
        const int nr_uniform = orbital.getNruniform();
        for (int ir = 0; ir < nr_uniform; ++ir)
        {
            const double r = ir * fine_dr;
            if (r >= cutoff)
            {
                EXPECT_DOUBLE_EQ(f[ir], 0.0);
                ++padding_points;
            }
        }
        EXPECT_GT(padding_points, 0);
    }
    const std::string error_text = ::testing::PrintToString(max_error_all);
    RecordProperty("max_derivative_error", error_text);
}
