#include "source_pw/module_proj/projgen.h"

#include "source_base/math_integral.h"

#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace
{

double radial_norm(const std::vector<double>& r, const std::vector<double>& radial)
{
    std::vector<double> dr(radial.size());
    std::vector<double> integrand(radial.size());
    std::adjacent_difference(r.begin(), r.begin() + radial.size(), dr.begin());
    for (std::size_t i = 0; i < radial.size(); ++i)
    {
        integrand[i] = r[i] * r[i] * radial[i] * radial[i];
    }
    return ModuleBase::Integral::simpson(radial.size(), integrand.data(), &dr[1]);
}

} // namespace

TEST(ProjgenTest, GeneratesNormalizedTruncatedProjector)
{
    const int nr = 201;
    const double dr = 0.05;
    const double rcut = 5.0;
    std::vector<double> r(nr);
    std::vector<double> chi(nr);
    for (int i = 0; i < nr; ++i)
    {
        r[i] = i * dr;
        chi[i] = std::exp(-r[i]);
    }

    std::vector<double> alpha;
    projgen(0, nr, r.data(), chi.data(), rcut, 4, alpha);

    ASSERT_EQ(alpha.size(), 101U);
    EXPECT_TRUE(std::all_of(alpha.begin(), alpha.end(), [](const double value) { return std::isfinite(value); }));
    EXPECT_NEAR(radial_norm(r, alpha), 1.0, 1.0e-10);
}

TEST(ProjgenTest, MatchesReferenceProfile)
{
    // normalized test orbital r^2 * exp(-r), l = 2
    const int l = 2;
    const double dr = 0.01;
    const double rcut_nao = 10.0;
    const int nr = static_cast<int>(rcut_nao / dr) + 1;
    std::vector<double> r(nr);
    std::vector<double> chi(nr);
    for (int i = 0; i < nr; ++i)
    {
        r[i] = i * dr;
        chi[i] = r[i] * r[i] * std::exp(-r[i]);
    }
    std::vector<double> integrand(nr);
    for (int i = 0; i < nr; ++i)
    {
        integrand[i] = std::pow(r[i] * chi[i], 2);
    }
    const double norm = 1.0 / std::sqrt(ModuleBase::Integral::simpson(nr, integrand.data(), dr));
    for (int i = 0; i < nr; ++i)
    {
        chi[i] *= norm;
    }

    const double rcut_proj = 7.0;
    const int nbes = 7;
    std::vector<double> alpha;
    projgen(l, nr, r.data(), chi.data(), rcut_proj, nbes, alpha);

    // reference values from the python prototype
    const std::vector<double> ref = {
        0.000000000000e+00, 2.344902364599e-05, 9.378381332712e-05, 2.109675345121e-04, 3.749388271050e-04,
        5.856118515995e-04, 8.428763536364e-04, 1.146597746904e-03, 1.496617214310e-03, 1.892751827321e-03,
        2.334794683381e-03, 2.822515061259e-03, 3.355658594204e-03, 3.933947460740e-03, 4.557080592928e-03,
        5.224733901903e-03, 5.936560520491e-03, 6.692191062668e-03, 7.491233899644e-03, 8.333275452302e-03,
    };

    ASSERT_GE(alpha.size(), ref.size());
    for (std::size_t i = 0; i < ref.size(); ++i)
    {
        EXPECT_NEAR(alpha[i], ref[i], 1.0e-12);
    }
}

