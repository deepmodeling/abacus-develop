#include "source_basis/module_nao/smoothgen.h"

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

TEST(SmoothgenTest, SmoothsAndNormalizesAtCutoff)
{
    const int nr = 121;
    const double dr = 0.05;
    const double rcut = 4.0;
    std::vector<double> r(nr);
    std::vector<double> chi(nr);
    for (int i = 0; i < nr; ++i)
    {
        r[i] = i * dr;
        chi[i] = std::exp(-0.5 * r[i] * r[i]);
    }

    std::vector<double> alpha;
    smoothgen(nr, r.data(), chi.data(), rcut, alpha);

    ASSERT_EQ(alpha.size(), 81U);
    EXPECT_TRUE(std::all_of(alpha.begin(), alpha.end(), [](const double value) { return std::isfinite(value); }));
    EXPECT_NEAR(alpha.back(), 0.0, 1.0e-14);
    EXPECT_NEAR(radial_norm(r, alpha), 1.0, 1.0e-10);
}

TEST(SmoothgenTest, MatchesReferenceProfile)
{
    // normalized test orbital r^2 * exp(-r)
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
    std::vector<double> alpha;
    smoothgen(nr, r.data(), chi.data(), rcut_proj, alpha);

    // reference values from the python prototype
    const std::vector<double> ref = {
        0.0,
        4.3350439973614511e-05,
        0.00017167638355532129,
        0.00038242839374460959,
        0.0006731078535961104,
        0.0010412661227313883,
        0.0014845037064463464,
        0.0020004694372387291,
        0.0025868596685825239,
        0.0032414174807783905,
        0.0039619318987114734,
        0.0047462371213502306,
        0.0055922117628221264,
        0.006497778104904167,
        0.0074609013607684835,
        0.0084795889498252546,
        0.0095518897835073727,
        0.010675893561843302,
        0.01184973008066669,
        0.013071568549313281,
    };

    ASSERT_GE(alpha.size(), ref.size());
    for (std::size_t i = 0; i < ref.size(); ++i)
    {
        EXPECT_NEAR(alpha[i], ref[i], 1.0e-12);
    }
}
