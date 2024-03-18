#include "module_basis/module_nao/projgen.h"
#include "gtest/gtest.h"

#include "module_base/math_integral.h"
#include "module_base/math_sphbes.h"
#include "module_base/cubic_spline.h"

#include <numeric>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <cstdio>

/***********************************************************
 *      Unit test of functions in projgen.cpp
 *     "projgen" : generate the projector coefficients
 *     "smoothgen"  : smooth the projector coefficients and optimize the sigma
 * 
 ***********************************************************/

TEST(projgen_test, projgen)
{
    // test orbital r^2 * exp(-r)
    int l = 2;
    double dr = 0.01;
    double rcut_nao = 10;
    int nr_nao = int(rcut_nao / dr) + 1;
    std::vector<double> r(nr_nao);
    std::vector<double> orb(nr_nao);

    for (int i = 0; i < nr_nao; ++i) {
        r[i] = i * dr;
        orb[i] = r[i] * r[i] * std::exp(-r[i]);
    }

    // normalize the input orbital
    std::vector<double> integrand(nr_nao);
    std::transform(r.begin(), r.end(), orb.begin(), integrand.begin(),
            [](double r_i, double orb_i) { return std::pow(r_i * orb_i, 2); });
    double N = 1.0 / std::sqrt(ModuleBase::Integral::simpson(nr_nao, integrand.data(), dr));
    std::for_each(orb.begin(), orb.end(), [N](double& chi_i) { chi_i *= N; });

    // projector information
    double rcut_proj = 7.0;
    int nbes = 7;
    std::vector<double> alpha;

    projgen(l, nr_nao, r.data(), orb.data(), rcut_proj, nbes, alpha);

    // compare with python script result
    std::vector<double> ref = {
        0.000000000000e+00, 
        2.344902364599e-05,
        9.378381332712e-05,
        2.109675345121e-04,
        3.749388271050e-04,
        5.856118515995e-04,
        8.428763536364e-04,
        1.146597746904e-03,
        1.496617214310e-03,
        1.892751827321e-03,
        2.334794683381e-03,
        2.822515061259e-03,
        3.355658594204e-03,
        3.933947460740e-03,
        4.557080592928e-03,
        5.224733901903e-03,
        5.936560520491e-03,
        6.692191062668e-03,
        7.491233899644e-03,
        8.333275452302e-03,
    };

    for (int i = 0; i < 20; ++i) {
        EXPECT_NEAR(alpha[i], ref[i], 1e-12);
    }
}

TEST(smoothgen_test, smoothgen)
{
    // test orbital r^2 * exp(-r)
    int l = 2;
    double dr = 0.01;
    double rcut_nao = 10;
    int nr_nao = int(rcut_nao / dr) + 1;
    std::vector<double> r(nr_nao);
    std::vector<double> orb(nr_nao);

    for (int i = 0; i < nr_nao; ++i) {
        r[i] = i * dr;
        orb[i] = r[i] * r[i] * std::exp(-r[i]);
    }

    // normalize the input orbital
    std::vector<double> integrand(nr_nao);
    std::transform(r.begin(), r.end(), orb.begin(), integrand.begin(),
            [](double r_i, double orb_i) { return std::pow(r_i * orb_i, 2); });
    double N = 1.0 / std::sqrt(ModuleBase::Integral::simpson(nr_nao, integrand.data(), dr));
    std::for_each(orb.begin(), orb.end(), [N](double& chi_i) { chi_i *= N; });

    // projector information
    double rcut_proj = 7.0;
    int nbes = 7;
    std::vector<double> alpha;

    smoothgen(nr_nao, r.data(), orb.data(), rcut_proj, alpha);

    // compare with python script result
    std::vector<double> ref = {
        0,
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
        0.013071568549313281
    };

    for (int i = 0; i < 20; ++i) {
        //std::cout<<alpha[i]<<std::endl;
        EXPECT_NEAR(alpha[i], ref[i], 1e-12);
    }
}