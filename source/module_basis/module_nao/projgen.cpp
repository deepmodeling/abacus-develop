#include "module_base/math_integral.h"
#include "module_base/math_sphbes.h"
#include "projgen.h"

#include <numeric>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <cstdio>

using namespace ModuleBase;

void projgen(const int l, const int nr, const double* r, const double* chi, const double rcut, const int nbes, std::vector<double>& alpha) 
{
    assert(rcut < r[nr - 1]);
    assert(std::is_sorted(r, r + nr));

    std::vector<double> dr(nr - 1);
    std::adjacent_difference(r, r + nr, dr.begin());

    // lower_bound returns the first element that is equal or larger than rcut
    int nr_proj = std::distance(r, std::lower_bound(r, r + nr, rcut)) + 1;

    // zeros of spherical Bessel function
    std::vector<double> theta(nbes);
    Sphbes::sphbes_zeros(l, nbes, theta.data());

    // z & w vectors (see notes)
    std::vector<double> z(nbes);
    std::vector<double> w(nbes);

    std::transform(theta.begin(), theta.end(), z.begin(), [rcut, l](double theta_p)
            { return 0.5 * std::pow(rcut, 3) * std::pow(Sphbes::sphbesj(l+1, theta_p), 2); });

    // r^2 * chi (independent from p)
    std::vector<double> tmp(nr_proj);
    std::transform(r, r + nr_proj, chi, tmp.begin(), [](double r_i, double chi_i) { return r_i * r_i * chi_i; });

    // r^2 * chi * j_l(theta[p] * r / rcut) (dependent on p)
    std::vector<double> integrand(nr_proj);

    for (int p = 0; p < nbes; ++p) {
        std::transform(r, r + nr_proj, tmp.begin(), integrand.begin(), [theta, p, rcut, l](double r_i, double tmp_i)
                { return tmp_i * Sphbes::sphbesj(l, theta[p] * r_i / rcut); });
        w[p] = Integral::simpson(nr_proj, integrand.data(), &dr[1]);
    }

    // optimal coefficients
    std::vector<double> c(nbes, 0.0);
    std::transform(w.begin(), w.end(), z.begin(), c.begin(), [](double w_p, double z_p) { return w_p * w_p / z_p; });
    double prefac = 1.0 / std::sqrt(std::accumulate(c.begin(), c.end(), 0.0));
    std::transform(w.begin(), w.end(), z.begin(), c.begin(), [prefac](double w_p, double z_p)
            { return prefac * w_p / z_p; });

    // new radial function
    alpha.resize(nr_proj);
    std::fill(alpha.begin(), alpha.end(), 0.0);
    for (int i = 0; i < nr_proj; ++i) {
        for (int p = 0; p < nbes; ++p) {
            alpha[i] += c[p] * Sphbes::sphbesj(l, theta[p] * r[i] / rcut);
        }
    }
}

/*int main() {

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
    double N = 1.0 / std::sqrt(Integral::simpson(nr_nao, integrand.data(), dr));
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
        assert(std::abs(alpha[i] - ref[i]) < 1e-12);
    }
    printf("test passed!\n");

    return 0;
}*/

