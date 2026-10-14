#include "projgen.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <numeric>

#include "source_base/math_integral.h"
#include "source_base/math_sphbes.h"

using namespace ModuleBase;

void projgen(const int l,
             const int nr,
             const double* r,
             const double* chi,
             const double rcut,
             const int nbes,
             std::vector<double>& alpha)
{
    assert(rcut < r[nr - 1]);
    assert(std::is_sorted(r, r + nr));

    std::vector<double> dr(nr);
    std::adjacent_difference(r, r + nr, dr.begin());

    // lower_bound returns the first element that is equal or larger than rcut
    int nr_proj = std::distance(r, std::lower_bound(r, r + nr, rcut)) + 1;

    // zeros of spherical Bessel function
    std::vector<double> theta(nbes);
    Sphbes::sphbes_zeros(l, nbes, theta.data());

    // z & w vectors (see notes)
    std::vector<double> z(nbes);
    std::vector<double> w(nbes);

    std::transform(theta.begin(), theta.end(), z.begin(), [rcut, l](double theta_p) {
        return 0.5 * std::pow(rcut, 3) * std::pow(Sphbes::sphbesj(l + 1, theta_p), 2);
    });

    // r^2 * chi (independent from p)
    std::vector<double> tmp(nr_proj);
    std::transform(r, r + nr_proj, chi, tmp.begin(), [](double r_i, double chi_i) { return r_i * r_i * chi_i; });

    // r^2 * chi * j_l(theta[p] * r / rcut) (dependent on p)
    std::vector<double> integrand(nr_proj);

    for (int p = 0; p < nbes; ++p)
    {
        std::transform(r, r + nr_proj, tmp.begin(), integrand.begin(), [theta, p, rcut, l](double r_i, double tmp_i) {
            return tmp_i * Sphbes::sphbesj(l, theta[p] * r_i / rcut);
        });
        w[p] = Integral::simpson(nr_proj, integrand.data(), &dr[1]);
    }

    // optimal coefficients
    std::vector<double> c(nbes, 0.0);
    std::transform(w.begin(), w.end(), z.begin(), c.begin(), [](double w_p, double z_p) { return w_p * w_p / z_p; });
    double prefac = 1.0 / std::sqrt(std::accumulate(c.begin(), c.end(), 0.0));
    std::transform(w.begin(), w.end(), z.begin(), c.begin(), [prefac](double w_p, double z_p) {
        return prefac * w_p / z_p;
    });

    // new radial function
    alpha.resize(nr_proj);
    std::fill(alpha.begin(), alpha.end(), 0.0);
    for (int i = 0; i < nr_proj; ++i)
    {
        for (int p = 0; p < nbes; ++p)
        {
            alpha[i] += c[p] * Sphbes::sphbesj(l, theta[p] * r[i] / rcut);
        }
    }
}

