#include "projgen.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <numeric>

#include "module_base/cubic_spline.h"
#include "module_base/math_integral.h"
#include "module_base/math_sphbes.h"

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

void smoothgen(const int nr, const double* r, const double* chi, const double rcut, std::vector<double>& alpha)
{
    // lambda function for generate the new radial function
    assert(rcut < r[nr - 1]);
    assert(std::is_sorted(r, r + nr));

    std::vector<double> dr(nr - 1);
    std::adjacent_difference(r, r + nr, dr.begin());

    // lower_bound returns the first element that is equal or larger than rcut
    int nr_proj = std::distance(r, std::lower_bound(r, r + nr, rcut)) + 1;
    alpha.resize(nr_proj);
    auto smooth_sigma = [&](double sigma_in) {
        for (int i = 0; i < nr_proj; i++)
        {
            alpha[i] = chi[i] * (1 - std::exp(-std::pow((r[i] - rcut), 2) / 2 / sigma_in / sigma_in));
        }
        // r^2 * chi (independent from p)
        std::vector<double> tmp(nr_proj);
        std::transform(r, r + nr_proj, alpha.data(), tmp.begin(), [](double r_i, double chi_i) {
            return r_i * r_i * chi_i;
        });

        // r^2 * chi * chi
        std::vector<double> integrand(nr_proj);

        std::transform(alpha.data(),
                       alpha.data() + nr_proj,
                       tmp.begin(),
                       integrand.begin(),
                       [](double chi_i, double tmp_i) { return tmp_i * chi_i; });
        double overlap = ModuleBase::Integral::simpson(nr_proj, integrand.data(), &dr[1]);
        for (int i = 0; i < nr_proj; i++)
        {
            alpha[i] /= std::sqrt(overlap);
        }
        return;
    };

    // cubic spline interpolation
    ModuleBase::CubicSpline cubspl;
    cubspl.build(nr_proj, r, chi);
    std::vector<double> dchi(nr_proj);
    cubspl.eval(nr_proj, r, nullptr, dchi.data());

    // function for calculating the overlap between dalpha and dchi
    auto overlap_dalpha_dchi = [&]() {
        // calculate dalpha first
        ModuleBase::CubicSpline cubspl_alpha;
        cubspl_alpha.build(nr_proj, r, alpha.data());
        std::vector<double> dalpha(nr_proj);
        cubspl_alpha.eval(nr_proj, r, nullptr, dalpha.data());
        for (int i = 0; i < nr_proj; i++)
            dalpha[i] -= dchi[i];
        // r^2 * dchi (independent from p)
        std::vector<double> tmp(nr_proj);
        std::transform(r, r + nr_proj, dalpha.data(), tmp.begin(), [](double r_i, double dalpha_i) {
            return r_i * r_i * dalpha_i;
        });

        // r^2 * dalpha * dchi
        std::vector<double> integrand(nr_proj);

        std::transform(dalpha.data(),
                       dalpha.data() + nr_proj,
                       tmp.begin(),
                       integrand.begin(),
                       [](double dalpha_i, double tmp_i) { return tmp_i * dalpha_i; });
        return ModuleBase::Integral::simpson(nr_proj, integrand.data(), &dr[1]);
    };

    // optimize sigma
    double sigma_left = 0.1;
    smooth_sigma(sigma_left);
    double overlap_alpha_chi_left = overlap_dalpha_dchi();
    double sigma_right = 1.0;
    smooth_sigma(sigma_right);
    double overlap_alpha_chi_right = overlap_dalpha_dchi();
    double overlap_alpha_chi = 0.0;
    double sigma = 0.0;
    while (std::abs(overlap_alpha_chi_right - overlap_alpha_chi_left) > 1e-6)
    {
        sigma = (sigma_left + sigma_right) / 2;
        smooth_sigma(sigma);
        overlap_alpha_chi = overlap_dalpha_dchi();
        if (overlap_alpha_chi < overlap_alpha_chi_left && overlap_alpha_chi < overlap_alpha_chi_right)
        { // the minimum is in the middle
            if (overlap_alpha_chi_left > overlap_alpha_chi_right)
            {
                sigma_left = sigma;
                overlap_alpha_chi_left = overlap_alpha_chi;
            }
            else
            {
                sigma_right = sigma;
                overlap_alpha_chi_right = overlap_alpha_chi;
            }
        }
        else
        { // the minimum is on the left or right
            if (overlap_alpha_chi_left < overlap_alpha_chi_right)
            {
                sigma_right = sigma;
                overlap_alpha_chi_right = overlap_alpha_chi;
                sigma_left = sigma_left - (sigma_right - sigma_left) * 0.5;
                smooth_sigma(sigma_left);
                overlap_alpha_chi_left = overlap_dalpha_dchi();
            }
            else
            {
                sigma_left = sigma;
                overlap_alpha_chi_left = overlap_alpha_chi;
                sigma_right = sigma_right + (sigma_right - sigma_left) * 0.5;
                smooth_sigma(sigma_right);
                overlap_alpha_chi_right = overlap_dalpha_dchi();
            }
        }
    }
}
