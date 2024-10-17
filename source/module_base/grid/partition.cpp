#include "module_base/grid/partition.h"
#include "module_base/constants.h"

#include <cmath>
#include <functional>
#include <numeric>
#include <algorithm>

namespace Grid {
namespace Partition {

const double stratmann_a = 0.64;
const double stratmann_mod_b = 0.8;


double w_becke(
    int nR0,
    const double* drR,
    const double* dRR,
    int nR,
    const int* iR,
    int c
) {
    double P[nR];
    std::fill(P, P + nR, 1.0);
    for (int i = 0; i < nR; ++i) {
        int I = iR[i];
        for (int j = i + 1; j < nR; ++j) {
            int J = iR[j];
            double mu = (drR[I] - drR[J]) / dRR[I*nR0 + J];
            double s = s_becke(mu);
            P[I] *= s;
            P[J] *= (1.0 - s); // s(-mu) = 1 - s(mu)
        }
    }
    return P[c] / std::accumulate(P, P + nR, 0.0);
}


double s_becke(double mu) {
    /* 
     * Becke's iterated polynomials (3rd order)
     *
     * s(mu) = 0.5 * (1 - p(p(p(mu))))
     *
     * p(x) = 0.5 * x * (3 - x^2)
     *
     */
    double p = 0.5 * mu * (3.0 - mu*mu);
    p = 0.5 * p * (3.0 - p*p);
    p = 0.5 * p * (3.0 - p*p);
    return 0.5 * (1.0 - p);
}


double w_stratmann(
    int nR0,
    const double* drR,
    const double* dRR,
    const double* drR_thr,
    const int nR,
    int* iR,
    int c
) {
    int I = iR[c], J = 0;

    // If r falls within the exclusive zone of a center, return immediately.
    for (int j = 0; j < nR; ++j) {
        J = iR[j];
        if (drR[J] <= drR_thr[J]) {
            return static_cast<double>(I == J);
        }
    }

    // Even if the grid point does not fall within the exclusive zone of any
    // center, the normalized weight could still be 0 or 1, and this can be
    // figured out by examining the unnormalized weight alone.

    // Move the center to the first position for convenience. Swap back later.
    std::swap(iR[0], iR[c]);

    double P[nR];
    for (int j = 1; j < nR; ++j) {
        J = iR[j];
        double mu = (drR[I] - drR[J]) / dRR[I*nR0 + J];
        P[j] = s_stratmann(mu);
    }
    P[0] = std::accumulate(P + 1, P + nR, 1.0, std::multiplies<double>());

    if (P[0] == 0.0 || P[0] == 1.0) {
        std::swap(iR[0], iR[c]);
        return P[0];
    }

    // If it passes all the screening, all unnormalized weights have to be
    // calculated in order to get the normalized weight.

    std::for_each(P + 1, P + nR, [](double& s) { s = 1.0 - s; });
    for (int i = 1; i < nR; ++i) {
        I = iR[i];
        for (int j = i + 1; j < nR; ++j) {
            J = iR[j];
            double mu = (drR[I] - drR[J]) / dRR[I*nR0 + J];
            double s = s_stratmann(mu);
            P[i] *= s;
            P[j] *= (1.0 - s); // s(-mu) = 1 - s(mu)
        }
    }

    std::swap(iR[0], iR[c]);
    return P[0] / std::accumulate(P, P + nR, 0.0);
}


double s_stratmann(double mu) {
    /*
     * Stratmann's piecewise cell function
     *
     * s(mu) = 0.5 * (1 - g(mu/a))
     *
     *        /             -1                          x <= -1
     *        |
     * g(x) = | (35x - 35x^3 + 21x^5 - 5x^7) / 16       |x| < 1
     *        |
     *        \             +1                          x >= +1
     *
     */
    double x = mu / stratmann_a;
    double x2 = x * x;
    double h = 0.0625 * x * (35 + x2 * (-35 + x2 * (21 - 5 * x2)));

    bool mid = std::abs(x) < 1;
    double g = !mid * (1 - 2 * std::signbit(x)) + mid * h;
    return 0.5 * (1 - g);
}


double w_stratmann_mod(
    int nR0,
    const double* drR,
    const double* dRR,
    const double* drR_thr,
    const double* Rcut,
    int nR,
    int* iR,
    int c
) {
    int I = iR[c], J = 0;

    // If r falls within the exclusive zone of a center, return immediately.
    for (int j = 0; j < nR; ++j) {
        J = iR[j];
        if (drR[J] <= drR_thr[J]) {
            return static_cast<double>(I == J);
        }
    }

    // Even if the grid point does not fall within the exclusive zone of any
    // center, the normalized weight could still be 0 or 1, and this can be
    // figured out by examining the unnormalized weight alone.

    // Move the center to the first position for convenience. Swap back later.
    std::swap(iR[0], iR[c]);

    double P[nR];
    for (int j = 1; j < nR; ++j) {
        J = iR[j];
        double mu = (drR[I] - drR[J]) / dRR[I*nR0 + J];
        P[j] = s_stratmann(mu);
    }
    P[0] = std::accumulate(P + 1, P + nR, 1.0, std::multiplies<double>());

    if (P[0] == 0.0 || P[0] == 1.0) {
        std::swap(iR[0], iR[c]);
        return P[0];
    }

    // If it passes all the screening, all unnormalized weights have to be
    // calculated in order to get the normalized weight.

    std::for_each(P + 1, P + nR, [](double& s) { s = 1.0 - s; });
    for (int i = 1; i < nR; ++i) {
        I = iR[i];
        for (int j = i + 1; j < nR; ++j) {
            J = iR[j];
            double mu = (drR[I] - drR[J]) / dRR[I*nR0 + J];
            double s = s_stratmann(mu);
            P[i] *= s;
            P[j] *= (1.0 - s); // s(-mu) = 1 - s(mu)
        }
    }

    std::swap(iR[0], iR[c]);
    return P[0] / std::accumulate(P, P + nR, 0.0);
}


double s_stratmann_mod(double mu, double y) {
    using ModuleBase::PI;
    bool core = y <= stratmann_mod_b;
    bool edge = !core && y < 1.0;
    double u = core + edge * 0.5 * (1.0 +
            std::cos(PI * (y - stratmann_mod_b) / (1.0 - stratmann_mod_b)));
    return 1.0 + u * (s_stratmann(mu) - 1.0);
}



} // end of namespace Partition
} // end of namespace Grid
