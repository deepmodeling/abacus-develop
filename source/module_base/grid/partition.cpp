#include "module_base/grid/partition.h"
#include "module_base/constants.h"

#include <cmath>
#include <functional>
#include <numeric>
#include <algorithm>
#include <vector>
#include <cassert>

namespace Grid {
namespace Partition {

const double stratmann_a = 0.64;

double w_becke(
    int nR0,
    const double* drR,
    const double* dRR,
    int nR,
    const int* iR,
    int c
) {
    assert(nR > 0 && nR0 >= nR);
    std::vector<double> P(nR, 1.0);
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
    return P[c] / std::accumulate(P.begin(), P.end(), 0.0);
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
    int nR,
    int* iR,
    int c
) {
    assert(nR > 0 && nR0 >= nR);
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

    // Swap the grid center to the first position in iteration for convenience.
    // Restore the original order before return.
    std::swap(iR[0], iR[c]);

    std::vector<double> P(nR);
    for (int j = 1; j < nR; ++j) {
        J = iR[j];
        double mu = (drR[I] - drR[J]) / dRR[I*nR0 + J];
        P[j] = s_stratmann(mu);
    }
    P[0] = std::accumulate(P.begin() + 1, P.end(), 1.0,
                           std::multiplies<double>());

    if (P[0] == 0.0 || P[0] == 1.0) {
        std::swap(iR[0], iR[c]); // restore the original order
        return P[0];
    }

    // If it passes all the screening above, all unnormalized weights
    // have to be calculated in order to get the normalized weight.

    std::for_each(P.begin() + 1, P.end(), [](double& s) { s = 1.0 - s; });
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

    std::swap(iR[0], iR[c]); // restore the original order
    return P[0] / std::accumulate(P.begin(), P.end(), 0.0);
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
    return 0.5 * (1.0 - g);
}


} // end of namespace Partition
} // end of namespace Grid
