#include "module_base/grid/partition.h"
#include "module_base/constants.h"

#include <cmath>
#include <numeric>
#include <vector>

namespace Grid {
namespace Partition {

double w_becke(
    int nR0,
    double* drR,
    double* dRR,
    int nR,
    int* iR,
    int c
) {
    std::vector<double> P(nR, 1.0);
    for (int i = 0; i < nR; ++i) {
        int I = iR[i];
        for (int j = i + 1; j < nR; ++j) {
            int J = iR[j];
            double mu = (drR[I] - drR[J]) / dRR[I*nR0 + J];
            double s = s_becke(mu);
            P[I] *= s;
            P[J] *= (1.0 - s); // s_becke(-mu) = 1 - s_becke(mu)
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


double s_stratmann(double mu, double a) {
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
    double x = mu / a;
    double x2 = x * x;
    double h = 0.625 * x * (35 + x2 * (-35 + x2 * (21 - 5 * x2)));

    bool mid = std::abs(x) < 1;
    double g = !mid * (1 - 2 * std::signbit(x)) + mid * h;
    return 0.5 * (1 - g);
}


double u_stratmann_mod(double y, double b) {
    using ModuleBase::PI;
    bool core = y <= b;
    bool edge = !core && y < 1.0;
    return core + edge * 0.5 * (std::cos(PI * (y - b) / (1.0 - b)) + 1.0);
}


double s_stratmann_mod(double mu, double y, double a, double b) {
    // Modified Stratmann's cell function by Knuth et al.
    // y = |r-R(J)| / Rcut(J)
    return 1.0 + u_stratmann_mod(y, b) * (s_stratmann(mu, a) - 1.0);
}


//double weight(
//    int nRtot,
//    double* drR,
//    double* dRR,
//    int nR,
//    int* iR,
//    int c,
//    CellFuncType type,
//    double* Rcut
//) {
//    // unnormalized weight
//    std::vector<double> P(nR, 1.0);
//
//    // confocal ellipsoidal coordinates
//    std::vector<double> mu(nR*nR, 0.0);
//    for (int i = 0; i < nR; ++i) {
//        int I = iR[i];
//        for (int j = i + 1; j < nR; ++j) {
//            int J = iR[j];
//            mu[I*nR + J] = (drR[I] - drR[J]) / dRR[I*nRtot + J];
//            mu[J*nR + I] = -mu[I*nR + J];
//        }
//    }
//
//    for (int i = 0; i < nR; ++i) {
//        if (i == c) {
//            continue;
//        }
//        int J = iR[i];
//        double mu = (drR[I] - drR[J]) / dRR[I*nRtot + J];
//        P[i] *= s_becke(mu);
//    }
//
//    return P[c] / std::accumulate(P.begin(), P.end(), 0.0);
//}

} // end of namespace Partition
} // end of namespace Grid
