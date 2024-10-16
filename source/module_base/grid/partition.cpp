#include "module_base/grid/partition.h"
#include "module_base/constants.h"

#include <cmath>
#include <numeric>
#include <vector>

namespace Grid {
namespace Partition {

const double stratmann_a = 0.64;
const double stratmann_mod_b = 0.8;

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


double w_stratmann(
    int nR0,
    double* drR,
    double* dRR,
    double* drR_thr,
    int nR,
    int* iR,
    int c
) {

    // if r falls within the exclusive zone of the center
    // whom this grid point belongs to
    int I = iR[c];
    if (drR[I] <= drR_thr[I]) {
        return 1.0;
    }

    // if r falls within the exclusive zone of a center
    // other than the one whom this grid point belongs to
    for (int J = 0; J < nR; ++J) {
        // no need to exclude J == c because it was checked before
        if (drR[iR[J]] <= drR_thr[iR[J]]) {
            return 0.0;
        }
    }

    double Pc = 1.0;
    for (int i = 0; i < c; ++i) {
        int J = iR[i];
        double mu = (drR[I] - drR[J]) / dRR[I*nR0 + J];
        Pc *= s_stratmann(mu);
    }
    for (int i = c + 1; i < nR; ++i) {
        int J = iR[i];
        double mu = (drR[I] - drR[J]) / dRR[I*nR0 + J];
        Pc *= s_stratmann(mu);
    }
    if (Pc == 0.0 || Pc == 1.0) {
        return Pc;
    }

    std::vector<double> P(nR, 1.0);
    for (int i = 0; i < nR; ++i) {
        int I = iR[i];
        for (int j = i + 1; j < nR; ++j) {
            int J = iR[j];
            double mu = (drR[I] - drR[J]) / dRR[I*nR0 + J];
            double s = s_stratmann(mu);
            P[I] *= s;
            P[J] *= (1.0 - s); // s(-mu) = 1 - s(mu)
        }
    }
    return P[c] / std::accumulate(P.begin(), P.end(), 0.0);
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


double u_stratmann_mod(double y) {
    using ModuleBase::PI;
    bool core = y <= stratmann_mod_b;
    bool edge = !core && y < 1.0;
    return core + edge * 0.5 * (1.0 +
            std::cos(PI * (y - stratmann_mod_b) / (1.0 - stratmann_mod_b)));
}


double s_stratmann_mod(double mu, double y) {
    // Modified Stratmann's cell function by Knuth et al.
    // y = |r-R(J)| / Rcut(J)
    return 1.0 + u_stratmann_mod(y) * (s_stratmann(mu) - 1.0);
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
