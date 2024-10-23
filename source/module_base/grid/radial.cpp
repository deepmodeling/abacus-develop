#include "module_base/grid/radial.h"

#include <cmath>

namespace Grid {
namespace Radial {

void baker(int nbase, double R, double* r, double* w, int mult) {
    int n = (nbase+1) * mult - 1;
    double r0 = -R / std::log((2.0 * nbase + 1.0) / ((nbase+1)*(nbase+1)));
    for (int i = 1; i <= n; ++i) {
        r[i-1] = -r0 * std::log(1.0 - static_cast<double>(i)*i/((n+1)*(n+1)));
        w[i-1] = 2.0 * i * r0 * r[i-1] * r[i-1] / ((n+1+i)*(n+1-i));
    }
}


void baker(int nbase, double R, std::vector<double>& r,
           std::vector<double>& w, int mult) {
    int n = (nbase+1) * mult - 1;
    r.resize(n);
    w.resize(n);
    baker(nbase, R, r.data(), w.data(), mult);
}


void murray(int n, double R, double* r, double* w) {
    for (int i = 1; i <= n; ++i) {
        double x = static_cast<double>(i) / (n + 1);
        r[i-1] = std::pow(x / (1.0 - x), 2) * R;
        w[i-1] = 2.0 / (n + 1) * std::pow(R, 3) * std::pow(x, 5)
                 / std::pow(1.0 - x, 7);
    }
}


void treutler_m4(int n, double R, double* r, double* w, double alpha) {
    const double pi = std::acos(-1.0);
    const double inv_ln2 = 1.0 / std::log(2.0);

    for (int i = 1; i <= n; ++i) {
        double x = std::cos(i * pi / (n + 1));
        double beta = std::sqrt((1.0 + x) / (1.0 - x));
        double gamma = std::log((1.0 - x) / 2.0);
        double delta = std::pow(1.0 + x, alpha);
        r[i-1] = -R * inv_ln2 * delta * gamma;
        w[i-1] = pi / (n + 1) * std::pow(delta * R * inv_ln2, 3)
                 * gamma * gamma * (beta - alpha / beta * gamma);
    }
}


void mura(int n, double R, double* r, double* w) {
    for (int i = 1; i <= n; ++i) {
        double x = static_cast<double>(i) / (n + 1);
        double alpha = 1.0 - x * x * x;
        r[i-1] = -R * std::log(alpha);
        w[i-1] = 3.0 * R * std::pow(x * r[i-1], 2) / ((n+1) * alpha);
    }
}


} // end of namespace Radial
} // end of namespace Grid
