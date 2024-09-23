#include "module_base/grid/radial.h"

#include <cmath>

namespace {
const double pi = std::acos(-1.0);
const double inv_ln2 = 1.0 / std::log(2.0);
}

namespace Grid {
namespace Radial {

void baker(int nbase, double R, int mult, double* grid, double* weight) {
    int ngrid = (nbase+1) * mult - 1;
    double r0 = -R / std::log(1.0 - nbase*nbase/((nbase+1)*(nbase+1)));
    for (int i = 1; i <= ngrid; ++i) {
        grid[i-1] = -r0 * std::log(1.0 - i*i/((ngrid+1)*(ngrid+1)));
        weight[i-1] = 2.0 * i * r0 * grid[i-1] * grid[i-1]
                        / ((ngrid + 1 + i) * (ngrid + 1 - i));
    }
}

void baker(int nbase, double R, int mult, std::vector<double>& grid,
           std::vector<double>& weight) {
    int ngrid = (nbase+1) * mult - 1;
    grid.resize(ngrid);
    weight.resize(ngrid);
    baker(nbase, R, mult, grid.data(), weight.data());
}


void murray(int ngrid, double R, double* grid, double* weight) {
    for (int i = 1; i <= ngrid; ++i) {
        double x = static_cast<double>(i) / (ngrid + 1);
        grid[i-1] = std::pow(x / (1.0 - x), 2) * R;
        weight[i-1] = 2.0 * std::pow(R, 3) * std::pow(x, 5)
                        / (std::pow(1.0 - x, 7) * (ngrid + 1));
    }
}


void ta3(int ngrid, double R, double* grid, double* weight) {
    for (int i = 1; i <= ngrid; ++i) {
        double x = std::cos(i * pi / (ngrid + 1));
        double beta = std::sqrt((1.0 + x) / (1.0 - x));
        double gamma = std::log((1.0 - x) / 2.0);
        grid[i-1] = -R * inv_ln2 * gamma;
        weight[i-1] = pi / (ngrid + 1) * std::pow(R * inv_ln2, 3)
                      * gamma * gamma * beta;
    }
}


void ta4(int ngrid, double R, double* grid, double* weight, double alpha) {
    for (int i = 1; i <= ngrid; ++i) {
        double x = std::cos(i * pi / (ngrid + 1));
        double beta = std::sqrt((1.0 + x) / (1.0 - x));
        double gamma = std::log((1.0 - x) / 2.0);
        double delta = std::pow(1.0 + x, alpha);
        grid[i-1] = -R * inv_ln2 * delta * gamma;
        weight[i-1] = pi / (ngrid + 1) * std::pow(delta * R * inv_ln2, 3)
                      * gamma * gamma * (beta - alpha / beta * gamma);
    }
}

}
}
