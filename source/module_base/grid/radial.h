#ifndef GRID_RADIAL_H
#define GRID_RADIAL_H

#include <vector>

namespace Grid {
namespace Radial {

void baker(int nbase, double R, int mult, double* grid, double* weight);

void baker(int nbase, double R, int mult, std::vector<double>& grid,
           std::vector<double>& weight);

void murray(int ngrid, double R, double* grid, double* weight);

void ta3(int ngrid, double R, double* grid, double* weight);

void ta4(int ngrid, double R, double* grid, double* weight, double alpha = 0.6);

}
}

#endif
