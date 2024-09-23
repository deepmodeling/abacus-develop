#ifndef GRID_ANGULAR_DELLEY_H
#define GRID_ANGULAR_DELLEY_H

#include <vector>

/**
 * @brief Delley's grid for quadrature on the unit sphere.
 *
 * Reference:
 * Delley, B. (1996). High order integration schemes on the unit sphere.
 * Journal of computational chemistry, 17(9), 1152-1155.
 *
 */
namespace Grid {
namespace Angular {
namespace Delley {

/**
 * @brief Number of Delley's grid points for a certain order of accuracy.
 *
 * This function finds the minimum Delley's grid with an accuracy
 * of order "lmax". On exit, lmax is set to the order of this grid, and
 * its corresponding number of grid points is returned. If no such grid
 * is available, lmax is left unchanged and the function will return -1.
 *
 * For example, if lmax = 20 on input, the function will return 194 and
 * lmax will be set to 23.
 *
 */
int ngrid(int& lmax);


/**
 * @brief Delley's quadrature grid and weights.
 *
 * This function retrieves the minimum Delley's grid with an accuray
 * of order "lmax". On exit, lmax is set to the order of this grid, and
 * the coordinates & weights are returned in "grid" & "weight".
 *
 * Coordinates are stored in the following order:
 *
 *      x0, y0, z0, x1, y1, z1, x2, y2, z2, ...
 *
 * "grid" and "weight" must be pre-allocated to hold 3*ngrid(lmax) and
 * ngrid(lmax) elements, respectively. The function will return 0 if
 * successful, or -1 if the requested order cannot be fulfilled.
 *
 */
int get(int& lmax, double* grid, double* weight);

// a handy wrapper doing the same as above
int get(int& lmax, std::vector<double>& grid, std::vector<double>& weight);

} // end of namespace Delley
} // end of namespace Angular
} // end of namespace Grid

#endif
