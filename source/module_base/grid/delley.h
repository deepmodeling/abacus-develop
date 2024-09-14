#ifndef GRID_DELLEY_H
#define GRID_DELLEY_H

#include <vector>

class Delley {

public:

    /**
     * @brief Number of grid points required to satisfy a certain order. 
     *
     * This function finds the minimum order of the Delley grid that is
     * equal or above lmax. On exit, lmax is set to this minimum order,
     * and its corresponding number of grid points is returned.
     *
     * This function will return -1 and leave lmax unchanged if the requested
     * lmax is not supported (i.e., greater than 59).
     *
     */
    static int ngrid(int& lmax);

    /**
     * @brief Maximum order below a certain number of grid points
     *
     * This function finds the maximum order (lmax) of the Delley grid whose
     * number of grid points is less than or equal to ngrid. On exit, this
     * order is returned, and ngrid is updated to its corresponding number
     * of grid points.
     *
     * This function will return -1 and leave ngrid unchanged if ngrid is
     * less than 110, which is the size of the smallest Delley grid.
     *
     */
    static int lmax(int& ngrid);

    /**
     * @brief Delley grid and weights.
     *
     * This function generates a Delley grid of order lmax. The grid
     * coordinates are returned in the array "grid" as
     *
     *      x0, y0, z0, x1, y1, z1, x2, y2, z2, ...
     *
     * and the corresponding weights are returned in the array "weight".
     * Their memory must be pre-allocated!
     *
     */
    static void get(
        int lmax,
        double* grid,
        double* weight
    );

    // convenient wrapper doing the same as above
    static void get(
        int lmax,
        std::vector<double>& grid,
        std::vector<double>& weight
    );

private:

};

#endif
