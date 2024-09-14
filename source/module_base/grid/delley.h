#ifndef GRID_DELLEY_H
#define GRID_DELLEY_H

#include <vector>
#include <functional>

/**
 * @brief Delley's grid for quadrature on the unit sphere.
 *
 * Reference:
 * Delley, B. (1996). High order integration schemes on the unit sphere.
 * Journal of computational chemistry, 17(9), 1152-1155.
 *
 */
class Delley {

public:

    /**
     * @brief Number of grid points required for a certain order of accuracy.
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
     * @brief Delley grid and weights.
     *
     * This function retrieves the Delley grid of order lmax. The grid
     * coordinates are returned in the array "grid" as
     *
     *      x0, y0, z0, x1, y1, z1, x2, y2, z2, ...
     *
     * "grid" and "weight" must be pre-allocated to hold 3*ngrid(lmax) and
     * ngrid(lmax) elements, respectively. The function will abort if the
     * requested order is not supported.
     *
     */
    static void get(
        int lmax,
        double* grid,
        double* weight
    );

    // a handy wrapper function doing the same as above
    static void get(
        int lmax,
        std::vector<double>& grid,
        std::vector<double>& weight
    );


private:

    struct DelleyTable {
        const int lmax_;
        const int ngrid_;
        const int ntype_[6];
        const std::vector<double> data_;
    };

    using FillFunc = std::function<void(double, double, double*)>;

    static constexpr int group_size_[] = {6, 8, 12, 24, 24, 48};
    static const std::vector<DelleyTable> table_;
    static const std::vector<FillFunc> fill_;

    static const DelleyTable* _retrieve(int lmax);
    static void _get(const DelleyTable* tab, double* grid, double* weight);
};

#endif
