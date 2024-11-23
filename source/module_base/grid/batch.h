#ifndef GRID_BATCH_H
#define GRID_BATCH_H

#include <vector>

namespace Grid {
namespace Batch {

/**
 * @brief Divide a set of points into batches by the "MaxMin" algorithm.
 *
 * This function recursively uses a cut plane to divide a set of grid
 * points into two subsets using the "MaxMin" algorithm, until the
 * number of points in each subset (batch) is no more than m_thr.
 *
 * @param[in]       grid        Coordinates of all grid points.
 * @param[in,out]   idx         Indices of the initial set within grid.
 *                              On return, idx will be rearranged such
 *                              that points belonging to the same subset
 *                              are grouped together.
 * @param[in]       m           Number of points in the initial set.
 * @param[in]       m_thr       Size limit of subset.
 *
 * @return          Indices (within idx) of the first point in each batch.
 *
 */
std::vector<int> maxmin(const double* grid, int* idx, int m, int m_thr);

} // end of namespace Batch
} // end of namespace Grid

#endif
