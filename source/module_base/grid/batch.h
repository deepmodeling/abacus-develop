#ifndef GRID_BATCH_H
#define GRID_BATCH_H

#include <vector>

namespace Grid {
namespace Batch {

/**
 * @brief Divide a set of points into batches by the "MaxMin" algorithm.
 *
 * This function recursively divides a given set of grid points into two
 * subsets by a cut plane using the "MaxMin" algorithm, until the number
 * of points in each subset is no more than n_max.
 *
 */
std::vector<int> maxmin(int n_max, int n, const double* grid, int* idx);

} // end of namespace Batch
} // end of namespace Grid

#endif
