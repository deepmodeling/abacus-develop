#ifndef GRID_BATCH_H
#define GRID_BATCH_H

#include <vector>

namespace Grid {
namespace Batch {

/**
 * @brief Divide a set of points into batches by the "MaxMin" algorithm.
 *
 * This function recursively uses cut planes to divide grid points into
 * two subsets using the "MaxMin" algorithm, until the number of points
 * in each subset (batch) is no more than m_thr.
 *
 * @param[in]       grid        Coordinates of all grid points.
 *                              grid[3*j], grid[3*j+1], grid[3*j+2] are
 *                              the x, y, z coordinates of the j-th point.
 * @param[in,out]   idx         Indices of the initial set within grid.
 *                              On return, idx will be rearranged such
 *                              that points belonging to the same batch
 *                              have their indices placed together.
 * @param[in]       m           Number of points in the initial set.
 *                              (length of idx)
 * @param[in]       m_thr       Size limit of a batch.
 *
 * @return          Indices (for idx) of the first point in each batch.
 *
 * For example, given grid (illustrated by their indices) located as follows:
 *
 *      0  1  16          2  3  18
 *      4  5  17            6  7
 *
 *
 *      8  9             20 10 11
 *     12 13 19           14 15
 *
 * a possible outcome with m_thr = 4 and idx(in) = {0, 1, 2, ..., 15}
 * (idx may correspond to a subset of grid and does not have to be sorted,
 * but it must not contain duplicates) is:
 *
 * idx(out): 0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15
 * return  : {0, 4, 8, 12}
 *
 * which means the selected set (labeled 0-15) is divided into 4 batches:
 * {0, 1, 4, 5}, {8, 9, 12, 13}, {2, 3, 6, 7}, {10, 11, 14, 15}.
 *
 */
std::vector<int> maxmin(const double* grid, int* idx, int m, int m_thr);

} // end of namespace Batch
} // end of namespace Grid

#endif
