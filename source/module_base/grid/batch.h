#ifndef GRID_BATCH_H
#define GRID_BATCH_H

#include <vector>

namespace Grid {
namespace Batch {

std::vector<int> maxmin(int n_max, int n, const double* grid, int* idx);

} // end of namespace Batch
} // end of namespace Grid

#endif
