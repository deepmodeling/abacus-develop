#ifndef GRID_BATCH_H
#define GRID_BATCH_H

#include <vector>

#include "module_base/blas_connector.h"
#include "module_base/lapack_connector.h"
#include <algorithm>
#include <cassert>

namespace Grid {
namespace Batch {

int _maxmin_core(int n, const double* grid, int* idx) {
    assert(n > 1);
    if (n == 2) {
        return 1;
    }

    // center of the grid
    std::vector<double> center(3, 0.0);
    for (int i = 0; i < n; ++i) {
        int ig = idx[i];
        center[0] += grid[3*ig    ];
        center[1] += grid[3*ig + 1];
        center[2] += grid[3*ig + 2];
    }
    center[0] /= n;
    center[1] /= n;
    center[2] /= n;

    // positions w.r.t. the center
    std::vector<double> R(3*n, 0.0);
    for (int i = 0; i < n; ++i) {
        int ig = idx[i];
        R[3*i    ] = grid[3*ig    ] - center[0];
        R[3*i + 1] = grid[3*ig + 1] - center[1];
        R[3*i + 2] = grid[3*ig + 2] - center[2];
    }

    // R*R^T
    std::vector<double> A(9, 0.0);
    int three_i = 3;
    double one_d = 1.0, zero_d = 0.0;
    dgemm_("N", "T", &three_i, &three_i, &n, &one_d, R.data(), &three_i, R.data(), &n, &zero_d, A.data(), &three_i);


	//void dsyev_(const char* jobz,const char* uplo,const int* n,double *a,
    //            const int* lda,double* w,double* work,const int* lwork, int* info);
}


std::vector<int> maxmin(int n_max, int n, const double* grid, int* idx) {
    if (n <= n_max) {
        return std::vector<int>{};
    }

    int n_left = _maxmin_core(n, grid, idx);

    std::vector<int> delim_left = maxmin(n_max, n_left, grid, idx);
    std::vector<int> delim_right = maxmin(n_max, n - n_left, grid + n_left, idx + n_left);;
    std::for_each(delim_right.begin(), delim_right.end(), [n_left](int& x) { x += n_left; });

    // merge all delimiters into delim_left
    delim_left.resize(delim_left.size() + delim_right.size() + 1);
    delim_left[delim_left.size()] = n_left;
    std::copy(delim_right.begin(), delim_right.end(), delim_left.begin() + delim_left.size() + 1);
    return delim_left;
}

} // end of namespace Batch
} // end of namespace Grid

#endif
