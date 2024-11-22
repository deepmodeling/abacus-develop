#include "module_base/grid/batch.h"

#include "module_base/blas_connector.h"
#include "module_base/lapack_connector.h"
#include <algorithm>
#include <cassert>
#include <iterator>

namespace {

/**
 * @brief Bisect a set of points by the "MaxMin" algorithm.
 *
 * Given a selected set of grid points specified by the indices `idx`,
 * bisect this set by a cut plane {x|n^T*(x-c) = 0} where the normal
 * vector n and the point c are determined by the "MaxMin" problem:
 *
 *      max min sum_{i=1}^{m} [n^T*(r[idx[i]] - c)]^2
 *       n   c
 *
 * here r[j] = (grid[3*j], grid[3*j+1], grid[3*j+2]) is the position of
 * the j-th point.
 *
 * It can be shown that the optimal c is the centroid of the points, and
 * the optimal n is the eigenvector corresponding to the largest eigenvalue
 * of the matrix R*R^T, where R is the matrix whose i-th column is
 * r[idx[i]] - c.
 *
 * param[in]        m       Number of the selected points (size of idx).
 * param[in]        grid    Coordinates of all grid points.
 *                          grid[3*j], grid[3*j+1], grid[3*j+2] are the
 *                          x, y, z coordinates of the j-th point.
 *                          The size of grid is at least 3*m.
 * param[in,out]    idx     Indices of the selected points within grid.
 *
 */
int _maxmin_bisect(int m, const double* grid, int* idx) {
    assert(m > 1);
    if (m == 2) {
        return 1;
    }

    std::vector<double> centroid(3, 0.0);
    for (int i = 0; i < m; ++i) {
        int ig = idx[i];
        centroid[0] += grid[3*ig    ];
        centroid[1] += grid[3*ig + 1];
        centroid[2] += grid[3*ig + 2];
    }
    centroid[0] /= m;
    centroid[1] /= m;
    centroid[2] /= m;

    // positions w.r.t. the centroid
    std::vector<double> R(3*m, 0.0);
    for (int i = 0; i < m; ++i) {
        int j = idx[i];
        R[3*i    ] = grid[3*j    ] - centroid[0];
        R[3*i + 1] = grid[3*j + 1] - centroid[1];
        R[3*i + 2] = grid[3*j + 2] - centroid[2];
    }

    // A = R*R^T is a 3-by-3 matrix
    std::vector<double> A(9, 0.0);
    int i3 = 3, i1 = 1;
    double d0 = 0.0, d1 = 1.0;
    dsyrk_("U", "N", &i3, &m, &d1, R.data(), &i3, &d0, A.data(), &i3);

    // eigenpairs of A
    int info = 0, lwork = 102 /* determined by a work space query */;
    std::vector<double> e(3), work(lwork);
    dsyev_("V", "U", &i3, A.data(), &i3, e.data(), work.data(), &lwork, &info);

    // normal vector of the cut plane
    // (eigenvector corresponding to the largest eigenvalue)
    double* n = A.data() + 6;

    // (signed) distance w.r.t. the cut plane
    std::vector<double> dist(m);
    for (int i = 0; i < m; ++i) {
        dist[i] = ddot_(&i3, R.data() + 3*i, &i1, n, &i1);
    }

    // 
    int *head = idx;
    std::reverse_iterator<int*> tail(idx + m), rend(idx);
    auto is_negative = [&dist](int j) { return dist[j] < 0; };
    auto is_nonnegative = [&dist](int j) { return dist[j] >= 0; };
    while ( (head = std::find(head, idx + m, is_negative)) <
            (tail = std::find(tail, rend, is_nonnegative)).base()) {
        std::swap(*head, *tail);
        std::swap(dist[head - idx], dist[tail.base() - idx]);
        ++head;
        ++tail;
    }

    return std::find(idx, idx + m, is_nonnegative) - idx;
}

} // end of anonymous namespace

std::vector<int> Grid::Batch::maxmin(
    int n_max,
    int n,
    const double* grid,
    int* idx
) {
    if (n <= n_max) {
        return std::vector<int>{};
    }

    int n_left = _maxmin_bisect(n, grid, idx);

    std::vector<int> delim_left = maxmin(n_max, n_left, grid, idx);
    std::vector<int> delim_right = maxmin(n_max, n - n_left, grid + n_left, idx + n_left);
    std::for_each(delim_right.begin(), delim_right.end(),
        [n_left](int& x) { x += n_left; }
    );

    // merge all delimiters into delim_left
    delim_left.resize(delim_left.size() + delim_right.size() + 1);
    delim_left[delim_left.size()] = n_left;
    std::copy(delim_right.begin(), delim_right.end(), delim_left.begin() + delim_left.size() + 1);
    return delim_left;
}


