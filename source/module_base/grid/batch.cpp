#include "module_base/grid/batch.h"

#include <algorithm>
#include <cassert>
#include <iterator>

#include "module_base/blas_connector.h"
#include "module_base/lapack_connector.h"

namespace {

/**
 * @brief Divide a set of points into two subsets by the "MaxMin" algorithm.
 *
 * This function divides a given set of grid points by a cut plane
 * {x|n^T*(x-c) = 0} where the normal vector n and the point c are
 * determined by the "MaxMin" problem:
 *
 *      max min sum_{i=1}^{m} [n^T * (r[idx[i]] - c)]^2
 *       n   c
 *
 * here r[j] = (grid[3*j], grid[3*j+1], grid[3*j+2]) is the position of
 * the j-th point.
 *
 * It can be shown that the optimal c is the centroid of the points, and
 * the optimal n is the eigenvector corresponding to the largest eigenvalue
 * of the matrix R*R^T, where the i-th column of R is r[idx[i]] - c.
 *
 * @param[in]       grid    Coordinates of all grid points.
 *                          grid[3*j], grid[3*j+1], grid[3*j+2] are the
 *                          x, y, z coordinates of the j-th point.
 * @param[in,out]   idx     Indices of the selected points within grid.
 *                          On return, idx will be rearranged such that
 *                          points belonging to the same subset have their
 *                          indices placed together.
 * @param[in]       m       Number of selected points (length of idx).
 *
 * @return The number of points in the first subset within idx.
 *
 */
int _maxmin_divide(const double* grid, int* idx, int m) {
    assert(m > 1);
    if (m == 2) {
        return 1;
    }

    std::vector<double> centroid(3, 0.0);
    for (int i = 0; i < m; ++i) {
        int j = idx[i];
        centroid[0] += grid[3*j    ];
        centroid[1] += grid[3*j + 1];
        centroid[2] += grid[3*j + 2];
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

    // The normal vector of the cut plane is taken to be the eigenvector
    // corresponding to the largest eigenvalue of the 3x3 matrix A = R*R^T.
    std::vector<double> A(9, 0.0);
    int i3 = 3, i1 = 1;
    double d0 = 0.0, d1 = 1.0;
    dsyrk_("U", "N", &i3, &m, &d1, R.data(), &i3, &d0, A.data(), &i3);

    int info = 0, lwork = 102 /* determined by a work space query */;
    std::vector<double> e(3), work(lwork);
    dsyev_("V", "U", &i3, A.data(), &i3, e.data(), work.data(), &lwork, &info);
    double* n = A.data() + 6; // normal vector of the cut plane

    // Rearrange the indices to put points in each subset together by
    // examining the signed distances of points to the cut plane (R^T*n).
    std::vector<double> dist(m);
    dgemv_("T", &i3, &m, &d1, R.data(), &i3, n, &i1, &d0, dist.data(), &i1);

    int *head = idx;
    std::reverse_iterator<int*> tail(idx + m), rend(idx);
    auto is_negative = [&dist, &idx](int& j) { return dist[&j - idx] < 0; };
    while ( ( head = std::find_if(head, idx + m, is_negative) ) <
            ( tail = std::find_if_not(tail, rend, is_negative) ).base() ) {
        std::swap(*head, *tail);
        std::swap(dist[head - idx], dist[tail.base() - idx - 1]);
        ++head;
        ++tail;
    }

    return head - idx;
}

} // end of anonymous namespace


std::vector<int> Grid::Batch::maxmin(
    const double* grid,
    int* idx,
    int m,
    int m_thr
) {
    if (m <= m_thr) {
        return std::vector<int>{0};
    }

    int m_left = _maxmin_divide(grid, idx, m);

    std::vector<int> left = maxmin(grid, idx, m_left, m_thr);
    std::vector<int> right = maxmin(grid, idx + m_left, m - m_left, m_thr);
    std::for_each(right.begin(), right.end(),
        [m_left](int& x) { x += m_left; }
    );

    left.insert(left.end(), right.begin(), right.end());
    return left;
}


