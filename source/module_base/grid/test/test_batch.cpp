#include "module_base/grid/batch.h"

#include "gtest/gtest.h"
#include <algorithm>
#include <random>
//#include <cstdio>

#ifdef __MPI
#include <mpi.h>
#endif

using namespace Grid::Batch;


class BatchTest: public ::testing::Test
{
protected:

    std::vector<double> grid_;
    std::vector<int> idx_;

    // parameters for 8-octant clusters
    const int n_batch_oct_ = 10;
    const double width_oct_ = 1.0;
    const double offset_x_ = 7.0;
    const double offset_y_ = 8.0;
    const double offset_z_ = 9.0;
    // NOTE: These offsets should be different from each other as maxmin
    // might fail for highly symmetric, well-separated clusters.
    // Consider the case where the 8 clusters as a whole have octahedral
    // symmetry. In this case, R*R^T must be proprotional to the identity,
    // and eigenvalues are three-fold degenerate, because xy, yz and zx
    // plane are equivalent in terms of the maxmin optimization problem.
    // This means eigenvectors are arbitrary in this case.


    // parameters for a random cluster
    const int n_grid_rand_ = 1000;
    const int n_batch_rand_ = 200;
    const double width_rand_ = 10.0;
    const double xc_ = 1.0;
    const double yc_ = 1.0;
    const double zc_ = 2.0;
};


void gen_random(
    int ngrid,
    double xc,
    double yc,
    double zc,
    double width,
    std::vector<double>& grid,
    std::vector<int>& idx
) {

    // Generates a set of points centered around (xc, yc, zc).

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-width, width);

    grid.resize(3 * ngrid);
    for (int i = 0; i < ngrid; ++i) {
        grid[3*i    ] = xc + dis(gen);
        grid[3*i + 1] = yc + dis(gen);
        grid[3*i + 2] = zc + dis(gen);
    }

    idx.resize(ngrid);
    std::iota(idx.begin(), idx.end(), 0);
    std::shuffle(idx.begin(), idx.end(), gen);
}


void gen_octant(
    int n_each,
    double offset_x,
    double offset_y,
    double offset_z,
    double width,
    std::vector<double>& grid,
    std::vector<int>& idx
) {

    // Generates a set of points consisting of 8 well-separated, equal-sized
    // clusters located in individual octants.

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-width, width);

    int ngrid = 8 * n_each;
    grid.resize(3 * ngrid);
    int I = 0;
    for (int sign_x : {-1, 1}) {
        for (int sign_y : {-1, 1}) {
            for (int sign_z : {-1, 1}) {
                for (int i = 0; i < n_each; ++i, ++I) {
                    grid[3*I    ] = sign_x * offset_x + dis(gen);
                    grid[3*I + 1] = sign_y * offset_y + dis(gen);
                    grid[3*I + 2] = sign_z * offset_z + dis(gen);
                }
            }
        }
    }

    idx.resize(ngrid);
    std::iota(idx.begin(), idx.end(), 0);
    std::shuffle(idx.begin(), idx.end(), gen);
}


bool is_same_octant(int ngrid, const double* grid) {
    if (ngrid == 0) {
        return true;
    }
    const bool is_positive_x = grid[0] > 0;
    const bool is_positive_y = grid[1] > 0;
    const bool is_positive_z = grid[2] > 0;
    const double* end = grid + 3 * ngrid;
    for (; grid != end; grid += 3) {
        if ( is_positive_x != (grid[0] > 0) ||
             is_positive_y != (grid[1] > 0) ||
             is_positive_z != (grid[2] > 0) ) {
            return false;
        }
    }
    return true;
}


bool good_batch_size(
    const std::vector<int>& idx,
    const std::vector<int>& delim,
    int n_batch_thr
) {
    // checks if the sizes of batches are within the specified limit

    bool flag = (delim[0] == 0);

    size_t i = 1;
    while (flag && i < delim.size()) {
        int sz_batch = delim[i] - delim[i-1];
        flag = flag && (sz_batch > 0) && (sz_batch <= n_batch_thr);
        ++i;
    }

    return flag && ( ((int)idx.size() - delim.back()) < n_batch_thr );
}


TEST_F(BatchTest, MaxMinRandom)
{
    // This test verifies that the sizes of batches produced by maxmin
    // do not exceed the specified limit.

    gen_random(n_grid_rand_, xc_, yc_, zc_, width_rand_, grid_, idx_);

    std::vector<int> delim = 
        maxmin(grid_.data(), idx_.data(), idx_.size(), n_batch_rand_);

    EXPECT_TRUE(good_batch_size(idx_, delim, n_batch_rand_));

    // write grid, idx & delim to file
    //FILE* fp = fopen("grid.dat", "w");
    //for (size_t i = 0; i < grid_.size() / 3; ++i) {
    //    std::fprintf(fp, "% 12.6f % 12.6f % 12.6f\n",
    //        grid_[3*i], grid_[3*i + 1], grid_[3*i + 2]);
    //}
    //fclose(fp);

    //fp = fopen("idx.dat", "w");
    //for (size_t i = 0; i < idx_.size(); ++i) {
    //    std::fprintf(fp, "%d\n", idx_[i]);
    //}
    //fclose(fp);

    //fp = fopen("delim.dat", "w");
    //for (size_t i = 0; i < delim.size(); ++i) {
    //    std::fprintf(fp, "%d\n", delim[i]);
    //}
    //fclose(fp);
}


TEST_F(BatchTest, MaxMinOctant)
{
    // This test applies maxmin to a set of points consisting of 8
    // well-separated, equal-sized clusters located in individual octants.
    // The resulting batches should be able to recover this structure.

    gen_octant(n_batch_oct_, offset_x_, offset_y_, offset_z_, width_oct_,
               grid_, idx_);

    std::vector<int> delim = 
        maxmin(grid_.data(), idx_.data(), idx_.size(), n_batch_oct_);

    EXPECT_EQ(delim.size(), 8);

    std::vector<double> grid_batch(3 * n_batch_oct_);
    for (int i = 0; i < 8; ++i) {

        EXPECT_EQ(delim[i], i * n_batch_oct_);

        // collect points within the present batch
        for (int j = 0; j < n_batch_oct_; ++j) {
            int ig = idx_[delim[i] + j];
            grid_batch[3*j    ] = grid_[3*ig    ];
            grid_batch[3*j + 1] = grid_[3*ig + 1];
            grid_batch[3*j + 2] = grid_[3*ig + 2];
        }

        // verify that points in a batch reside in the same octant
        EXPECT_TRUE(is_same_octant(n_batch_oct_, grid_batch.data()));
    }
}


int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    return result;
}
