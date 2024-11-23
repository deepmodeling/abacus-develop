#include "module_base/grid/batch.h"

#include "gtest/gtest.h"
#include <algorithm>
#include <random>

#ifdef __MPI
#include <mpi.h>
#endif

using namespace Grid::Batch;

class BatchTest: public ::testing::Test
{
protected:
    void SetUp();

    std::vector<double> grid_;
    std::vector<int> idx_;

    // parameters for cluster generation
    int n_each_ = 10;
    double width_ = 1.0;

    // These offsets should be different from each other as maxmin might
    // fail for highly symmetric, well-separated clusters.
    // Consider the case where the 8 clusters as a whole have octahedral
    // symmetry. In this case, R*R^T must be proprotional to the identity,
    // and eigenvalues are three-fold degenerate, because xy, yz and zx
    // plane are equivalent in terms of the maxmin optimization problem.
    // This means eigenvectors are arbitrary in this case.
    double offset_x_ = 7.0;
    double offset_y_ = 8.0;
    double offset_z_ = 9.0;
};

std::vector<double> gen_octant_cluster(int n_each, double offset_x, double offset_y, double offset_z, double width) {

    // Generates a set of points consisting of 8 well-separated, equal-sized
    // clusters located in individual octants.

    std::vector<double> grid(n_each * 8 * 3);
    int I = 0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-width, width);

    for (int sign_x : {-1, 1}) {
        for (int sign_y : {-1, 1}) {
            for (int sign_z : {-1, 1}) {
                for (int i = 0; i < n_each; ++i) {
                    grid[3*I    ] = sign_x * offset_x + dis(gen);
                    grid[3*I + 1] = sign_y * offset_y + dis(gen);
                    grid[3*I + 2] = sign_z * offset_z + dis(gen);
                    ++I;
                }
            }
        }
    }
    return grid;
}

bool is_same_octant(int ngrid, const double* grid) {
    if (ngrid == 0) {
        return true;
    }
    bool is_positive_x = grid[0] > 0;
    bool is_positive_y = grid[1] > 0;
    bool is_positive_z = grid[2] > 0;
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


void BatchTest::SetUp()
{
    grid_ = gen_octant_cluster(n_each_, offset_x_, offset_y_, offset_z_, width_);

    idx_.resize(grid_.size() / 3);
    std::iota(idx_.begin(), idx_.end(), 0);

    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(idx_.begin(), idx_.end(), g);
}


TEST_F(BatchTest, MaxMinOctantCluster)
{
    // This test applies maxmin to a set of points consisting of 8
    // well-separated, equal-sized clusters located in individual octants.
    // The resulting batches should be able to recover this structure.

    std::vector<int> delim = 
        maxmin(grid_.data(), idx_.data(), grid_.size() / 3, n_each_);

    EXPECT_EQ(delim.size(), 8);

    std::vector<double> grid_batch(3 * n_each_);
    for (int i = 0; i < 8; ++i) {

        EXPECT_EQ(delim[i], i * n_each_);

        // collect points within the present batch
        for (int j = 0; j < n_each_; ++j) {
            int ig = idx_[delim[i] + j];
            grid_batch[3*j    ] = grid_[3*ig    ];
            grid_batch[3*j + 1] = grid_[3*ig + 1];
            grid_batch[3*j + 2] = grid_[3*ig + 2];
        }

        // verify that points in a batch reside in the same octant
        EXPECT_TRUE(is_same_octant(n_each_, grid_batch.data()));
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
