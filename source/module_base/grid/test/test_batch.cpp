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

    int n_each_ = 10;
    double offset_ = 10.0;
    double width_ = 1.0;
};

std::vector<double> gen_octant_cluster(int n_each, double offset, double width) {

    // Generates a set of points consisting of 8 well-separated, equal-sized
    // clusters located in individual octants.

    std::vector<double> grid(n_each * 8);
    int I = 0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-width, width);

    for (int sign_x : {-1, 1}) {
        for (int sign_y : {-1, 1}) {
            for (int sign_z : {-1, 1}) {
                for (int i = 0; i < n_each; ++i) {
                    grid[3*I    ] = sign_x * offset + dis(gen);
                    grid[3*I + 1] = sign_y * offset + dis(gen);
                    grid[3*I + 2] = sign_z * offset + dis(gen);
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
    grid_ = gen_octant_cluster(n_each_, offset_, width_);

    idx_.resize(grid_.size());
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
        maxmin(n_each_, grid_.size(), grid_.data(), idx_.data());

    EXPECT_EQ(delim.size(), 7);
    for (int i = 0; i < 7; ++i) {
        // check number of points in each batch via index delimiters
        EXPECT_EQ(delim[i], (i+1) * n_each_);

        // verify that points in each batch is in the same octant
        std::vector<double> batch(3 * n_each_);
        for (int j = 0; j < n_each_; ++j) {
            for (int k = 0; k < 3; ++k) {
                batch[3*j + k] = grid_[3*(i*n_each_ + j) + k];
            }
        }
        EXPECT_TRUE(is_same_octant(n_each_, batch.data()));
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
