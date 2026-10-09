#include "gtest/gtest.h"

#include "../esolver_dp.h"

#include <vector>

class ESolverDPNeighborListTest : public ::testing::Test
{
  public:
    static std::vector<int> sort_neighbors(const std::vector<double>& coord,
                                           int central_atom,
                                           const int* neighbor_indices,
                                           int neighbor_count)
    {
        return ModuleESolver::ESolver_DP::sort_neighbor_indices_by_distance(coord,
                                                                              central_atom,
                                                                              neighbor_indices,
                                                                              neighbor_count);
    }
};

TEST_F(ESolverDPNeighborListTest, SortsNeighborsByDistanceWithStableIndexTieBreak)
{
    const std::vector<double> coord = {
        0.0, 0.0, 0.0,
        2.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 3.0, 0.0,
        -1.0, 0.0, 0.0};
    const int neighbors[] = {3, 1, 4, 2};

    const std::vector<int> sorted = ESolverDPNeighborListTest::sort_neighbors(coord, 0, neighbors, 4);

    const std::vector<int> expected = {2, 1, 4, 3};
    EXPECT_EQ(sorted, expected);
}
