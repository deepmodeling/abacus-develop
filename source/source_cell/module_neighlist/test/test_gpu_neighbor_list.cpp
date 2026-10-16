#include <gtest/gtest.h>

#include "source_cell/module_neighlist/gpu_neighbor_list.h"

#include <string>
#include <vector>

using ModuleNeighList::GPU_NeighborList;

TEST(GPUNeighborListTest, RejectsInvalidBuildInput)
{
    GPU_NeighborList neighbor_list;
    std::vector<int> neighbor_count;
    std::vector<int> neighbor_indices;
    std::string error;
    int max_neighbors = -1;

    EXPECT_FALSE(neighbor_list.build(0,
                                     1.0,
                                     std::vector<double>(),
                                     neighbor_count,
                                     neighbor_indices,
                                     max_neighbors,
                                     error));
    EXPECT_EQ(error, "invalid GPU neighbor-list input");
    EXPECT_TRUE(neighbor_count.empty());
    EXPECT_TRUE(neighbor_indices.empty());
}

TEST(GPUNeighborListTest, RejectsInvalidFilterInput)
{
    GPU_NeighborList neighbor_list;
    std::vector<int> neighbor_count;
    std::vector<int> neighbor_indices;
    std::string error;
    int max_neighbors = -1;
    const std::vector<double> position(3, 0.0);
    const std::vector<int> candidate_count(1, 0);

    EXPECT_FALSE(neighbor_list.filter(1,
                                      1.0,
                                      position,
                                      candidate_count,
                                      std::vector<int>(),
                                      1,
                                      neighbor_count,
                                      neighbor_indices,
                                      max_neighbors,
                                      error));
    EXPECT_EQ(error, "invalid GPU neighbor-list filter input");
}

TEST(GPUNeighborListTest, RejectsInvalidDeviceInput)
{
    GPU_NeighborList neighbor_list;
    std::string error;
    int max_neighbors = -1;

    EXPECT_FALSE(neighbor_list.build_device(0, 1.0, std::vector<double>(), max_neighbors, error));
    EXPECT_FALSE(neighbor_list.filter_device(0, 1.0, std::vector<double>(), 0, max_neighbors, error));
    EXPECT_EQ(neighbor_list.device_nall(), 0);
    EXPECT_EQ(neighbor_list.device_max_neighbors(), 0);
    EXPECT_EQ(neighbor_list.device_neighbor_count(), nullptr);
    EXPECT_EQ(neighbor_list.device_neighbor_indices(), nullptr);
}

TEST(GPUNeighborListTest, RejectsInvalidUploadInput)
{
    GPU_NeighborList neighbor_list;
    std::string error;

    EXPECT_FALSE(neighbor_list.upload(1, 1, std::vector<int>(), std::vector<int>(), error));
    EXPECT_EQ(neighbor_list.device_nall(), 0);
    EXPECT_EQ(neighbor_list.device_max_neighbors(), 0);
}
