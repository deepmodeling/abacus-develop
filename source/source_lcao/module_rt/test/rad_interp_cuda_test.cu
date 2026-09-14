#include "source_lcao/module_rt/radial_interpolation.h"

#include <cuda_runtime.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace
{

__global__ void interpolate_queries(const double* radial_grid,
                                    const double* radial_values,
                                    const int mesh,
                                    const module_rt::RadialGridInfo grid_info,
                                    const double* queries,
                                    const int query_count,
                                    double* results)
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < query_count)
    {
        results[index]
            = module_rt::interpolate_radial(radial_grid, radial_values, mesh, grid_info, queries[index]);
    }
}

void compare_cpu_and_gpu(const std::vector<double>& grid,
                         const std::vector<double>& values,
                         const std::vector<double>& queries)
{
    module_rt::RadialGridInfo grid_info;
    ASSERT_TRUE(module_rt::analyze_radial_grid(grid.data(),
                                               values.data(),
                                               static_cast<int>(grid.size()),
                                               grid_info));

    double* grid_device = nullptr;
    double* values_device = nullptr;
    double* queries_device = nullptr;
    double* results_device = nullptr;
    ASSERT_EQ(cudaMalloc(&grid_device, grid.size() * sizeof(double)), cudaSuccess);
    ASSERT_EQ(cudaMalloc(&values_device, values.size() * sizeof(double)), cudaSuccess);
    ASSERT_EQ(cudaMalloc(&queries_device, queries.size() * sizeof(double)), cudaSuccess);
    ASSERT_EQ(cudaMalloc(&results_device, queries.size() * sizeof(double)), cudaSuccess);

    ASSERT_EQ(cudaMemcpy(grid_device, grid.data(), grid.size() * sizeof(double), cudaMemcpyHostToDevice), cudaSuccess);
    ASSERT_EQ(cudaMemcpy(values_device,
                         values.data(),
                         values.size() * sizeof(double),
                         cudaMemcpyHostToDevice),
              cudaSuccess);
    ASSERT_EQ(cudaMemcpy(queries_device,
                         queries.data(),
                         queries.size() * sizeof(double),
                         cudaMemcpyHostToDevice),
              cudaSuccess);

    constexpr int block_size = 128;
    const int block_count = (static_cast<int>(queries.size()) + block_size - 1) / block_size;
    interpolate_queries<<<block_count, block_size>>>(grid_device,
                                                     values_device,
                                                     static_cast<int>(grid.size()),
                                                     grid_info,
                                                     queries_device,
                                                     static_cast<int>(queries.size()),
                                                     results_device);
    ASSERT_EQ(cudaGetLastError(), cudaSuccess);
    ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);

    std::vector<double> gpu_results(queries.size());
    ASSERT_EQ(cudaMemcpy(gpu_results.data(),
                         results_device,
                         gpu_results.size() * sizeof(double),
                         cudaMemcpyDeviceToHost),
              cudaSuccess);

    double max_difference = 0.0;
    for (size_t i = 0; i < queries.size(); ++i)
    {
        const double cpu_result = module_rt::interpolate_radial(grid.data(),
                                                                values.data(),
                                                                static_cast<int>(grid.size()),
                                                                grid_info,
                                                                queries[i]);
        const double difference = std::abs(gpu_results[i] - cpu_result);
        max_difference = std::max(max_difference, difference);
        EXPECT_LE(difference, 5.0e-15) << "query = " << queries[i];
    }
    std::cout << std::scientific << std::setprecision(12) << "radial interpolation CPU/GPU mesh=" << grid.size()
              << ": max difference = " << max_difference << std::defaultfloat << std::endl;

    EXPECT_EQ(cudaFree(grid_device), cudaSuccess);
    EXPECT_EQ(cudaFree(values_device), cudaSuccess);
    EXPECT_EQ(cudaFree(queries_device), cudaSuccess);
    EXPECT_EQ(cudaFree(results_device), cudaSuccess);
}

} // namespace

TEST(RadialInterpolationCuda, MatchesCpuOnUniformAndNonuniformGrids)
{
    int device_count = 0;
    const cudaError_t device_status = cudaGetDeviceCount(&device_count);
    if (device_status != cudaSuccess)
    {
        GTEST_SKIP() << "cudaGetDeviceCount failed: " << cudaGetErrorString(device_status);
    }
    if (device_count == 0)
    {
        GTEST_SKIP() << "cudaGetDeviceCount reported zero CUDA devices";
    }

    const std::vector<double> uniform_grid = {0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25};
    std::vector<double> uniform_values;
    for (const double radius: uniform_grid)
    {
        uniform_values.push_back(std::sin(radius) + radius * radius);
    }
    compare_cpu_and_gpu(uniform_grid,
                        uniform_values,
                        {0.25,
                         0.5,
                         2.1,
                         2.6,
                         3.1,
                         3.25,
                         3.5,
                         std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::infinity()});

    const std::vector<double> nonuniform_grid = {0.2, 0.201, 0.35, 1.4, 1.45, 4.0, 9.0};
    std::vector<double> nonuniform_values;
    for (const double radius: nonuniform_grid)
    {
        nonuniform_values.push_back(std::cos(radius) - 0.25 * radius);
    }
    compare_cpu_and_gpu(nonuniform_grid, nonuniform_values, {0.2, 0.2005, 1.43, 3.2, 8.5, 9.0, 9.1});
}
