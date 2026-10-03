#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

namespace gpu_neighlist_detail
{

int set_message(char* message, int message_size, const char* text)
{
    if (message != nullptr && message_size > 0)
    {
        std::snprintf(message, message_size, "%s", text);
    }
    return 1;
}

int check_cuda(cudaError_t status, char* message, int message_size, const char* where)
{
    if (status == cudaSuccess)
    {
        return 0;
    }
    if (message != nullptr && message_size > 0)
    {
        std::snprintf(message, message_size, "%s: %s", where, cudaGetErrorString(status));
    }
    return 1;
}

__global__ void count_bins_kernel(int nall,
                                  const double* position,
                                  double x_min,
                                  double y_min,
                                  double z_min,
                                  double bin_size,
                                  int nbinx,
                                  int nbiny,
                                  int nbinz,
                                  int* bin_counts)
{
    const int atom = blockIdx.x * blockDim.x + threadIdx.x;
    if (atom >= nall)
    {
        return;
    }
    const int ix = min(max(static_cast<int>((position[atom] - x_min) / bin_size), 0), nbinx - 1);
    const int iy = min(max(static_cast<int>((position[nall + atom] - y_min) / bin_size), 0), nbiny - 1);
    const int iz = min(max(static_cast<int>((position[2 * nall + atom] - z_min) / bin_size), 0), nbinz - 1);
    const int bin = ix * nbiny * nbinz + iy * nbinz + iz;
    atomicAdd(&bin_counts[bin], 1);
}

__global__ void fill_bins_kernel(int nall,
                                 const double* position,
                                 double x_min,
                                 double y_min,
                                 double z_min,
                                 double bin_size,
                                 int nbinx,
                                 int nbiny,
                                 int nbinz,
                                 int* bin_cursor,
                                 int* bin_atoms)
{
    const int atom = blockIdx.x * blockDim.x + threadIdx.x;
    if (atom >= nall)
    {
        return;
    }
    const int ix = min(max(static_cast<int>((position[atom] - x_min) / bin_size), 0), nbinx - 1);
    const int iy = min(max(static_cast<int>((position[nall + atom] - y_min) / bin_size), 0), nbiny - 1);
    const int iz = min(max(static_cast<int>((position[2 * nall + atom] - z_min) / bin_size), 0), nbinz - 1);
    const int bin = ix * nbiny * nbinz + iy * nbinz + iz;
    const int slot = atomicAdd(&bin_cursor[bin], 1);
    bin_atoms[slot] = atom;
}

__global__ void count_neighbors_kernel(int nall,
                                       const double* position,
                                       double x_min,
                                       double y_min,
                                       double z_min,
                                       double bin_size,
                                       double cutoff2,
                                       int nbinx,
                                       int nbiny,
                                       int nbinz,
                                       const int* bin_offsets,
                                       const int* bin_atoms,
                                       int* neighbor_count)
{
    const int center = blockIdx.x * blockDim.x + threadIdx.x;
    if (center >= nall)
    {
        return;
    }
    const int ix = min(max(static_cast<int>((position[center] - x_min) / bin_size), 0), nbinx - 1);
    const int iy = min(max(static_cast<int>((position[nall + center] - y_min) / bin_size), 0), nbiny - 1);
    const int iz = min(max(static_cast<int>((position[2 * nall + center] - z_min) / bin_size), 0), nbinz - 1);
    const double x = position[center];
    const double y = position[nall + center];
    const double z = position[2 * nall + center];
    int count = 0;
    for (int dx = -1; dx <= 1; ++dx)
    {
        for (int dy = -1; dy <= 1; ++dy)
        {
            for (int dz = -1; dz <= 1; ++dz)
            {
                const int jx = ix + dx;
                const int jy = iy + dy;
                const int jz = iz + dz;
                if (jx < 0 || jx >= nbinx || jy < 0 || jy >= nbiny || jz < 0 || jz >= nbinz)
                {
                    continue;
                }
                const int bin = jx * nbiny * nbinz + jy * nbinz + jz;
                for (int slot = bin_offsets[bin]; slot < bin_offsets[bin + 1]; ++slot)
                {
                    const int neighbor = bin_atoms[slot];
                    if (neighbor == center)
                    {
                        continue;
                    }
                    const double dx12 = x - position[neighbor];
                    const double dy12 = y - position[nall + neighbor];
                    const double dz12 = z - position[2 * nall + neighbor];
                    if (dx12 * dx12 + dy12 * dy12 + dz12 * dz12 <= cutoff2)
                    {
                        ++count;
                    }
                }
            }
        }
    }
    neighbor_count[center] = count;
}

__global__ void fill_neighbors_kernel(int nall,
                                      int max_neighbors,
                                      const double* position,
                                      double x_min,
                                      double y_min,
                                      double z_min,
                                      double bin_size,
                                      double cutoff2,
                                      int nbinx,
                                      int nbiny,
                                      int nbinz,
                                      const int* bin_offsets,
                                      const int* bin_atoms,
                                      const int* neighbor_count,
                                      int* neighbor_indices)
{
    const int center = blockIdx.x * blockDim.x + threadIdx.x;
    if (center >= nall)
    {
        return;
    }
    const int ix = min(max(static_cast<int>((position[center] - x_min) / bin_size), 0), nbinx - 1);
    const int iy = min(max(static_cast<int>((position[nall + center] - y_min) / bin_size), 0), nbiny - 1);
    const int iz = min(max(static_cast<int>((position[2 * nall + center] - z_min) / bin_size), 0), nbinz - 1);
    const double x = position[center];
    const double y = position[nall + center];
    const double z = position[2 * nall + center];
    int count = 0;
    for (int dx = -1; dx <= 1; ++dx)
    {
        for (int dy = -1; dy <= 1; ++dy)
        {
            for (int dz = -1; dz <= 1; ++dz)
            {
                const int jx = ix + dx;
                const int jy = iy + dy;
                const int jz = iz + dz;
                if (jx < 0 || jx >= nbinx || jy < 0 || jy >= nbiny || jz < 0 || jz >= nbinz)
                {
                    continue;
                }
                const int bin = jx * nbiny * nbinz + jy * nbinz + jz;
                for (int slot = bin_offsets[bin]; slot < bin_offsets[bin + 1]; ++slot)
                {
                    const int neighbor = bin_atoms[slot];
                    if (neighbor == center)
                    {
                        continue;
                    }
                    const double dx12 = x - position[neighbor];
                    const double dy12 = y - position[nall + neighbor];
                    const double dz12 = z - position[2 * nall + neighbor];
                    if (dx12 * dx12 + dy12 * dy12 + dz12 * dz12 <= cutoff2)
                    {
                        neighbor_indices[center + nall * count] = neighbor;
                        ++count;
                    }
                }
            }
        }
    }
    static_cast<void>(neighbor_count);
}

__global__ void filter_neighbors_kernel(int nall,
                                        int candidate_max_neighbors,
                                        const double* position,
                                        double cutoff2,
                                        const int* candidate_count,
                                        const int* candidate_indices,
                                        int* neighbor_count,
                                        int* neighbor_indices)
{
    const int center = blockIdx.x * blockDim.x + threadIdx.x;
    if (center >= nall)
    {
        return;
    }
    const double x = position[center];
    const double y = position[nall + center];
    const double z = position[2 * nall + center];
    int count = 0;
    for (int slot = 0; slot < candidate_count[center]; ++slot)
    {
        const int neighbor = candidate_indices[center + nall * slot];
        const double dx = x - position[neighbor];
        const double dy = y - position[nall + neighbor];
        const double dz = z - position[2 * nall + neighbor];
        if (dx * dx + dy * dy + dz * dz <= cutoff2)
        {
            neighbor_indices[center + nall * count] = neighbor;
            ++count;
        }
    }
    neighbor_count[center] = count;
    static_cast<void>(candidate_max_neighbors);
}

__global__ void max_neighbor_count_kernel(int nall,
                                          const int* neighbor_count,
                                          int* max_neighbors)
{
    const int atom = blockIdx.x * blockDim.x + threadIdx.x;
    if (atom < nall)
    {
        atomicMax(max_neighbors, neighbor_count[atom]);
    }
}


} // namespace gpu_neighlist_detail

