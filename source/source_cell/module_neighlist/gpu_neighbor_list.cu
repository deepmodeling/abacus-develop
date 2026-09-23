#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

namespace
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

int build_impl(int nall,
               double cutoff,
               const double* position,
               int* neighbor_count,
               int* neighbor_indices,
               int* max_neighbors,
               char* message,
               int message_size)
{
    if (nall <= 0 || cutoff <= 0.0 || position == nullptr || neighbor_count == nullptr
        || max_neighbors == nullptr)
    {
        return set_message(message, message_size, "invalid GPU neighbor-list arguments");
    }

    double x_min = position[0];
    double x_max = position[0];
    double y_min = position[nall];
    double y_max = position[nall];
    double z_min = position[2 * nall];
    double z_max = position[2 * nall];
    for (int i = 1; i < nall; ++i)
    {
        x_min = std::min(x_min, position[i]);
        x_max = std::max(x_max, position[i]);
        y_min = std::min(y_min, position[nall + i]);
        y_max = std::max(y_max, position[nall + i]);
        z_min = std::min(z_min, position[2 * nall + i]);
        z_max = std::max(z_max, position[2 * nall + i]);
    }

    const double bin_size = cutoff;
    const int nbinx = std::max(1, static_cast<int>(std::ceil((x_max - x_min) / bin_size)) + 1);
    const int nbiny = std::max(1, static_cast<int>(std::ceil((y_max - y_min) / bin_size)) + 1);
    const int nbinz = std::max(1, static_cast<int>(std::ceil((z_max - z_min) / bin_size)) + 1);
    const long long total_bins_ll = static_cast<long long>(nbinx) * nbiny * nbinz;
    if (total_bins_ll > std::numeric_limits<int>::max())
    {
        return set_message(message, message_size, "GPU neighbor-list bin count exceeds int range");
    }
    const int total_bins = static_cast<int>(total_bins_ll);
    const int threads = 256;
    const int blocks = (nall + threads - 1) / threads;

    double* d_position = nullptr;
    int* d_bin_counts = nullptr;
    int* d_bin_offsets = nullptr;
    int* d_bin_cursor = nullptr;
    int* d_bin_atoms = nullptr;
    int* d_neighbor_count = nullptr;
    int* d_neighbor_indices = nullptr;
    std::vector<int> bin_counts(static_cast<std::size_t>(total_bins), 0);
    std::vector<int> bin_offsets(static_cast<std::size_t>(total_bins + 1), 0);
    int status = 0;
    status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_position), sizeof(double) * 3 * nall),
                        message, message_size, "cudaMalloc position");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_bin_counts), sizeof(int) * total_bins),
                        message, message_size, "cudaMalloc bin counts");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_bin_offsets), sizeof(int) * (total_bins + 1)),
                        message, message_size, "cudaMalloc bin offsets");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_bin_cursor), sizeof(int) * total_bins),
                        message, message_size, "cudaMalloc bin cursor");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_bin_atoms), sizeof(int) * nall),
                        message, message_size, "cudaMalloc bin atoms");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_neighbor_count), sizeof(int) * nall),
                        message, message_size, "cudaMalloc neighbor count");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMemcpy(d_position, position, sizeof(double) * 3 * nall, cudaMemcpyHostToDevice),
                        message, message_size, "cudaMemcpy position");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMemset(d_bin_counts, 0, sizeof(int) * total_bins),
                        message, message_size, "cudaMemset bin counts");
    if (status != 0) goto cleanup;
    count_bins_kernel<<<blocks, threads>>>(nall, d_position, x_min, y_min, z_min, bin_size,
                                           nbinx, nbiny, nbinz, d_bin_counts);
    status = check_cuda(cudaGetLastError(), message, message_size, "count bins kernel");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaDeviceSynchronize(), message, message_size, "count bins synchronize");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMemcpy(bin_counts.data(), d_bin_counts, sizeof(int) * total_bins,
                                   cudaMemcpyDeviceToHost),
                        message, message_size, "copy bin counts");
    if (status != 0) goto cleanup;
    for (int bin = 0; bin < total_bins; ++bin)
    {
        bin_offsets[static_cast<std::size_t>(bin + 1)]
            = bin_offsets[static_cast<std::size_t>(bin)] + bin_counts[static_cast<std::size_t>(bin)];
    }
    status = check_cuda(cudaMemcpy(d_bin_offsets, bin_offsets.data(), sizeof(int) * (total_bins + 1),
                                   cudaMemcpyHostToDevice),
                        message, message_size, "copy bin offsets");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMemcpy(d_bin_cursor, bin_offsets.data(), sizeof(int) * total_bins,
                                   cudaMemcpyHostToDevice),
                        message, message_size, "copy bin cursor");
    if (status != 0) goto cleanup;
    fill_bins_kernel<<<blocks, threads>>>(nall, d_position, x_min, y_min, z_min, bin_size,
                                          nbinx, nbiny, nbinz, d_bin_cursor, d_bin_atoms);
    status = check_cuda(cudaGetLastError(), message, message_size, "fill bins kernel");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaDeviceSynchronize(), message, message_size, "fill bins synchronize");
    if (status != 0) goto cleanup;
    count_neighbors_kernel<<<blocks, threads>>>(nall, d_position, x_min, y_min, z_min, bin_size,
                                                cutoff * cutoff, nbinx, nbiny, nbinz,
                                                d_bin_offsets, d_bin_atoms, d_neighbor_count);
    status = check_cuda(cudaGetLastError(), message, message_size, "count neighbors kernel");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaDeviceSynchronize(), message, message_size, "count neighbors synchronize");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMemcpy(neighbor_count, d_neighbor_count, sizeof(int) * nall,
                                   cudaMemcpyDeviceToHost),
                        message, message_size, "copy neighbor counts");
    if (status != 0) goto cleanup;
    *max_neighbors = 0;
    for (int i = 0; i < nall; ++i)
    {
        *max_neighbors = std::max(*max_neighbors, neighbor_count[i]);
    }
    if (neighbor_indices != nullptr && *max_neighbors > 0)
    {
        const std::size_t list_size = static_cast<std::size_t>(nall)
                                      * static_cast<std::size_t>(*max_neighbors);
        status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_neighbor_indices),
                                       sizeof(int) * list_size),
                            message, message_size, "cudaMalloc neighbor indices");
        if (status != 0) goto cleanup;
        fill_neighbors_kernel<<<blocks, threads>>>(nall, *max_neighbors, d_position,
                                                   x_min, y_min, z_min, bin_size,
                                                   cutoff * cutoff, nbinx, nbiny, nbinz,
                                                   d_bin_offsets, d_bin_atoms, d_neighbor_count,
                                                   d_neighbor_indices);
        status = check_cuda(cudaGetLastError(), message, message_size, "fill neighbors kernel");
        if (status != 0) goto cleanup;
        status = check_cuda(cudaDeviceSynchronize(), message, message_size, "fill neighbors synchronize");
        if (status != 0) goto cleanup;
        status = check_cuda(cudaMemcpy(neighbor_indices, d_neighbor_indices, sizeof(int) * list_size,
                                       cudaMemcpyDeviceToHost),
                            message, message_size, "copy neighbor indices");
    }

cleanup:
    cudaFree(d_position);
    cudaFree(d_bin_counts);
    cudaFree(d_bin_offsets);
    cudaFree(d_bin_cursor);
    cudaFree(d_bin_atoms);
    cudaFree(d_neighbor_count);
    cudaFree(d_neighbor_indices);
    return status;
}

int filter_impl(int nall,
                double cutoff,
                const double* position,
                const int* candidate_count,
                const int* candidate_indices,
                int candidate_max_neighbors,
                int* neighbor_count,
                int* neighbor_indices,
                int* max_neighbors,
                char* message,
                int message_size)
{
    if (nall <= 0 || cutoff <= 0.0 || position == nullptr || candidate_count == nullptr
        || max_neighbors == nullptr || candidate_max_neighbors < 0
        || (candidate_max_neighbors > 0 && (candidate_indices == nullptr || neighbor_indices == nullptr)))
    {
        return set_message(message, message_size, "invalid GPU neighbor-list filter arguments");
    }

    const int threads = 256;
    const int blocks = (nall + threads - 1) / threads;
    double* d_position = nullptr;
    int* d_candidate_count = nullptr;
    int* d_candidate_indices = nullptr;
    int* d_neighbor_count = nullptr;
    int* d_neighbor_indices = nullptr;
    int status = 0;
    status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_position), sizeof(double) * 3 * nall),
                        message, message_size, "cudaMalloc filter position");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_candidate_count), sizeof(int) * nall),
                        message, message_size, "cudaMalloc candidate count");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_neighbor_count), sizeof(int) * nall),
                        message, message_size, "cudaMalloc filtered count");
    if (status != 0) goto cleanup;
    if (candidate_max_neighbors > 0)
    {
        const std::size_t list_size = static_cast<std::size_t>(nall)
                                      * static_cast<std::size_t>(candidate_max_neighbors);
        status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_candidate_indices),
                                       sizeof(int) * list_size),
                            message, message_size, "cudaMalloc candidate indices");
        if (status != 0) goto cleanup;
        status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_neighbor_indices),
                                       sizeof(int) * list_size),
                            message, message_size, "cudaMalloc filtered indices");
        if (status != 0) goto cleanup;
        status = check_cuda(cudaMemcpy(d_candidate_indices,
                                       candidate_indices,
                                       sizeof(int) * list_size,
                                       cudaMemcpyHostToDevice),
                            message, message_size, "copy candidate indices");
        if (status != 0) goto cleanup;
    }
    status = check_cuda(cudaMemcpy(d_position, position, sizeof(double) * 3 * nall,
                                   cudaMemcpyHostToDevice),
                        message, message_size, "copy filter position");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMemcpy(d_candidate_count, candidate_count, sizeof(int) * nall,
                                   cudaMemcpyHostToDevice),
                        message, message_size, "copy candidate count");
    if (status != 0) goto cleanup;
    filter_neighbors_kernel<<<blocks, threads>>>(nall,
                                                 candidate_max_neighbors,
                                                 d_position,
                                                 cutoff * cutoff,
                                                 d_candidate_count,
                                                 d_candidate_indices,
                                                 d_neighbor_count,
                                                 d_neighbor_indices);
    status = check_cuda(cudaGetLastError(), message, message_size, "filter neighbors kernel");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaDeviceSynchronize(), message, message_size, "filter neighbors synchronize");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMemcpy(neighbor_count, d_neighbor_count, sizeof(int) * nall,
                                   cudaMemcpyDeviceToHost),
                        message, message_size, "copy filtered count");
    if (status != 0) goto cleanup;
    *max_neighbors = 0;
    for (int i = 0; i < nall; ++i)
    {
        *max_neighbors = std::max(*max_neighbors, neighbor_count[i]);
    }
    if (candidate_max_neighbors > 0)
    {
        const std::size_t list_size = static_cast<std::size_t>(nall)
                                      * static_cast<std::size_t>(candidate_max_neighbors);
        status = check_cuda(cudaMemcpy(neighbor_indices,
                                       d_neighbor_indices,
                                       sizeof(int) * list_size,
                                       cudaMemcpyDeviceToHost),
                            message, message_size, "copy filtered indices");
    }

cleanup:
    cudaFree(d_position);
    cudaFree(d_candidate_count);
    cudaFree(d_candidate_indices);
    cudaFree(d_neighbor_count);
    cudaFree(d_neighbor_indices);
    return status;
}

int build_device_impl(int nall,
                      double cutoff,
                      const double* position,
                      void** device_neighbor_count,
                      void** device_neighbor_indices,
                      int* max_neighbors,
                      char* message,
                      int message_size)
{
    if (nall <= 0 || cutoff <= 0.0 || position == nullptr
        || device_neighbor_count == nullptr || device_neighbor_indices == nullptr
        || max_neighbors == nullptr)
    {
        return set_message(message, message_size, "invalid GPU device neighbor-list arguments");
    }
    *device_neighbor_count = nullptr;
    *device_neighbor_indices = nullptr;
    *max_neighbors = 0;

    double x_min = position[0];
    double x_max = position[0];
    double y_min = position[nall];
    double y_max = position[nall];
    double z_min = position[2 * nall];
    double z_max = position[2 * nall];
    for (int i = 1; i < nall; ++i)
    {
        x_min = std::min(x_min, position[i]);
        x_max = std::max(x_max, position[i]);
        y_min = std::min(y_min, position[nall + i]);
        y_max = std::max(y_max, position[nall + i]);
        z_min = std::min(z_min, position[2 * nall + i]);
        z_max = std::max(z_max, position[2 * nall + i]);
    }
    const double bin_size = cutoff;
    const int nbinx = std::max(1, static_cast<int>(std::ceil((x_max - x_min) / bin_size)) + 1);
    const int nbiny = std::max(1, static_cast<int>(std::ceil((y_max - y_min) / bin_size)) + 1);
    const int nbinz = std::max(1, static_cast<int>(std::ceil((z_max - z_min) / bin_size)) + 1);
    const long long total_bins_ll = static_cast<long long>(nbinx) * nbiny * nbinz;
    if (total_bins_ll > std::numeric_limits<int>::max())
    {
        return set_message(message, message_size, "GPU neighbor-list bin count exceeds int range");
    }
    const int total_bins = static_cast<int>(total_bins_ll);
    const int threads = 256;
    const int blocks = (nall + threads - 1) / threads;
    double* d_position = nullptr;
    int* d_bin_counts = nullptr;
    int* d_bin_offsets = nullptr;
    int* d_bin_cursor = nullptr;
    int* d_bin_atoms = nullptr;
    int* d_neighbor_count = nullptr;
    int* d_neighbor_indices = nullptr;
    int* d_max_neighbors = nullptr;
    std::vector<int> bin_counts(static_cast<std::size_t>(total_bins), 0);
    std::vector<int> bin_offsets(static_cast<std::size_t>(total_bins + 1), 0);
    int status = 0;
    status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_position), sizeof(double) * 3 * nall),
                        message, message_size, "cudaMalloc device neighbor position");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_bin_counts), sizeof(int) * total_bins),
                        message, message_size, "cudaMalloc device bin counts");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_bin_offsets),
                                   sizeof(int) * (total_bins + 1)),
                        message, message_size, "cudaMalloc device bin offsets");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_bin_cursor), sizeof(int) * total_bins),
                        message, message_size, "cudaMalloc device bin cursor");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_bin_atoms), sizeof(int) * nall),
                        message, message_size, "cudaMalloc device bin atoms");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_neighbor_count), sizeof(int) * nall),
                        message, message_size, "cudaMalloc device neighbor count");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_max_neighbors), sizeof(int)),
                        message, message_size, "cudaMalloc device max neighbor count");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMemcpy(d_position, position, sizeof(double) * 3 * nall,
                                   cudaMemcpyHostToDevice),
                        message, message_size, "copy device neighbor position");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMemset(d_bin_counts, 0, sizeof(int) * total_bins),
                        message, message_size, "clear device bin counts");
    if (status != 0) goto cleanup;
    count_bins_kernel<<<blocks, threads>>>(nall, d_position, x_min, y_min, z_min, bin_size,
                                           nbinx, nbiny, nbinz, d_bin_counts);
    status = check_cuda(cudaGetLastError(), message, message_size, "count device bins");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaDeviceSynchronize(), message, message_size, "synchronize device bins");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMemcpy(bin_counts.data(), d_bin_counts, sizeof(int) * total_bins,
                                   cudaMemcpyDeviceToHost),
                        message, message_size, "copy device bin counts");
    if (status != 0) goto cleanup;
    for (int bin = 0; bin < total_bins; ++bin)
    {
        bin_offsets[static_cast<std::size_t>(bin + 1)]
            = bin_offsets[static_cast<std::size_t>(bin)] + bin_counts[static_cast<std::size_t>(bin)];
    }
    status = check_cuda(cudaMemcpy(d_bin_offsets, bin_offsets.data(),
                                   sizeof(int) * (total_bins + 1), cudaMemcpyHostToDevice),
                        message, message_size, "copy device bin offsets");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMemcpy(d_bin_cursor, bin_offsets.data(), sizeof(int) * total_bins,
                                   cudaMemcpyHostToDevice),
                        message, message_size, "copy device bin cursor");
    if (status != 0) goto cleanup;
    fill_bins_kernel<<<blocks, threads>>>(nall, d_position, x_min, y_min, z_min, bin_size,
                                          nbinx, nbiny, nbinz, d_bin_cursor, d_bin_atoms);
    status = check_cuda(cudaGetLastError(), message, message_size, "fill device bins");
    if (status != 0) goto cleanup;
    count_neighbors_kernel<<<blocks, threads>>>(nall, d_position, x_min, y_min, z_min, bin_size,
                                                cutoff * cutoff, nbinx, nbiny, nbinz,
                                                d_bin_offsets, d_bin_atoms, d_neighbor_count);
    status = check_cuda(cudaGetLastError(), message, message_size, "count device neighbors");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMemset(d_max_neighbors, 0, sizeof(int)),
                        message, message_size, "clear device max neighbors");
    if (status != 0) goto cleanup;
    max_neighbor_count_kernel<<<blocks, threads>>>(nall, d_neighbor_count, d_max_neighbors);
    status = check_cuda(cudaGetLastError(), message, message_size, "reduce device max neighbors");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMemcpy(max_neighbors, d_max_neighbors, sizeof(int), cudaMemcpyDeviceToHost),
                        message, message_size, "copy device max neighbors");
    if (status != 0) goto cleanup;
    if (*max_neighbors > 0)
    {
        const std::size_t list_size = static_cast<std::size_t>(nall)
                                      * static_cast<std::size_t>(*max_neighbors);
        status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_neighbor_indices),
                                       sizeof(int) * list_size),
                            message, message_size, "cudaMalloc device neighbor indices");
        if (status != 0) goto cleanup;
        fill_neighbors_kernel<<<blocks, threads>>>(nall, *max_neighbors, d_position,
                                                   x_min, y_min, z_min, bin_size,
                                                   cutoff * cutoff, nbinx, nbiny, nbinz,
                                                   d_bin_offsets, d_bin_atoms, d_neighbor_count,
                                                   d_neighbor_indices);
        status = check_cuda(cudaGetLastError(), message, message_size, "fill device neighbors");
        if (status != 0) goto cleanup;
    }
    *device_neighbor_count = d_neighbor_count;
    *device_neighbor_indices = d_neighbor_indices;
    d_neighbor_count = nullptr;
    d_neighbor_indices = nullptr;

cleanup:
    cudaFree(d_position);
    cudaFree(d_bin_counts);
    cudaFree(d_bin_offsets);
    cudaFree(d_bin_cursor);
    cudaFree(d_bin_atoms);
    cudaFree(d_neighbor_count);
    cudaFree(d_neighbor_indices);
    cudaFree(d_max_neighbors);
    return status;
}

int filter_device_impl(int nall,
                       double cutoff,
                       const double* position,
                       int candidate_max_neighbors,
                       const int* device_candidate_count,
                       const int* device_candidate_indices,
                       void** device_neighbor_count,
                       void** device_neighbor_indices,
                       int* max_neighbors,
                       char* message,
                       int message_size)
{
    if (nall <= 0 || cutoff <= 0.0 || position == nullptr
        || candidate_max_neighbors < 0 || device_candidate_count == nullptr
        || device_neighbor_count == nullptr || device_neighbor_indices == nullptr
        || (candidate_max_neighbors > 0 && device_candidate_indices == nullptr)
        || max_neighbors == nullptr)
    {
        return set_message(message, message_size, "invalid GPU device filter arguments");
    }
    *device_neighbor_count = nullptr;
    *device_neighbor_indices = nullptr;
    *max_neighbors = 0;
    const int threads = 256;
    const int blocks = (nall + threads - 1) / threads;
    double* d_position = nullptr;
    int* d_neighbor_count = nullptr;
    int* d_neighbor_indices = nullptr;
    int* d_max_neighbors = nullptr;
    int status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_position), sizeof(double) * 3 * nall),
                            message, message_size, "cudaMalloc filter device position");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_neighbor_count), sizeof(int) * nall),
                        message, message_size, "cudaMalloc filtered device count");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_max_neighbors), sizeof(int)),
                        message, message_size, "cudaMalloc filtered max neighbors");
    if (status != 0) goto cleanup;
    if (candidate_max_neighbors > 0)
    {
        const std::size_t list_size = static_cast<std::size_t>(nall)
                                      * static_cast<std::size_t>(candidate_max_neighbors);
        status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_neighbor_indices),
                                       sizeof(int) * list_size),
                            message, message_size, "cudaMalloc filtered device indices");
        if (status != 0) goto cleanup;
    }
    status = check_cuda(cudaMemcpy(d_position, position, sizeof(double) * 3 * nall,
                                   cudaMemcpyHostToDevice),
                        message, message_size, "copy filter device position");
    if (status != 0) goto cleanup;
    filter_neighbors_kernel<<<blocks, threads>>>(nall,
                                                 candidate_max_neighbors,
                                                 d_position,
                                                 cutoff * cutoff,
                                                 device_candidate_count,
                                                 device_candidate_indices,
                                                 d_neighbor_count,
                                                 d_neighbor_indices);
    status = check_cuda(cudaGetLastError(), message, message_size, "filter device neighbors");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMemset(d_max_neighbors, 0, sizeof(int)),
                        message, message_size, "clear filtered max neighbors");
    if (status != 0) goto cleanup;
    max_neighbor_count_kernel<<<blocks, threads>>>(nall, d_neighbor_count, d_max_neighbors);
    status = check_cuda(cudaGetLastError(), message, message_size, "reduce filtered max neighbors");
    if (status != 0) goto cleanup;
    status = check_cuda(cudaMemcpy(max_neighbors, d_max_neighbors, sizeof(int), cudaMemcpyDeviceToHost),
                        message, message_size, "copy filtered max neighbors");
    if (status != 0) goto cleanup;
    // The filter kernel writes with candidate_max_neighbors as its row stride.
    // Keep that stride for the NEP kernels; the actual count remains in d_neighbor_count.
    *max_neighbors = candidate_max_neighbors;
    *device_neighbor_count = d_neighbor_count;
    *device_neighbor_indices = d_neighbor_indices;
    d_neighbor_count = nullptr;
    d_neighbor_indices = nullptr;

cleanup:
    cudaFree(d_position);
    cudaFree(d_neighbor_count);
    cudaFree(d_neighbor_indices);
    cudaFree(d_max_neighbors);
    return status;
}

} // namespace

int gpu_build_neighbor_list(int nall,
                                double cutoff,
                                const double* position,
                                int* neighbor_count,
                                int* neighbor_indices,
                                int* max_neighbors,
                                char* message,
                                int message_size)
{
    return build_impl(nall, cutoff, position, neighbor_count, neighbor_indices,
                      max_neighbors, message, message_size);
}

int gpu_filter_neighbor_list(int nall,
                                 double cutoff,
                                 const double* position,
                                 const int* candidate_count,
                                 const int* candidate_indices,
                                 int candidate_max_neighbors,
                                 int* neighbor_count,
                                 int* neighbor_indices,
                                 int* max_neighbors,
                                 char* message,
                                 int message_size)
{
    return filter_impl(nall, cutoff, position, candidate_count, candidate_indices,
                       candidate_max_neighbors, neighbor_count, neighbor_indices,
                       max_neighbors, message, message_size);
}

int gpu_upload_neighbor_list(int nall,
                                 int max_neighbors,
                                 const int* neighbor_count,
                                 const int* neighbor_indices,
                                 void** device_neighbor_count,
                                 void** device_neighbor_indices,
                                 char* message,
                                 int message_size)
{
    if (nall <= 0 || max_neighbors < 0 || neighbor_count == nullptr
        || device_neighbor_count == nullptr || device_neighbor_indices == nullptr
        || (max_neighbors > 0 && neighbor_indices == nullptr))
    {
        return set_message(message, message_size, "invalid GPU neighbor-list upload arguments");
    }
    *device_neighbor_count = nullptr;
    *device_neighbor_indices = nullptr;
    int* d_count = nullptr;
    int* d_indices = nullptr;
    int status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_count), sizeof(int) * nall),
                            message, message_size, "cudaMalloc uploaded neighbor count");
    if (status != 0)
    {
        return status;
    }
    status = check_cuda(cudaMemcpy(d_count, neighbor_count, sizeof(int) * nall,
                                   cudaMemcpyHostToDevice),
                        message, message_size, "copy uploaded neighbor count");
    if (status != 0)
    {
        cudaFree(d_count);
        return status;
    }
    if (max_neighbors > 0)
    {
        const std::size_t list_size = static_cast<std::size_t>(nall)
                                      * static_cast<std::size_t>(max_neighbors);
        status = check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_indices),
                                       sizeof(int) * list_size),
                            message, message_size, "cudaMalloc uploaded neighbor indices");
        if (status != 0)
        {
            cudaFree(d_count);
            return status;
        }
        status = check_cuda(cudaMemcpy(d_indices, neighbor_indices, sizeof(int) * list_size,
                                       cudaMemcpyHostToDevice),
                            message, message_size, "copy uploaded neighbor indices");
        if (status != 0)
        {
            cudaFree(d_count);
            cudaFree(d_indices);
            return status;
        }
    }
    *device_neighbor_count = d_count;
    *device_neighbor_indices = d_indices;
    return 0;
}

void gpu_release_neighbor_list(void* device_neighbor_count,
                                   void* device_neighbor_indices)
{
    cudaFree(device_neighbor_count);
    cudaFree(device_neighbor_indices);
}

int gpu_build_neighbor_list_device(int nall,
                                       double cutoff,
                                       const double* position,
                                       void** device_neighbor_count,
                                       void** device_neighbor_indices,
                                       int* max_neighbors,
                                       char* message,
                                       int message_size)
{
    return build_device_impl(nall, cutoff, position, device_neighbor_count,
                             device_neighbor_indices, max_neighbors, message, message_size);
}

int gpu_filter_neighbor_list_device(int nall,
                                        double cutoff,
                                        const double* position,
                                        int candidate_max_neighbors,
                                        const void* device_candidate_count,
                                        const void* device_candidate_indices,
                                        void** device_neighbor_count,
                                        void** device_neighbor_indices,
                                        int* max_neighbors,
                                        char* message,
                                        int message_size)
{
    return filter_device_impl(nall, cutoff, position, candidate_max_neighbors,
                              static_cast<const int*>(device_candidate_count),
                              static_cast<const int*>(device_candidate_indices),
                              device_neighbor_count, device_neighbor_indices,
                              max_neighbors, message, message_size);
}
