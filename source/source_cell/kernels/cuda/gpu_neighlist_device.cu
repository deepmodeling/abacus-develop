#include "gpu_neighlist_cuda.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace gpu_neighlist_detail
{

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


} // namespace gpu_neighlist_detail

int gpu_build_neighbor_list_device(int nall,
                                       double cutoff,
                                       const double* position,
                                       void** device_neighbor_count,
                                       void** device_neighbor_indices,
                                       int* max_neighbors,
                                       char* message,
                                       int message_size)
{
    return gpu_neighlist_detail::build_device_impl(nall, cutoff, position, device_neighbor_count,
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
    return gpu_neighlist_detail::filter_device_impl(nall, cutoff, position, candidate_max_neighbors,
                              static_cast<const int*>(device_candidate_count),
                              static_cast<const int*>(device_candidate_indices),
                              device_neighbor_count, device_neighbor_indices,
                              max_neighbors, message, message_size);
}

