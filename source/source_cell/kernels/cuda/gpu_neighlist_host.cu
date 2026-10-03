#include "gpu_neighlist_cuda.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

namespace gpu_neighlist_detail
{

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


} // namespace gpu_neighlist_detail

int gpu_build_neighbor_list(int nall,
                                double cutoff,
                                const double* position,
                                int* neighbor_count,
                                int* neighbor_indices,
                                int* max_neighbors,
                                char* message,
                                int message_size)
{
    return gpu_neighlist_detail::build_impl(nall, cutoff, position, neighbor_count, neighbor_indices,
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
    return gpu_neighlist_detail::filter_impl(nall, cutoff, position, candidate_count, candidate_indices,
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
        return gpu_neighlist_detail::set_message(message, message_size, "invalid GPU neighbor-list upload arguments");
    }
    *device_neighbor_count = nullptr;
    *device_neighbor_indices = nullptr;
    int* d_count = nullptr;
    int* d_indices = nullptr;
    int status = gpu_neighlist_detail::check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_count), sizeof(int) * nall),
                            message, message_size, "cudaMalloc uploaded neighbor count");
    if (status != 0)
    {
        return status;
    }
    status = gpu_neighlist_detail::check_cuda(cudaMemcpy(d_count, neighbor_count, sizeof(int) * nall,
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
        status = gpu_neighlist_detail::check_cuda(cudaMalloc(reinterpret_cast<void**>(&d_indices),
                                       sizeof(int) * list_size),
                            message, message_size, "cudaMalloc uploaded neighbor indices");
        if (status != 0)
        {
            cudaFree(d_count);
            return status;
        }
        status = gpu_neighlist_detail::check_cuda(cudaMemcpy(d_indices, neighbor_indices, sizeof(int) * list_size,
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
