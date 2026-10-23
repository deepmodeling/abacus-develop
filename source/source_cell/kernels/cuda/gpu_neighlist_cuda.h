#ifndef SOURCE_CELL_KERNELS_CUDA_GPU_NEIGHLIST_CUDA_H
#define SOURCE_CELL_KERNELS_CUDA_GPU_NEIGHLIST_CUDA_H

#include <cuda_runtime.h>

namespace gpu_neighlist_detail
{

int set_message(char* message, int message_size, const char* text);
int check_cuda(cudaError_t status, char* message, int message_size, const char* where);

__global__ void count_bins_kernel(int nall, const double* position, double x_min, double y_min,
                                  double z_min, double bin_size, int nbinx, int nbiny,
                                  int nbinz, int* bin_counts);
__global__ void fill_bins_kernel(int nall, const double* position, double x_min, double y_min,
                                 double z_min, double bin_size, int nbinx, int nbiny,
                                 int nbinz, int* bin_cursor, int* bin_atoms);
__global__ void count_neighbors_kernel(int nall, const double* position, double x_min,
                                       double y_min, double z_min, double bin_size,
                                       double cutoff2, int nbinx, int nbiny, int nbinz,
                                       const int* bin_offsets, const int* bin_atoms,
                                       int* neighbor_count);
__global__ void fill_neighbors_kernel(int nall, int max_neighbors, const double* position,
                                      double x_min, double y_min, double z_min, double bin_size,
                                      double cutoff2, int nbinx, int nbiny, int nbinz,
                                      const int* bin_offsets, const int* bin_atoms,
                                      const int* neighbor_count, int* neighbor_indices);
__global__ void filter_neighbors_kernel(int nall, int candidate_max_neighbors,
                                        const double* position, double cutoff2,
                                        const int* candidate_count, const int* candidate_indices,
                                        int* neighbor_count, int* neighbor_indices);
__global__ void max_neighbor_count_kernel(int nall, const int* neighbor_count, int* max_neighbors);

} // namespace gpu_neighlist_detail

#endif

