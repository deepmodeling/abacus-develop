#include "norm_psi_kernel.cuh"

#include <math_constants.h>

#ifdef __CUBLASMP
namespace module_rt
{
namespace gpu
{

// Device function for global index mapping
__device__ inline int get_global_index_dev(int local_idx, int block_size, int num_procs, int proc_coord)
{
    return (local_idx / block_size) * (num_procs * block_size) + proc_coord * block_size + (local_idx % block_size);
}

// CUDA kernel to normalize the Cij matrix directly on the GPU
__global__ void normalize_cij_kernel(cuDoubleComplex* d_Cij,
                                     int lld,
                                     int local_elems,
                                     int nb,
                                     int dim0,
                                     int dim1,
                                     int my_prow,
                                     int my_pcol,
                                     int nband)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // We iterate over the entire allocated buffer (local_elems)
    if (idx >= local_elems)
    {
        return;
    }

    // Column-major indexing using the formal Leading Dimension (LLD)
    int i = idx % lld;
    int j = idx / lld;

    int grow = get_global_index_dev(i, nb, dim0, my_prow);
    int gcol = get_global_index_dev(j, nb, dim1, my_pcol);

    // Filter out the "empty" spaces that are outside the nband x nband logic
    if (grow >= nband || gcol >= nband)
    {
        return;
    }

    if (grow == gcol)
    {
        double val = cuCreal(d_Cij[idx]);
        if (val < 1e-12)
            val = 1e-12;
        d_Cij[idx] = make_cuDoubleComplex(1.0 / sqrt(val), 0.0);
    }
    else
    {
        d_Cij[idx] = make_cuDoubleComplex(0.0, 0.0);
    }
}

// Wrapper function implementation
void launch_normalize_cij_kernel(cuDoubleComplex* d_Cij,
                                 int nrow,
                                 int ncol,
                                 int nb,
                                 int dim0,
                                 int dim1,
                                 int my_prow,
                                 int my_pcol,
                                 int nband,
                                 cudaStream_t stream)
{
    int total_elems = nrow * ncol;
    if (total_elems > 0)
    {
        int threads_per_block = 256;
        int blocks_per_grid = (total_elems + threads_per_block - 1) / threads_per_block;

        normalize_cij_kernel<<<blocks_per_grid, threads_per_block, 0, stream>>>(d_Cij,
                                                                                nrow,
                                                                                ncol,
                                                                                nb,
                                                                                dim0,
                                                                                dim1,
                                                                                my_prow,
                                                                                my_pcol,
                                                                                nband);
    }
}

} // namespace gpu
} // namespace module_rt
#endif // __CUBLASMP
