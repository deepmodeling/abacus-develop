#include "band_energy_kernel.cuh"

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

// Kernel to extract the real part of the diagonal elements of Eij
__global__ void extract_ekb_kernel(const cuDoubleComplex* d_Eij,
                                   double* d_eii,
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

    // Iterate over the allocated local buffer
    if (idx >= local_elems)
    {
        return;
    }

    // Column-major indexing using the formal Leading Dimension (LLD)
    int i = idx % lld;
    int j = idx / lld;

    int grow = get_global_index_dev(i, nb, dim0, my_prow);
    int gcol = get_global_index_dev(j, nb, dim1, my_pcol);

    // Filter out invalid blocks
    if (grow >= nband || gcol >= nband)
    {
        return;
    }

    // Extract the diagonal elements
    if (grow == gcol)
    {
        d_eii[grow] = cuCreal(d_Eij[idx]);
    }
}

// Wrapper implementation
void launch_extract_ekb_kernel(const cuDoubleComplex* d_Eij,
                               double* d_eii,
                               int lld,
                               int local_elems,
                               int nb,
                               int dim0,
                               int dim1,
                               int my_prow,
                               int my_pcol,
                               int nband,
                               cudaStream_t stream)
{
    if (local_elems > 0)
    {
        int threads_per_block = 256;
        int blocks_per_grid = (local_elems + threads_per_block - 1) / threads_per_block;

        extract_ekb_kernel<<<blocks_per_grid, threads_per_block, 0, stream>>>(d_Eij,
                                                                              d_eii,
                                                                              lld,
                                                                              local_elems,
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
