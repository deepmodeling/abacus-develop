#ifndef NORM_PSI_KERNEL_CUH
#define NORM_PSI_KERNEL_CUH

#include <cuda_runtime.h>

#ifdef __CUBLASMP
#include <cuComplex.h>

namespace module_rt
{
namespace gpu
{

// Standard C++ wrapper to launch the normalization kernel
void launch_normalize_cij_kernel(cuDoubleComplex* d_Cij,
                                 int nrow,
                                 int ncol,
                                 int nb,
                                 int dim0,
                                 int dim1,
                                 int my_prow,
                                 int my_pcol,
                                 int nband,
                                 cudaStream_t stream);

} // namespace gpu
} // namespace module_rt
#endif // __CUBLASMP

#endif // NORM_PSI_KERNEL_CUH
