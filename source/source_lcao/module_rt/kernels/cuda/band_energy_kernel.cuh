#ifndef BAND_ENERGY_KERNEL_CUH
#define BAND_ENERGY_KERNEL_CUH

#include <cuda_runtime.h>

#ifdef __CUBLASMP
#include <cuComplex.h>

namespace module_rt
{
namespace gpu
{

// Standard C++ wrapper to launch the diagonal extraction kernel
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
                               cudaStream_t stream);

} // namespace gpu
} // namespace module_rt
#endif // __CUBLASMP

#endif // BAND_ENERGY_KERNEL_CUH
