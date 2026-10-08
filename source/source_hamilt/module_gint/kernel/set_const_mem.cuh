#ifndef ABACUS_SOURCE_HAMILT_MODULE_GINT_KERNEL_SET_CONST_MEM_CUH
#define ABACUS_SOURCE_HAMILT_MODULE_GINT_KERNEL_SET_CONST_MEM_CUH

#include <cuda_runtime.h>

namespace ModuleGint
{
__host__ void set_ylmcoe_d(const double* ylmcoe_h, double** ylmcoe_d_addr);
}

#endif // ABACUS_SOURCE_HAMILT_MODULE_GINT_KERNEL_SET_CONST_MEM_CUH
