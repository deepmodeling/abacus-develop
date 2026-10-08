#ifndef ABACUS_SOURCE_HAMILT_MODULE_GINT_KERNEL_SPH_CUH
#define ABACUS_SOURCE_HAMILT_MODULE_GINT_KERNEL_SPH_CUH


#include "source_base/kernels/cuda/sph_harm_gpu.cuh"

namespace ModuleGint
{
    // Import unified GPU sph_harm functions from ModuleBase
    using ModuleBase::sph_harm;
    using ModuleBase::grad_rl_sph_harm;
}

#endif // ABACUS_SOURCE_HAMILT_MODULE_GINT_KERNEL_SPH_CUH
