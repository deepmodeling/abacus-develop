#pragma once

#include "source_base/kernels/cuda/sph_harm_gpu.cuh"

namespace ModuleGint
{
    // Import unified GPU sph_harm functions from ModuleBase
    using ModuleBase::sph_harm;
    using ModuleBase::grad_rl_sph_harm;
}
