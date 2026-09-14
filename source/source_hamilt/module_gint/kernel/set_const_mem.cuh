#pragma once
#include <cuda_runtime.h>

namespace ModuleGint
{
__host__ void set_ylmcoe_d(const double* ylmcoe_h, double** ylmcoe_d_addr);
}