#pragma once
#include <cuda_runtime.h>

__constant__ double ylmcoe_d[100];

namespace ModuleGint
{
__host__ void set_ylmcoe_d(const double* ylmcoe_h, double** ylmcoe_d_addr);
}