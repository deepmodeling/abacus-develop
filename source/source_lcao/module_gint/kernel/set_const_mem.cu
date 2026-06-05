#include "set_const_mem.cuh"
#include "source_base/module_device/device_check.h"

__constant__ double ylmcoe_d[100];

namespace ModuleGint
{
    __host__ void set_ylmcoe_d(const double* ylmcoe_h, double** ylmcoe_d_addr)
    {
        CHECK_CUDA(cudaMemcpyToSymbol(ylmcoe_d, ylmcoe_h, sizeof(double) * 100));
        CHECK_CUDA(cudaGetSymbolAddress((void**)ylmcoe_d_addr, ylmcoe_d));
    }
}