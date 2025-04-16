#include "set_const_mem.cuh"
#include "module_gint/kernels/cuda/cuda_tools.cuh"

namespace ModuleGint
{
    __host__ void set_ylmcoe_d(const double* ylmcoe_h, double** ylmcoe_d_addr)
    {
        checkCuda(cudaMemcpyToSymbol(ylmcoe_d, ylmcoe_h, sizeof(double) * 100));
        checkCuda(cudaGetSymbolAddress((void**)ylmcoe_d_addr, ylmcoe_d));
    }
}