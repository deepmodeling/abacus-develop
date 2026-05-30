#include "esolver_nep_postprocess.h"

#include "source_base/module_device/device_check.h"
#include "source_base/module_device/kernel_compat.h"

#include <cuda_runtime.h>

namespace ModuleESolver
{
namespace
{

__global__ void nep_postprocess_kernel(const int nat,
                                       const double* atomic_energy,
                                       const double* raw_force,
                                       const double* raw_virial,
                                       const double fact_e,
                                       const double fact_f,
                                       const double fact_v,
                                       double* potential,
                                       double* force,
                                       double* virial)
{
    const int stride = blockDim.x * gridDim.x;
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < nat; i += stride)
    {
        atomicAdd(potential, atomic_energy[i] * fact_e);
        force[3 * i] = raw_force[i] * fact_f;
        force[3 * i + 1] = raw_force[i + nat] * fact_f;
        force[3 * i + 2] = raw_force[i + 2 * nat] * fact_f;

        for (int j = 0; j < 9; ++j)
        {
            atomicAdd(&virial[j], raw_virial[j * nat + i] * fact_v);
        }
    }
}

} // namespace

void postprocess_nep_cuda(const int nat,
                          const double* atomic_energy,
                          const double* raw_force,
                          const double* raw_virial,
                          const double fact_e,
                          const double fact_f,
                          const double fact_v,
                          double& potential,
                          ModuleBase::matrix& force,
                          ModuleBase::matrix& virial)
{
    double* d_energy = nullptr;
    double* d_raw_force = nullptr;
    double* d_raw_virial = nullptr;
    double* d_potential = nullptr;
    double* d_force = nullptr;
    double* d_virial = nullptr;

    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_energy), sizeof(double) * nat));
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_raw_force), sizeof(double) * 3 * nat));
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_raw_virial), sizeof(double) * 9 * nat));
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_potential), sizeof(double)));
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_force), sizeof(double) * 3 * nat));
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_virial), sizeof(double) * 9));

    CHECK_CUDA(cudaMemcpy(d_energy, atomic_energy, sizeof(double) * nat, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_raw_force, raw_force, sizeof(double) * 3 * nat, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_raw_virial, raw_virial, sizeof(double) * 9 * nat, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemset(d_potential, 0, sizeof(double)));
    CHECK_CUDA(cudaMemset(d_virial, 0, sizeof(double) * 9));

    const int block_size = 256;
    const int grid_size = (nat + block_size - 1) / block_size;
    nep_postprocess_kernel<<<grid_size, block_size>>>(nat,
                                                      d_energy,
                                                      d_raw_force,
                                                      d_raw_virial,
                                                      fact_e,
                                                      fact_f,
                                                      fact_v,
                                                      d_potential,
                                                      d_force,
                                                      d_virial);
    CHECK_LAST_CUDA_ERROR("nep_postprocess_kernel");
    CHECK_CUDA_SYNC();

    CHECK_CUDA(cudaMemcpy(&potential, d_potential, sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(force.c, d_force, sizeof(double) * 3 * nat, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(virial.c, d_virial, sizeof(double) * 9, cudaMemcpyDeviceToHost));

    CHECK_CUDA(cudaFree(d_energy));
    CHECK_CUDA(cudaFree(d_raw_force));
    CHECK_CUDA(cudaFree(d_raw_virial));
    CHECK_CUDA(cudaFree(d_potential));
    CHECK_CUDA(cudaFree(d_force));
    CHECK_CUDA(cudaFree(d_virial));
}

} // namespace ModuleESolver
