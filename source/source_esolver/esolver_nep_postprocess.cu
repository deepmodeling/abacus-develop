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

void init_nep_cuda_postprocess_workspace(NepCudaPostprocessWorkspace& workspace, const int nat)
{
    if (workspace.capacity >= nat)
    {
        return;
    }

    release_nep_cuda_postprocess_workspace(workspace);

    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&workspace.energy), sizeof(double) * nat));
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&workspace.raw_force), sizeof(double) * 3 * nat));
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&workspace.raw_virial), sizeof(double) * 9 * nat));
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&workspace.potential), sizeof(double)));
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&workspace.force), sizeof(double) * 3 * nat));
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&workspace.virial), sizeof(double) * 9));
    workspace.capacity = nat;
}

void release_nep_cuda_postprocess_workspace(NepCudaPostprocessWorkspace& workspace)
{
    if (workspace.energy != nullptr)
    {
        CHECK_CUDA(cudaFree(workspace.energy));
    }
    if (workspace.raw_force != nullptr)
    {
        CHECK_CUDA(cudaFree(workspace.raw_force));
    }
    if (workspace.raw_virial != nullptr)
    {
        CHECK_CUDA(cudaFree(workspace.raw_virial));
    }
    if (workspace.potential != nullptr)
    {
        CHECK_CUDA(cudaFree(workspace.potential));
    }
    if (workspace.force != nullptr)
    {
        CHECK_CUDA(cudaFree(workspace.force));
    }
    if (workspace.virial != nullptr)
    {
        CHECK_CUDA(cudaFree(workspace.virial));
    }

    workspace = NepCudaPostprocessWorkspace{};
}

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
    NepCudaPostprocessWorkspace workspace;
    postprocess_nep_cuda(nat,
                         atomic_energy,
                         raw_force,
                         raw_virial,
                         fact_e,
                         fact_f,
                         fact_v,
                         potential,
                         force,
                         virial,
                         workspace);
    release_nep_cuda_postprocess_workspace(workspace);
}

void postprocess_nep_cuda(const int nat,
                          const double* atomic_energy,
                          const double* raw_force,
                          const double* raw_virial,
                          const double fact_e,
                          const double fact_f,
                          const double fact_v,
                          double& potential,
                          ModuleBase::matrix& force,
                          ModuleBase::matrix& virial,
                          NepCudaPostprocessWorkspace& workspace)
{
    init_nep_cuda_postprocess_workspace(workspace, nat);

    CHECK_CUDA(cudaMemcpy(workspace.energy, atomic_energy, sizeof(double) * nat, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(workspace.raw_force, raw_force, sizeof(double) * 3 * nat, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(workspace.raw_virial, raw_virial, sizeof(double) * 9 * nat, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemset(workspace.potential, 0, sizeof(double)));
    CHECK_CUDA(cudaMemset(workspace.virial, 0, sizeof(double) * 9));

    const int block_size = 256;
    const int grid_size = (nat + block_size - 1) / block_size;
    nep_postprocess_kernel<<<grid_size, block_size>>>(nat,
                                                      workspace.energy,
                                                      workspace.raw_force,
                                                      workspace.raw_virial,
                                                      fact_e,
                                                      fact_f,
                                                      fact_v,
                                                      workspace.potential,
                                                      workspace.force,
                                                      workspace.virial);
    CHECK_LAST_CUDA_ERROR("nep_postprocess_kernel");
    CHECK_CUDA_SYNC();

    CHECK_CUDA(cudaMemcpy(&potential, workspace.potential, sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(force.c, workspace.force, sizeof(double) * 3 * nat, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(virial.c, workspace.virial, sizeof(double) * 9, cudaMemcpyDeviceToHost));
}

} // namespace ModuleESolver
