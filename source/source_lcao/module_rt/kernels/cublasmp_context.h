#ifndef CUBLASMP_CONTEXT_H
#define CUBLASMP_CONTEXT_H

#ifdef __MPI
#include <mpi.h>
#endif

#ifdef __CUDA
#include <cuda_runtime.h>
#endif

#ifdef __CUBLASMP
#include "source_base/global_variable.h"
#include "source_base/module_device/device.h"

#include <cublasmp.h>
#include <cusolverMp.h>
#include <iostream>
#include <nccl.h>

extern "C"
{
#include "source_hsolver/module_genelpa/Cblacs.h"
}

#define LOG_DEBUG(msg)                                                                                                 \
    do                                                                                                                 \
    {                                                                                                                  \
        if (g_EnableDebugLog)                                                                                          \
        {                                                                                                              \
            std::cerr << "[DEBUG] " << msg << " (at " << __func__ << ")" << std::endl;                                 \
        }                                                                                                              \
    } while (0)
#endif // __CUBLASMP

// The struct is ALWAYS available.
struct CublasMpResources
{
    bool is_initialized = false;

#ifdef __MPI
    MPI_Comm mpi_comm = MPI_COMM_NULL;
#endif

#ifdef __CUDA
    cudaStream_t stream = nullptr;
#endif

#ifdef __CUBLASMP
    ncclComm_t nccl_comm = nullptr;

    cublasMpHandle_t cublasmp_handle = nullptr;
    cublasMpGrid_t cublasmp_grid = nullptr;

    cusolverMpHandle_t cusolvermp_handle = nullptr;
    cusolverMpGrid_t cusolvermp_grid = nullptr;
#endif
};

// API functions are only visible when cuBLASMp is enabled.
#ifdef __CUBLASMP

inline void init_cublasmp_resources(CublasMpResources& res, MPI_Comm mpi_comm, const int* desc)
{
    if (res.is_initialized)
    {
        return;
    }

    res.mpi_comm = mpi_comm;
    MPI_Barrier(res.mpi_comm);

    // 1. Get BLACS topology info
    int cblacs_ctxt = desc[1];
    int nprows, npcols, myprow, mypcol;
    Cblacs_gridinfo(cblacs_ctxt, &nprows, &npcols, &myprow, &mypcol);

    GlobalV::ofs_running << "nprows = " << nprows << std::endl;
    GlobalV::ofs_running << "npcols = " << npcols << std::endl;
    GlobalV::ofs_running << "myprow = " << myprow << std::endl;
    GlobalV::ofs_running << "mypcol = " << mypcol << std::endl;
    GlobalV::ofs_running << "device = " << base_device::DeviceContext::instance().get_device_id() << std::endl;

    int rank, size;
    MPI_Comm_rank(res.mpi_comm, &rank);
    MPI_Comm_size(res.mpi_comm, &size);

    int device_id = base_device::DeviceContext::instance().get_device_id();
    cudaSetDevice(device_id);
    cudaStreamCreate(&res.stream);

    // 2. Initialize NCCL communicator
    ncclUniqueId id;
    if (rank == 0)
    {
        ncclGetUniqueId(&id);
    }
    // Broadcast the unique NCCL ID to all ranks
    MPI_Bcast((void*)&id, sizeof(id), MPI_BYTE, 0, res.mpi_comm);
    // Initialize NCCL with the generated ID
    ncclCommInitRank(&res.nccl_comm, size, id, rank);

    // 3. Initialize cuBLASMp specific resources
    cublasMpCreate(&res.cublasmp_handle, res.stream);
    cublasMpGridCreate(nprows, npcols, CUBLASMP_GRID_LAYOUT_ROW_MAJOR, res.nccl_comm, &res.cublasmp_grid);

    // 4. Initialize cuSOLVERMp specific resources
    cusolverMpCreate(&res.cusolvermp_handle, device_id, res.stream);
    cusolverMpCreateDeviceGrid(res.cusolvermp_handle,
                               &res.cusolvermp_grid,
                               res.nccl_comm,
                               nprows,
                               npcols,
                               CUSOLVERMP_GRID_MAPPING_ROW_MAJOR);

    res.is_initialized = true;
}

inline void finalize_cublasmp_resources(CublasMpResources& res)
{
    if (!res.is_initialized)
    {
        return;
    }

    if (res.stream)
    {
        cudaStreamSynchronize(res.stream);
    }

    // Destroy cuBLASMp resources
    if (res.cublasmp_grid)
    {
        cublasMpGridDestroy(res.cublasmp_grid);
    }
    if (res.cublasmp_handle)
    {
        cublasMpDestroy(res.cublasmp_handle);
    }

    // Destroy cuSOLVERMp resources
    if (res.cusolvermp_grid)
    {
        cusolverMpDestroyGrid(res.cusolvermp_grid);
    }
    if (res.cusolvermp_handle)
    {
        cusolverMpDestroy(res.cusolvermp_handle);
    }

    // Destroy NCCL communicator
    if (res.nccl_comm)
    {
        ncclCommDestroy(res.nccl_comm);
    }

    if (res.stream)
    {
        cudaStreamDestroy(res.stream);
    }

    res.is_initialized = false;
}

#endif // __CUBLASMP

#endif // CUBLASMP_CONTEXT_H
