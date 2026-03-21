#include "device.h"

#include "source_base/tool_quit.h"

#include <base/macros/macros.h>
#include <cstring>
#include <iostream>
#include <string>
#ifdef __MPI
#include "mpi.h"
#endif

#if defined(__CUDA)
#include <cuda_runtime.h>
#endif

#if defined(__ROCM)
#include <hip/hip_runtime.h>
#endif

namespace base_device {

namespace information {

#if __MPI
int get_node_rank_with_mpi_shared(const MPI_Comm mpi_comm) {
  // 20240530 zhanghaochong
  // The main difference between this function and the above is that it does not
  // use hostname, but uses MPI's built-in function to achieve similar
  // functions.
  MPI_Comm localComm;
  int localMpiRank;
  MPI_Comm_split_type(mpi_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                      &localComm);
  MPI_Comm_rank(localComm, &localMpiRank);
  MPI_Comm_free(&localComm);
  return localMpiRank;
}
#endif

bool probe_gpu_availability() {
#if defined(__CUDA)
    int device_count = 0;
    // Directly call cudaGetDeviceCount without CHECK_CUDA to prevent program exit
    cudaError_t error_id = cudaGetDeviceCount(&device_count);
    if (error_id == cudaSuccess && device_count > 0) {
        return true;
    }
    return false;
#elif defined(__ROCM)
    int device_count = 0;
    hipError_t error_id = hipGetDeviceCount(&device_count);
    if (error_id == hipSuccess && device_count > 0) {
        return true;
    }
    return false;
#else
    // If not compiled with GPU support, GPU is not available
    return false;
#endif
}

std::string get_device_flag(const std::string &device,
                            const std::string &basis_type) {
    // 1. Validate input string
    if (device != "cpu" && device != "gpu" && device != "auto") {
        ModuleBase::WARNING_QUIT("device", "Parameter \"device\" can only be set to \"cpu\", \"gpu\", or \"auto\"!");
    }
    
    // NOTE: This function is called only on rank 0 during input parsing.
    // The result will be broadcast to other ranks via the standard bcast mechanism.
    // DO NOT use MPI_Bcast here as other ranks are not in this code path.
    
    std::string result = "cpu";
    
    if (device == "gpu") {
        if (probe_gpu_availability()) {
            result = "gpu";
            // std::cout << " INFO: 'device=gpu' specified. GPU will be used." << std::endl;
        } else {
            ModuleBase::WARNING_QUIT("device", "Device is set to 'gpu', but no available GPU was found. Please check your hardware/drivers or set 'device=cpu'.");
        }
    } else if (device == "auto") {
        if (probe_gpu_availability()) {
            result = "gpu";
            // std::cout << " INFO: 'device=auto' specified. GPU detected and will be used." << std::endl;
        } else {
            result = "cpu";
            // std::cout << " WARNING: 'device=auto' specified, but no GPU was found. Falling back to CPU." << std::endl;
            // std::cout << "          To suppress this warning, please explicitly set 'device=cpu' in your input." << std::endl;
        }
    } else { // device == "cpu"
        result = "cpu";
        // std::cout << " INFO: 'device=cpu' specified. CPU will be used." << std::endl;
    }

    // 2. Final check for incompatible basis type
    if (result == "gpu" && basis_type == "lcao_in_pw") {
        ModuleBase::WARNING_QUIT("device", "The GPU currently does not support the basis type \"lcao_in_pw\"!");
    }

    // 3. Return the final decision
    return result;
}

} // end of namespace information

// ============================================================================
// DeviceContext singleton implementation
// ============================================================================

DeviceContext& DeviceContext::instance() {
    static DeviceContext instance;
    return instance;
}

void DeviceContext::init() {
    // Thread-safe initialization using mutex
    std::lock_guard<std::mutex> lock(init_mutex_);

    // If already initialized, do nothing (idempotent)
    if (initialized_) {
        return;
    }

#if defined(__CUDA) || defined(__ROCM)

#ifdef __MPI
    // Get local rank within the node using MPI_COMM_TYPE_SHARED
    // This is the modern and recommended way to get node-local rank
    // Use MPI_COMM_WORLD as the default communicator
    MPI_Comm local_comm;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &local_comm);
    MPI_Comm_rank(local_comm, &local_rank_);
    MPI_Comm_free(&local_comm);
#else
    local_rank_ = 0;
#endif

    // Get the number of available GPU devices
#if defined(__CUDA)
    cudaError_t err = cudaGetDeviceCount(&device_count_);
    if (err != cudaSuccess || device_count_ <= 0) {
        ModuleBase::WARNING_QUIT("DeviceContext::init",
            "No CUDA-capable GPU device found! Please check your hardware/drivers.");
        return;
    }

    // Bind to GPU device based on local rank
    device_id_ = local_rank_ % device_count_;
    err = cudaSetDevice(device_id_);
    if (err != cudaSuccess) {
        ModuleBase::WARNING_QUIT("DeviceContext::init",
            "cudaSetDevice failed! Device ID: " + std::to_string(device_id_));
        return;
    }
#elif defined(__ROCM)
    hipError_t err = hipGetDeviceCount(&device_count_);
    if (err != hipSuccess || device_count_ <= 0) {
        ModuleBase::WARNING_QUIT("DeviceContext::init",
            "No ROCm-capable GPU device found! Please check your hardware/drivers.");
        return;
    }

    // Bind to GPU device based on local rank
    device_id_ = local_rank_ % device_count_;
    err = hipSetDevice(device_id_);
    if (err != hipSuccess) {
        ModuleBase::WARNING_QUIT("DeviceContext::init",
            "hipSetDevice failed! Device ID: " + std::to_string(device_id_));
        return;
    }
#endif

    gpu_enabled_ = true;
    initialized_ = true;

#else
    // No GPU support compiled in
    initialized_ = true;
    gpu_enabled_ = false;
    device_id_ = -1;
    device_count_ = 0;
#endif
}

} // end of namespace base_device
