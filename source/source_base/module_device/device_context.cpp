#include "device_context.h"
#include "source_base/tool_quit.h"

#include <string>

#if defined(__CUDA)
#include <cuda_runtime.h>
#endif

#if defined(__ROCM)
#include <hip/hip_runtime.h>
#endif

namespace base_device {

DeviceContext& DeviceContext::instance() {
    static DeviceContext instance;
    return instance;
}

#ifdef __MPI
void DeviceContext::init(MPI_Comm comm) {
#else
void DeviceContext::init() {
#endif
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
    MPI_Comm local_comm;
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &local_comm);
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

} // namespace base_device
