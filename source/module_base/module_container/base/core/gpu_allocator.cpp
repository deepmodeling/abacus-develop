#include <base/core/gpu_allocator.h>

#if defined(__CUDA)
#include <cuda_runtime.h> // for CUDA APIs
#define device_malloc cudaMalloc
#define device_free cudaFree
#define device_result_t cudaError_t
#define device_success cudaSuccess
#elif defined(__ROCM)
#include <hip/hip_runtime.h> // for ROCm APIs
#define device_malloc hipMalloc
#define device_free hipFree
#define device_result_t hipError_t
#define device_success hipSuccess
#endif 

namespace base {
namespace core {
// Allocate a block of memory with the given size and default alignment on GPU.
void *GPUAllocator::allocate(size_t size) {
    void *ptr;
    device_result_t result = device_malloc(&ptr, size);
    if (result != device_success) {
        return nullptr;
    }
    this->allocated_size_ = size;
    return ptr;
}

// Allocate a block of CPU memory with the given size and alignment.
void *GPUAllocator::allocate(size_t size, size_t alignment) {
    void *ptr;
    device_result_t result = device_malloc(&ptr, size);
    if (result != device_success) {
        return nullptr;
    }
    this->allocated_size_ = size;
    return ptr;
}

// Free a block of CPU memory that was previously allocated by this allocator.
void GPUAllocator::free(void *ptr) {
    device_free(ptr);
    this->allocated_size_ = 0;
}

// Get the type of device used by the TensorBuffer.
container::DeviceType GPUAllocator::GetDeviceType() {
    return container::DeviceType::GpuDevice;
}

} // namespace core
} // namespace base
