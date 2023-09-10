#if defined(__CUDA) || defined(__ROCM)

#include <cuda_runtime.h> // for CUDA APIs
#include <base/core/gpu_allocator.h>

namespace container {
namespace base {
// Allocate a block of memory with the given size and default alignment on GPU.
void* GPUAllocator::allocate(size_t size) {
    void* ptr;
    cudaError_t result = cudaMalloc(&ptr, size);
    if (result != cudaSuccess) {
        return nullptr;
    }
    this->allocated_size_ = size;
    return ptr;
}

// Allocate a block of CPU memory with the given size and alignment.
void* GPUAllocator::allocate(size_t size, size_t alignment) {
    void* ptr;
    cudaError_t result = cudaMalloc(&ptr, size);
    if (result != cudaSuccess) {
        return nullptr;
    }
    this->allocated_size_ = size;
    return ptr;
}

// Free a block of CPU memory that was previously allocated by this allocator.
void GPUAllocator::free(void* ptr) {
    cudaFree(ptr);
    this->allocated_size_ = 0;
}

// Get the type of device used by the TensorBuffer.
DeviceType GPUAllocator::GetDeviceType() {
    return DeviceType::GpuDevice;
}

} // namespace base
} // namespace container

#endif // __CUDA || __ROCM
