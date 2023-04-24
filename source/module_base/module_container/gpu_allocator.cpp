#include <cassert> // for assert
#include <cuda_runtime.h> // for CUDA APIs

#include "gpu_allocator.h"

namespace container {

// Allocate a block of memory with the given size and default alignment on GPU.
void* GPUAllocator::allocate(size_t size) {
    void* ptr;
    cudaError_t result = cudaMalloc(&ptr, size);
    if (result != cudaSuccess) {
        return nullptr;
    }
    return ptr;
}

// Allocate a block of CPU memory with the given size and alignment.
void* GPUAllocator::allocate(size_t size, size_t alignment) {
    void* ptr;
    cudaError_t result = cudaMalloc(&ptr, size);
    if (result != cudaSuccess) {
        return nullptr;
    }
    return ptr;
}

// Free a block of CPU memory that was previously allocated by this allocator.
void GPUAllocator::free(void* ptr) {
    cudaFree(ptr);
}

// Get the allocated size of a given pointer.
size_t GPUAllocator::AllocatedSize(void* ptr) {
    assert(false && "not implemented");
    return 0;
}

// Get the type of device used by the TensorBuffer.
DeviceType GPUAllocator::GetDeviceType() {
    return DeviceType::GpuDevice;
}

} // namespace container
