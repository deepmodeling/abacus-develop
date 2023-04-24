#include <cstdlib> // for ::operator new, ::operator delete
#include <cassert> // for assert

#include "cpu_allocator.h"

namespace container {

// Allocate a block of CPU memory with the given size and default alignment.
void* CPUAllocator::allocate(size_t size) {
    return ::operator new(size);
}

// Allocate a block of CPU memory with the given size and alignment.
void* CPUAllocator::allocate(size_t size, size_t alignment) {
    void* ptr = nullptr;
    if (posix_memalign(&ptr, alignment, size) != 0) {
        ptr = nullptr;
    }
    return ptr;
}

// Free a block of CPU memory that was previously allocated by this allocator.
void CPUAllocator::free(void* ptr) {
    ::operator delete(ptr);
}

// Get the allocated size of a given pointer.
size_t CPUAllocator::AllocatedSize(void* ptr) {
    std::cerr << "Method `AllocatedSize` not implemented!" << std::endl;
    exit(EXIT_FAILURE);
}

//  Get the type of device used by the TensorBuffer.
DeviceType CPUAllocator::GetDeviceType() {
    return DeviceType::CpuDevice;
}

} // namespace container
