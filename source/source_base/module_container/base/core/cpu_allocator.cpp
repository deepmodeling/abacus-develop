#include <base/core/cpu_allocator.h>

#include <cstddef>
#ifdef _WIN32
#include <malloc.h> // _aligned_malloc / _aligned_free
#endif

namespace base {
namespace core {

// Allocate a block of CPU memory with the given size and default alignment.
// On Windows the aligned allocator family is used consistently so that every
// pointer handed out by this class can be released through free() below
// (_aligned_malloc memory must not be released with ::operator delete).
void *CPUAllocator::allocate(size_t size) {
    this->allocated_size_ = size;
#ifdef _WIN32
    return _aligned_malloc(size, alignof(std::max_align_t));
#else
    return ::operator new(size);
#endif
}

// Allocate a block of CPU memory with the given size and alignment.
void *CPUAllocator::allocate(size_t size, size_t alignment) {
    this->allocated_size_ = size;
    void *ptr = nullptr;
#ifdef _WIN32
    ptr = _aligned_malloc(size, alignment); // posix_memalign has no Windows CRT equivalent
#else
    if (posix_memalign(&ptr, alignment, size) != 0) {
        ptr = nullptr;
    }
#endif
    return ptr;
}

// Free a block of CPU memory that was previously allocated by this allocator.
void CPUAllocator::free(void *ptr) {
    this->allocated_size_ = 0;
#ifdef _WIN32
    _aligned_free(ptr);
#else
    ::operator delete(ptr);
#endif
}

//  Get the type of device used by the TensorBuffer.
container::DeviceType CPUAllocator::GetDeviceType() {
    return container::DeviceType::CpuDevice;
}

} // namespace core
} // namespace base
