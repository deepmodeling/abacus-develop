#ifndef ATEN_CORE_CPU_ALLOCATOR_H_
#define ATEN_CORE_CPU_ALLOCATOR_H_

#include <ATen/core/allocator.h>

namespace container {

/**
 * @brief An Allocator subclass for CPU memory.
 *
 * This class provides an implementation of the Allocator interface for CPU memory. It
 * uses the standard library functions std::malloc, std::free, and std::aligned_alloc
 * to allocate and deallocate memory blocks.
 */
class CPUAllocator : public Allocator {
  public:

    /**
     * @brief Allocate a block of CPU memory with the given size and default alignment.
     *
     * @param size The size of the memory block to allocate.
     *
     * @return A pointer to the allocated memory block, or nullptr if the allocation fails.
     */
    void* allocate(size_t size) override;

    /**
     * @brief Allocate a block of CPU memory with the given size and alignment.
     *
     * @param size The size of the memory block to allocate.
     * @param alignment The alignment of the memory block to allocate.
     *
     * @return A pointer to the allocated memory block, or nullptr if the allocation fails.
     */
    void* allocate(size_t size, size_t alignment) override;

    /**
     * @brief Free a block of CPU memory that was previously allocated by this allocator.
     *
     * @param ptr A pointer to the memory block to free.
     */
    void free(void* ptr) override;

    /**
     * @brief Get the type of device used by the TensorBuffer.
     *
     * @return MemoryType The type of memory used by the TensorBuffer.
     */
    DeviceType GetDeviceType() override;

};

} // namespace container

#endif // ATEN_CORE_CPU_ALLOCATOR_H_
