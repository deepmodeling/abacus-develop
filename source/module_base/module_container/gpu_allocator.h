#ifndef CONTAINER_GPU_ALLOCATOR_H
#define CONTAINER_GPU_ALLOCATOR_H

#include "allocator.h"

namespace container {

/**
 * @brief An allocator that allocates memory on a GPU device.
 *
 * This class provides an implementation of the Allocator interface that allocates memory
 * on a GPU device using CUDA APIs.
 */
class GPUAllocator : public Allocator {
public:
    /**
     * @brief Allocate a block of memory with the given size and default alignment on GPU.
     *
     * @param size The size of the memory block to allocate.
     *
     * @return A pointer to the allocated memory block, or nullptr if the allocation fails.
     */
    void* allocate(size_t size) override;

    /**
     * @brief Allocate a block of memory with the given size and alignment on GPU.
     *
     * @param size The size of the memory block to allocate.
     * @param alignment The alignment of the memory block to allocate.
     *
     * @return A pointer to the allocated memory block, or nullptr if the allocation fails.
     */
    void* allocate(size_t size, size_t alignment) override;

    /**
     * @brief Free a block of GPU memory that was previously allocated by this allocator.
     *
     * @param ptr A pointer to the memory block to free.
     */
    void free(void* ptr) override;

    /**
     * @brief Get the allocated size of a given pointer.
     *
     * @param ptr The pointer to get the allocated size of.
     * @return size_t The size of the allocated block of memory, in bytes.
     *
     * @note This function is not implemented for CPUAllocator and always returns 0.
     */
    size_t AllocatedSize(void* ptr) override;

    /**
     * @brief Get the type of memory used by the TensorBuffer.
     *
     * @return MemoryType The type of memory used by the TensorBuffer.
     */
    DeviceType GetDeviceType() override;
};

} // namespace ABACUS

#endif // CONTAINER_GPU_ALLOCATOR_H
