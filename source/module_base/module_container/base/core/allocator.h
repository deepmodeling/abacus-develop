#ifndef BASE_CORE_ALLOCATOR_H_
#define BASE_CORE_ALLOCATOR_H_

#include <ATen/core/tensor_types.h>

namespace container {

namespace base {
/**
 * @brief An abstract base class for memory allocators.
 *
 * This class defines an interface for memory allocators. Subclasses of this class
 * can provide different implementations of memory allocation/deallocation strategies.
 *
 * All memory allocated by an Allocator must be freed using the same allocator that
 * allocated it.
 */
class Allocator {
  public:
    /**
     * @brief Allocate a block of memory with the given size and default alignment.
     *
     * @param size The size of the memory block to allocate.
     *
     * @return A pointer to the allocated memory block, or nullptr if the allocation fails.
     */
    virtual void* allocate(size_t size) = 0;

    /**
     * @brief Allocate a block of memory with the given size and alignment.
     *
     * @param size The size of the memory block to allocate.
     * @param alignment The alignment of the memory block to allocate.
     *
     * @return A pointer to the allocated memory block, or nullptr if the allocation fails.
     */
    virtual void* allocate(size_t size, size_t alignment) = 0;

    /**
     * @brief Free a block of memory that was previously allocated by this allocator.
     *
     * @param ptr A pointer to the memory block to free.
     */
    virtual void free(void* ptr) = 0;

    /**
     * @brief Get the allocated size of a given pointer.
     *
     * @param ptr The pointer to get the allocated size of.
     * @return size_t The size of the allocated block of memory, in bytes.
     */
    virtual size_t AllocatedSize(void* ptr) {
        return allocated_size_;
    }

    /**
     * @brief Get the type of memory used by the TensorBuffer.
     *
     * @return MemoryType The type of memory used by the TensorBuffer.
     */
    virtual DeviceType GetDeviceType() = 0;

    virtual ~Allocator() = default;

  protected:
    /**
     * @brief The total number of bytes allocated by this allocator.
     */
    size_t allocated_size_ = 0;
};

} // namespace base
} // namespace container

#endif // BASE_CORE_ALLOCATOR_H_
