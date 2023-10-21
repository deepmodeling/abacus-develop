#ifndef BASE_CORE_BFC_ALLOCATOR_H_
#define BASE_CORE_BFC_ALLOCATOR_H_

#include <base/core/allocator.h>

namespace container {
namespace base {
/**
 * @brief An allocator that allocates memory on a GPU device.
 *
 * This class provides an implementation of the Allocator interface that allocates memory
 * on a GPU device using CUDA APIs.
 */
class BFCAllocator : public Allocator {
public:

    struct Options {
        bool allow_growth = true;
        double fragment_fraction = 0.0;
    };

    BFCAllocator(std::unique_ptr<Allocator> sub_alloc, const size_t& total_memory, const Options& options = Options());

    ~BFCAllocator();
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
     * @brief Get the type of memory used by the TensorBuffer.
     *
     * @return MemoryType The type of memory used by the TensorBuffer.
     */
    DeviceType GetDeviceType() override;

    private:

    // The sub allocator to use for extending the BFC's memory pool.
    std::unique_ptr<Allocator> sub_alloc_;

    struct bin;
    mutable std::mutex mtx_;

    // A chunk_handle is an index into the chunks_ vector in BFCAllocator
    // kInvalidChunkHandle means an invalid chunk index.
    typedef size_t chunk_handle_t;
    static constexpr chunk_handle_t kInvalidChunkHandle = UINT64_MAX;

    typedef int bin_index_t;
    static constexpr int kInvalidBinNum = -1;
    // The following means that the largest bin'd chunk size is 256 << 21 = 512MB.
    static constexpr int kNumBins = 21;

    struct chunk {
        // The size of the chunk in bytes.
        size_t size = 0;
        // The bin index of the chunk.
        bin_index_t bin_index = kInvalidBinNum;
        // We sometimes give chunks that are larger than needed to reduce
        // fragmentation.  requested_size keeps track of what the client
        // actually wanted so we can understand whether our splitting
        // strategy is efficient.
        size_t requested_size = 0;
        // allocation_id is set to -1 when the chunk is not in use. It is assigned a
        // value greater than zero before the chunk is returned from
        // AllocateRaw, and this value is unique among values assigned by
        // the parent allocator.
        int64_t allocation_id = -1;
        // pointer to granted subbuffer.
        void* ptr = nullptr;  
        chunk_handle_t next_chunk_handle = kInvalidChunkHandle;
        // The handle of the previous chunk in the bin.
        chunk_handle_t prev_chunk_handle = kInvalidChunkHandle;
        // Whether the chunk is allocated.
        bool allocated() const { return allocation_id > 0; }
    };

    struct bin {
        // The size of the chunks in this bin.
        size_t bin_size = 0;
        // The number of chunks in this bin.
        size_t num_chunks = 0;
        // The handle of the first chunk in the bin.
        chunk_handle_t first_chunk_handle = kInvalidChunkHandle;
        // The handle of the last chunk in the bin.
        chunk_handle_t last_chunk_handle = kInvalidChunkHandle;

        class chunk_comparator {
          public:
            explicit chunk_comparator(BFCAllocator* allocator) : allocator_(allocator) {}
            // Sort first by size and then use pointer address as a tie breaker.
            bool operator()(const chunk_handle_t ha,
                            const chunk_handle_t hb) const {
                const chunk* a = allocator_->chunk_from_handle(ha);
                const chunk* b = allocator_->chunk_from_handle(hb);
                if (a->size != b->size) {
                    return a->size < b->size;
                }
                return a->ptr < b->ptr;
            }

          private:
            BFCAllocator* allocator_;  // The parent allocator
        };

        using free_chunk_set_t = std::set<ChunkHandle, ChunkComparator>;
    };

    
};

} // namespace base
} // namespace container

#endif // BASE_CORE_BFC_ALLOCATOR_H_
