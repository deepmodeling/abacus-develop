#include <gtest/gtest.h>

#include "../tensor_buffer.h"
#include "../cpu_allocator.h"

// Test the GetAllocatedBytes() method.
TEST(TensorBuffer, GetAllocatedBytes) {
    // Create an allocator and allocate memory for a TensorBuffer.
    container::CPUAllocator alloc;
    const size_t buffer_size = 100;
    void* buffer = alloc.allocate(buffer_size);

    // Create a TensorBuffer.
    container::TensorBuffer tensor_buffer(&alloc, buffer);

    // Check the allocator 
    EXPECT_EQ(&alloc, tensor_buffer.allocator());

    // Check the DeviceType
    EXPECT_EQ(container::DeviceType::CpuDevice, tensor_buffer.GetDeviceType());

    // Check that the allocated bytes are correct.
    ASSERT_EXIT(
        tensor_buffer.GetAllocatedBytes(),
        ::testing::ExitedWithCode(EXIT_FAILURE),
        "Method `AllocatedSize` not implemented!"
    );
}

// Test the resize() method.
TEST(TensorBuffer, resize) {
    // Create an allocator and allocate memory for a TensorBuffer.
    container::CPUAllocator alloc;
    const size_t initial_buffer_size = 100;
    void* buffer = alloc.allocate(initial_buffer_size);

    // Create a TensorBuffer.
    container::TensorBuffer tensor_buffer(&alloc, buffer);

    // Check the allocator 
    EXPECT_EQ(&alloc, tensor_buffer.allocator());

    // Check the DeviceType
    EXPECT_EQ(container::DeviceType::CpuDevice, tensor_buffer.GetDeviceType());

    // Resize the buffer.
    const size_t new_buffer_size = 200;
    tensor_buffer.resize(new_buffer_size);

    // Check that the allocated bytes are correct.
    ASSERT_EXIT(
        tensor_buffer.GetAllocatedBytes(),
        ::testing::ExitedWithCode(EXIT_FAILURE),
        "Method `AllocatedSize` not implemented!"
    );

    // Free the memory.
    // auto free by the destructor
    // alloc.free(buffer);
}

// Test the root_buffer() method.
TEST(TensorBuffer, root_buffer) {
    // Create an allocator and allocate memory for a TensorBuffer.
    container::CPUAllocator alloc;
    const size_t buffer_size = 100;
    void* buffer = alloc.allocate(buffer_size);

    // Create a root TensorBuffer.
    container::TensorBuffer root_buffer(&alloc, buffer);

    // Check the allocator 
    EXPECT_EQ(&alloc, root_buffer.allocator());

    // Check the DeviceType
    EXPECT_EQ(container::DeviceType::CpuDevice, root_buffer.GetDeviceType());

    // Check that the sub-buffer's root buffer is correct.
    EXPECT_EQ(&root_buffer, root_buffer.root_buffer());

    // Free the memory.
    // alloc.free(buffer);
}

TEST(TensorBuffer, empty_allocator) {
    // Create an allocator and allocate memory for a TensorBuffer.
    container::CPUAllocator alloc;
    const size_t buffer_size = 100;
    void* buffer = alloc.allocate(buffer_size);

    // Create a root TensorBuffer.
    container::TensorBuffer root_buffer(buffer);

    // Check the DeviceType
    EXPECT_EQ(container::DeviceType::UnKnown, root_buffer.GetDeviceType());

    // Free the memory.
    alloc.free(buffer);
}