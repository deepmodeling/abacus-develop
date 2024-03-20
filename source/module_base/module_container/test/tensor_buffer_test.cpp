#include <gtest/gtest.h>

#include <ATen/core/tensor_buffer.h>
#include <base/core/cpu_allocator.h>

// Test the GetAllocatedBytes() method.
TEST(TensorBuffer, GetAllocatedBytes) {
    // Create an allocator and allocate memory for a TensorBuffer.
    base::core::Allocator* alloc = new base::core::CPUAllocator();
    const size_t buffer_size = 100;

    // Create a TensorBuffer.
    container::TensorBuffer tensor_buffer(alloc, 100);

    // Check the allocator 
    EXPECT_EQ(alloc, tensor_buffer.allocator());

    // Check the DeviceType
    EXPECT_EQ(container::DeviceType::CpuDevice, tensor_buffer.GetDeviceType());

    // Check the size of the buffer.
    EXPECT_EQ(buffer_size, tensor_buffer.GetAllocatedBytes());
}

// Test the resize() method.
TEST(TensorBuffer, resize) {
    // Create an allocator and allocate memory for a TensorBuffer.
    base::core::Allocator* alloc = new base::core::CPUAllocator();
    const size_t initial_buffer_size = 100;

    // Create a TensorBuffer.
    container::TensorBuffer tensor_buffer(alloc, initial_buffer_size);

    // Check the allocator 
    EXPECT_EQ(alloc, tensor_buffer.allocator());

    // Check the DeviceType
    EXPECT_EQ(container::DeviceType::CpuDevice, tensor_buffer.GetDeviceType());

    // Resize the buffer.
    const size_t new_buffer_size = 200;
    tensor_buffer.resize(new_buffer_size);

    // Free the memory.
    // auto free by the destructor
    // alloc.free(buffer);
}

// Test the root_buffer() method.
TEST(TensorBuffer, root_buffer) {
    // Create an allocator and allocate memory for a TensorBuffer.
    base::core::Allocator* alloc = new base::core::CPUAllocator();
    const size_t buffer_size = 100;

    // Create a root TensorBuffer.
    container::TensorBuffer root_buffer(alloc, buffer_size);

    // Check the allocator 
    EXPECT_EQ(alloc, root_buffer.allocator());

    // Check the DeviceType
    EXPECT_EQ(container::DeviceType::CpuDevice, root_buffer.GetDeviceType());

    // Check that the sub-buffer's root buffer is correct.
    EXPECT_EQ(&root_buffer, root_buffer.root_buffer());

    // Free the memory.
    // alloc.free(buffer);
}

TEST(TensorBuffer, empty_allocator) {
    // Create an allocator and allocate memory for a TensorBuffer.
    base::core::CPUAllocator alloc;
    const size_t buffer_size = 100;
    void* buffer = alloc.allocate(buffer_size);

    // Create a root TensorBuffer.
    container::TensorBuffer root_buffer(buffer);

    // Check the DeviceType
    EXPECT_EQ(container::DeviceType::UnKnown, root_buffer.GetDeviceType());

    // Free the memory.
    alloc.free(buffer);
}