#include <vector>
#include <gtest/gtest.h>

#include <ATen/core/tensor.h>
#include <base/core/allocator.h>
#include <base/core/cpu_allocator.h>


TEST(CPUAllocator, AllocateAndFree) {
  base::core::CPUAllocator alloc;
  // Allocate memory of size 100.
  void* ptr = alloc.allocate(100);
  EXPECT_NE(nullptr, ptr);
  alloc.free(ptr);

  // Allocate memory of size 200 with alignment 16.
  ptr = alloc.allocate(200, 16);
  EXPECT_NE(nullptr, ptr);
  alloc.free(ptr);

  // Allocate memory of size 200 with alignment 16.
  ptr = alloc.allocate(0, 0);
  EXPECT_EQ(nullptr, ptr);
}

TEST(CPUAllocator, AllocatedSize) {
  base::core::CPUAllocator alloc;
  // Allocate memory of size 100 and check its size.
  void* ptr = alloc.allocate(100);
  EXPECT_NE(nullptr, ptr);
  alloc.free(ptr);
}

TEST(CPUAllocator, GetDeviceType) {
  base::core::CPUAllocator alloc;
  EXPECT_EQ(container::DeviceType::CpuDevice,
            alloc.GetDeviceType());
}