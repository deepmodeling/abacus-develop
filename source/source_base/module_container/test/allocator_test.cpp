#include <vector>
#include <gtest/gtest.h>

#include <ATen/core/tensor.h>
#include <base/core/allocator.h>
#include <base/core/cpu_allocator.h>
#include <base/utils/logging.h>


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
  base::core::Allocator* alloc = new base::core::CPUAllocator();
  // Allocate memory of size 100 and check its size.
  void* ptr = alloc->allocate(100);
  EXPECT_NE(nullptr, ptr);
  EXPECT_EQ(alloc->AllocatedSize(ptr), 100);
  alloc->free(ptr);
  EXPECT_EQ(alloc->AllocatedSize(ptr), 0);
  delete alloc;
}

TEST(CPUAllocator, GetDeviceType) {
  base::core::CPUAllocator alloc;
  EXPECT_EQ(container::DeviceType::CpuDevice,
            alloc.GetDeviceType());
}

TEST(Logging, ReturnsMessage)
{
  const char* message = "allocation failed";
  EXPECT_STREQ(base::utils::check_msg_impl(message), message);
}

TEST(LoggingDeathTest, AbortsWithContext)
{
  EXPECT_DEATH(base::utils::check_exit_impl("allocate", "allocator_test.cpp", 42, "allocation failed"),
               "Fatal error.*allocation failed");
}

namespace
{
class TestCounted : public base::core::counted_base
{
  public:
    explicit TestCounted(bool* destroyed) : destroyed_(destroyed) {}

    ~TestCounted() override
    {
        *destroyed_ = true;
    }

  private:
    bool* destroyed_;
};
} // namespace

TEST(RefCount, TracksReferencesAndDeletes)
{
    bool destroyed = false;
    TestCounted* counted = new TestCounted(&destroyed);
    EXPECT_TRUE(counted->ref_count_is_one());
    counted->ref();
    EXPECT_EQ(counted->ref_count(), 2);
    EXPECT_FALSE(counted->unref());
    EXPECT_TRUE(counted->unref());
    EXPECT_TRUE(destroyed);
}
