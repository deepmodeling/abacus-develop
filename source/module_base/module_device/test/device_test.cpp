#include "../device.h"

#include <complex>
#include <gtest/gtest.h>
#include <iostream>

class TestModulePsiDevice : public ::testing::Test
{
  protected:
    const device::CPU* cpu_ctx = {};
    const device::GPU* gpu_ctx = {};

    void SetUp() override
    {
    }
    void TearDown() override
    {
    }
};

TEST_F(TestModulePsiDevice, get_device_type_cpu)
{
    device::AbacusDevice_t device = device::get_device_type<device::CPU>(cpu_ctx);
    EXPECT_EQ(device, device::CpuDevice);
}

#if __UT_USE_CUDA || __UT_USE_ROCM
TEST_F(TestModulePsiDevice, get_device_type_gpu)
{
    device::AbacusDevice_t device = device::get_device_type<device::GPU>(gpu_ctx);
    EXPECT_EQ(device, device::GpuDevice);
}
#endif // __UT_USE_CUDA || __UT_USE_ROCM
