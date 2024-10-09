#include "module_base/module_device/device.h"

#include <complex>
#include <gtest/gtest.h>
#include <iostream>

class TestModulePsiDevice : public ::testing::Test
{
  protected:
    const base_device::DEVICE_CPU* cpu_ctx = {};
    const base_device::DEVICE_GPU* gpu_ctx = {};

    void SetUp() override
    {
    }
    void TearDown() override
    {
    }
};

TEST_F(TestModulePsiDevice, get_device_type_cpu)
{
    base_device::AbacusDevice_t device = base_device::get_device_type<base_device::DEVICE_CPU>(cpu_ctx);
    EXPECT_EQ(device, base_device::CpuDevice);
}

#if __UT_USE_CUDA || __UT_USE_ROCM
TEST_F(TestModulePsiDevice, get_device_type_gpu)
{
    base_device::AbacusDevice_t device = base_device::get_device_type<base_device::DEVICE_GPU>(gpu_ctx);
    EXPECT_EQ(device, base_device::GpuDevice);
}
#endif // __UT_USE_CUDA || __UT_USE_ROCM
