#include "source_base/module_device/device.h"

#include <complex>
#include <gmock/gmock.h>
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
    base_device::AbacusDevice_t device = base_device::get_device_type(cpu_ctx);
    EXPECT_EQ(device, base_device::CpuDevice);
}

TEST_F(TestModulePsiDevice, resolve_device_flag)
{
    using base_device::information::resolve_device_flag;

    EXPECT_EQ(resolve_device_flag("auto", "pw", true), "gpu");
    EXPECT_EQ(resolve_device_flag("auto", "pw", false), "cpu");
    EXPECT_EQ(resolve_device_flag("cpu", "pw", true), "cpu");
    EXPECT_EQ(resolve_device_flag("gpu", "pw", true), "gpu");

    testing::internal::CaptureStdout();
    EXPECT_EXIT(resolve_device_flag("gpu", "pw", false), ::testing::ExitedWithCode(1), "");
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("no available GPU"));
}

#if __UT_USE_CUDA || __UT_USE_ROCM
TEST_F(TestModulePsiDevice, get_device_type_gpu)
{
    base_device::AbacusDevice_t device = base_device::get_device_type(gpu_ctx);
    EXPECT_EQ(device, base_device::GpuDevice);
}
#endif // __UT_USE_CUDA || __UT_USE_ROCM
