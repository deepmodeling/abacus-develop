#include <complex>
#include <iostream>
#include <gtest/gtest.h>
#include "module_psi/kernels/types.h"
#include "module_psi/kernels/device.h"

class TestModulePsiDevice : public ::testing::Test
{
  protected:

    const psi::DEVICE_CPU * cpu_ctx = {};
    const psi::DEVICE_GPU * gpu_ctx = {};

    void SetUp() override {
    }
    void TearDown() override {
    }
};

TEST_F(TestModulePsiDevice, get_device_type_cpu)
{
  psi::AbacusDevice_t device = psi::device::get_device_type<psi::DEVICE_CPU>(cpu_ctx);
  EXPECT_EQ(device, psi::CpuDevice);
}

#if __UT_USE_CUDA || __UT_USE_ROCM
TEST_F(TestModulePsiDevice, get_device_type_gpu)
{
  psi::AbacusDevice_t device = psi::device::get_device_type<psi::DEVICE_GPU>(gpu_ctx);
  EXPECT_EQ(device, psi::GpuDevice);
}
#endif // __UT_USE_CUDA || __UT_USE_ROCM

