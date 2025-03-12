#include "module_base/kernels/math_op.h"

#include "module_base/module_device/memory_op.h"

#include <complex>
#include <gtest/gtest.h>
#include <vector>

class TestModuleBaseMathMultiDevice : public ::testing::Test
{
protected:
    // xx = tf.random.uniform([100], minval=-4, maxval=4, dtype = tf.float64)

  const base_device::DEVICE_CPU* cpu_ctx = {};
  const base_device::DEVICE_GPU* gpu_ctx = {};

  int ng = 59, lmax = 1;

  double SQRT2 = 1.4142135623730951, PI = 3.1415926535897931, PI_HALF = 1.5707963267948966,
         FOUR_PI = 12.566370614359172, SQRT_INVERSE_FOUR_PI = 0.28209479177387814;

  std::vector<double> g
      = {2,  -2, -2, 1,  -1, -1, 0,  0,  0,  -1, 1,  1,  -2, 2,  2,  2,  -2, 0,  1,  -1, 1,  0,  0,  2,  -1, 1,
         3,  2,  -2, 2,  1,  -1, 3,  -1, 1,  -3, -2, 2,  -2, 1,  -1, -3, 0,  0,  -2, -1, 1,  -1, -2, 2,  0,  2,
         0,  -2, 1,  1,  -1, 0,  2,  0,  -1, 3,  1,  3,  -1, -1, 2,  0,  0,  1,  1,  1,  0,  2,  2,  3,  -1, 1,
         2,  0,  2,  1,  1,  3,  1,  1,  -3, 0,  2,  -2, -1, 3,  -1, 2,  2,  -2, 1,  3,  -1, 3,  1,  -1, 2,  2,
         0,  1,  3,  1,  3,  1,  1,  2,  2,  2,  -1, -3, 1,  -2, -2, 2,  -2, -2, -2, -3, -1, -1, -1, -3, -1, -2,
         -2, 0,  -3, -1, 1,  1,  -3, -1, 0,  -2, 0,  -1, -1, 1,  -2, 0,  2,  1,  -3, 1,  0,  -2, 2,  -1, -1, 3,
         -1, -1, -3, -2, 0,  -2, -3, 1,  -1, 0,  -2, -2, -1, -1, -1, -2, 0,  0,  -3, 1,  1};
  std::vector<double> expected_ylm = {0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      0.282095,
                                      -0.282095,
                                      -0.282095,
                                      0,
                                      0.282095,
                                      0.282095,
                                      0,
                                      0.282095,
                                      0.488603,
                                      0.441958,
                                      0.282095,
                                      0.441958,
                                      -0.441958,
                                      -0.282095,
                                      -0.441958,
                                      -0.488603,
                                      -0.282095,
                                      0,
                                      -0.345494,
                                      -0.282095,
                                      0,
                                      0.147319,
                                      -0.147319,
                                      0,
                                      0.282095,
                                      0.345494,
                                      0.147319,
                                      0.345494,
                                      0.441958,
                                      -0.441958,
                                      -0.345494,
                                      -0.147319,
                                      -0.282095,
                                      -0.147319,
                                      -0.147319,
                                      0,
                                      0.147319,
                                      0.147319,
                                      0.282095,
                                      0.147319,
                                      0.282095,
                                      -0.282095,
                                      -0.147319,
                                      -0.147319,
                                      0,
                                      0.147319,
                                      -0.147319,
                                      0,
                                      0.282095,
                                      0.345494,
                                      0.147319,
                                      0.345494,
                                      0.441958,
                                      -0.441958,
                                      -0.345494,
                                      -0.147319,
                                      -0.345494,
                                      -0.282095,
                                      0,
                                      0.147319,
                                      -0.282095,
                                      -0.282095,
                                      -2.99183e-17,
                                      0.282095,
                                      0.282095,
                                      -0.345494,
                                      -0.282095,
                                      -0,
                                      0.147319,
                                      -0.282095,
                                      -0.147319,
                                      0.147319,
                                      0.282095,
                                      -0.147319,
                                      -0,
                                      0.282095,
                                      0.345494,
                                      -0.345494,
                                      -0.282095,
                                      -2.99183e-17,
                                      0.147319,
                                      -0.441958,
                                      -0.488603,
                                      -0.282095,
                                      -2.11554e-17,
                                      -0.441958,
                                      -0.345494,
                                      -0.147319,
                                      -0.147319,
                                      -2.11554e-17,
                                      0.147319,
                                      -0.282095,
                                      -0.147319,
                                      -0.441958,
                                      -0.345494,
                                      -0.147319,
                                      -0.441958,
                                      -0.282095,
                                      0.147319,
                                      0.282095,
                                      0.282095,
                                      0.441958,
                                      0.147319,
                                      0.345494,
                                      0.441958,
                                      -0.147319,
                                      -2.99183e-17,
                                      0.282095,
                                      0.345494,
                                      -0.147319,
                                      -2.11554e-17,
                                      0.147319,
                                      0.147319,
                                      0.345494,
                                      0.441958,
                                      -2.11554e-17,
                                      0.282095,
                                      0.488603,
                                      0.441958,
                                      0.282095,
                                      0.282095,
                                      -0.488603,
                                      -0.282095,
                                      -0.282095,
                                      0.345494,
                                      0.282095,
                                      -0,
                                      -0.147319,
                                      0.282095,
                                      0.147319,
                                      -0.147319,
                                      -0.282095,
                                      0.147319,
                                      -0,
                                      -0.282095,
                                      -0.345494,
                                      -0,
                                      -0.282095,
                                      -0.488603,
                                      -0.441958,
                                      0.147319,
                                      -0,
                                      -0.282095,
                                      -0.345494,
                                      0.147319,
                                      -0,
                                      -0.147319,
                                      -0.147319,
                                      -0.345494,
                                      -0.441958,
                                      -0.282095,
                                      -0.441958,
                                      -0.147319,
                                      -0.345494,
                                      -0.441958,
                                      -0.147319,
                                      -0.282095,
                                      0.441958,
                                      0.282095,
                                      0.282095,
                                      0.147319,
                                      0.441958,
                                      0.345494,
                                      0.147319,
                                      0.441958,
                                      0.488603,
                                      0.282095,
                                      -4.23108e-17,
                                      0.441958,
                                      0.345494,
                                      0.147319,
                                      0.147319,
                                      -4.23108e-17,
                                      -0.147319,
                                      0.345494,
                                      0.282095,
                                      -5.98366e-17,
                                      -0.147319};

  using delmem_var_op = base_device::memory::delete_memory_op<double, base_device::DEVICE_GPU>;
  using resmem_var_op = base_device::memory::resize_memory_op<double, base_device::DEVICE_GPU>;
  using syncmem_var_h2d_op
      = base_device::memory::synchronize_memory_op<double, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
  using syncmem_var_d2h_op
      = base_device::memory::synchronize_memory_op<double, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;

  void SetUp() override
  {
  }
    void TearDown() override {
    }
};

TEST_F(TestModuleBaseMathMultiDevice, cal_ylm_real_op_cpu)
{
    std::vector<double> p((lmax + 1) * (lmax + 1) * ng, 0.0);
    std::vector<double> ylm(expected_ylm.size(), 0.0);
    ModuleBase::cal_ylm_real_op<double, base_device::DEVICE_CPU>()(cpu_ctx,
                                                                   ng,
                                                                   lmax,
                                                                   SQRT2,
                                                                   PI,
                                                                   PI_HALF,
                                                                   FOUR_PI,
                                                                   SQRT_INVERSE_FOUR_PI,
                                                                   g.data(),
                                                                   p.data(),
                                                                   ylm.data());

    for (int ii = 0; ii < ylm.size(); ii++) {
        EXPECT_LT(fabs(ylm[ii] - expected_ylm[ii]), 6e-5);
    }
}

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
TEST_F(TestModuleBaseMathMultiDevice, cal_ylm_real_op_gpu)
{
    std::vector<double> p((lmax + 1) * (lmax + 1) * ng, 0.0);
    std::vector<double> ylm(expected_ylm.size(), 0.0);
    double * d_ylm = nullptr, * d_g = nullptr, * d_p = nullptr;

    resmem_var_op()(d_g, g.size());
    resmem_var_op()(d_p, p.size());
    resmem_var_op()(d_ylm, ylm.size());

    syncmem_var_h2d_op()(d_g, g.data(), g.size());
    syncmem_var_h2d_op()(d_p, p.data(), p.size());
    syncmem_var_h2d_op()(d_ylm, ylm.data(), ylm.size());

    ModuleBase::cal_ylm_real_op<double, base_device::DEVICE_GPU>()(gpu_ctx,
                                                                   ng,
                                                                   lmax,
                                                                   SQRT2,
                                                                   PI,
                                                                   PI_HALF,
                                                                   FOUR_PI,
                                                                   SQRT_INVERSE_FOUR_PI,
                                                                   d_g,
                                                                   d_p,
                                                                   d_ylm);

    syncmem_var_d2h_op()(ylm.data(), d_ylm, ylm.size());

    for (int ii = 0; ii < ylm.size(); ii++) {
        EXPECT_LT(fabs(ylm[ii] - expected_ylm[ii]), 6e-5);
    }

    delmem_var_op()(d_g);
    delmem_var_op()(d_p);
    delmem_var_op()(d_ylm);
}

#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM