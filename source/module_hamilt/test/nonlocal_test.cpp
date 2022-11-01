#include <vector>
#include <complex>
#include <iostream>
#include <gtest/gtest.h>
#include "module_psi/include/memory.h"
#include "module_hamilt/include/nonlocal.h"

class TestModuleHamiltNonlocal : public ::testing::Test
{
  protected:
    // xx = tf.random.uniform([100], minval=-4, maxval=4, dtype = tf.float64)
    const int l1 = 2;
    const int l2 = 8;
    const int l3 = 4;

    int sum = 0, expected_sum = 8;
    int iat = 0, expected_iat = 2;

    const int spin = 0;
    const int nkb = 8;
    const int deeq_x = 2;
    const int deeq_y = 4;
    const int deeq_z = 4;

    std::vector<double> deeq = {
      1.52389, 0, 0, 0, 0, 3.6833, 0, 0, 0, 0, 3.6833, 0, 0, 0, 0, 3.6833, 1.52389, 0, 0, 0, 0, 3.6833, 0, 0, 0, 0, 3.6833, 0, 0, 0, 0, 3.6833
    };

    std::vector<std::complex<double>> expected_ps = {
      std::complex<double>(1.08811,0), std::complex<double>(0,2.11482e-17), std::complex<double>(0,4.22963e-17), std::complex<double>(0,-6.34445e-17), std::complex<double>(0.473301,2.11482e-17), std::complex<double>(-3.80667e-16,-5.28704e-18), std::complex<double>(-2.11482e-17,-2.11482e-17), std::complex<double>(4.22963e-17,2.64352e-17), std::complex<double>(1.50743e-33,1.75712e-17), std::complex<double>(0.93285,-5.71206e-17), std::complex<double>(0,-2.83751e-33), std::complex<double>(-3.83371e-17,1.41876e-33), 
      std::complex<double>(6.38951e-18,2.79541e-18), std::complex<double>(-0.0996506,1.91685e-17), std::complex<double>(1.2779e-17,5.11161e-17), std::complex<double>(2.87528e-17,3.83371e-17), std::complex<double>(-2.48283e-33,7.02847e-17), std::complex<double>(1.2779e-17,-2.83751e-33), std::complex<double>(0.93285,-5.71206e-17), std::complex<double>(-1.02232e-16,3.54689e-34), std::complex<double>(0,6.38951e-18), std::complex<double>(9.58427e-18,5.11161e-17), std::complex<double>(-0.0996506,-6.38951e-18), 
      std::complex<double>(-7.83714e-18,0), std::complex<double>(1.41876e-33,1.24596e-16), std::complex<double>(-5.11161e-17,1.41876e-33), std::complex<double>(-1.91685e-17,3.54689e-34), std::complex<double>(0.93285,-5.71206e-17), std::complex<double>(0,5.09164e-18), std::complex<double>(2.87528e-17,2.55581e-17), std::complex<double>(-7.78722e-18,0), std::complex<double>(-0.0996506,6.38951e-18), std::complex<double>(0.473301,-2.11482e-17), std::complex<double>(3.80667e-16,-5.28704e-18), 
      std::complex<double>(2.11482e-17,-2.11482e-17), std::complex<double>(-4.22963e-17,2.64352e-17), std::complex<double>(1.08811,-6.92094e-34), std::complex<double>(-2.95657e-33,-2.11482e-17), std::complex<double>(-1.31476e-33,4.22963e-17), std::complex<double>(3.35836e-34,1.26889e-16), std::complex<double>(-6.38951e-18,-3.99345e-18), std::complex<double>(-0.0996506,6.38951e-18), std::complex<double>(1.2779e-17,-4.47266e-17), std::complex<double>(2.87528e-17,-4.47266e-17), 
      std::complex<double>(-5.32034e-34,-6.38951e-18), std::complex<double>(0.93285,-5.71206e-17), std::complex<double>(-1.2779e-17,1.41876e-33), std::complex<double>(2.55581e-17,-5.67503e-33), std::complex<double>(0,-6.38951e-18), std::complex<double>(1.2779e-17,-4.47266e-17), std::complex<double>(-0.0996506,1.91685e-17), std::complex<double>(9.98362e-19,0), std::complex<double>(1.41876e-33,-1.91685e-17), std::complex<double>(0,-1.41876e-33), std::complex<double>(0.93285,-5.71206e-17), 
      std::complex<double>(1.02232e-16,7.80316e-33), std::complex<double>(-6.38951e-18,-3.5941e-18), std::complex<double>(2.87528e-17,-1.2779e-17), std::complex<double>(-5.34123e-18,0), std::complex<double>(-0.0996506,0), std::complex<double>(0,-3.51423e-17), std::complex<double>(1.2779e-17,-4.96565e-33), std::complex<double>(8.30637e-17,-2.12814e-33), std::complex<double>(0.93285,-5.71206e-17)
    };

    std::vector<std::complex<double> > becp = {
      std::complex<double>(0.714037,0), std::complex<double>(4.0926e-34,4.77049e-18), std::complex<double>(-6.74075e-34,1.9082e-17), std::complex<double>(3.85186e-34,3.38271e-17), std::complex<double>(0.310588,-1.38778e-17), std::complex<double>(-1.73472e-18,-1.0842e-18), std::complex<double>(0,-1.73472e-18), std::complex<double>(-1.73472e-18,-9.75782e-19), std::complex<double>(0,1.38778e-17), std::complex<double>(0.253264,-1.5508e-17), std::complex<double>(3.46945e-18,-7.70372e-34), std::complex<double>(-1.38778e-17,3.85186e-34),
      std::complex<double>(2.498e-16,-3.46945e-18), std::complex<double>(-0.0270547,1.73472e-18), std::complex<double>(3.46945e-18,-1.21431e-17), std::complex<double>(7.80626e-18,-3.46945e-18), std::complex<double>(0,2.77556e-17), std::complex<double>(0,-7.70372e-34), std::complex<double>(0.253264,-1.5508e-17), std::complex<double>(-5.20417e-18,9.62965e-35), std::complex<double>(1.38778e-17,-1.38778e-17), std::complex<double>(3.46945e-18,-1.21431e-17), std::complex<double>(-0.0270547,5.20417e-18),
      std::complex<double>(-1.45012e-18,0), std::complex<double>(0,-4.16334e-17), std::complex<double>(-1.04083e-17,3.85186e-34), std::complex<double>(-2.77556e-17,9.62965e-35), std::complex<double>(0.253264,-1.5508e-17), std::complex<double>(-2.77556e-17,1.73472e-17), std::complex<double>(7.80626e-18,-1.21431e-17), std::complex<double>(2.71051e-19,0), std::complex<double>(-0.0270547,0), std::complex<double>(0.310588,1.38778e-17), std::complex<double>(1.73472e-18,7.58942e-19), std::complex<double>(0,1.73472e-18),
      std::complex<double>(0,1.38236e-18), std::complex<double>(0.714037,-4.54164e-34), std::complex<double>(-1.44445e-34,-1.73472e-18), std::complex<double>(3.85186e-34,-5.20417e-18), std::complex<double>(0,-9.54098e-18), std::complex<double>(-2.498e-16,-3.46945e-18), std::complex<double>(-0.0270547,5.20417e-18), std::complex<double>(2.60209e-18,1.38778e-17), std::complex<double>(7.80626e-18,6.93889e-18), std::complex<double>(-1.94015e-33,-1.38778e-17), std::complex<double>(0.253264,-1.5508e-17),
      std::complex<double>(0,-3.85186e-34), std::complex<double>(3.46945e-18,-1.34815e-33), std::complex<double>(-1.38778e-17,-1.38778e-17), std::complex<double>(3.46945e-18,1.38778e-17), std::complex<double>(-0.0270547,-1.73472e-18), std::complex<double>(-2.11419e-18,0), std::complex<double>(-8.62772e-34,2.77556e-17), std::complex<double>(-3.46945e-18,3.85186e-34), std::complex<double>(0.253264,-1.5508e-17), std::complex<double>(2.25514e-17,-5.77779e-34), std::complex<double>(2.77556e-17,1.73472e-17),
      std::complex<double>(7.80626e-18,1.04083e-17), std::complex<double>(-2.12775e-18,0), std::complex<double>(-0.0270547,1.73472e-18), std::complex<double>(2.20382e-34,8.32667e-17), std::complex<double>(6.93889e-18,-1.54074e-33), std::complex<double>(2.77556e-17,2.11852e-33), std::complex<double>(0.253264,-1.5508e-17) 
    };

    const psi::DEVICE_CPU * cpu_ctx = {};
    const psi::DEVICE_GPU * gpu_ctx = {};

    void SetUp() override {
    }
    void TearDown() override {
    }

    using nonlocal_cpu_op = hamilt::nonlocal_pw_op<double, psi::DEVICE_CPU>;
    using nonlocal_gpu_op = hamilt::nonlocal_pw_op<double, psi::DEVICE_GPU>;
    using set_memory_complex_double_op = psi::memory::set_memory_op<std::complex<double>, psi::DEVICE_GPU>;
    using delete_memory_double_op = psi::memory::delete_memory_op<double, psi::DEVICE_GPU>;
    using delete_memory_complex_double_op = psi::memory::delete_memory_op<std::complex<double>, psi::DEVICE_GPU>;
    using resize_memory_double_op = psi::memory::resize_memory_op<double, psi::DEVICE_GPU>;
    using resize_memory_complex_double_op = psi::memory::resize_memory_op<std::complex<double>, psi::DEVICE_GPU>;
    using syncmem_d_h2d_op = psi::memory::synchronize_memory_op<double, psi::DEVICE_GPU, psi::DEVICE_CPU>;
    using syncmem_cd_h2d_op = psi::memory::synchronize_memory_op<std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
    using syncmem_cd_d2h_op = psi::memory::synchronize_memory_op<std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_GPU>;
};

TEST_F(TestModuleHamiltNonlocal, nonlocal_pw_op_cpu)
{
  sum = 0; iat = 0;
  std::vector<std::complex<double>> ps(expected_ps.size(), std::complex<double>(0.0, 0.0));
  nonlocal_cpu_op()(
      cpu_ctx, 
      l1, l2, l3, 
      sum, iat, spin, nkb, 
      deeq_x, deeq_y, deeq_z,
      deeq.data(),
      ps.data(), becp.data());
  for (int ii = 0; ii < ps.size(); ii++) {
    EXPECT_LT(fabs(ps[ii] - expected_ps[ii]), 5 * 1e-6);
  }
  EXPECT_EQ(sum, expected_sum);
  EXPECT_EQ(iat, expected_iat);
}

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
TEST_F(TestModuleHamiltNonlocal, nonlocal_pw_op_gpu)
{
  sum = 0; iat = 0;
  double* deeq_dev = NULL;
  std::complex<double>* ps_dev = NULL, * becp_dev = NULL;
  std::vector<std::complex<double>> ps(expected_ps.size(), std::complex<double>(0.0, 0.0));
  resize_memory_double_op()(gpu_ctx, deeq_dev, deeq.size());
  resize_memory_complex_double_op()(gpu_ctx, ps_dev, ps.size());
  resize_memory_complex_double_op()(gpu_ctx, becp_dev, becp.size());
  syncmem_d_h2d_op()(gpu_ctx, cpu_ctx, deeq_dev, deeq.data(), deeq.size());
  syncmem_cd_h2d_op()(gpu_ctx, cpu_ctx, ps_dev, ps.data(), ps.size());
  syncmem_cd_h2d_op()(gpu_ctx, cpu_ctx, becp_dev, becp.data(), becp.size());
  nonlocal_gpu_op()(
      gpu_ctx, 
      l1, l2, l3, 
      sum, iat, spin, nkb, 
      deeq_x, deeq_y, deeq_z,
      deeq_dev,
      ps_dev, becp_dev);

  syncmem_cd_d2h_op()(cpu_ctx, gpu_ctx, ps.data(), ps_dev, ps.size());
  for (int ii = 0; ii < ps.size(); ii++) {
    EXPECT_LT(fabs(ps[ii] - expected_ps[ii]), 5 * 1e-6);
  }
  EXPECT_EQ(sum, expected_sum);
  EXPECT_EQ(iat, expected_iat);
  delete_memory_double_op()(gpu_ctx, deeq_dev);
  delete_memory_complex_double_op()(gpu_ctx, ps_dev);
  delete_memory_complex_double_op()(gpu_ctx, becp_dev);
}
#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM