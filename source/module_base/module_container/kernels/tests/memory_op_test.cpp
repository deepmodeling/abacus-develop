#include <complex>
#include <iostream>
#include <gtest/gtest.h>
#include "../memory_op.h"
#if __UT_USE_CUDA || __UT_USE_ROCM
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/complex.h>
#include <thrust/device_ptr.h>
#include <thrust/device_free.h>
#include <thrust/host_vector.h>
#include <thrust/device_malloc.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#endif

class TestModulePsiMemory : public ::testing::Test
{
  protected:
    const std::vector<double> xx = {
        -0.65412617, -0.74208893, -2.21731157,  0.42540039, 
        3.36373004,  -2.51647562, -2.985111  , -0.53251562, 
        0.37908265,   0.81605825,  1.66281318,  2.71761869, 
        2.2010268 ,   0.65498149,  1.51153638,  0.71501482, 
        0.53546578,   1.4564317 , -2.36701143,  1.23009056, 
        3.41302551,  -2.3175205 , -0.27628221, -1.35701656
    };

    const std::vector<std::complex<double> > z_xx = {
        {-0.65412617, -0.74208893}, {-2.21731157,  0.42540039},
        {3.36373004,  -2.51647562}, {-2.985111  , -0.53251562},
        {0.37908265,   0.81605825}, { 1.66281318,  2.71761869},
        {2.2010268 ,   0.65498149}, { 1.51153638,  0.71501482},
        {0.53546578,   1.4564317 }, {-2.36701143,  1.23009056},
        {3.41302551,  -2.3175205 }, {-0.27628221, -1.35701656}
    };

    const int z_dim = z_xx.size();

    const psi::DEVICE_CPU * cpu_ctx = {};
    const psi::DEVICE_GPU * gpu_ctx = {};

    void SetUp() override {
    }
    void TearDown() override {
    }

    using set_memory_double_cpu_op = psi::memory::set_memory_op<double, psi::DEVICE_CPU>;
    using set_memory_complex_double_cpu_op = psi::memory::set_memory_op<std::complex<double>, psi::DEVICE_CPU>;
    using resize_memory_double_cpu_op = psi::memory::resize_memory_op<double, psi::DEVICE_CPU>;
    using resize_memory_comlex_double_cpu_op = psi::memory::resize_memory_op<std::complex<double>, psi::DEVICE_CPU>;
    using synchronize_memory_double_cpu_to_cpu_op = psi::memory::synchronize_memory_op<double, psi::DEVICE_CPU, psi::DEVICE_CPU>;
    using synchronize_memory_complex_double_cpu_to_cpu_op = psi::memory::synchronize_memory_op<std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_CPU>;
    using delete_memory_double_cpu_op = psi::memory::delete_memory_op<double, psi::DEVICE_CPU>;
    using delete_memory_complex_double_cpu_op = psi::memory::delete_memory_op<std::complex<double>, psi::DEVICE_CPU>;

#if __UT_USE_CUDA || __UT_USE_ROCM
    using set_memory_double_gpu_op = psi::memory::set_memory_op<double, psi::DEVICE_GPU>;
    using set_memory_complex_double_gpu_op = psi::memory::set_memory_op<std::complex<double>, psi::DEVICE_GPU>;
    using resize_memory_double_gpu_op = psi::memory::resize_memory_op<double, psi::DEVICE_GPU>;
    using resize_memory_comlex_double_gpu_op = psi::memory::resize_memory_op<std::complex<double>, psi::DEVICE_GPU>;
    using synchronize_memory_double_cpu_to_gpu_op = psi::memory::synchronize_memory_op<double, psi::DEVICE_GPU, psi::DEVICE_CPU>;
    using synchronize_memory_double_gpu_to_cpu_op = psi::memory::synchronize_memory_op<double, psi::DEVICE_CPU, psi::DEVICE_GPU>;
    using synchronize_memory_double_gpu_to_gpu_op = psi::memory::synchronize_memory_op<double, psi::DEVICE_GPU, psi::DEVICE_GPU>;
    using synchronize_memory_complex_double_cpu_to_gpu_op = psi::memory::synchronize_memory_op<std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
    using synchronize_memory_complex_double_gpu_to_cpu_op = psi::memory::synchronize_memory_op<std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_GPU>;
    using synchronize_memory_complex_double_gpu_to_gpu_op = psi::memory::synchronize_memory_op<std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_GPU>;
    using delete_memory_double_gpu_op = psi::memory::delete_memory_op<double, psi::DEVICE_GPU>;
    using delete_memory_complex_double_gpu_op = psi::memory::delete_memory_op<std::complex<double>, psi::DEVICE_GPU>;
#endif // __UT_USE_CUDA || __UT_USE_ROCM
};

TEST_F(TestModulePsiMemory, set_memory_op_double_cpu)
{
  std::vector<double> v_xx = xx;
  set_memory_double_cpu_op()(cpu_ctx, v_xx.data(), 0, xx.size());
  for (int ii = 0; ii < xx.size(); ii++) {
    EXPECT_EQ(v_xx[ii], 0.0);
  }
}

TEST_F(TestModulePsiMemory, set_memory_op_complex_double_cpu)
{
  std::vector<std::complex<double>> vz_xx = z_xx;
  set_memory_complex_double_cpu_op()(cpu_ctx, vz_xx.data(), 0, z_xx.size());
  for (int ii = 0; ii < z_xx.size(); ii++) {
    EXPECT_EQ(vz_xx[ii], std::complex<double>(0.0, 0.0));
  }
}

TEST_F(TestModulePsiMemory, resize_memory_op_double_cpu)
{
  double* xx_tmp = NULL;
  resize_memory_double_cpu_op()(cpu_ctx, xx_tmp, xx.size());
  for (int ii = 0; ii < xx.size(); ii++) {
    xx_tmp[ii] = xx[ii];
  }
  for (int ii = 0; ii < xx.size(); ii++) {
    EXPECT_EQ(xx_tmp[ii], xx[ii]);
  }
  free(xx_tmp);
}

TEST_F(TestModulePsiMemory, resize_memory_op_comlex_double_cpu)
{
  std::complex<double>* z_xx_tmp = NULL;
  resize_memory_comlex_double_cpu_op()(cpu_ctx, z_xx_tmp, z_xx.size());
  for (int ii = 0; ii < z_xx.size(); ii++) {
    z_xx_tmp[ii] = z_xx[ii];
  }
  for (int ii = 0; ii < z_xx.size(); ii++) {
    EXPECT_EQ(z_xx_tmp[ii], z_xx[ii]);
  }
  free(z_xx_tmp);
}

TEST_F(TestModulePsiMemory, synchronize_memory_op_double_cpu_to_cpu)
{
  std::vector<double> h_xx(xx.size(), 0);
  synchronize_memory_double_cpu_to_cpu_op()(cpu_ctx, cpu_ctx, h_xx.data(), xx.data(), xx.size());
  for (int ii = 0; ii < z_xx.size(); ii++) {
    EXPECT_EQ(h_xx[ii], xx[ii]);
  }
}

TEST_F(TestModulePsiMemory, synchronize_memory_op_complex_double_cpu_to_cpu)
{
  std::vector<std::complex<double>> hz_xx(z_xx.size(), std::complex<double>(0, 0));
  synchronize_memory_complex_double_cpu_to_cpu_op()(cpu_ctx, cpu_ctx, hz_xx.data(), z_xx.data(), z_xx.size());
  for (int ii = 0; ii < z_xx.size(); ii++) {
    EXPECT_EQ(hz_xx[ii], z_xx[ii]);
  }
}

TEST_F(TestModulePsiMemory, delete_memory_op_double_cpu)
{
  double * h_xx = (double*)malloc(sizeof(double) * xx.size());
  delete_memory_double_cpu_op()(cpu_ctx, h_xx);
}

TEST_F(TestModulePsiMemory, delete_memory_op_complex_double_cpu)
{
  std::complex<double> * hz_xx = (std::complex<double>*)malloc(sizeof(std::complex<double>) * z_xx.size());
  delete_memory_complex_double_cpu_op()(cpu_ctx, hz_xx);
}


#if __UT_USE_CUDA || __UT_USE_ROCM
TEST_F(TestModulePsiMemory, set_memory_op_double_gpu)
{
  thrust::device_ptr<double> d_xx = thrust::device_malloc<double>(xx.size());
  thrust::copy(xx.begin(), xx.end(), d_xx);
  set_memory_double_gpu_op()(gpu_ctx, thrust::raw_pointer_cast(d_xx), 0, xx.size());
  thrust::host_vector<double> h_xx(xx.size());
  thrust::copy(d_xx, d_xx + xx.size(), h_xx.begin());
  for (int ii = 0; ii < xx.size(); ii++) {
    EXPECT_EQ(h_xx[ii], 0.0);
  }
}

TEST_F(TestModulePsiMemory, set_memory_op_complex_double_gpu)
{
  thrust::device_ptr<std::complex<double>> dz_xx = thrust::device_malloc<std::complex<double>>(z_xx.size());
  thrust::copy(z_xx.begin(), z_xx.end(), dz_xx);
  set_memory_complex_double_gpu_op()(gpu_ctx, thrust::raw_pointer_cast(dz_xx), 0, z_xx.size());
  thrust::host_vector<std::complex<double>> h_xx(z_xx.size());
  thrust::copy(dz_xx, dz_xx + z_xx.size(), h_xx.begin());
  for (int ii = 0; ii < z_xx.size(); ii++) {
    EXPECT_EQ(h_xx[ii], std::complex<double>(0.0, 0.0));
  }
}

TEST_F(TestModulePsiMemory, resize_memory_op_double_gpu)
{
  double* xx_tmp = NULL;
  resize_memory_double_gpu_op()(gpu_ctx, xx_tmp, xx.size());

  thrust::device_ptr<double> d_xx(xx_tmp);
  thrust::copy(xx.begin(), xx.end(), d_xx);

  thrust::host_vector<double> h_xx(xx.size());
  thrust::copy(d_xx, d_xx + xx.size(), h_xx.begin());
  for (int ii = 0; ii < xx.size(); ii++) {
    EXPECT_EQ(h_xx[ii], xx[ii]);
  }
  thrust::device_free(d_xx);
}

TEST_F(TestModulePsiMemory, resize_memory_op_complex_double_gpu)
{
  std::complex<double>* z_xx_tmp = NULL;
  resize_memory_comlex_double_gpu_op()(gpu_ctx, z_xx_tmp, z_xx.size());

  thrust::device_ptr<std::complex<double>> dz_xx(z_xx_tmp);
  thrust::copy(z_xx.begin(), z_xx.end(), dz_xx);

  thrust::host_vector<std::complex<double>> h_z_xx(z_xx.size());
  thrust::copy(dz_xx, dz_xx + z_xx.size(), h_z_xx.begin());
  for (int ii = 0; ii < z_xx.size(); ii++) {
    EXPECT_EQ(h_z_xx[ii], z_xx[ii]);
  }
  thrust::device_free(dz_xx);
}

TEST_F(TestModulePsiMemory, synchronize_memory_op_double_cpu_to_gpu)
{
  thrust::device_ptr<double> d_xx = thrust::device_malloc<double>(xx.size());
  std::vector<double> hv_xx(xx.size(), 0);
  thrust::copy(hv_xx.begin(), hv_xx.end(), d_xx);
  synchronize_memory_double_cpu_to_gpu_op()(
    gpu_ctx, 
    cpu_ctx, 
    thrust::raw_pointer_cast(d_xx), 
    xx.data(), 
    xx.size());

  thrust::host_vector<double> h_xx(xx.size());
  thrust::copy(d_xx, d_xx + xx.size(), h_xx.begin());
  for (int ii = 0; ii < xx.size(); ii++) {
    EXPECT_EQ(h_xx[ii], xx[ii]);
  }
  thrust::device_free(d_xx);
}

TEST_F(TestModulePsiMemory, synchronize_memory_op_double_gpu_to_cpu)
{
  thrust::device_ptr<double> d_xx = thrust::device_malloc<double>(xx.size());
  thrust::copy(xx.begin(), xx.end(), d_xx);
  thrust::host_vector<double> h_xx(xx.size());
  synchronize_memory_double_gpu_to_cpu_op()(
    cpu_ctx, 
    gpu_ctx, 
    thrust::raw_pointer_cast(h_xx.data()), 
    thrust::raw_pointer_cast(d_xx), 
    xx.size());

  for (int ii = 0; ii < xx.size(); ii++) {
    EXPECT_EQ(h_xx[ii], xx[ii]);
  }
  thrust::device_free(d_xx);
}

TEST_F(TestModulePsiMemory, synchronize_memory_op_double_gpu_to_gpu)
{
  thrust::device_ptr<double> d1_xx = thrust::device_malloc<double>(xx.size());
  thrust::device_ptr<double> d2_xx = thrust::device_malloc<double>(xx.size());
  thrust::copy(xx.begin(), xx.end(), d1_xx);
  synchronize_memory_double_gpu_to_gpu_op()(
    gpu_ctx, 
    gpu_ctx, 
    thrust::raw_pointer_cast(d2_xx), 
    thrust::raw_pointer_cast(d1_xx), 
    xx.size());

  thrust::host_vector<double> h_xx(xx.size());
  thrust::copy(d2_xx, d2_xx + xx.size(), h_xx.begin());
  for (int ii = 0; ii < xx.size(); ii++) {
    EXPECT_EQ(h_xx[ii], xx[ii]);
  }
  thrust::device_free(thrust::device_ptr<double>(d1_xx));
  thrust::device_free(thrust::device_ptr<double>(d2_xx));
}

TEST_F(TestModulePsiMemory, synchronize_memory_op_complex_double_cpu_to_gpu)
{
  thrust::device_ptr<std::complex<double>> dz_xx = thrust::device_malloc<std::complex<double>>(z_xx.size());
  std::vector<std::complex<double>> hvz_xx(z_xx.size(), 0);
  thrust::copy(hvz_xx.begin(), hvz_xx.end(), dz_xx);
  synchronize_memory_complex_double_cpu_to_gpu_op()(
    gpu_ctx, 
    cpu_ctx, 
    thrust::raw_pointer_cast(dz_xx), 
    z_xx.data(), 
    z_xx.size());

  thrust::host_vector<std::complex<double>> hz_xx(z_xx.size());
  thrust::copy(dz_xx, dz_xx + z_xx.size(), hz_xx.begin());
  for (int ii = 0; ii < z_xx.size(); ii++) {
    EXPECT_EQ(hz_xx[ii], z_xx[ii]);
  }
  thrust::device_free(dz_xx);
}

TEST_F(TestModulePsiMemory, synchronize_memory_op_complex_double_gpu_to_cpu)
{
  thrust::device_ptr<std::complex<double>> dz_xx = thrust::device_malloc<std::complex<double>>(z_xx.size());
  thrust::copy(z_xx.begin(), z_xx.end(), dz_xx);
  thrust::host_vector<std::complex<double>> hz_xx(z_xx.size());
  synchronize_memory_complex_double_gpu_to_cpu_op()(
    cpu_ctx, 
    gpu_ctx, 
    thrust::raw_pointer_cast(hz_xx.data()), 
    thrust::raw_pointer_cast(dz_xx), 
    z_xx.size());

  for (int ii = 0; ii < z_xx.size(); ii++) {
    EXPECT_EQ(hz_xx[ii], z_xx[ii]);
  }
  thrust::device_free(dz_xx);
}

TEST_F(TestModulePsiMemory, synchronize_memory_op_complex_double_gpu_to_gpu)
{
  thrust::device_ptr<std::complex<double>> dz1_xx = thrust::device_malloc<std::complex<double>>(z_xx.size());
  thrust::device_ptr<std::complex<double>> dz2_xx = thrust::device_malloc<std::complex<double>>(z_xx.size());
  thrust::copy(z_xx.begin(), z_xx.end(), dz1_xx);
  synchronize_memory_complex_double_gpu_to_gpu_op()(
    gpu_ctx, 
    gpu_ctx, 
    thrust::raw_pointer_cast(dz2_xx), 
    thrust::raw_pointer_cast(dz1_xx), 
    z_xx.size());

  thrust::host_vector<std::complex<double>> h_xx(z_xx.size());
  thrust::copy(dz2_xx, dz2_xx + z_xx.size(), h_xx.begin());
  for (int ii = 0; ii < z_xx.size(); ii++) {
    EXPECT_EQ(h_xx[ii], z_xx[ii]);
  }
  thrust::device_free(thrust::device_ptr<std::complex<double>>(dz1_xx));
  thrust::device_free(thrust::device_ptr<std::complex<double>>(dz2_xx));
}

TEST_F(TestModulePsiMemory, delete_memory_op_double_gpu)
{
  thrust::device_ptr<double> d_xx = thrust::device_malloc<double>(xx.size());
  delete_memory_double_gpu_op()(gpu_ctx, thrust::raw_pointer_cast(d_xx));
}

TEST_F(TestModulePsiMemory, delete_memory_op_complex_double_gpu)
{
  thrust::device_ptr<std::complex<double>> dz_xx = thrust::device_malloc<std::complex<double>>(z_xx.size());
  delete_memory_complex_double_gpu_op()(gpu_ctx, thrust::raw_pointer_cast(dz_xx));
}

#endif // __UT_USE_CUDA || __UT_USE_ROCM
