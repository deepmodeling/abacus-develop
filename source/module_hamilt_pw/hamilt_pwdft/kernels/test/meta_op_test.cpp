#include <vector>
#include <complex>
#include <gtest/gtest.h>
#include "module_psi/kernels/memory_op.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/meta_op.h"

class TestModuleHamiltMeta : public ::testing::Test
{
protected:
    // xx = tf.random.uniform([100], minval=-4, maxval=4, dtype = tf.float64)
    const int ik = 0, pol = 0, npw = 15, npwx = 15;
    const double tpiba = 0.61599855952741045;

    const std::vector<double> kvec_c = {0, 0, 0};
    const std::vector<double> gcar = {1, -1, -1, 0, 0, 0, -1, 1, 1, 1, -1, 1, 0, 0, 2, 0, 0, -2, -1, 1, -1, 1, 1, -1, 0, 2, 0, 2, 0, 0, 1, 1, 1, 0, -2, 0, -1, -1, 1, -1, -1, -1, -2, 0, 0};

    const std::vector<std::complex<double> > in = {{0.073755, 0.215648}, {-0.661877, 0.517527}, {0.0586699, -0.186777}, {0.00957812, -0.0832565}, {-0.0523913, -0.0184799}, {-0.0259417, -0.124468}, {-0.090873, -0.00766968}, {-0.0277014, -0.0220417}, {0.000145658, 0.00325686}, {-0.0254601, 0.018268}, {0.0251554, 0.0204952}, {-0.126029, -0.0193828}, {0.0681057, -0.108652}, {0.0164339, -0.130657}, {0.0972758, -0.0168042}};
    const std::vector<std::complex<double> > expected_out = {{-0.132839, 0.045433}, {-0, 0}, {-0.115054, -0.0361406}, {0.0512859, 0.00590011}, {0, -0}, {0, -0}, {-0.00472451, 0.0559777}, {0.0135776, -0.0170641}, {0, 0}, {-0.0225061, -0.0313667}, {-0.012625, 0.0154957}, {0, -0}, {-0.0669297, -0.041953}, {-0.0804844, -0.0101233}, {-0.0207027, -0.119844}};

    const psi::DEVICE_CPU * cpu_ctx = {};
    const psi::DEVICE_GPU * gpu_ctx = {};

    using meta_cpu_op = hamilt::meta_pw_op<double, psi::DEVICE_CPU>;
    using meta_gpu_op = hamilt::meta_pw_op<double, psi::DEVICE_GPU>;

    using delmem_complex_op = psi::memory::delete_memory_op<std::complex<double>, psi::DEVICE_GPU>;
    using resmem_complex_op = psi::memory::resize_memory_op<std::complex<double>, psi::DEVICE_GPU>;
    using syncmem_complex_h2d_op = psi::memory::synchronize_memory_op<std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
    using syncmem_complex_d2h_op = psi::memory::synchronize_memory_op<std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_GPU>;

    using delmem_var_op = psi::memory::delete_memory_op<double, psi::DEVICE_GPU>;
    using resmem_var_op = psi::memory::resize_memory_op<double, psi::DEVICE_GPU>;
    using syncmem_var_h2d_op = psi::memory::synchronize_memory_op<double, psi::DEVICE_GPU, psi::DEVICE_CPU>;

    void SetUp() override {
    }
    void TearDown() override {
    }
};

TEST_F(TestModuleHamiltMeta, meta_pw_op_cpu)
{
    std::vector<std::complex<double>> res(expected_out.size(), std::complex<double> {0, 0});
    meta_cpu_op()(cpu_ctx, ik, pol, npw, npwx, tpiba, gcar.data(), kvec_c.data(), in.data(), res.data());

    for (int ii = 0; ii < res.size(); ii++) {
        EXPECT_LT(std::abs(res[ii] - expected_out[ii]), 6e-5);
    }
}


#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
TEST_F(TestModuleHamiltMeta, meta_pw_op_gpu)
{
    std::vector<std::complex<double>> res(expected_out.size(), std::complex<double> {0, 0});
    double * d_gcar = nullptr, * d_kvec_c = nullptr;
    std::complex<double>* d_in = nullptr, * d_res = nullptr;
    resmem_var_op()(gpu_ctx, d_gcar, gcar.size());
    resmem_var_op()(gpu_ctx, d_kvec_c, kvec_c.size());
    resmem_complex_op()(gpu_ctx, d_in, in.size());
    resmem_complex_op()(gpu_ctx, d_res, res.size());
    syncmem_var_h2d_op()(gpu_ctx, cpu_ctx, d_gcar, gcar.data(), gcar.size());
    syncmem_var_h2d_op()(gpu_ctx, cpu_ctx, d_kvec_c, kvec_c.data(), kvec_c.size());
    syncmem_complex_h2d_op()(gpu_ctx, cpu_ctx, d_in, in.data(), in.size());
    syncmem_complex_h2d_op()(gpu_ctx, cpu_ctx, d_res, res.data(), res.size());

    meta_gpu_op()(gpu_ctx, ik, pol, npw, npwx, tpiba, d_gcar, d_kvec_c, d_in, d_res);

    syncmem_complex_d2h_op()(cpu_ctx, gpu_ctx, res.data(), d_res, res.size());
    for (int ii = 0; ii < res.size(); ii++) {
        EXPECT_LT(fabs(res[ii] - expected_out[ii]), 6e-5);
    }
    delmem_var_op()(gpu_ctx, d_gcar);
    delmem_var_op()(gpu_ctx, d_kvec_c);
    delmem_complex_op()(gpu_ctx, d_in);
    delmem_complex_op()(gpu_ctx, d_res);
}
#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM