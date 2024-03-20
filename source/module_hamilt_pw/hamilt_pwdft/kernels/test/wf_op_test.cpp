#include <vector>
#include <complex>
#include <gtest/gtest.h>
#include "module_psi/kernels/memory_op.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/wf_op.h"

class TestSrcPWWfMultiDevice : public ::testing::Test
{
protected:
    // xx = tf.random.uniform([100], minval=-4, maxval=4, dtype = tf.float64)

    const psi::DEVICE_CPU * cpu_ctx = {};
    const psi::DEVICE_GPU * gpu_ctx = {};

    int ik = 0, ntype = 1, nx = 12, ny = 12, nz = 12, rho_nx = 12, rho_ny = 12, rho_nz = 12, npw = 59, npwx = 70,
        fftny = 12, eigts1_nc = 25, eigts2_nc = 25, eigts3_nc = 25;

    double TWO_PI = 6.2831853071795862;

    std::vector<int> atom_na {2};
    std::vector<int> igl2isz {10, 11, 0, 1, 2, 23, 12, 13, 14, 24, 25, 71, 60, 82, 83, 72, 73, 95, 84, 85, 86, 107, 96, 97, 98, 108, 109, 110, 155, 144, 145, 156, 157, 168, 169, 170, 181, 182, 323, 312, 358, 359, 370, 371, 360, 382, 383, 372, 373, 395, 384, 385, 430, 431, 420, 442, 443, 432, 433, -1074790400, 0, 0, 0, 1072693248, 0, 0, 0, -1074790400, 0, 1072693248, 10, 11, 0, 1, 23, 12, 13, 70, 71, 60, 82, 83, 72, 73, 95, 84, 85, 107, 96, 97, 98, 109, 155, 144, 145, 169, 285, 286, 297, 298, 299, 310, 311, 322, 323, 312, 345, 346, 347, 357, 358, 359, 348, 369, 370, 371, 360, 382, 383, 372, 373, 395, 384, 385, 418, 419, 429, 430, 431, 420, 441, 442, 443, 432, 433, -1074790400, 0, 0, 0, -1074790400, 10, 11, 0, 1, 2, 23, 12, 13, 14, 24, 25, 70, 71, 60, 82, 83, 72, 73, 95, 84, 85, 86, 107, 96, 97, 98, 108, 109, 155, 144, 145, 156, 157, 168, 169, 298, 299, 310, 311, 322, 323, 312, 346, 347, 357, 358, 359, 348, 369, 370, 371, 360, 382, 383, 372, 373, 395, 384, 385, 418, 419, 429, 430, 431, 420, 441, 442, 443, 432, 433, 10, 11, 0, 1, 23, 12, 13, 58, 59, 69, 70, 71, 60, 81, 82, 83, 72, 73, 94, 95, 84, 85, 107, 96, 97, 142, 143, 132, 154, 155, 144, 145, 167, 156, 157, 168, 169, 215, 204, 345, 346, 357, 358, 359, 370, 371, 382, 383, 372, 417, 418, 419, 429, 430, 431, 420, 441, 442, 443, 432, 0, 1072693248, 0, 0, 0, 0, 0, 0, 0, 1072693248, 11, 0, 1, 2, 23, 12, 13, 14, 70, 71, 60, 61, 82, 83, 72, 73, 74, 95, 84, 85, 86, 96, 97, 98, 109, 155, 144, 145, 169, 286, 287, 298, 299, 310, 311, 300, 322, 323, 312, 313, 346, 347, 357, 358, 359, 348, 370, 371, 360, 361, 382, 383, 372, 373, 374, 395, 384, 385, 418, 419, 408, 430, 431, 420, 421, 442, 443, 432, 433, 0, 10, 11, 0, 1, 2, 23, 12, 13, 14, 70, 71, 60, 82, 83, 72, 73, 95, 84, 85, 86, 96, 97, 98, 109, 155, 144, 145, 169, 286, 297, 298, 299, 310, 311, 322, 323, 312, 345, 346, 347, 357, 358, 359, 348, 369, 370, 371, 360, 361, 382, 383, 372, 373, 395, 384, 385, 418, 419, 429, 430, 431, 420, 421, 442, 443, 432, 433, 0, 0, 0, 10, 11, 0, 1, 23, 12, 13, 70, 71, 60, 82, 83, 72, 73, 94, 95, 84, 85, 86, 107, 96, 97, 98, 108, 109, 143, 132, 154, 155, 144, 145, 167, 156, 157, 168, 169, 170, 181, 204, 323, 346, 357, 358, 359, 370, 371, 360, 382, 383, 372, 373, 395, 384, 418, 419, 429, 430, 431, 420, 441, 442, 443, 432, 433, 0, 1072693248, 0, 1072693248, 0, 0, 10, 11, 0, 1, 2, 23, 12, 13, 14, 24, 25, 71, 60, 82, 83, 72, 73, 95, 84, 85, 86, 96, 97, 98, 109, 110, 155, 144, 145, 169, 286, 298, 299, 310, 311, 300, 322, 323, 312, 313, 324, 346, 347, 357, 358, 359, 348, 370, 371, 360, 361, 382, 383, 372, 373, 374, 395, 384, 385, 419, 430, 431, 420, 421, 442, 443, 432, 433, 0, 1072693248};
    std::vector<int> is2fftixy {0, 1, 2, 3, 9, 10, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, 27, 35, 36, 37, 38, 39, 108, 117, 118, 119, 120, 121, 129, 130, 131, 132, 133, 134, 141, 142, 143};

    std::vector<double> kvec_c = {0, 0, 0, 0.75, 0.75, 0.75, 0.5, 0.5};
    std::vector<double> atom_tau = {0, 0, 0, 0.25, 0.25, 0.25};

    std::vector<std::complex<double>> eigts1 = {{1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -7.34788e-16}, {-2.44991e-15, -1}, {-1, 6.12323e-16}, {5.51091e-16, 1}, {1, -4.89859e-16}, {-4.28626e-16, -1}, {-1, 3.67394e-16}, {3.06162e-16, 1}, {1, -2.44929e-16}, {-1.83697e-16, -1}, {-1, 1.22465e-16}, {6.12323e-17, 1}, {1, -0}, {6.12323e-17, -1}, {-1, -1.22465e-16}, {-1.83697e-16, 1}, {1, 2.44929e-16}, {3.06162e-16, -1}, {-1, -3.67394e-16}, {-4.28626e-16, 1}, {1, 4.89859e-16}, {5.51091e-16, -1}, {-1, -6.12323e-16}, {-2.44991e-15, 1}, {1, 7.34788e-16}};
    std::vector<std::complex<double>> eigts2 = {{1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -7.34788e-16}, {-2.44991e-15, -1}, {-1, 6.12323e-16}, {5.51091e-16, 1}, {1, -4.89859e-16}, {-4.28626e-16, -1}, {-1, 3.67394e-16}, {3.06162e-16, 1}, {1, -2.44929e-16}, {-1.83697e-16, -1}, {-1, 1.22465e-16}, {6.12323e-17, 1}, {1, -0}, {6.12323e-17, -1}, {-1, -1.22465e-16}, {-1.83697e-16, 1}, {1, 2.44929e-16}, {3.06162e-16, -1}, {-1, -3.67394e-16}, {-4.28626e-16, 1}, {1, 4.89859e-16}, {5.51091e-16, -1}, {-1, -6.12323e-16}, {-2.44991e-15, 1}, {1, 7.34788e-16}};
    std::vector<std::complex<double>> eigts3 = {{1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -7.34788e-16}, {-2.44991e-15, -1}, {-1, 6.12323e-16}, {5.51091e-16, 1}, {1, -4.89859e-16}, {-4.28626e-16, -1}, {-1, 3.67394e-16}, {3.06162e-16, 1}, {1, -2.44929e-16}, {-1.83697e-16, -1}, {-1, 1.22465e-16}, {6.12323e-17, 1}, {1, -0}, {6.12323e-17, -1}, {-1, -1.22465e-16}, {-1.83697e-16, 1}, {1, 2.44929e-16}, {3.06162e-16, -1}, {-1, -3.67394e-16}, {-4.28626e-16, 1}, {1, 4.89859e-16}, {5.51091e-16, -1}, {-1, -6.12323e-16}, {-2.44991e-15, 1}, {1, 7.34788e-16}};
    std::vector<std::complex<double>> expected_sk = {{1, 0}, {1, 0}, {1, -0}, {1, -0}, {1, -0}, {1, 0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, -0}, {1, -0}, {1, -0}, {1, 0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, 0}, {1, 0}, {1, 0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, -0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {-1, 1.22465e-16}, {6.12323e-17, 1}, {1, -0}, {6.12323e-17, -1}, {-1, -1.22465e-16}, {1, 0}, {6.12323e-17, -1}, {-1, -1.22465e-16}, {-1.83697e-16, 1}, {-1, -1.22465e-16}, {-1.83697e-16, 1}, {-1.83697e-16, -1}, {-1, 1.22465e-16}, {-1.83697e-16, -1}, {-1, 1.22465e-16}, {6.12323e-17, 1}, {1, 0}, {1, 0}, {6.12323e-17, -1}, {-1, -1.22465e-16}, {-1.83697e-16, 1}, {6.12323e-17, -1}, {-1, -1.22465e-16}, {-1.83697e-16, 1}, {1, 2.44929e-16}, {-1.83697e-16, 1}, {1, 2.44929e-16}, {3.06162e-16, -1}, {6.12323e-17, 1}, {1, 0}, {6.12323e-17, -1}, {-1, -1.22465e-16}, {-1.83697e-16, 1}, {-1.83697e-16, 1}, {1, 2.44929e-16}, {3.06162e-16, -1}, {3.06162e-16, -1}, {-1, -3.67394e-16}, {-1.83697e-16, -1}, {-1, 1.22465e-16}, {-1, 3.67394e-16}, {3.06162e-16, 1}, {3.06162e-16, 1}, {1, -2.44929e-16}, {-1.83697e-16, -1}, {-1.83697e-16, -1}, {-1, 1.22465e-16}, {6.12323e-17, 1}, {1, 0}, {6.12323e-17, 1}, {1, 0}, {6.12323e-17, -1}, {3.06162e-16, 1}, {1, -2.44929e-16}, {-1.83697e-16, -1}, {1, -2.44929e-16}, {-1.83697e-16, -1}, {-1, 1.22465e-16}, {6.12323e-17, 1}};

    using delmem_complex_op = psi::memory::delete_memory_op<std::complex<double>, psi::DEVICE_GPU>;
    using resmem_complex_op = psi::memory::resize_memory_op<std::complex<double>, psi::DEVICE_GPU>;
    using syncmem_complex_h2d_op = psi::memory::synchronize_memory_op<std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
    using syncmem_complex_d2h_op = psi::memory::synchronize_memory_op<std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_GPU>;

    using delmem_var_op = psi::memory::delete_memory_op<double, psi::DEVICE_GPU>;
    using resmem_var_op = psi::memory::resize_memory_op<double, psi::DEVICE_GPU>;
    using syncmem_var_h2d_op = psi::memory::synchronize_memory_op<double, psi::DEVICE_GPU, psi::DEVICE_CPU>;
    using syncmem_var_d2h_op = psi::memory::synchronize_memory_op<double, psi::DEVICE_CPU, psi::DEVICE_GPU>;

    using delmem_int_op = psi::memory::delete_memory_op<int, psi::DEVICE_GPU>;
    using resmem_int_op = psi::memory::resize_memory_op<int, psi::DEVICE_GPU>;
    using syncmem_int_h2d_op = psi::memory::synchronize_memory_op<int, psi::DEVICE_GPU, psi::DEVICE_CPU>;

    void SetUp() override {
    }
    void TearDown() override {
    }
};

TEST_F(TestSrcPWWfMultiDevice, cal_sk_op_cpu)
{
    std::vector<std::complex<double>> sk(expected_sk.size(), 0);
    hamilt::cal_sk_op<double, psi::DEVICE_CPU>()(cpu_ctx,
                                                 ik,
                                                 ntype,
                                                 nx,
                                                 ny,
                                                 nz,
                                                 rho_nx,
                                                 rho_ny,
                                                 rho_nz,
                                                 npw,
                                                 npwx,
                                                 fftny,
                                                 eigts1_nc,
                                                 eigts2_nc,
                                                 eigts3_nc,
                                                 atom_na.data(),
                                                 igl2isz.data(),
                                                 is2fftixy.data(),
                                                 TWO_PI,
                                                 kvec_c.data(),
                                                 atom_tau.data(),
                                                 eigts1.data(),
                                                 eigts2.data(),
                                                 eigts3.data(),
                                                 sk.data());

    for (int ii = 0; ii < sk.size(); ii++) {
        EXPECT_LT(fabs(sk[ii] - expected_sk[ii]), 6e-5);
    }
}


#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
TEST_F(TestSrcPWWfMultiDevice, cal_sk_op_gpu)
{
    std::vector<std::complex<double>> sk(expected_sk.size(), 0);
    int * d_atom_na = nullptr, * d_igl2isz = nullptr, * d_is2fftixy = nullptr;
    double * d_kvec_c = nullptr, * d_atom_tau = nullptr;
    std::complex<double> * d_sk = nullptr, * d_eigts1 = nullptr, * d_eigts2 = nullptr, * d_eigts3 = nullptr;

    resmem_int_op()(gpu_ctx, d_atom_na, atom_na.size());
    resmem_int_op()(gpu_ctx, d_igl2isz, igl2isz.size());
    resmem_int_op()(gpu_ctx, d_is2fftixy, is2fftixy.size());
    syncmem_int_h2d_op()(gpu_ctx, cpu_ctx, d_atom_na, atom_na.data(), atom_na.size());
    syncmem_int_h2d_op()(gpu_ctx, cpu_ctx, d_igl2isz, igl2isz.data(), igl2isz.size());
    syncmem_int_h2d_op()(gpu_ctx, cpu_ctx, d_is2fftixy, is2fftixy.data(), is2fftixy.size());

    resmem_var_op()(gpu_ctx, d_kvec_c, kvec_c.size());
    resmem_var_op()(gpu_ctx, d_atom_tau, atom_tau.size());
    syncmem_var_h2d_op()(gpu_ctx, cpu_ctx, d_kvec_c, kvec_c.data(), kvec_c.size());
    syncmem_var_h2d_op()(gpu_ctx, cpu_ctx, d_atom_tau, atom_tau.data(), atom_tau.size());

    resmem_complex_op()(gpu_ctx, d_sk, sk.size());
    resmem_complex_op()(gpu_ctx, d_eigts1, eigts1.size());
    resmem_complex_op()(gpu_ctx, d_eigts2, eigts2.size());
    resmem_complex_op()(gpu_ctx, d_eigts3, eigts3.size());
    syncmem_complex_h2d_op()(gpu_ctx, cpu_ctx, d_sk, sk.data(), sk.size());
    syncmem_complex_h2d_op()(gpu_ctx, cpu_ctx, d_eigts1, eigts1.data(), eigts1.size());
    syncmem_complex_h2d_op()(gpu_ctx, cpu_ctx, d_eigts2, eigts2.data(), eigts2.size());
    syncmem_complex_h2d_op()(gpu_ctx, cpu_ctx, d_eigts3, eigts3.data(), eigts3.size());

    hamilt::cal_sk_op<double, psi::DEVICE_GPU>()(gpu_ctx,
                                                 ik,
                                                 ntype,
                                                 nx,
                                                 ny,
                                                 nz,
                                                 rho_nx,
                                                 rho_ny,
                                                 rho_nz,
                                                 npw,
                                                 npwx,
                                                 fftny,
                                                 eigts1_nc,
                                                 eigts2_nc,
                                                 eigts3_nc,
                                                 d_atom_na,
                                                 d_igl2isz,
                                                 d_is2fftixy,
                                                 TWO_PI,
                                                 d_kvec_c,
                                                 d_atom_tau,
                                                 d_eigts1,
                                                 d_eigts2,
                                                 d_eigts3,
                                                 d_sk);

    syncmem_complex_d2h_op()(cpu_ctx, gpu_ctx, sk.data(), d_sk, sk.size());

    for (int ii = 0; ii < sk.size(); ii++) {
        EXPECT_LT(fabs(sk[ii] - expected_sk[ii]), 6e-5);
    }

    delmem_int_op()(gpu_ctx, d_atom_na);
    delmem_int_op()(gpu_ctx, d_igl2isz);
    delmem_int_op()(gpu_ctx, d_is2fftixy);

    delmem_var_op()(gpu_ctx, d_kvec_c);
    delmem_var_op()(gpu_ctx, d_atom_tau);

    delmem_complex_op()(gpu_ctx, d_sk);
    delmem_complex_op()(gpu_ctx, d_eigts1);
    delmem_complex_op()(gpu_ctx, d_eigts2);
    delmem_complex_op()(gpu_ctx, d_eigts3);
}
#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM