#include "cuda_runtime.h"
#include "fftw3.h"
#include "module_base/module_device/device.h"
#include "module_base/vector3.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_basis/module_pw/module_fft/fft_bundle.h"
#include "pw_test.h"
#include <complex>
#include <vector>
#include <gtest/gtest.h>
#include <typeinfo>

using namespace std;

class PW_BASIS_K_BATCH_GPU_TEST : public ::testing::Test
{
    public:
        const int batch = 2;  // Number of batches
        const int npwk = 30;   // Number of planewaves
        const int nxyz = 1000; // Size of the 3D grid
        std::vector<int> box_index;  // Index mapping for 3D grid
        int* d_box_index=nullptr; // Device memory for box_index
        std::vector<std::complex<double>> psig; // psig(K space) data for the test,
        std::complex<double>* d_psig = nullptr; // Device memory for psig data
        std::complex<double>* d_psig_batch = nullptr; // Device memory for psig output data
        std::vector<std::complex<double>> psir; // Device memory for output psir(R space) data
        std::complex<double>* d_psir = nullptr; // Device memory for psir data
        std::complex<double>* d_psir_batch = nullptr; // Device memory for psir output batch data
        ModulePW::FFT_Bundle ft_gpu;            // FFT bundle for 3D FFT operations on GPU
        ModulePW::FFT_Bundle ft_gpu_batch;      // FFT bundle for 3D FFT operations on batch-GPU
    void SetUp() override
    {
        box_index.resize(npwk);
        psig.resize(npwk * batch);
        psir.resize(nxyz * batch);

        resize_memory_int_gpu_op()(d_box_index, npwk);
        resize_memory_complex_gpu_op()(d_psig, npwk * batch);
        resize_memory_complex_gpu_op()(d_psir, nxyz * batch);
        resize_memory_complex_gpu_op()(d_psig_batch, npwk * batch);
        resize_memory_complex_gpu_op()(d_psir_batch, nxyz * batch);

        // Initialize the box_index and input with some values
        int idx = 0;
        std::generate_n(box_index.begin(), npwk, [&idx] { return idx * idx++; });
        idx =0;
        int npwk = box_index.size();
        // Initialize psig with some complex values,it generates a complex number 
        // with real part as sqrt(idx) and imaginary part as 1/(idx+1),
        // thus in different batches the values of psig will be different.
        std::generate_n(psig.begin(), npwk * batch, [&idx,npwk] 
        {   
            idx ++;
            return std::complex<double>(std::sqrt(idx), 1.0/(idx+1));
        });
        synchronize_memory_int_h2d_op()(d_box_index, box_index.data(), npwk);
        synchronize_memory_complex_h2d_op()(d_psig, psig.data(), npwk * batch);
        synchronize_memory_complex_h2d_op()(d_psig_batch, psig.data(), npwk * batch);
        // Initialize the box_index with some values
        ft_gpu.setfft("gpu", "double");
        ft_gpu.initfft(10, 10, 10 , 1, 1, 1, 1, 1, 1);
        ft_gpu.setupFFT();
        ft_gpu_batch.setfft("gpu_batch", "double");
        ft_gpu_batch.initfft(10, 10, 10 , 1, 1, 1, 1, 1, 1);
        ft_gpu_batch.setupFFT();
    }
    void TearDown() override
    {
        box_index.clear();
        psig.clear();
        psir.clear();
        delete_memory_int_gpu_op()(d_box_index);
        delete_memory_complex_gpu_op()(d_psig);
        delete_memory_complex_gpu_op()(d_psir);
        ft_gpu.clear();
        ft_gpu_batch.clear();
    }
};

TEST_F(PW_BASIS_K_BATCH_GPU_TEST,convulution)
{
    // STEP 1 set the 3D FFT box operation for CPU
    for (int i = 0; i < npwk; ++i)
    {
        EXPECT_EQ(box_index[i], i * i);
    }

    // STEP 2 check the input psig has been
    // correctly mapped to the 3D grid
    std::vector<std::complex<double>> compute_psir(nxyz * batch);
    std::vector<std::complex<double>> compute_psir_batch(nxyz * batch);
    for (int i = 0; i< batch; i++)
    {
        ModulePW::set_3d_fft_box_op<double, 
            base_device::DEVICE_GPU>()
        (
            npwk,
            d_box_index,
            d_psig + i * npwk,
            d_psir + i * nxyz
        );
        synchronize_memory_complex_d2h_op()(compute_psir.data()+i * nxyz, d_psir + i *nxyz, nxyz);
    }
    ModulePW::set_3d_fft_box_op<double, 
        base_device::DEVICE_GPU>()
    (
        npwk,
        nxyz,
        d_box_index,
        d_psig_batch,
        d_psir_batch,
        batch
    );
    
    synchronize_memory_complex_d2h_op()(compute_psir_batch.data(), d_psir_batch,nxyz * batch);
    for (int i = 0; i < nxyz*batch ; ++i)
    {
        EXPECT_NEAR(compute_psir[i].real(), compute_psir_batch[i].real(), 1e-7);
        EXPECT_NEAR(compute_psir[i].imag(), compute_psir_batch[i].imag(), 1e-7);
    }

    // STEP 3 perform the 3D FFT forward operation

    for (int i=0;i<batch;i++)
    {
        ft_gpu.fft3D_forward(d_psir + i *nxyz, d_psir + i *nxyz);
        
    }
    ft_gpu_batch.fft3D_forward(d_psir_batch, d_psir_batch );
    synchronize_memory_complex_d2h_op()(compute_psir.data(),d_psir , nxyz * batch);
    synchronize_memory_complex_d2h_op()(compute_psir_batch.data(), d_psir_batch,nxyz * batch);

    for (int i = 0; i < nxyz *batch ; ++i)
    {
        EXPECT_NEAR(compute_psir[i].real(), compute_psir_batch[i].real(), 1e-4);
        EXPECT_NEAR(compute_psir[i].imag(), compute_psir_batch[i].imag(), 1e-4);
    }
    
}