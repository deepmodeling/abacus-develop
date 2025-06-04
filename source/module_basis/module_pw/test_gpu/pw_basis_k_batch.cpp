#include "cuda_runtime.h"
#include "fftw3.h"
#include "module_base/module_device/device.h"
#include "module_base/vector3.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "pw_test.h"
#include <complex>
#include <vector>
#include <gtest/gtest.h>
#include <typeinfo>

using namespace std;

class PW_BASIS_K_BATCH_GPU_TEST : public ::testing::Test
{
    public:
        const int batch = 10; // Number of batches
        const int npwk = 30;   // Number of planewaves
        const int nxyz = 1000; // Size of the 3D grid
        std::vector<int> box_index;  // Index mapping for 3D grid
        int* d_box_index=nullptr; // Device memory for box_index
        std::vector<std::complex<double>> rhog;   // Input data for the test,
        std::complex<double>* d_rhog = nullptr; // Device memory for rhoG data
        std::vector<std::complex<double>> rhor = nullptr; // Device memory for output rhoG data
        std::complex<double>* d_rhor = nullptr; // Device memory for output data
    void SetUp() override
    {
        box_index.resize(npwk);
        rhog.resize(npwk);
        resize_memory_int_gpu_op()(d_box_index, npwk);
        resize_memory_complex_gpu_op()(d_rhog, npwk);
        // Initialize the box_index and input with some values
        int idx = 0;
        std::generate_n(box_index.begin(), npwk, [&idx] { return idx * idx++; });
        idx =0;
        std::generate_n(rhog.begin(), npwk, [&idx] { return std::complex<double>(sqrt(idx), 1/(idx+1)); });
        synchronize_memory_int_h2d_op()(d_box_index, box_index.data(), npwk);
        synchronize_memory_complex_h2d_op()(d_rhog, rhog.data(), npwk);
        // Initialize the box_index with some values
        
        // resize_memory_int_gpu_op
    }
    void TearDown() override
    {
        box_index.clear();
        rhog.clear();
        delete_memory_int_gpu_op()(d_box_index);
        delete_memory_complex_gpu_op()(d_rhog);
    }
};

TEST_F(PW_BASIS_K_BATCH_GPU_TEST,convulution)
{
    for (int i = 0; i < npwk; ++i)
    {
        EXPECT_EQ(box_index[i], i * i);
    }
}