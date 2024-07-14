#include <bits/stdc++.h>
#include "../kernels/cuda/sph.cuh"

#include "float.h"
#include "cuda_runtime.h"
#include "device_functions.h"
#include "device_launch_parameters.h"
#include "gtest/gtest.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "test_sph.h"
#include "module_base/array_pool.h"
using namespace std;

class gintTest : public ::testing::Test
{
  public:
};

__global__ void cuda_test(double* dr, int nwl, double* ylma_g, double* ylmcoef)
{
    double ylma[49] = {0.0};
    GintKernel::spherical_harmonics(dr, nwl, ylma, ylmcoef);
    for (int i = 0; i < 49; i++)
    {
        ylma_g[i] = ylma[i];
    }
}
__global__ void cuda_test2(double* dr, double distance, int nwl, double* dylma_g, double* ylmcoef)
{
    double ylma[49] = {0.0};
    double grly[49][3] = {0.0};
    GintKernel::spherical_harmonics_d(dr, distance, grly, nwl, ylma, ylmcoef);
    for (int i = 0; i < 49; i++)
    {
        dylma_g[i] = ylma[i];
    }
}

void get_random_double(int min, int max, double* result, int length)
{
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<double> distribution(0, 10);
    for (int i = 0; i < 3; i++)
    {
        result[i] = distribution(eng);
    }
}
void get_random_int(int min, int max, int& result)
{
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_int_distribution<int> distribution(min, max);
    result = distribution(eng);
}
// __global__ void cuda_test
TEST_F(gintTest, test)
{
    int nwl;
    double distance;

    double* dr = new double[3];
    double* dr_g;

    double ylma[49];
    double dylma[49];
    double ylma_ans[49];

    double* ylmcoef_g;
    double* ylma_g;
    double* dylma_g;
    double* ylmcoef = new double[100];

    std::vector<double> ylma_cpu(49, 0.0);
    std::vector<double> ylma_cpu_dpsir(49, 0.0);
    ModuleBase::Array_Pool<double> ylma_cpu_ddpsir(49, 3);
    
    nwl=3;
    for (int i=0;i<3;i++){
        dr[i]=i*1.0;
        distance += dr[i] * dr[i];
    }
    for (int i=0;i<100;i++)
    {
        ylmcoef[i]=i*0.1;
    }

    cudaMalloc((void**)&ylmcoef_g, 100 * sizeof(double));
    cudaMalloc((void**)&dr_g, 3 * sizeof(double));
    cudaMalloc((void**)&ylma_g, 49 * sizeof(double));
    cudaMalloc((void**)&dylma_g, 49 * 3 * sizeof(double));

    cudaMemcpy(ylmcoef_g, ylmcoef, 100 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dr_g, dr, 3 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemset(ylma_g, 0, 49 * sizeof(double));
    cudaMemset(dylma_g, 0, 49 * sizeof(double));

    cuda_test<<<1, 1>>>(dr_g, nwl, ylma_g, ylmcoef_g);
    cuda_test2<<<1, 1>>>(dr_g, distance, nwl, dylma_g, ylmcoef_g);
    sph_harm(nwl, dr[0], dr[1], dr[2], ylma_cpu, ylmcoef);
    grad_rl_sph_harm(nwl, dr[0], dr[1], dr[2], ylma_cpu_dpsir.data(), ylma_cpu_ddpsir.get_ptr_2D(), ylmcoef);
    cudaMemcpy(ylma, ylma_g, 49 * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(dylma, dylma_g, 49 * sizeof(double), cudaMemcpyDeviceToHost);
    cudaDeviceReset();

    for (int i = 0; i < 49; i++)
    {
        ylma_ans[i] = ylma_cpu[i];
        if ((abs(ylma[i])!= 0) && (ylma_ans[i]==ylma_ans[i]) && (ylma[i]==ylma[i]))
        {
            EXPECT_LT(abs(ylma_ans[i] - ylma[i]) / abs(ylma[i]), 1e-15);
        }
        ylma_ans[i] = ylma_cpu_dpsir[i];
        if ((abs(dylma[i]) != 0) &&(ylma_ans[i]==ylma_ans[i]) && (dylma[i]==dylma[i]))
        {
            EXPECT_LT(abs(ylma_ans[i] - dylma[i]) / abs(dylma[i]), 1e-15);
        }
    }
    delete[] dr;
    delete[] ylmcoef;
    
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
    int Nproc=0;
    int Rank=0;
    MPI_Comm_size(MPI_COMM_WORLD, &Nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
#endif
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    return result;
}