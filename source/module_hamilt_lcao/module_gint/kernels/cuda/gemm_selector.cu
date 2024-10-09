#include <iostream>

#include "gemm_selector.cuh"
#include "vbatch_matrix_mul.cuh"
#include "cuda_tools.cuh"
#include "module_base/blas_connector.h"
#include "code_gen.cuh"

/*
 * Here we have utilized a very straightforward and brute-force method to select
 * the optimal matrix multiplication kernel for a given scale of computation: we
 * compute with all scales of kernels under the current computational task to
 * find the fastest parameter combination. This approach can lead to an increase
 * in compilation time.
 */
void gemm_algo_selector(int matrix_k, matrix_multiple_func_type& fastest_algo,const UnitCell& ucell)
{
    int batchCount_per_type = 32;
    int batchCount
        = batchCount_per_type * ucell.ntype * ucell.ntype;

    Cuda_Mem_Wrapper<int> m(batchCount);
    Cuda_Mem_Wrapper<int> n(batchCount);
    Cuda_Mem_Wrapper<int> k(batchCount);

    int max_m = ucell.nwmax, max_n = ucell.nwmax;

    Cuda_Mem_Wrapper<double> A(batchCount * max_m * matrix_k);
    Cuda_Mem_Wrapper<double> B(batchCount * max_n * matrix_k);
    Cuda_Mem_Wrapper<double> C(batchCount * max_m * max_n);

    Cuda_Mem_Wrapper<int> lda(batchCount);
    Cuda_Mem_Wrapper<int> ldb(batchCount);
    Cuda_Mem_Wrapper<int> ldc(batchCount);

    Cuda_Mem_Wrapper<double*> A_array(batchCount);
    Cuda_Mem_Wrapper<double*> B_array(batchCount);
    Cuda_Mem_Wrapper<double*> C_array(batchCount);

    for (int i = 0; i < batchCount * max_m * matrix_k; ++i)
    {
        A.get_host_pointer()[i] = i * 0.001;
    }
    for (int i = 0; i < batchCount * max_n * matrix_k; ++i)
    {
        B.get_host_pointer()[i] = i * 0.002;
    }

    double* cpu_result = new double[batchCount * max_m * max_n];
    memset(cpu_result, 0, batchCount * max_m * max_n * sizeof(double));
    int index = 0;
    for (int i = 0; i < batchCount_per_type; ++i)
    {
        for (int j = 0; j < ucell.ntype; j++)
        {
            for (int l = 0; l < ucell.ntype; l++)
            {
                m.get_host_pointer()[index] = ucell.atoms[j].nw;
                n.get_host_pointer()[index] = ucell.atoms[l].nw;
                k.get_host_pointer()[index] = matrix_k;

                lda.get_host_pointer()[index] = matrix_k;
                ldb.get_host_pointer()[index] = matrix_k;
                ldc.get_host_pointer()[index] = ucell.atoms[l].nw;

                A_array.get_host_pointer()[index]
                    = &A.get_device_pointer()[index * max_m * matrix_k];
                B_array.get_host_pointer()[index]
                    = &B.get_device_pointer()[index * max_n * matrix_k];
                C_array.get_host_pointer()[index]
                    = &C.get_device_pointer()[index * max_n
                                              * max_m]; // test atom add
                BlasConnector::gemm(
                    'N',
                    'T',
                    m.get_host_pointer()[index],
                    n.get_host_pointer()[index],
                    matrix_k,
                    1.0,
                    &A.get_host_pointer()[index * max_m * matrix_k],
                    matrix_k,
                    &B.get_host_pointer()[index * max_n * matrix_k],
                    matrix_k,
                    1.0,
                    &cpu_result[index * max_m * max_n],
                    n.get_host_pointer()[index]);
                index++;
            }
        }
    }

    m.copy_host_to_device_sync();
    n.copy_host_to_device_sync();
    k.copy_host_to_device_sync();

    lda.copy_host_to_device_sync();
    ldb.copy_host_to_device_sync();
    ldc.copy_host_to_device_sync();

    A.copy_host_to_device_sync();
    B.copy_host_to_device_sync();
    A_array.copy_host_to_device_sync();
    B_array.copy_host_to_device_sync();
    C_array.copy_host_to_device_sync();

    cudaStream_t temp_stream;
    checkCuda(cudaStreamCreate(&temp_stream));

    float fastest_time = 1000000;
    fastest_algo = vbatched_gemm_impl<double, 16, 4, 32, 16, 16, 16, 4, 16, 4>;

    int* d_m = m.get_device_pointer();
    int* d_n = n.get_device_pointer();
    int* d_k = k.get_device_pointer();

    double** d_global_A_array = A_array.get_device_pointer();
    double** d_global_B_array = B_array.get_device_pointer();
    double** d_global_C_array = C_array.get_device_pointer();

    double* h_global_C = C.get_host_pointer();
    double* d_global_C = C.get_device_pointer();

    int* d_global_lda = lda.get_device_pointer();
    int* d_global_ldb = ldb.get_device_pointer();
    int* d_global_ldc = ldc.get_device_pointer();

/*
 * Please do not manually modify the code in the following file;
 * it should simply be generated through a loop using a short Python program.
 */
#include "code_gen.cpp"
    checkCuda(cudaStreamDestroy(temp_stream));
    std::cout << " gemm_algo_selector::Fastest time: " << fastest_time << " ms"
              << std::endl;
    // fastest_algo = vbatched_gemm_impl<double, 16, 4, 32, 16, 16, 16, 4, 16,
    // 4>;
    delete[] cpu_result;
}