#ifndef VBATCH_MATRIX_MUL_H
#define VBATCH_MATRIX_MUL_H
#include <cuda.h>    // for CUDA_VERSION
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <functional>
#include <stdio.h>    // for fprintf and stderr
#include <assert.h>   // for assert


cudaError_t checkCuda(cudaError_t result);
cudaError_t checkCudaLastError();

template <typename T, int DIM_X, int DIM_Y, int BLK_M, int BLK_N, int BLK_K,
         int DIM_XA, int DIM_YA, int DIM_XB, int DIM_YB>
void vbatched_gemm_impl(int max_m, int max_n,
                 int* m, int* n, int* k,
                 T ** global_A_array, int* global_lda,
                 T ** global_B_array, int* global_ldb,
                 T ** global_C_array, int* global_ldc,
                 int batchCount, cudaStream_t stream);

typedef std::function<void(int, int,
                    int *, int *, int*,
                    double **, int *,
                    double **, int *,
                    double **, int *,
                    int, cudaStream_t)> matrix_multiple_func_type;

void gemm_algo_selector(int k, matrix_multiple_func_type &func);
#endif // VBATCH_MATRIX_MUL_H