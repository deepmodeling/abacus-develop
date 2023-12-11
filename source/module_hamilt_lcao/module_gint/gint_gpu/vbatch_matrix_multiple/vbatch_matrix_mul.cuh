#ifndef VBATCH_MATRIX_MUL_H
#define VBATCH_MATRIX_MUL_H
#include <cuda.h>    // for CUDA_VERSION
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <functional>

template<typename T, int DIM_X, int DIM_Y, int BLK_M, int BLK_N, int BLK_K,
         int DIM_XA, int DIM_YA, int DIM_XB, int DIM_YB,
         int THR_M, int THR_N>
static __device__
void vbatched_gemm_device(
    int M, int N, int K,
    T* __restrict__ A, int LDA,
    T* __restrict__ B, int LDB,
    T*       __restrict__ C, int LDC,
    T* sA, int slda,
    T* sB, int sldb);

template <typename T, int DIM_X, int DIM_Y, int BLK_M, int BLK_N, int BLK_K,
         int DIM_XA, int DIM_YA, int DIM_XB, int DIM_YB>
static __global__
void vbatched_gemm_kernel(
    int* M, int* N, int K,
    T * * Aarray, int* LDA,
    T * * Barray, int* LDB,
    T              ** Carray, int* LDC);

template <typename T, int DIM_X, int DIM_Y, int BLK_M, int BLK_N, int BLK_K,
         int DIM_XA, int DIM_YA, int DIM_XB, int DIM_YB>
void vbatched_gemm_impl(int max_m, int max_n,
                 int* m, int* n, int k,
                 T  * * global_A_array, int* global_lda,
                 T * * global_B_array, int* global_ldb,
                 T ** global_C_array, int* global_ldc,
                 int batchCount, cudaStream_t stream);
typedef std::function<void(int, int,
                    int *, int *, int,
                    double **, int *,
                    double **, int *,
                    double **, int *,
                    int, cudaStream_t)> func_type;
void gemm_algo_selector(int k, func_type &func);
#endif // VBATCH_MATRIX_MUL_H