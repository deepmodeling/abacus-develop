#ifndef VBATCH_MATRIX_MUL_H
#define VBATCH_MATRIX_MUL_H
#include <cuda.h>    // for CUDA_VERSION
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <functional>
#include <stdio.h>    // for fprintf and stderr
#include <assert.h>   // for assert


/**
 * Performs a batched matrix multiplication using the vbatched_gemm_impl function.
 *
 * C = alpha * A * B + C
 * @tparam T The data type of the matrices.
 * @tparam DIM_X The number of threads in the x-dimension of each block.
 * @tparam DIM_Y The number of threads in the y-dimension of each block.
 * @tparam BLK_M The number of rows processed by each thread block.
 * @tparam BLK_N The number of columns processed by each thread block.
 * @tparam BLK_K The number of elements processed by each thread block along the K dimension.
 * @tparam DIM_XA The number of threads in the x-dimension used for loading matrix A.
 * @tparam DIM_YA The number of threads in the y-dimension used for loading matrix A.
 * @tparam DIM_XB The number of threads in the x-dimension used for loading matrix B.
 * @tparam DIM_YB The number of threads in the y-dimension used for loading matrix B.
 * @param max_m The maximum number of rows in the matrices.
 * @param max_n The maximum number of columns in the matrices.
 * @param m An array of batch sizes for the number of rows in each matrix.
 * @param n An array of batch sizes for the number of columns in each matrix.
 * @param k An array of batch sizes for the number of elements in each matrix along the K dimension.
 * @param global_A_array An array of pointers to the input matrices A.
 * @param global_lda An array of leading dimensions for the input matrices A.
 * @param global_B_array An array of pointers to the input matrices B.
 * @param global_ldb An array of leading dimensions for the input matrices B.
 * @param global_C_array An array of pointers to the output matrices C.
 * @param global_ldc An array of leading dimensions for the output matrices C.
 * @param batchCount The number of matrices in the batch.
 * @param stream The CUDA stream to use for the computation.
 * @param alpha The scalar value to multiply the matrices by (optional, default is nullptr).
 * generate by copilot
 */
template <typename T, int DIM_X, int DIM_Y, int BLK_M, int BLK_N, int BLK_K,
         int DIM_XA, int DIM_YA, int DIM_XB, int DIM_YB>
void vbatched_gemm_impl(int max_m, int max_n,
                 int* m, int* n, int* k,
                 T ** global_A_array, int* global_lda,
                 T ** global_B_array, int* global_ldb,
                 T ** global_C_array, int* global_ldc,
                 int batchCount, cudaStream_t stream, T* alpha = nullptr);

typedef std::function<void(int, int,
                    int *, int *, int*,
                    double **, int *,
                    double **, int *,
                    double **, int *,
                    int, cudaStream_t,  double* alpha)> matrix_multiple_func_type;

void gemm_algo_selector(int k, matrix_multiple_func_type &func);
#endif // VBATCH_MATRIX_MUL_H