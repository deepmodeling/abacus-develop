#ifndef VBATCH_MATRIX_MUL_CUH
#define VBATCH_MATRIX_MUL_CUH
#include <assert.h> // for assert
#include <cublas_v2.h>
#include <cuda.h> // for CUDA_VERSION
#include <cuda_runtime.h>
#include <stdio.h> // for fprintf and stderr

#include "cuda_tools.cuh"
#include <functional>
#include "../../grid_technique.h"
#include "module_cell/unitcell.h"
#include "module_base/module_device/device.h"


#define sA(i, j) sA[(j)*slda + (i)]
#define sB(i, j) sB[(j)*sldb + (i)]
#define fetch(A, m, n, bound) offs_d##A[min(n * LD##A + m, bound)]

template <typename T,
          int DIM_X,
          int DIM_Y,
          int BLK_M,
          int BLK_N,
          int BLK_K,
          int DIM_XA,
          int DIM_YA,
          int DIM_XB,
          int DIM_YB,
          int THR_M,
          int THR_N>
static __device__ void vbatched_gemm_device(int M,
                                            int N,
                                            int K,
                                            T* __restrict__ A,
                                            int LDA,
                                            T* __restrict__ B,
                                            int LDB,
                                            T* __restrict__ C,
                                            int LDC,
                                            T* sA,
                                            int slda,
                                            T* sB,
                                            int sldb,
                                            T alpha)
{
    int idx = threadIdx.x; // thread's m dimension
    int idy = threadIdx.y; // thread's n dimension

    int idt = DIM_X * idy + idx; // thread's global number

    int idxA = idt % DIM_XA; // idx within A
    int idyA = idt / DIM_XA; // idy within A

    int idxB = idt % DIM_XB; // idx within B
    int idyB = idt / DIM_XB; // idy within B

    int blx = blockIdx.x; // block's m dimension
    int bly = blockIdx.y; // block's n dimension

    // Registers for the innermost loop
    T rC[THR_N][THR_M];
    T rA[THR_M];
    T rB[THR_N];

    // Registers for the dev->shmem copy
    T ra[BLK_M / DIM_YA][BLK_K / DIM_XA];
    T rb[BLK_N / DIM_YB][BLK_K / DIM_XB];

    // bound is the correction to offs_d in order to not get out of memory bound
    // so bound could be negative value since offs_d could be out of bound
    T* offs_dA = A + blx * BLK_M * LDA + idyA * LDA + idxA;
    int boundA
        = (LDA * (M - 1) + K) - (blx * BLK_M * LDA + idyA * LDA + idxA) - 1;

    T* offs_dB = B + bly * BLK_N * LDB + idyB * LDB + idxB;
    int boundB
        = (LDB * (N - 1) + K) - (bly * BLK_N * LDB + idyB * LDB + idxB) - 1;

    int m, n, k, kk;

// Zero C
#pragma unroll
    for (n = 0; n < THR_N; n++)
    {
#pragma unroll
        for (m = 0; m < THR_M; m++)
        {
            rC[n][m] = 0.0;
        }
    }

// Load A dev->shmem
#pragma unroll
    for (n = 0; n < BLK_M; n += DIM_YA)
    {
#pragma unroll
        for (m = 0; m < BLK_K; m += DIM_XA)
        {
            sA(n + idyA, m + idxA) = fetch(A, m, n, boundA);
        }
    }

#pragma unroll
    for (n = 0; n < BLK_N; n += DIM_YB)
    {
#pragma unroll
        for (m = 0; m < BLK_K; m += DIM_XB)
        {
            sB(m + idxB, n + idyB) = fetch(B, m, n, boundB);
        }
    }

    __syncthreads();

    for (kk = 0; kk < K - BLK_K; kk += BLK_K)
    {
        offs_dA += BLK_K;
        boundA -= BLK_K;

        offs_dB += BLK_K;
        boundB -= BLK_K;

// Load A dev->regs
#pragma unroll
        for (n = 0; n < BLK_M / DIM_YA; n++)
        {
#pragma unroll
            for (m = 0; m < BLK_K / DIM_XA; m++)
            {
                ra[n][m] = fetch(A, m * DIM_XA, n * DIM_YA, boundA);
            }
        }

// Load B dev->regs
#pragma unroll
        for (n = 0; n < BLK_N / DIM_YB; n++)
        {
#pragma unroll
            for (m = 0; m < BLK_K / DIM_XB; m++)
            {
                rb[n][m] = fetch(B, m * DIM_XB, n * DIM_YB, boundB);
            }
        }

// Multiply
#pragma unroll
        for (k = 0; k < BLK_K; k++)
        {
// Load A shmem->regs
#pragma unroll
            for (m = 0; m < THR_M; m++)
            {
                rA[m] = sA(m * DIM_X + idx, k);
            }

// Load B shmem->regs
#pragma unroll
            for (n = 0; n < THR_N; n++)
            {
                rB[n] = sB(k, n * DIM_Y + idy);
            }

// Compute
#pragma unroll
            for (n = 0; n < THR_N; n++)
            {
#pragma unroll
                for (m = 0; m < THR_M; m++)
                {
                    rC[n][m] += rA[m] * rB[n];
                }
            }
        }

        __syncthreads();

// Load A regs->shmem
#pragma unroll
        for (n = 0; n < BLK_M / DIM_YA; n++)
        {
#pragma unroll
            for (m = 0; m < BLK_K / DIM_XA; m++)
            {
                sA(n * DIM_YA + idyA, m * DIM_XA + idxA) = ra[n][m];
            }
        }

// Load B regs->shmem
#pragma unroll
        for (n = 0; n < BLK_N / DIM_YB; n++)
        {
#pragma unroll
            for (m = 0; m < BLK_K / DIM_XB; m++)
            {
                sB(m * DIM_XB + idxB, n * DIM_YB + idyB) = rb[n][m];
            }
        }
        __syncthreads();
    }

    // Multiply last full (BLK_K) or partial block of
    // columns of op(A) and rows of op(B).
    // It's okay that m,n exceed matrix bounds as all work is in registers
    // or shared memory, and out-of-bounds rC[n][m] will not be saved later.
    kk = K - kk;
#pragma unroll
    for (k = 0; k < kk; k++)
    {
// Load A shmem->regs
#pragma unroll
        for (m = 0; m < THR_M; m++)
        {
            rA[m] = sA(m * DIM_X + idx, k);
        }

// Load B shmem->regs
#pragma unroll
        for (n = 0; n < THR_N; n++)
        {
            rB[n] = sB(k, n * DIM_Y + idy);
        }

// Compute
#pragma unroll
        for (n = 0; n < THR_N; n++)
        {
#pragma unroll
            for (m = 0; m < THR_M; m++)
            {
                rC[n][m] += rA[m] * rB[n];
            }
        }
    }

// Store C regs->dev
#pragma unroll
    for (n = 0; n < THR_N; n++)
    {
        int coord_dCn = bly * BLK_N + n * DIM_Y + idy;
#pragma unroll
        for (m = 0; m < THR_M; m++)
        {
            int coord_dCm = blx * BLK_M + m * DIM_X + idx;
            if (coord_dCm < M && coord_dCn < N)
            {
                int offsC = coord_dCn * LDC + coord_dCm;

                atomicAdd(C + offsC, rC[n][m] * alpha);
            }
        }
    }
}

/******************************************************************************/
template <typename T,
          int DIM_X,
          int DIM_Y,
          int BLK_M,
          int BLK_N,
          int BLK_K,
          int DIM_XA,
          int DIM_YA,
          int DIM_XB,
          int DIM_YB>
static __global__ void vbatched_gemm_kernel(int* M,
                                            int* N,
                                            int* K,
                                            T** global_A_array,
                                            int* global_lda,
                                            T** global_B_array,
                                            int* global_ldb,
                                            T** global_C_array,
                                            int* global_ldc,
                                            T* alpha)
{
    extern __shared__ __align__(sizeof(T)) unsigned char smem[];
    T* shared_mem = reinterpret_cast<T*>(smem);

    int batchid = blockIdx.z;
    int local_M = (int)M[batchid];
    int local_N = (int)N[batchid];
    int local_K = (int)K[batchid];

    if (blockIdx.x >= (local_M + BLK_M - 1) / BLK_M)
        return;
    if (blockIdx.y >= (local_N + BLK_N - 1) / BLK_N)
        return;

    int shared_lda = BLK_M + 1;
    int shared_ldb = BLK_K + 1;
    T* shared_A = (T*)shared_mem;
    T* shared_B = shared_A + shared_lda * BLK_K;
    double alpha_tmp = 1.0;
    if (alpha != nullptr)
    {
        alpha_tmp = alpha[batchid];
    }
    vbatched_gemm_device<T,
                         DIM_X,
                         DIM_Y,
                         BLK_M,
                         BLK_N,
                         BLK_K,
                         DIM_XA,
                         DIM_YA,
                         DIM_XB,
                         DIM_YB,
                         (BLK_M / DIM_X),
                         (BLK_N / DIM_Y)>(local_M,
                                          local_N,
                                          local_K,
                                          global_A_array[batchid],
                                          (int)global_lda[batchid],
                                          global_B_array[batchid],
                                          (int)global_ldb[batchid],
                                          global_C_array[batchid],
                                          (int)global_ldc[batchid],
                                          shared_A,
                                          shared_lda,
                                          shared_B,
                                          shared_ldb,
                                          alpha_tmp);
}

static inline int ceildiv(int x, int y)
{
    return (x + y - 1) / y;
}

/**
 * Performs a batched matrix multiplication using the vbatched_gemm_impl
 * function.
 *
 * C = alpha * A * B + C
 * @tparam T The data type of the matrices.
 * @tparam DIM_X The number of threads in the x-dimension of each block.
 * @tparam DIM_Y The number of threads in the y-dimension of each block.
 * @tparam BLK_M The number of rows processed by each thread block.
 * @tparam BLK_N The number of columns processed by each thread block.
 * @tparam BLK_K The number of elements processed by each thread block along the
 * K dimension.
 * @tparam DIM_XA The number of threads in the x-dimension used for loading
 * matrix A.
 * @tparam DIM_YA The number of threads in the y-dimension used for loading
 * matrix A.
 * @tparam DIM_XB The number of threads in the x-dimension used for loading
 * matrix B.
 * @tparam DIM_YB The number of threads in the y-dimension used for loading
 * matrix B.
 * @param max_m The maximum number of rows in the matrices.
 * @param max_n The maximum number of columns in the matrices.
 * @param m An array of batch sizes for the number of rows in each matrix.
 * @param n An array of batch sizes for the number of columns in each matrix.
 * @param k An array of batch sizes for the number of elements in each matrix
 * along the K dimension.
 * @param global_A_array An array of pointers to the input matrices A.
 * @param global_lda An array of leading dimensions for the input matrices A.
 * @param global_B_array An array of pointers to the input matrices B.
 * @param global_ldb An array of leading dimensions for the input matrices B.
 * @param global_C_array An array of pointers to the output matrices C.
 * @param global_ldc An array of leading dimensions for the output matrices C.
 * @param batchCount The number of matrices in the batch.
 * @param stream The CUDA stream to use for the computation.
 * @param alpha The scalar value to multiply the matrices by (optional, default
 * is nullptr). generate by copilot
 */

/*
 * Why do we need to implement our own matrix multiplication based on the magma
 * code? There are two main reasons. First is when we are doing batch matrix
 * multiplication, since we need to accumulate the results of the
 * multiplications, it is necessary to pass the same memory address of matrix C
 * to different multiplications. This way, the accumulation can be done directly
 * through atomic operations during the matrix multiplication, avoiding the
 * reduction operations after the multiplication. Secondly, when calculating the
 * charge density, where C = alpha * A * B + C, the value of alpha might be
 * different for the same batch of matrices. Using the standard matrix
 * multiplication interface would require breaking down the batch matrix
 * multiplication into smaller batches. In practice, it is difficult to
 * accumulate a batch.
 *
 * Moreover, taking into account the specific requirements of our application,
 * especially the fact that we can relatively easily control the arrangement of
 * the matrix elements, we have only implemented one type of requirement for
 * matrix transposition. That is, we have implemented the operation C = alpha *
 * trans(A) * B + C under the constraint of column-major order.
 *
 * Finally, we would like to thank Magma for its contributions to the field of
 * scientific computing.
 */

template <typename T,
          int DIM_X,
          int DIM_Y,
          int BLK_M,
          int BLK_N,
          int BLK_K,
          int DIM_XA,
          int DIM_YA,
          int DIM_XB,
          int DIM_YB>
void vbatched_gemm_impl(int max_m,
                               int max_n,
                               int* m,
                               int* n,
                               int* k,
                               T** global_A_array,
                               int* global_lda,
                               T** global_B_array,
                               int* global_ldb,
                               T** global_C_array,
                               int* global_ldc,
                               int batchCount,
                               cudaStream_t stream,
                               T* alpha = nullptr)
{
    // The positions of A and B have been swapped here.
    // This is because the original code is for column-major matrices.
    // We use row-major matrices, so we need to swap A and B.
    // The vbatched_gemm_impl is for C = trans(A) * B + C, but we need trans(C).
    // Which means: trans(C) = trans(trans(A)*B + C) = trans(B) * A + trans(C)
    // Then, ldc should be N, lda and ldb should be K

    size_t shared_mem_size = 0;
    shared_mem_size += (BLK_M + 1) * BLK_K * sizeof(T);
    shared_mem_size += (BLK_K + 1) * BLK_N * sizeof(T);
    dim3 dimBlock(DIM_X, DIM_Y);
    const int max_batch_count = 32768;
    const int loop_num = batchCount / max_batch_count;
    const int remain_num = batchCount % max_batch_count;

    for (int i = 0; i < loop_num; ++i)
    {
        dim3 dimGrid(ceildiv(max_n, BLK_M),
                     ceildiv(max_m, BLK_N),
                     max_batch_count);
        T* alpha_tmp = nullptr;
        if (alpha != nullptr)
        {
            alpha_tmp = alpha + i * max_batch_count;
        }

        vbatched_gemm_kernel<T,
                             DIM_X,
                             DIM_Y,
                             BLK_M,
                             BLK_N,
                             BLK_K,
                             DIM_XA,
                             DIM_YA,
                             DIM_XB,
                             DIM_YB>
            <<<dimGrid, dimBlock, shared_mem_size, stream>>>(
                n + i * max_batch_count,
                m + i * max_batch_count,
                k + i * max_batch_count,
                global_B_array + i * max_batch_count,
                global_ldb + i * max_batch_count,
                global_A_array + i * max_batch_count,
                global_lda + i * max_batch_count,
                global_C_array + i * max_batch_count,
                global_ldc + i * max_batch_count,
                alpha_tmp);
        checkCudaLastError();
    }
    if (remain_num > 0)
    {
        dim3 dimGrid(ceildiv(max_n, BLK_M), ceildiv(max_m, BLK_N), remain_num);
        T* alpha_tmp = nullptr;
        if (alpha != nullptr)
        {
            alpha_tmp = alpha + loop_num * max_batch_count;
        }
        vbatched_gemm_kernel<T,
                             DIM_X,
                             DIM_Y,
                             BLK_M,
                             BLK_N,
                             BLK_K,
                             DIM_XA,
                             DIM_YA,
                             DIM_XB,
                             DIM_YB>
            <<<dimGrid, dimBlock, shared_mem_size, stream>>>(
                n + loop_num * max_batch_count,
                m + loop_num * max_batch_count,
                k + loop_num * max_batch_count,
                global_B_array + loop_num * max_batch_count,
                global_ldb + loop_num * max_batch_count,
                global_A_array + loop_num * max_batch_count,
                global_lda + loop_num * max_batch_count,
                global_C_array + loop_num * max_batch_count,
                global_ldc + loop_num * max_batch_count,
                alpha_tmp);
        checkCudaLastError();
    }
}

template <typename T,
          int DIM_X,
          int DIM_Y,
          int BLK_M,
          int BLK_N,
          int BLK_K,
          int DIM_XA,
          int DIM_YA,
          int DIM_XB,
          int DIM_YB>
void gemm_time_measure(int max_m,
                       int max_n,
                       int* m,
                       int* n,
                       int* k,
                       T** global_A_array,
                       int* global_lda,
                       T** global_B_array,
                       int* global_ldb,
                       T** global_C_array,
                       int* global_ldc,
                       int batchCount,
                       cudaStream_t stream,
                       float& fast_time,
                       matrix_multiple_func_type& fastest_algo,
                       double* cpu_result,
                       double* h_global_C,
                       double* d_global_C)
{
    cudaEvent_t start, stop;
    checkCuda(
        cudaMemset(d_global_C, 0, batchCount * max_m * max_n * sizeof(double)));
    checkCuda(cudaEventCreate(&start));
    checkCuda(cudaEventCreate(&stop));
    checkCuda(cudaEventRecord(start, stream));
    vbatched_gemm_impl<T,
                       DIM_X,
                       DIM_Y,
                       BLK_M,
                       BLK_N,
                       BLK_K,
                       DIM_XA,
                       DIM_YA,
                       DIM_XB,
                       DIM_YB>(max_m,
                               max_n,
                               m,
                               n,
                               k,
                               global_A_array,
                               global_lda,
                               global_B_array,
                               global_ldb,
                               global_C_array,
                               global_ldc,
                               batchCount,
                               stream);
    checkCuda(cudaEventRecord(stop, stream));
    cudaError_t cuda_status = cudaGetLastError();
    checkCuda(cudaStreamSynchronize(stream));
    float milliseconds = 0;
    checkCuda(cudaEventElapsedTime(&milliseconds, start, stop));

    // WARNING !!!!! Here we assume that all m and n are the same
    checkCuda(cudaMemcpy(h_global_C,
                         d_global_C,
                         batchCount * max_m * max_n * sizeof(double),
                         cudaMemcpyDeviceToHost));
    bool check_result = true;
    for (int i = 0; i < batchCount * max_m * max_n; ++i)
    {
        if (abs(cpu_result[i] - h_global_C[i]) > 0.001)
        {
            check_result = false;
            break;
        }
    }
    if (milliseconds < fast_time && cuda_status == cudaSuccess && check_result)
    {
        fast_time = milliseconds;
        fastest_algo = vbatched_gemm_impl<T,
                                          DIM_X,
                                          DIM_Y,
                                          BLK_M,
                                          BLK_N,
                                          BLK_K,
                                          DIM_XA,
                                          DIM_YA,
                                          DIM_XB,
                                          DIM_YB>;
#ifdef __DEBUG
        std::cout << "found! fastest time: " << fast_time << std::endl;
        std::cout << DIM_X << "," << DIM_Y << "," << BLK_M << "," << BLK_N
                  << "," << BLK_K << "," << DIM_XA << "," << DIM_YA << ","
                  << DIM_XB << "," << DIM_YB << std::endl;
#endif
    }
}
#endif // VBATCH_MATRIX_MUL_CUH