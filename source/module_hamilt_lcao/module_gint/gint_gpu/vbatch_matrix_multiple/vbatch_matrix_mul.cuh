#ifndef VBATCH_MATRIX_MUL_H
#define VBATCH_MATRIX_MUL_H

#define sA(i,j)    sA[(j)*slda + (i)]
#define sB(i,j)    sB[(j)*sldb + (i)]
#define sC(i,j)    sC[(j)*sldc + (i)]

template<typename T, const int DIM_X, const int DIM_Y, const int BLK_M, const int BLK_N, const int BLK_K,
         const int DIM_XA, const int DIM_YA, const int DIM_XB, const int DIM_YB,
         const int THR_M, const int THR_N>
static __device__
void vbatched_gemm_device(
    int M, int N, int K,
    const T* __restrict__ A, int LDA,
    const T* __restrict__ B, int LDB,
    T*       __restrict__ C, int LDC,
    T* sA, int slda,
    T* sB, int sldb);

template <typename T, const int DIM_X, const int DIM_Y, const int BLK_M, const int BLK_N, const int BLK_K,
         const int DIM_XA, const int DIM_YA, const int DIM_XB, const int DIM_YB>
static __global__
void vbatched_gemm_kernel(
    int* M, int* N, int K,
    T const * const * Aarray, int* LDA,
    T const * const * Barray, int* LDB,
    T              ** Carray, int* LDC);

template <typename T, const int DIM_X, const int DIM_Y, const int BLK_M, const int BLK_N, const int BLK_K, const int dim_vec,
         const int DIM_XA, const int DIM_YA, const int DIM_XB, const int DIM_YB>
void vbatched_gemm(
    int max_m, int max_n,
    int* m, int* n, int k,
    T const * const * dA_array, int* ldda,
    T const * const * dB_array, int* lddb,
    T ** dC_array, int* lddc,
    int batchCount, cudaStream_t stream);

#endif // VBATCH_MATRIX_MUL_H