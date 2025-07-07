#include "gemm_tn_vbatch.cuh"
#include "gemm_nn_vbatch.cuh"
#include "dgemm_vbatch.h"

void dgemm_nn_vbatch(
    int max_m, int max_n, int max_k,
    const int* m_d, const int* n_d, const int* k_d,
    const double* const* A_array_d, const int* lda_d,
    const double* const* B_array_d, const int* ldb_d,
    double** C_array_d, const int* ldc_d,
    int batchCount, cudaStream_t stream,
    const double* alpha)
{
    vbatched_gemm_nn_impl<double, 8, 4, 16, 16, 8, 8, 4, 8, 4>
    (max_m, max_n, m_d, n_d, k_d,
    A_array_d, lda_d,
    B_array_d, ldb_d,
    C_array_d, ldc_d,
    batchCount, stream, alpha);

}

void dgemm_tn_vbatch(
    int max_m, int max_n, int max_k,
    const int* m_d, const int* n_d, const int* k_d,
    const double* const* A_array_d, const int* lda_d,
    const double* const* B_array_d, const int* ldb_d,
    double** C_array_d, const int* ldc_d,
    int batchCount, cudaStream_t stream,
    const double* alpha)
{
    vbatched_gemm_tn_impl<double, 8,4,16,16,4,8,4,8,4>
        (max_m, max_n, m_d, n_d, k_d,
        A_array_d, lda_d,
        B_array_d, ldb_d,
        C_array_d, ldc_d,
        batchCount, stream, alpha);
}
