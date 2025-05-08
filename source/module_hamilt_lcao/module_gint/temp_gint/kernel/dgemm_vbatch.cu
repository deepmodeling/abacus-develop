#include "gemm_tn_vbatch.cuh"
#include "gemm_nn_vbatch.cuh"
#include "dgemm_vbatch.h"

// The template parameter settings for the function are based on the MAGMA source code settings.
// Specifically, they refer to the settings for the "nn" shape in dgemm_vbatched_core.
void dgemm_nn_vbatch(
    int max_m, int max_n, int max_k,
    const int* m_d, const int* n_d, const int* k_d,
    const double* const* A_array_d, const int* lda_d,
    const double* const* B_array_d, const int* ldb_d,
    double** C_array_d, const int* ldc_d,
    int batchCount, cudaStream_t stream,
    const double* alpha)
{
    if (max_k < 32)
    {
        if(max_k == 8 && max_m ==24)
        {
            vbatched_gemm_nn_impl<double, 8, 8, 16, 24, 8, 8, 8, 8, 8>
            (max_m, max_n, m_d, n_d, k_d,
            A_array_d, lda_d,
            B_array_d, ldb_d,
            C_array_d, ldc_d,
            batchCount, stream, alpha);
        }
        else if (max_m < 32)
        {
            vbatched_gemm_nn_impl<double, 8, 8, 32, 16, 8, 8, 8, 8, 8>
            (max_m, max_n, m_d, n_d, k_d,
            A_array_d, lda_d,
            B_array_d, ldb_d,
            C_array_d, ldc_d,
            batchCount, stream, alpha);
        }
        else
        {
            vbatched_gemm_nn_impl<double, 16, 16, 48, 32, 16, 16, 16, 16, 16>
            (max_m, max_n, m_d, n_d, k_d,
            A_array_d, lda_d,
            B_array_d, ldb_d,
            C_array_d, ldc_d,
            batchCount, stream, alpha);
        }
    }
    else
    {
        if (max_n < 80)
        {
            vbatched_gemm_nn_impl<double, 16, 8, 32, 24, 16, 16, 8, 16, 8>
            (max_m, max_n, m_d, n_d, k_d,
            A_array_d, lda_d,
            B_array_d, ldb_d,
            C_array_d, ldc_d,
            batchCount, stream, alpha);
        }
        else
        {
            vbatched_gemm_nn_impl<double, 16, 16, 48, 32, 16, 16, 16, 16, 16>
            (max_m, max_n, m_d, n_d, k_d,
            A_array_d, lda_d,
            B_array_d, ldb_d,
            C_array_d, ldc_d,
            batchCount, stream, alpha);
        }
    }
}

// the template parameters refer to the settings for the "nt" shape in dgemm_vbatched_core.
void dgemm_tn_vbatch(
    int max_m, int max_n, int max_k,
    const int* m_d, const int* n_d, const int* k_d,
    const double* const* A_array_d, const int* lda_d,
    const double* const* B_array_d, const int* ldb_d,
    double** C_array_d, const int* ldc_d,
    int batchCount, cudaStream_t stream,
    const double* alpha)
{
    if (max_k < 128)
    {
        vbatched_gemm_tn_impl<double, 16, 8, 32, 32, 8, 16, 8, 16, 8>
        (max_m, max_n, m_d, n_d, k_d,
        A_array_d, lda_d,
        B_array_d, ldb_d,
        C_array_d, ldc_d,
        batchCount, stream, alpha);
    }
    else
    {
        if (max_n < 256)
        {
            vbatched_gemm_tn_impl<double, 16, 8, 32, 32, 8, 16, 8, 16, 8>
            (max_m, max_n, m_d, n_d, k_d,
            A_array_d, lda_d,
            B_array_d, ldb_d,
            C_array_d, ldc_d,
            batchCount, stream, alpha);
        }
        else
        {
            vbatched_gemm_tn_impl<double, 16, 16, 48, 48, 16, 16, 16, 16, 16>
            (max_m, max_n, m_d, n_d, k_d,
            A_array_d, lda_d,
            B_array_d, ldb_d,
            C_array_d, ldc_d,
            batchCount, stream, alpha);
        }
    }
}
