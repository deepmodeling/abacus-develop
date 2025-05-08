#pragma once

#include <cuda_runtime.h>

// C(batch_id) = alpha * A(batch_id) * B(batch_id) + C(batch_id)
void dgemm_nn_vbatch(
    int max_m, int max_n, int max_k,
    const int* m_d, const int* n_d, const int* k_d,
    const double* const* A_array_d, const int* lda_d,
    const double* const* B_array_d, const int* ldb_d,
    double** C_array_d, const int* ldc_d,
    int batchCount, cudaStream_t stream,
    const double* alpha = nullptr);

// C(batch_id) = alpha * A(batch_id)^T * B(batch_id) + C(batch_id)
void dgemm_tn_vbatch(
    int max_m, int max_n, int max_k,
    const int* m_d, const int* n_d, const int* k_d,
    const double* const* A_array_d, const int* lda_d,
    const double* const* B_array_d, const int* ldb_d,
    double** C_array_d, const int* ldc_d,
    int batchCount, cudaStream_t stream,
    const double* alpha = nullptr);