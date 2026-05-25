#pragma once

#include <cuda_runtime.h>

// Template version: C(batch_id) = alpha * A(batch_id) * B(batch_id) + C(batch_id)
// As with gemm_tn_vbatch, the C accumulator is always double regardless of the
// input type T so the per-block reduction and device-side atomicAdd run in fp64.
template<typename T>
void gemm_nn_vbatch(
    int max_m, int max_n, int max_k,
    const int* m_d, const int* n_d, const int* k_d,
    const T* const* A_array_d, const int* lda_d,
    const T* const* B_array_d, const int* ldb_d,
    double** C_array_d, const int* ldc_d,
    int batchCount, cudaStream_t stream,
    const T* alpha = nullptr);

// Template version: C(batch_id) = alpha * A(batch_id)^T * B(batch_id) + C(batch_id)
// The C accumulator is always double regardless of input type T: a fp32 GEMM
// path (T=float) feeds fp32 multiplies into fp64 accumulators (registers and
// device-side atomicAdds) to avoid catastrophic precision loss across many
// atom-pair contributions to the same hr_gint element.
template<typename T>
void gemm_tn_vbatch(
    int max_m, int max_n, int max_k,
    const int* m_d, const int* n_d, const int* k_d,
    const T* const* A_array_d, const int* lda_d,
    const T* const* B_array_d, const int* ldb_d,
    double** C_array_d, const int* ldc_d,
    int batchCount, cudaStream_t stream,
    const T* alpha = nullptr);

// Legacy double-only aliases for backward compatibility
inline void dgemm_nn_vbatch(
    int max_m, int max_n, int max_k,
    const int* m_d, const int* n_d, const int* k_d,
    const double* const* A_array_d, const int* lda_d,
    const double* const* B_array_d, const int* ldb_d,
    double** C_array_d, const int* ldc_d,
    int batchCount, cudaStream_t stream,
    const double* alpha = nullptr)
{
    gemm_nn_vbatch<double>(max_m, max_n, max_k,
        m_d, n_d, k_d, A_array_d, lda_d, B_array_d, ldb_d,
        C_array_d, ldc_d, batchCount, stream, alpha);
}

inline void dgemm_tn_vbatch(
    int max_m, int max_n, int max_k,
    const int* m_d, const int* n_d, const int* k_d,
    const double* const* A_array_d, const int* lda_d,
    const double* const* B_array_d, const int* ldb_d,
    double** C_array_d, const int* ldc_d,
    int batchCount, cudaStream_t stream,
    const double* alpha = nullptr)
{
    // T=double path: A, B, and C are all double — the C-channel double-fix
    // matches the legacy signature here.
    gemm_tn_vbatch<double>(max_m, max_n, max_k,
        m_d, n_d, k_d, A_array_d, lda_d, B_array_d, ldb_d,
        C_array_d, ldc_d, batchCount, stream, alpha);
}