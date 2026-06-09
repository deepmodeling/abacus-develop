#pragma once

#include <cuda_runtime.h>

// Shape-exact batched GEMM dispatchers.
//
// Every (A_i, B_i, C_i) in the batch has exactly the same (m, n, k); the
// caller (phi_operator_gpu.cu) guarantees this by bucketing atom pairs on
// (nw1, nw2) before calling in. The scalar m/n/k drive tile-ladder selection,
// grid sizing, and the kernel itself -- there is no per-batch-id M/N/K
// indirection.
//
// The C output is always double, independent of T. For T=float the per-item
// inner products accumulate in fp32, but the cross-item accumulation into a
// shared C element is done with a device-side fp64 atomicAdd (see the kernels'
// store loop), so summing many atom-pair contributions into the same
// hr_gint / phi_dm element does not drift. For T=double, A, B and C are all
// double.

// C(batch) = alpha * A(batch) * B(batch) + C(batch)
template<typename T>
void gemm_nn_vbatch(
    int m, int n, int k,
    const T* const* A_array_d, const int* lda_d,
    const T* const* B_array_d, const int* ldb_d,
    double** C_array_d, const int* ldc_d,
    int batchCount, cudaStream_t stream,
    const T* alpha = nullptr);

// C(batch) = alpha * A(batch)^T * B(batch) + C(batch)
template<typename T>
void gemm_tn_vbatch(
    int m, int n, int k,
    const T* const* A_array_d, const int* lda_d,
    const T* const* B_array_d, const int* ldb_d,
    double** C_array_d, const int* ldc_d,
    int batchCount, cudaStream_t stream,
    const T* alpha = nullptr);
