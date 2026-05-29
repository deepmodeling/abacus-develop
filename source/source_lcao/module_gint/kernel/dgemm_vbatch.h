#pragma once

#include <cuda_runtime.h>

// Shape-exact batched GEMM dispatchers.
//
// Every (A_i, B_i, C_i) in the batch has exactly the same (m, n, k); the
// caller (phi_operator_gpu.cu) enforces this by bucketing atom pairs on
// (nw1, nw2). The scalars drive tile-ladder selection, grid sizing, and
// flow all the way through the kernel -- there is no per-batchid M/N/K
// indirection left.
//
// The C accumulator is always double regardless of the input type T: a fp32
// GEMM path (T=float) feeds fp32 multiplies into fp64 registers and a
// device-side fp64 atomicAdd, so summing many atom-pair contributions into the
// same hr_gint / phi_dm element does not drift. For T=double, A, B and C are
// all double and this matches the legacy signature.

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
