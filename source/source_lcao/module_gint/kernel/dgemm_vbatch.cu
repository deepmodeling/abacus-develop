#include "gemm_tn_vbatch.cuh"
#include "gemm_nn_vbatch.cuh"
#include "dgemm_vbatch.h"
#include "source_base/module_device/device.h"

// Tile ladder
// -----------
// The caller splits each batch into buckets of identical (m, n, k) and calls
// in once per bucket. The dispatchers below pick, for each bucket, the kernel
// instantiation whose (BLK_M, BLK_N) tile is the smallest rung that still
// covers the bucket's output shape, so boundary blocks don't spend most of
// their work on masked-off padding.
//
// Each thread owns a THR_M x THR_N register accumulator tile, i.e. it computes
//   THR = THR_M * THR_N = (BLK_M / DIM_X) * (BLK_N / DIM_Y)
// output elements. We aim to keep THR in roughly [16, 36]: below that the inner
// FMAs don't amortize the shared-memory traffic and there's too little ILP;
// above it register pressure starts cutting occupancy. The "(in band)" /
// "(under)" notes on each case below mark where that rung lands.

template<typename T>
void gemm_nn_vbatch(
    int m, int n, int k,
    const T* const* A_array_d, const int* lda_d,
    const T* const* B_array_d, const int* ldb_d,
    double** C_array_d, const int* ldc_d,
    int batchCount, cudaStream_t stream,
    const T* alpha)
{
    // 4 (nw2 bracket) x 2 (bxyz bracket) = 8 instantiations.
    //
    // Mapping into the impl's parameter list is:
    //   <T, DIM_X, DIM_Y, BLK_M, BLK_N, BLK_K=16,
    //    DIM_XA=DIM_X, DIM_YA=DIM_Y, DIM_XB=DIM_X, DIM_YB=DIM_Y>
    // which satisfies the kernel's tile-divisibility asserts because every
    // (BLK_M, BLK_N, BLK_K=16) chosen below is a multiple of the matching
    // (DIM_X, DIM_Y) pair.
    #define NN_DISPATCH(DX, DY, BM, BN)                                 \
        vbatched_gemm_nn_impl<T, DX, DY, BM, BN, 16, DX, DY, DX, DY>(   \
            m, n, k,                                                    \
            A_array_d, lda_d, B_array_d, ldb_d,                         \
            C_array_d, ldc_d, batchCount, stream, alpha)

    // BLK_M bracket -- smallest tile in {8,16,32,48} covering nw2.
    const int blk_m_tag = (n <=  8) ? 0
                        : (n <= 16) ? 1
                        : (n <= 32) ? 2
                        :             3;

    // BLK_N bracket -- tiles the bxyz (mesh-grid) axis. Use 32 when bxyz<=32 so
    // a partial final block-row isn't mostly masked padding (e.g. bxyz=27 in a
    // 64-row tile leaves ~58% of the rows idle); use 64 above that, where the
    // larger tile gives better shared-memory reuse.
    const int blk_n_tag = (m <= 32) ? 0 : 1;

    switch (blk_m_tag * 2 + blk_n_tag)
    {
        // BLK_M=8  (nw2 <=8 ).  DIM=4x8  -> THR_M=2.
        case 0: NN_DISPATCH( 4, 8,  8, 32); break;  // THR=2*4=8   (under)
        case 1: NN_DISPATCH( 4, 8,  8, 64); break;  // THR=2*8=16  (in band)
        // BLK_M=16 (nw2<=16).  DIM=4x8  -> THR_M=4.
        case 2: NN_DISPATCH( 4, 8, 16, 32); break;  // THR=4*4=16  (in band)
        case 3: NN_DISPATCH( 4, 8, 16, 64); break;  // THR=4*8=32  (in band)
        // BLK_M=32 (nw2<=32).  DIM=8x8  -> THR_M=4.
        case 4: NN_DISPATCH( 8, 8, 32, 32); break;  // THR=4*4=16  (in band)
        case 5: NN_DISPATCH( 8, 8, 32, 64); break;  // THR=4*8=32  (in band)
        // BLK_M=48 (nw2<=48).  DIM=16x8 -> THR_M=3 (cap at 3 to keep
        // register pressure room for the BLK_N=64 sibling).
        case 6: NN_DISPATCH(16, 8, 48, 32); break;  // THR=3*4=12  (just under)
        case 7: NN_DISPATCH(16, 8, 48, 64); break;  // THR=3*8=24  (in band)
    }

    #undef NN_DISPATCH
}

template<typename T>
void gemm_tn_vbatch(
    int m, int n, int k,
    const T* const* A_array_d, const int* lda_d,
    const T* const* B_array_d, const int* ldb_d,
    double** C_array_d, const int* ldc_d,
    int batchCount, cudaStream_t stream,
    const T* alpha)
{
    // 4 (nw2 bracket) x 4 (nw1 bracket) = 16 instantiations.
    //
    // Both output axes here are the small nw axis, so we use the same
    // {8,16,32,48} ladder on both. BLK_K = 32 (the bxyz axis -- large).
    #define TN_DISPATCH(DX, DY, BM, BN)                                 \
        vbatched_gemm_tn_impl<T, DX, DY, BM, BN, 32, DX, DY, DX, DY>(   \
            m, n, k,                                                    \
            A_array_d, lda_d, B_array_d, ldb_d,                         \
            C_array_d, ldc_d, batchCount, stream, alpha)

    auto bracket = [](int x) {
        return (x <=  8) ? 0
             : (x <= 16) ? 1
             : (x <= 32) ? 2
             :             3;
    };
    const int blk_m_tag = bracket(n);  // BLK_M <- nw2
    const int blk_n_tag = bracket(m);  // BLK_N <- nw1

    switch (blk_m_tag * 4 + blk_n_tag)
    {
        // BLK_M=8  rungs (nw2<=8).  DIM_X=4, THR_M=2.
        case  0: TN_DISPATCH(4, 8,  8,  8); break;  // THR=2*1=2  (well under band)
        case  1: TN_DISPATCH(4, 8,  8, 16); break;  // THR=2*2=4
        case  2: TN_DISPATCH(4, 8,  8, 32); break;  // THR=2*4=8
        case  3: TN_DISPATCH(4, 8,  8, 48); break;  // THR=2*6=12
        // BLK_M=16 rungs (nw2<=16). DIM_X=4, THR_M=4.
        case  4: TN_DISPATCH(4, 8, 16,  8); break;  // THR=4*1=4
        case  5: TN_DISPATCH(4, 8, 16, 16); break;  // THR=4*2=8
        case  6: TN_DISPATCH(4, 8, 16, 32); break;  // THR=4*4=16  (in band)
        case  7: TN_DISPATCH(4, 8, 16, 48); break;  // THR=4*6=24  (in band)
        // BLK_M=32 rungs (nw2<=32). DIM_X=8, THR_M=4.
        case  8: TN_DISPATCH(8, 4, 32,  8); break;  // THR=4*2=8
        case  9: TN_DISPATCH(8, 4, 32, 16); break;  // THR=4*4=16  (in band)
        case 10: TN_DISPATCH(8, 8, 32, 32); break;  // THR=4*4=16  (in band)
        case 11: TN_DISPATCH(8, 8, 32, 48); break;  // THR=4*6=24  (in band)
        // BLK_M=48 rungs (nw2<=48). DIM_X=8, THR_M=6.
        case 12: TN_DISPATCH(8, 4, 48,  8); break;  // THR=6*2=12
        case 13: TN_DISPATCH(8, 4, 48, 16); break;  // THR=6*4=24  (in band)
        case 14: TN_DISPATCH(8, 8, 48, 32); break;  // THR=6*4=24  (in band)
        case 15: TN_DISPATCH(8, 8, 48, 48); break;  // THR=6*6=36  (top of band)
    }

    #undef TN_DISPATCH
}

// Explicit instantiations
template void gemm_nn_vbatch<double>(
    int, int, int,
    const double* const*, const int*, const double* const*, const int*,
    double**, const int*, int, cudaStream_t, const double*);

template void gemm_nn_vbatch<float>(
    int, int, int,
    const float* const*, const int*, const float* const*, const int*,
    double**, const int*, int, cudaStream_t, const float*);

template void gemm_tn_vbatch<double>(
    int, int, int,
    const double* const*, const int*, const double* const*, const int*,
    double**, const int*, int, cudaStream_t, const double*);

template void gemm_tn_vbatch<float>(
    int, int, int,
    const float* const*, const int*, const float* const*, const int*,
    double**, const int*, int, cudaStream_t, const float*);
