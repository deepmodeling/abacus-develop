
#include "vbatch_matrix_mul.cuh"
#include "cuda_tools.cuh"


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
    T* sB, int sldb)
{
    int idx = threadIdx.x;  // thread's m dimension
    int idy = threadIdx.y;  // thread's n dimension

    int idt = DIM_X * idy + idx;    // thread's global number

    int idxA = idt % DIM_XA;    // idx within A
    int idyA = idt / DIM_XA;    // idy within A

    int idxB = idt % DIM_XB;    // idx within B
    int idyB = idt / DIM_XB;    // idy within B

    int blx = blockIdx.x;   // block's m dimension
    int bly = blockIdx.y;   // block's n dimension

    // Registers for the innermost loop
    T rC[THR_N][THR_M];
    T rA[THR_M];
    T rB[THR_N];

    T ra[BLK_K/DIM_YA][BLK_M/DIM_XA];
    T rb[BLK_K/DIM_YB][BLK_N/DIM_XB];

    const T *offs_dA = A + blx*BLK_M     + idyA*LDA + idxA;
    int boundA = (LDA*(K-1) + M) - ( blx*BLK_M  + idyA*LDA + idxA ) -1;

    const T *offs_dB = B + bly*BLK_N     + idyB*LDB + idxB;
    int boundB = (LDB*(K-1) + N) - ( bly*BLK_N     + idyB*LDB + idxB ) -1;

    int m, n, k, kk;


    // Zero C
    #pragma unroll
    for (n = 0; n < THR_N; n++)
        #pragma unroll
        for (m = 0; m < THR_M; m++)
            rC[n][m] = make_FloatingPoint(0.0, 0.0);

    #pragma unroll
    for (n = 0; n < BLK_K; n += DIM_YA)
        #pragma unroll
        for (m = 0; m < BLK_M; m += DIM_XA)
            sA(m+idxA, n+idyA) = fetch(A, m, n, boundA);

    // Load B dev->shmem
    #pragma unroll
    for (n = 0; n < BLK_K; n += DIM_YB)
        #pragma unroll
        for (m = 0; m < BLK_N; m += DIM_XB)
            sB(n+idyB, m+idxB) = fetch(B, m, n, boundB);

    __syncthreads();

    for (kk = 0; kk < K-BLK_K; kk += BLK_K) {
        offs_dA += BLK_K*LDA;
        boundA  -= BLK_K*LDA;

        offs_dB += BLK_K*LDB;
        boundB  -= BLK_K*LDB;

        // Load A dev->regs
        #pragma unroll
        for (n = 0; n < BLK_K/DIM_YA; n++)
            #pragma unroll
            for (m = 0; m < BLK_M/DIM_XA; m++)
                ra[n][m] = fetch(A, m*DIM_XA, n*DIM_YA, boundA);

        // Load B dev->regs
        #pragma unroll
        for (n = 0; n < BLK_K/DIM_YB; n++)
            #pragma unroll
            for (m = 0; m < BLK_N/DIM_XB; m++)
                rb[n][m] = fetch(B, m*DIM_XB, n*DIM_YB, boundB);

        // Multiply
        #pragma unroll
        for (k = 0; k < BLK_K; k++) {
            // Load A shmem->regs
            #pragma unroll
            for (m = 0; m < THR_M; m++)
                rA[m] = sA(m*DIM_X+idx, k);

            // Load B shmem->regs
            #pragma unroll
            for (n = 0; n < THR_N; n++)
                rB[n] = sB(k, n*DIM_Y+idy);

            // Compute
            #pragma unroll
            for (n = 0; n < THR_N; n++) {
                #pragma unroll
                for (m = 0; m < THR_M; m++) {
                    fma(rA[m], rB[n], rC[n][m]);
                }
            }
        }

        __syncthreads();

        // Load A regs->shmem
        #pragma unroll
        for (n = 0; n < BLK_K/DIM_YA; n++)
            #pragma unroll
            for (m = 0; m < BLK_M/DIM_XA; m++)
                sA(m*DIM_XA+idxA, n*DIM_YA+idyA) = ra[n][m];

        // Load B regs->shmem
        #pragma unroll
        for (n = 0; n < BLK_K/DIM_YB; n++)
            #pragma unroll
            for (m = 0; m < BLK_N/DIM_XB; m++)
                sB(n*DIM_YB+idyB, m*DIM_XB+idxB) = rb[n][m];
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
            rA[m] = sA(m*DIM_X+idx, k);

        // Load B shmem->regs
        #pragma unroll
        for (n = 0; n < THR_N; n++)
            rB[n] = sB(k, n*DIM_Y+idy);

        // Compute
        #pragma unroll
        for (n = 0; n < THR_N; n++) {
            #pragma unroll
            for (m = 0; m < THR_M; m++) {
                fma(rA[m], rB[n], rC[n][m]);
            }
        }
    }

    // Store C regs->dev
    #pragma unroll
    for (n = 0; n < THR_N; n++) {
        int coord_dCn = bly*BLK_N + n*DIM_Y + idy;
        #pragma unroll
        for (m = 0; m < THR_M; m++) {
            int coord_dCm = blx*BLK_M + m*DIM_X + idx;
            if (coord_dCm < M && coord_dCn < N) {
                int offsC = coord_dCn*LDC + coord_dCm;

                T &regC = rC[n][m];
                atomic_add(C+offsC, regC);
            }
        }
    }
}


/******************************************************************************/
template <typename T, const int DIM_X, const int DIM_Y, const int BLK_M, const int BLK_N, const int BLK_K,
         const int DIM_XA, const int DIM_YA, const int DIM_XB, const int DIM_YB>
static __global__
void vbatched_gemm_kernel(
    int* M, int* N, int K,
    T const * const * Aarray, int* LDA,
    T const * const * Barray, int* LDB,
    T              ** Carray, int* LDC)
{
    extern __shared__ T* sdata_nt[];

    const int batchid = blockIdx.z;
    int my_M = (int)M[batchid];
    int my_N = (int)N[batchid];

    if( blockIdx.x >= (my_M+BLK_M-1)/BLK_M ) return;
    if( blockIdx.y >= (my_N+BLK_N-1)/BLK_N ) return;

    const int slda = BLK_M+1;    // +1 only required if A is transposed
    const int sldb = BLK_K+1;    // +1 always required
    T* sA = (T*)sdata_nt;        // sA is (BLK_M+1) x (BLK_K)
    T* sB = sA + slda * BLK_K;   // sB is (BLK_K+1) x (BLK_N)

    vbatched_gemm_device<T, DIM_X, DIM_Y, BLK_M, BLK_N, BLK_K, DIM_XA, DIM_YA, DIM_XB, DIM_YB, (BLK_M/DIM_X), (BLK_N/DIM_Y)>
    ( my_M, my_N, K,
      Aarray[batchid], (int)LDA[batchid],
      Barray[batchid], (int)LDB[batchid],
      Carray[batchid], (int)LDC[batchid],
      sA, slda, sB, sldb);
}

static inline int ceildiv( int x, int y )
{
    return (x + y - 1)/y;
}

template <typename T, const int DIM_X, const int DIM_Y, const int BLK_M, const int BLK_N, const int BLK_K, const int dim_vec,
         const int DIM_XA, const int DIM_YA, const int DIM_XB, const int DIM_YB>
void vbatched_gemm(
    int max_m, int max_n,
    int* m, int* n, int k,
    T const * const * dA_array, int* ldda,
    T const * const * dB_array, int* lddb,
    T ** dC_array, int* lddc,
    int batchCount, cudaStream_t stream)
{
    size_t shmem = 0;
    shmem += (BLK_M+1) * BLK_K * sizeof(T);  // sA
    shmem += (BLK_K+1) * BLK_N * sizeof(T);  // sB
    dim3 dimBlock(DIM_X, DIM_Y);
    dim3 dimGrid( ceildiv( max_m, BLK_M ), ceildiv( max_n, BLK_N ), batchCount );

    vbatched_gemm_kernel<T, DIM_X, DIM_Y, BLK_M, BLK_N, BLK_K, DIM_XA, DIM_YA, DIM_XB, DIM_YB>
    <<<dimGrid, dimBlock, shmem, stream>>>
    (m, n, k, dA_array, ldda, dB_array, lddb, dC_array, lddc);
}
