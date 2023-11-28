
#include "vbatch_matrix_mul.cuh"
#include "cuda_tools.cuh"
#define sA(i,j)    sA[(j)*slda + (i)]
#define sB(i,j)    sB[(j)*sldb + (i)]
#define fetch(A, m, n, bound)  offs_d##A[min(n*LD##A+m, bound)]

#define shared_A(i,j)    shared_A[(j)*shared_lda + (i)]
#define shared_B(i,j)    shared_B[(j)*shared_ldb + (i)]
/*template<typename T>
inline __device__ T fetch(const T* offs_d,const int  m,const int  n, const int bound, const int global_ld)
{
    return offs_d[min(n*global_ld+m, bound)];
}*/
__device__ inline double atomic_add(double* address, double val)
{
    return atomicAdd(address, val);
}

__device__ inline float atomic_add(float* address, float val)
{
    return atomicAdd(address, val);
}
__device__ inline cuFloatComplex atomic_add(cuFloatComplex* address, cuFloatComplex val)
{
    float re = atomicAdd( (float*) (&(*address).x) ,val.x);
    float im = atomicAdd( (float*) (&(*address).y) ,val.y);
    return make_cuFloatComplex(re, im);
}

__device__ inline cuDoubleComplex atomic_add(cuDoubleComplex* address, cuDoubleComplex val)
{
    double re = atomicAdd( (double*) (&(*address).x) ,val.x);
    double im = atomicAdd( (double*) (&(*address).y) ,val.y);
    return make_cuDoubleComplex(re, im);
}



template<typename T, const int DIM_X, const int DIM_Y, const int BLK_M, const int BLK_N, const int BLK_K,
         const int DIM_XA, const int DIM_YA, const int DIM_XB, const int DIM_YB,
         const int THR_M, const int THR_N>
static __device__
void vbatched_gemm_device(
    const int M, const int N, const int K,
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

    // Registers for the dev->shmem copy
    T ra[BLK_M/DIM_YA][BLK_K/DIM_XA];
    T rb[BLK_N/DIM_YB][BLK_K/DIM_XB];

    // bound is the correction to offs_d in order to not get out of memory bound
    // so bound could be negative value since offs_d could be out of bound
    const T *offs_dA = A + blx*BLK_M*LDA + idyA*LDA + idxA;
    int boundA = (LDA*(M-1) + K) - ( blx*BLK_M*LDA + idyA*LDA + idxA ) -1;

    const T *offs_dB = B + bly*BLK_N*LDB + idyB*LDB + idxB;
    int boundB = (LDB*(N-1) + K) - ( bly*BLK_N*LDB + idyB*LDB + idxB ) -1;

    int m, n, k, kk;

    // Zero C
    #pragma unroll
    for (n = 0; n < THR_N; n++)
        #pragma unroll
        for (m = 0; m < THR_M; m++)
            rC[n][m] = 0.0;

    // Load A dev->shmem
    #pragma unroll
    for (n = 0; n < BLK_M; n += DIM_YA)
        #pragma unroll
        for (m = 0; m < BLK_K; m += DIM_XA)
            sA(n+idyA, m+idxA) = fetch(A, m, n, boundA);

    #pragma unroll
    for (n = 0; n < BLK_N; n += DIM_YB)
        #pragma unroll
        for (m = 0; m < BLK_K; m += DIM_XB)
            sB(m+idxB, n+idyB) = fetch(B, m, n, boundB);

    __syncthreads();

    for (kk = 0; kk < K-BLK_K; kk += BLK_K) {
        offs_dA += BLK_K;
        boundA  -= BLK_K;

        offs_dB += BLK_K;
        boundB  -= BLK_K;

        // Load A dev->regs
        #pragma unroll
        for (n = 0; n < BLK_M/DIM_YA; n++)
            #pragma unroll
            for (m = 0; m < BLK_K/DIM_XA; m++)
                ra[n][m] = fetch(A, m*DIM_XA, n*DIM_YA, boundA);

        // Load B dev->regs
        #pragma unroll
        for (n = 0; n < BLK_N/DIM_YB; n++)
            #pragma unroll
            for (m = 0; m < BLK_K/DIM_XB; m++)
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
                    rC[n][m] += rA[m] * rB[n];
                }
            }
        }

        __syncthreads();

        // Load A regs->shmem
        #pragma unroll
        for (n = 0; n < BLK_M/DIM_YA; n++)
            #pragma unroll
            for (m = 0; m < BLK_K/DIM_XA; m++)
                sA(n*DIM_YA+idyA, m*DIM_XA+idxA) = ra[n][m];

        // Load B regs->shmem
        #pragma unroll
        for (n = 0; n < BLK_N/DIM_YB; n++)
            #pragma unroll
            for (m = 0; m < BLK_K/DIM_XB; m++)
                sB(m*DIM_XB+idxB, n*DIM_YB+idyB) = rb[n][m];

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
                rC[n][m] += rA[m] * rB[n];
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
                atomic_add(C + offsC, rC[n][m]);
            }
        }
    }
}


/******************************************************************************/
template <typename T, const int DIM_X, const int DIM_Y,
         const int BLK_M, const int BLK_N, const int BLK_K,
         const int DIM_XA, const int DIM_YA,
         const int DIM_XB, const int DIM_YB>
static __global__
void vbatched_gemm_kernel(
    const int* M, const int* N, const int K,
    T const * const * global_A_array, const int* global_lda,
    T const * const * global_B_array, const int* global_ldb,
    T              ** global_C_array, const int* global_ldc)
{
    //extern __shared__ __align__(sizeof(T)) unsigned char smem[];
    //T *shared_mem = reinterpret_cast<T *>(smem);
    extern __shared__ T* shared_mem[];

    const int batchid = blockIdx.z;
    int local_M = (int)M[batchid];
    int local_N = (int)N[batchid];

    if( blockIdx.x >= (local_M+BLK_M-1)/BLK_M ) return;
    if( blockIdx.y >= (local_N+BLK_N-1)/BLK_N ) return;

    const int shared_lda = BLK_M+1;
    const int shared_ldb = BLK_K+1;
    T* shared_A = (T*)shared_mem;
    T* shared_B = shared_A + shared_lda * BLK_K;

    vbatched_gemm_device<T, DIM_X, DIM_Y, 
                         BLK_M, BLK_N, BLK_K,
                         DIM_XA, DIM_YA,
                         DIM_XB, DIM_YB, 
                         (BLK_M/DIM_X), (BLK_N/DIM_Y)>
                        (local_M, local_N, K,
                        global_A_array[batchid], (int)global_lda[batchid],
                        global_B_array[batchid], (int)global_ldb[batchid],
                        global_C_array[batchid], (int)global_ldc[batchid],
                        shared_A, shared_lda, shared_B, shared_ldb);
}

static inline int ceildiv( int x, int y )
{
    return (x + y - 1)/y;
}

template <typename T, const int DIM_X, const int DIM_Y,
         const int BLK_M, const int BLK_N, const int BLK_K,
         const int DIM_XA, const int DIM_YA,
         const int DIM_XB, const int DIM_YB>
void vbatched_gemm_impl(const int max_m, const int max_n,
    const int* m, const int* n, const int k,
    T const * const * global_A_array, const int* global_lda,
    T const * const * global_B_array, const int* global_ldb,
    T ** global_C_array, const int* global_ldc,
    const int batchCount, cudaStream_t stream)
{
    size_t shared_mem_size = 0;
    shared_mem_size += (BLK_M+1) * BLK_K * sizeof(T);
    shared_mem_size += (BLK_K+1) * BLK_N * sizeof(T);
    dim3 dimBlock(DIM_X, DIM_Y);
    dim3 dimGrid(ceildiv( max_m, BLK_M ), ceildiv( max_n, BLK_N ), batchCount);

    vbatched_gemm_kernel<T, DIM_X, DIM_Y,
                         BLK_M, BLK_N, BLK_K,
                         DIM_XA, DIM_YA,
                         DIM_XB, DIM_YB>
                         <<<dimGrid, dimBlock, shared_mem_size, stream>>>
                         (m, n, k,
                         global_A_array, global_lda,
                         global_B_array, global_ldb,
                         global_C_array, global_ldc);
}
/*
template <>
void vbatch_gemm<float>(const int max_m, const int max_n,
                 const int* m, int* n, const int k,
                 float const * const * global_A_array, const int* global_lda,
                 float const * const * global_B_array, const int* global_ldb,
                 float ** global_C_array, const int* global_ldc,
                 const int batchCount, cudaStream_t stream)
{
    vbatched_gemm_impl<float, 4, 8, 8, 24, 8, 4, 8, 4, 8>
                    (max_m, max_n, m, n, k,
                    global_A_array, global_lda,
                    global_B_array, global_ldb,
                    global_C_array, global_ldc,
                    batchCount, stream);
}
*/
template <>
void vbatch_gemm<double>(const int max_m, const int max_n,
                 const int* m, int* n, const int k,
                 double  const * const * global_A_array, const int* global_lda,
                 double const * const * global_B_array, const int* global_ldb,
                 double ** global_C_array, const int* global_ldc,
                 const int batchCount, cudaStream_t stream)
{
    // The positions of A and B have been swapped here.
    // This is because the original code is for column-major matrices.
    // We use row-major matrices, so we need to swap A and B.
    // The vbatched_gemm_impl is for C = trans(A) * B + C, but we need trans(C).
    // Which means: trans(C) = trans(trans(A)*B + C) = trans(B) * A + trans(C)
    // Then, ldc should be N, lda and ldb should be K

    vbatched_gemm_impl<double, 16, 16, 48, 32, 16, 16, 16, 16, 16>
                    (max_n, max_m, n, m, k,
                    global_B_array, global_ldb,
                    global_A_array, global_lda,
                    global_C_array, global_ldc,
                    batchCount, stream);
}

