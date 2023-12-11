#include <functional>
#include "vbatch_matrix_mul.cuh"
#include "cuda_tools.cuh"
#include "module_base/blas_connector.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

#define sA(i,j)    sA[(j)*slda + (i)]
#define sB(i,j)    sB[(j)*sldb + (i)]
#define fetch(A, m, n, bound)  offs_d##A[min(n*LD##A+m, bound)]

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

template<typename T, int DIM_X, int DIM_Y, int BLK_M, int BLK_N, int BLK_K,
         int DIM_XA, int DIM_YA, int DIM_XB, int DIM_YB,
         int THR_M, int THR_N>
static __device__
void vbatched_gemm_device(
    int M, int N, int K,
    T* __restrict__ A, int LDA,
    T* __restrict__ B, int LDB,
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
    T *offs_dA = A + blx*BLK_M*LDA + idyA*LDA + idxA;
    int boundA = (LDA*(M-1) + K) - ( blx*BLK_M*LDA + idyA*LDA + idxA ) -1;

    T *offs_dB = B + bly*BLK_N*LDB + idyB*LDB + idxB;
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
template <typename T, int DIM_X, int DIM_Y,
         int BLK_M, int BLK_N, int BLK_K,
         int DIM_XA, int DIM_YA,
         int DIM_XB, int DIM_YB>
static __global__
void vbatched_gemm_kernel(
    int* M, int* N, int K,
    T * * global_A_array, int* global_lda,
    T * * global_B_array, int* global_ldb,
    T              ** global_C_array, int* global_ldc)
{
    extern __shared__ __align__(sizeof(T)) unsigned char smem[];
    T *shared_mem = reinterpret_cast<T *>(smem);

    int batchid = blockIdx.z;
    int local_M = (int)M[batchid];
    int local_N = (int)N[batchid];

    if( blockIdx.x >= (local_M+BLK_M-1)/BLK_M ) return;
    if( blockIdx.y >= (local_N+BLK_N-1)/BLK_N ) return;

    int shared_lda = BLK_M+1;
    int shared_ldb = BLK_K+1;
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

template <typename T, int DIM_X, int DIM_Y,
         int BLK_M, int BLK_N, int BLK_K,
         int DIM_XA, int DIM_YA,
         int DIM_XB, int DIM_YB>
void vbatched_gemm_impl(int max_m, int max_n,
                 int* m, int* n, int k,
                 T  * * global_A_array, int* global_lda,
                 T * * global_B_array, int* global_ldb,
                 T ** global_C_array, int* global_ldc,
                 int batchCount, cudaStream_t stream)
{
    // The positions of A and B have been swapped here.
    // This is because the original code is for column-major matrices.
    // We use row-major matrices, so we need to swap A and B.
    // The vbatched_gemm_impl is for C = trans(A) * B + C, but we need trans(C).
    // Which means: trans(C) = trans(trans(A)*B + C) = trans(B) * A + trans(C)
    // Then, ldc should be N, lda and ldb should be K

    size_t shared_mem_size = 0;
    shared_mem_size += (BLK_M+1) * BLK_K * sizeof(T);
    shared_mem_size += (BLK_K+1) * BLK_N * sizeof(T);
    dim3 dimBlock(DIM_X, DIM_Y);
    dim3 dimGrid(ceildiv( max_n, BLK_M ), ceildiv( max_m, BLK_N ), batchCount);

    vbatched_gemm_kernel<T, DIM_X, DIM_Y,
                         BLK_M, BLK_N, BLK_K,
                         DIM_XA, DIM_YA,
                         DIM_XB, DIM_YB>
                         <<<dimGrid, dimBlock, shared_mem_size, stream>>>
                         (n, m, k,
                         global_B_array, global_ldb,
                         global_A_array, global_lda,
                         global_C_array, global_ldc);
}

template <typename T, int DIM_X, int DIM_Y, int BLK_M, int BLK_N, int BLK_K,
          int DIM_XA, int DIM_YA, int DIM_XB, int DIM_YB>
void gemm_time_measure(int max_m, int max_n,
                 int* m, int* n, int k,
                 T ** global_A_array, int* global_lda,
                 T ** global_B_array, int* global_ldb,
                 T ** global_C_array, int* global_ldc,
                 int batchCount, cudaStream_t stream, float &fast_time, func_type &fastest_algo,
                 double *cpu_result, double * h_global_C, double *d_global_C)
{
    cudaEvent_t start, stop;
    cudaMemset(d_global_C, 0, batchCount * max_m * max_n * sizeof(double));
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, stream);
    vbatched_gemm_impl<T, DIM_X, DIM_Y, BLK_M, BLK_N, BLK_K, DIM_XA, DIM_YA, DIM_XB, DIM_YB>
                    (max_m, max_n, m, n, k,
                    global_A_array, global_lda,
                    global_B_array, global_ldb,
                    global_C_array, global_ldc,
                    batchCount, stream);
    cudaEventRecord(stop, stream);
    cudaError_t cuda_status = cudaGetLastError();
    cudaStreamSynchronize(stream);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);

    // WARNING !!!!! Here we assume that all m and n are the same
    cudaMemcpy(h_global_C, d_global_C, batchCount * max_m * max_n * sizeof(double), cudaMemcpyDeviceToHost);  
    bool check_result = true;
    for (int i = 0; i < batchCount * max_m * max_n; ++i)
    {
        if (abs(cpu_result[i] - h_global_C[i]) > 0.001)
        {
            check_result = false;
            break;
        }
    }
    if (milliseconds < fast_time && cuda_status == cudaSuccess && check_result)
    {
        fast_time = milliseconds;
        fastest_algo = vbatched_gemm_impl<T, DIM_X, DIM_Y, BLK_M, BLK_N, BLK_K, DIM_XA, DIM_YA, DIM_XB, DIM_YB>;
        #define __DEBUG
        #ifdef __DEBUG
        std::cout << "found! fastest time: " << fast_time << std::endl;
        std::cout << DIM_X << ","<< DIM_Y<< ","<< BLK_M<< ","<< BLK_N<< ","<< BLK_K<< ","<< DIM_XA<< ","<< DIM_YA<< ","<< DIM_XB<< ","<< DIM_YB << std::endl;
        #endif
    }
}

void gemm_algo_selector(int matrix_k, func_type & fastest_algo)
{

    int batchCount_per_type = 32;
    int batchCount = batchCount_per_type * GlobalC::ucell.ntype * GlobalC::ucell.ntype;
    int *h_m = new int[batchCount];
    int *h_n = new int[batchCount];
    int *h_global_lda = new int[batchCount];
    int *h_global_ldb = new int[batchCount];
    int *h_global_ldc = new int[batchCount];
    int max_m = GlobalC::ucell.nwmax, max_n = GlobalC::ucell.nwmax;
    double **h_global_A_array = new double *[batchCount];
    double **h_global_B_array = new double *[batchCount];
    double **h_global_C_array = new double *[batchCount];

    double *h_global_A = new double[batchCount * max_m * matrix_k];
    double *h_global_B = new double[batchCount * max_n * matrix_k];
    double *h_global_C = new double[batchCount * max_m * max_n];

    for (int i = 0; i < batchCount * max_m * matrix_k; ++i)
    {
        h_global_A[i] = i * 0.001;
    }
    for (int i = 0; i < batchCount * max_n * matrix_k; ++i)
    {
        h_global_B[i] = i * 0.002;
    }
    memset(h_global_C, 0, batchCount * max_m * max_n * sizeof(double));

    // Allocate device memory
    int *d_m;
    int *d_n;
    int *d_global_lda;
    int *d_global_ldb;
    int *d_global_ldc;

    double **d_global_A_array;
    double **d_global_B_array;
    double **d_global_C_array;

    double *d_global_A;
    double *d_global_B;
    double *d_global_C;

    double *cpu_result = new double[batchCount * max_m * max_n];
    memset(cpu_result, 0, batchCount * max_m * max_n * sizeof(double));

    cudaMalloc(&d_m, batchCount * sizeof(int));
    cudaMalloc(&d_n, batchCount * sizeof(int));
    cudaMalloc(&d_global_lda, batchCount * sizeof(int));
    cudaMalloc(&d_global_ldb, batchCount * sizeof(int));
    cudaMalloc(&d_global_ldc, batchCount * sizeof(int));
    cudaMalloc(&d_global_A_array, batchCount * sizeof(double *));
    cudaMalloc(&d_global_B_array, batchCount * sizeof(double *));
    cudaMalloc(&d_global_C_array, batchCount * sizeof(double *));

    cudaMalloc(&d_global_A, batchCount * max_m * matrix_k * sizeof(double));
    cudaMalloc(&d_global_B, batchCount * max_n * matrix_k * sizeof(double));
    cudaMalloc(&d_global_C, batchCount * max_m * max_n * sizeof(double));
     
    cudaMemset(d_global_C, 0, batchCount * max_m * max_n * sizeof(double));
    int index = 0;
    for (int i = 0; i < batchCount_per_type; ++i)
    {
        for (int j = 0; j < GlobalC::ucell.ntype; j++)
        {
            for (int k = 0; k < GlobalC::ucell.ntype; k++)
            {
                h_m[index] = GlobalC::ucell.atoms[j].nw;
                h_n[index] = GlobalC::ucell.atoms[k].nw;
                h_global_lda[index] = matrix_k;
                h_global_ldb[index] = matrix_k;
                h_global_ldc[index] = GlobalC::ucell.atoms[k].nw;

                h_global_A_array[index] = &d_global_A[index * max_m * matrix_k];
                h_global_B_array[index] = &d_global_B[index * max_n * matrix_k];
                h_global_C_array[index] = &d_global_C[index * max_n * max_m]; // test atom add
                BlasConnector::gemm('N', 'T', h_m[index], h_n[index], matrix_k, 1.0, &h_global_A[index * max_m * matrix_k], matrix_k, &h_global_B[index * max_n * matrix_k], matrix_k, 1.0, &cpu_result[index * max_m * max_n], h_n[index]);
                index++;
            }
        }
    }

    cudaMemcpy(d_m, h_m, batchCount * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_n, h_n, batchCount * sizeof(int), cudaMemcpyHostToDevice);

    cudaMemcpy(d_global_lda, h_global_lda, batchCount * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_global_ldb, h_global_ldb, batchCount * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_global_ldc, h_global_ldc, batchCount * sizeof(int), cudaMemcpyHostToDevice);

    cudaMemcpy(d_global_A_array, h_global_A_array, batchCount * sizeof(double *), cudaMemcpyHostToDevice);
    cudaMemcpy(d_global_B_array, h_global_B_array, batchCount * sizeof(double *), cudaMemcpyHostToDevice);
    cudaMemcpy(d_global_C_array, h_global_C_array, batchCount * sizeof(double *), cudaMemcpyHostToDevice);

    cudaMemcpy(d_global_A, h_global_A, batchCount * max_m * matrix_k * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_global_B, h_global_B, batchCount * max_n * matrix_k * sizeof(double), cudaMemcpyHostToDevice);

    cudaStream_t temp_stream;
    cudaStreamCreate(&temp_stream);

    float fastest_time = 1000000;
    fastest_algo = vbatched_gemm_impl<double, 16, 4, 32, 16, 16, 16, 4, 16, 4>;
    #include"code_gen.cpp"
    cudaStreamDestroy(temp_stream);
    std::cout << " gemm_algo_selector::Fastest time: " << fastest_time << " ms" << std::endl;
    // fastest_algo = vbatched_gemm_impl<double, 16, 4, 32, 16, 16, 16, 4, 16, 4>;
    delete[] h_global_A_array;
    delete[] h_global_B_array;
    delete[] h_global_C_array;

    delete[] h_m;
    delete[] h_n;

    delete[] h_global_lda;
    delete[] h_global_ldb;
    delete[] h_global_ldc;

    delete[] h_global_A;
    delete[] h_global_B;
    delete[] h_global_C;

    delete[] cpu_result;

    // Cleanup
    cudaFree(d_global_A_array);
    cudaFree(d_global_B_array);
    cudaFree(d_global_C_array);

    cudaFree(d_m);
    cudaFree(d_n);

    cudaFree(d_global_lda);
    cudaFree(d_global_ldb);
    cudaFree(d_global_ldc);

    cudaFree(d_global_A_array);
    cudaFree(d_global_B_array);
    cudaFree(d_global_C_array);

    cudaFree(d_global_A);
    cudaFree(d_global_B);
    cudaFree(d_global_C);

}