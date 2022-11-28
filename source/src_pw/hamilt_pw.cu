//obsolete code
//please remove the globalc::hm

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "global.h"
#include "hamilt_pw.cuh"
#include "../module_base/blas_connector.h"
#include "myfunc.h"
#include "../module_base/timer.h"
#include "../src_parallel/parallel_reduce.h"
using namespace CudaCheck;

__global__ void cast_d2f(float *dst, double *src, int size)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i < size)
    {
        dst[i] = __double2float_rn(src[i]);
    }
}

__global__ void cast_f2d(double *dst, float *src, int size)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i < size)
    {
        dst[i] = (double)(src[i]);
    }
}

__global__ void cast_d2f(float2 *dst, double2 *src, int size)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i < size)
    {
        dst[i].x = __double2float_rn(src[i].x);
        dst[i].y = __double2float_rn(src[i].y);
    }
}

__global__ void cast_f2d(double2 *dst, float2 *src, int size)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i < size)
    {
        dst[i].x = (double)(src[i].x);
        dst[i].y = (double)(src[i].y);
    }
}


template<class T2>
__global__ void kernel_copy(int size, T2* dst, const T2 *src)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        dst[idx].x = src[idx].x;
        dst[idx].y = src[idx].y;
    }
}

template<class T, class T2>
__global__ void kernel_get_tmhpsi(int size, T2 *dst, const T2 *src, T *g2kin)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        dst[idx].x = src[idx].x * g2kin[idx];
        dst[idx].y = src[idx].y * g2kin[idx];
    }
}

template<class T2>
__global__ void kernel_add_tmhpsi(int size, T2 *dst, T2 *src, int *index)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int p = index[idx];
    if(idx < size)
    {
        dst[idx].x += src[p].x;
        dst[idx].y += src[p].y;
    }
}

template<class T, class T2>
__global__ void kernel_addpp(T2 *ps, T *deeq, const T2 *becp, int nproj, int nprojx, int sum, int m, int nkb)
{
    int ip2 = blockDim.x * blockIdx.x + threadIdx.x;
    int ib = blockDim.y * blockIdx.y + threadIdx.y;
    if(ip2<nproj && ib<m)
    {
        ps[(sum+ip2) * m + ib].x = ps[(sum+ip2) * m + ib].y = 0;

        for(int ip=0; ip<nproj; ip++)
        {
            ps[(sum+ip2) * m + ib].x += deeq[ip * nprojx + ip2] * becp[ib * nkb + sum + ip].x;
            ps[(sum+ip2) * m + ib].y += deeq[ip * nprojx + ip2] * becp[ib * nkb + sum + ip].y;
        }
        // __syncthreads();
    }
}

template<class T>
void print_test(T *data, int size)
{
    T *h_data = (T*)malloc(size * sizeof(T));
    CHECK_CUDA(cudaMemcpy(h_data, data, size*sizeof(T), cudaMemcpyDeviceToHost));
    cout<<sizeof(h_data[0])<<endl;
    for(int i=0;i<size;i++){
        cout<<h_data[i].x<<" "<<h_data[i].y<<endl;
    }
    delete [] h_data;
}

int Hamilt_PW::moved = 0;

Hamilt_PW::Hamilt_PW()
{
    // hpsi = nullptr;
    // spsi = nullptr;
    // GR_index = nullptr;
    // Bec = nullptr;
#ifdef __CUDA
	CHECK_CUBLAS(cublasCreate(&hpw_handle));
#endif
}

Hamilt_PW::~Hamilt_PW()
{
    // delete[] hpsi;
    // delete[] spsi;
    delete[] GR_index;
    // delete[] Bec;
#ifdef __CUDA
	CHECK_CUDA(cudaFree(GR_index_d));
	CHECK_CUBLAS(cublasDestroy(hpw_handle));
	CHECK_CUDA(cudaFree(d_becp));
	CHECK_CUDA(cudaFree(d_ps));
    // CHECK_CUDA(cudaFree(GR_index_d));
#endif
}