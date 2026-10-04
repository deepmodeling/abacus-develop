#include "source_hsolver/kernels/linear_op.h"

#include <cuda_runtime.h>
#include <thrust/complex.h>
#include <complex>
#include <stdexcept>

namespace hsolver
{
namespace
{
constexpr int linear_threads = 256;

void check_launch()
{
    const cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(error));
    }
}

template <typename Real>
__global__ void combine_kernel(const int n, thrust::complex<Real>* out,
                               const thrust::complex<Real>* x, const thrust::complex<Real>* y,
                               const thrust::complex<Real> a, const thrust::complex<Real> b)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n)
    {
        out[i] = a * x[i] + b * y[i];
    }
}

template <typename Real>
__global__ void product_kernel(const int n, thrust::complex<Real>* out,
                               const thrust::complex<Real>* x, const thrust::complex<Real>* y)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n)
    {
        out[i] = x[i] * y[i];
    }
}

template <typename Real>
__global__ void swap_kernel(const int n, thrust::complex<Real>* x, thrust::complex<Real>* y)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n)
    {
        const thrust::complex<Real> tmp = x[i];
        x[i] = y[i];
        y[i] = tmp;
    }
}

template <typename Real>
__global__ void dot_kernel(const int ld, const int dim, const thrust::complex<Real>* x,
                           const thrust::complex<Real>* y, thrust::complex<Real>* out)
{
    __shared__ Real re[linear_threads];
    __shared__ Real im[linear_threads];
    const int tid = threadIdx.x;
    const int offset = blockIdx.x * ld;
    thrust::complex<Real> sum(0, 0);
    for (int i = tid; i < dim; i += blockDim.x)
    {
        sum += thrust::conj(x[offset + i]) * y[offset + i];
    }
    re[tid] = sum.real();
    im[tid] = sum.imag();
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride /= 2)
    {
        if (tid < stride)
        {
            re[tid] += re[tid + stride];
            im[tid] += im[tid + stride];
        }
        __syncthreads();
    }
    if (tid == 0)
    {
        out[blockIdx.x] = thrust::complex<Real>(re[0], im[0]);
    }
}
} // namespace

template <typename T>
void linear_op<T, base_device::DEVICE_GPU>::combine(const int n, T* out, const T* x, const T* y, const T a, const T b) const
{
    if (n <= 0) { return; }
    using Real = typename T::value_type;
    using Complex = thrust::complex<Real>;
    combine_kernel<<<(n - 1) / linear_threads + 1, linear_threads>>>(
        n, reinterpret_cast<Complex*>(out), reinterpret_cast<const Complex*>(x),
        reinterpret_cast<const Complex*>(y), Complex(a.real(), a.imag()), Complex(b.real(), b.imag()));
    check_launch();
}

template <typename T>
void linear_op<T, base_device::DEVICE_GPU>::product(const int n, T* out, const T* x, const T* y) const
{
    if (n <= 0) { return; }
    using Complex = thrust::complex<typename T::value_type>;
    product_kernel<<<(n - 1) / linear_threads + 1, linear_threads>>>(
        n, reinterpret_cast<Complex*>(out), reinterpret_cast<const Complex*>(x), reinterpret_cast<const Complex*>(y));
    check_launch();
}

template <typename T>
void linear_op<T, base_device::DEVICE_GPU>::swap(const int n, T* x, T* y) const
{
    if (n <= 0) { return; }
    using Complex = thrust::complex<typename T::value_type>;
    swap_kernel<<<(n - 1) / linear_threads + 1, linear_threads>>>(
        n, reinterpret_cast<Complex*>(x), reinterpret_cast<Complex*>(y));
    check_launch();
}

template <typename T>
void linear_op<T, base_device::DEVICE_GPU>::dot(const int ld, const int dim, const int nvec,
                                              const T* x, const T* y, T* out) const
{
    if (nvec <= 0) { return; }
    using Complex = thrust::complex<typename T::value_type>;
    dot_kernel<<<nvec, linear_threads>>>(ld, dim, reinterpret_cast<const Complex*>(x),
        reinterpret_cast<const Complex*>(y), reinterpret_cast<Complex*>(out));
    check_launch();
}


template <typename C>
__global__ void batch_kernel(int ld,int dim,int nvec,C* out,const C* x,const C* y,C a,C b,const C* ca,const C* cb,const int* skip)
{
    const int i=blockIdx.x*blockDim.x+threadIdx.x;
    const int band=blockIdx.y;
    if(i<dim && (!skip || !skip[band]))
    {
        const int j=band*ld+i;
        out[j]=(ca ? a*ca[band] : a)*x[j]+(cb ? b*cb[band] : b)*y[j];
    }
}
template <typename Real>
__global__ void dots_kernel(int ld,int dim,int nvec,int tiles,
    const thrust::complex<Real>* x,const thrust::complex<Real>* y,
    const thrust::complex<Real>* z,const thrust::complex<Real>* w,thrust::complex<Real>* partial)
{
    __shared__ Real re[256];
    __shared__ Real im[256];
    const int tid=threadIdx.x;
    const int band=blockIdx.x;
    const int pair=blockIdx.z;
    const thrust::complex<Real>* a=pair ? z : x;
    const thrust::complex<Real>* b=pair ? w : y;
    thrust::complex<Real> sum(0,0);
    for(int i=blockIdx.y*256+tid;i<dim;i+=tiles*256) sum+=thrust::conj(a[band*ld+i])*b[band*ld+i];
    re[tid]=sum.real(); im[tid]=sum.imag();
    __syncthreads();
    for(int n=128;n>0;n/=2) { if(tid<n){re[tid]+=re[tid+n];im[tid]+=im[tid+n];} __syncthreads(); }
    if(tid==0) partial[(pair*nvec+band)*tiles+blockIdx.y]=thrust::complex<Real>(re[0],im[0]);
}
template <typename C>
__global__ void sum_tiles(int n,int tiles,const C* partial,C* out)
{
    const int i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i<n) { C sum(0,0); for(int t=0;t<tiles;++t)sum+=partial[i*tiles+t]; out[i]=sum; }
}
template <typename C>
__global__ void gather_kernel(int ld,int dim,int stride,const C* in,C* out,const int* map)
{
    const int i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i<dim) out[blockIdx.z*stride+blockIdx.y*ld+i]=in[blockIdx.z*stride+map[blockIdx.y]*ld+i];
}
template <typename T>
void linear_op<T,base_device::DEVICE_GPU>::batch(int ld,int dim,int nvec,T* out,const T* x,const T* y,
    T a,T b,const T* ca,const T* cb,const int* skip) const
{
    if(dim<=0 || nvec<=0)return;
    using C=thrust::complex<typename T::value_type>;
    batch_kernel<<<dim3((dim+255)/256,nvec),256>>>(ld,dim,nvec,reinterpret_cast<C*>(out),
        reinterpret_cast<const C*>(x),reinterpret_cast<const C*>(y),C(a.real(),a.imag()),C(b.real(),b.imag()),
        reinterpret_cast<const C*>(ca),reinterpret_cast<const C*>(cb),skip);
    check_launch();
}
template <typename T>
void linear_op<T,base_device::DEVICE_GPU>::dots(int ld,int dim,int nvec,const T* x,const T* y,const T* z,const T* w,
    T* out,T* partial,int tiles) const
{
    if(nvec<=0)return;
    using C=thrust::complex<typename T::value_type>;
    const int pairs=z ? 2 : 1;
    C* target=reinterpret_cast<C*>(tiles==1 ? out : partial);
    dots_kernel<<<dim3(nvec,tiles,pairs),256>>>(ld,dim,nvec,tiles,reinterpret_cast<const C*>(x),
        reinterpret_cast<const C*>(y),reinterpret_cast<const C*>(z),reinterpret_cast<const C*>(w),target);
    check_launch();
    if(tiles>1)
    {
        sum_tiles<<<(pairs*nvec+255)/256,256>>>(pairs*nvec,tiles,target,reinterpret_cast<C*>(out));
        check_launch();
    }
}
template <typename T>
void linear_op<T,base_device::DEVICE_GPU>::gather(int ld,int dim,int nvec,int slots,int stride,const T* in,T* out,const int* map) const
{
    if(dim<=0 || nvec<=0)return;
    using C=thrust::complex<typename T::value_type>;
    gather_kernel<<<dim3((dim+255)/256,nvec,slots),256>>>(ld,dim,stride,reinterpret_cast<const C*>(in),reinterpret_cast<C*>(out),map);
    check_launch();
}


template <typename C>
__global__ void swaps_kernel(int ld,int dim,int stride,int count,C* vectors,const int* pairs)
{
    const int i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i>=dim)return;
    const int offset=blockIdx.y*stride+i;
    for(int k=0;k<count;++k)
    {
        C& a=vectors[offset+pairs[2*k]*ld];
        C& b=vectors[offset+pairs[2*k+1]*ld];
        const C value=a; a=b; b=value;
    }
}
template <typename T>
void linear_op<T,base_device::DEVICE_GPU>::swaps(int ld,int dim,int slots,int stride,int count,T* vectors,const int* pairs) const
{
    if(dim==0 || count==0)return;
    using C=thrust::complex<typename T::value_type>;
    swaps_kernel<<<dim3((dim+255)/256,slots),256>>>(ld,dim,stride,count,reinterpret_cast<C*>(vectors),pairs);
    check_launch();
}


template <typename C>
__global__ void bicg_update_kernel(int ld,int dim,bool direction,C* out,const C* x,const C* y,const C* a,const C* b)
{
    const int i=blockIdx.x*blockDim.x+threadIdx.x;
    const int band=blockIdx.y;
    if(i<dim)
    {
        const int j=band*ld+i;
        out[j]=direction ? x[j]+a[band]*(out[j]-b[band]*y[j]) : (out[j]+a[band]*x[j])+b[band]*y[j];
    }
}
template <typename T>
void linear_op<T,base_device::DEVICE_GPU>::bicg_update(int ld,int dim,int nvec,bool direction,T* out,const T* x,const T* y,const T* a,const T* b) const
{
    if(dim==0 || nvec==0)return;
    using C=thrust::complex<typename T::value_type>;
    bicg_update_kernel<<<dim3((dim+255)/256,nvec),256>>>(ld,dim,direction,reinterpret_cast<C*>(out),
        reinterpret_cast<const C*>(x),reinterpret_cast<const C*>(y),reinterpret_cast<const C*>(a),reinterpret_cast<const C*>(b));
    check_launch();
}


template <typename C>
__global__ void diagonal_kernel(int ld, int dim, const C* inverse, const C* x, C* y)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < dim) y[blockIdx.y * ld + i] = inverse[i] * x[blockIdx.y * ld + i];
}

template <typename C>
__global__ void cgs_direction_kernel(int ld, int dim, C* p, C* u, const C* r, const C* q, const C* beta)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int band = blockIdx.y;
    if (i < dim)
    {
        const int j = band * ld + i;
        u[j] = r[j] + beta[band] * q[j];
        p[j] = u[j] + beta[band] * (q[j] + beta[band] * p[j]);
    }
}

template <typename C>
__global__ void cgs_alpha_kernel(int ld, int dim, C* q, C* direction, const C* u, const C* v, const C* alpha)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int band = blockIdx.y;
    if (i < dim)
    {
        const int j = band * ld + i;
        q[j] = u[j] - alpha[band] * v[j];
        direction[j] = u[j] + q[j];
    }
}

template <typename C>
__global__ void cgs_finish_kernel(int ld, int dim, C* x, C* r, const C* direction, const C* ad, const C* alpha)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int band = blockIdx.y;
    if (i < dim)
    {
        const int j = band * ld + i;
        x[j] += alpha[band] * direction[j];
        r[j] -= alpha[band] * ad[j];
    }
}

template <typename T>
void linear_op<T, base_device::DEVICE_GPU>::diagonal(int ld, int dim, int nvec, const T* inverse, const T* x, T* y) const
{
    if (dim == 0 || nvec == 0) return;
    using C = thrust::complex<typename T::value_type>;
    diagonal_kernel<<<dim3((dim + 255) / 256, nvec), 256>>>(ld, dim,
        reinterpret_cast<const C*>(inverse),
        reinterpret_cast<const C*>(x),
        reinterpret_cast<C*>(y));
    check_launch();
}

template <typename T>
void linear_op<T, base_device::DEVICE_GPU>::cgs_direction(int ld, int dim, int nvec, T* p, T* u, const T* r, const T* q, const T* beta) const
{
    if (dim == 0 || nvec == 0) return;
    using C = thrust::complex<typename T::value_type>;
    cgs_direction_kernel<<<dim3((dim + 255) / 256, nvec), 256>>>(ld, dim,
        reinterpret_cast<C*>(p),
        reinterpret_cast<C*>(u),
        reinterpret_cast<const C*>(r),
        reinterpret_cast<const C*>(q),
        reinterpret_cast<const C*>(beta));
    check_launch();
}

template <typename T>
void linear_op<T, base_device::DEVICE_GPU>::cgs_alpha(int ld, int dim, int nvec, T* q, T* direction, const T* u, const T* v, const T* alpha) const
{
    if (dim == 0 || nvec == 0) return;
    using C = thrust::complex<typename T::value_type>;
    cgs_alpha_kernel<<<dim3((dim + 255) / 256, nvec), 256>>>(ld, dim,
        reinterpret_cast<C*>(q),
        reinterpret_cast<C*>(direction),
        reinterpret_cast<const C*>(u),
        reinterpret_cast<const C*>(v),
        reinterpret_cast<const C*>(alpha));
    check_launch();
}

template <typename T>
void linear_op<T, base_device::DEVICE_GPU>::cgs_finish(int ld, int dim, int nvec, T* x, T* r, const T* direction, const T* ad, const T* alpha) const
{
    if (dim == 0 || nvec == 0) return;
    using C = thrust::complex<typename T::value_type>;
    cgs_finish_kernel<<<dim3((dim + 255) / 256, nvec), 256>>>(ld, dim,
        reinterpret_cast<C*>(x),
        reinterpret_cast<C*>(r),
        reinterpret_cast<const C*>(direction),
        reinterpret_cast<const C*>(ad),
        reinterpret_cast<const C*>(alpha));
    check_launch();
}
template class linear_op<std::complex<float>, base_device::DEVICE_GPU>;
template class linear_op<std::complex<double>, base_device::DEVICE_GPU>;

} // namespace hsolver
