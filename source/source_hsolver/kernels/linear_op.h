#ifndef HSOLVER_KERNELS_LINEAR_OP_H
#define HSOLVER_KERNELS_LINEAR_OP_H

#include "source_base/module_device/types.h"

namespace hsolver
{

/** @brief Elementwise operations; pointers and dot outputs belong to Device. */
template <typename T, typename Device>
class linear_op
{
  public:

    /** @brief Apply the same inverse diagonal to every column. */
    void diagonal(int ld, int dim, int nvec, const T* inverse, const T* x, T* y) const;
    /** @brief Fused CGS search direction, alpha step, and final updates. */
    void cgs_direction(int ld, int dim, int nvec, T* p, T* u, const T* r, const T* q, const T* beta) const;
    void cgs_alpha(int ld, int dim, int nvec, T* q, T* direction, const T* u, const T* v, const T* alpha) const;
    void cgs_finish(int ld, int dim, int nvec, T* x, T* r, const T* direction, const T* ad, const T* alpha) const;

    /** @brief Batched updates with optional device coefficients and skip mask. */
    void batch(int ld, int dim, int nvec, T* out, const T* x, const T* y,
               T a, T b, const T* ca, const T* cb, const int* skip) const;
    /** @brief One or two inner products using caller-owned partial sums. */
    void dots(int ld, int dim, int nvec, const T* x, const T* y, const T* z, const T* w,
              T* out, T* partial, int tiles) const;
    /** @brief Apply an ordered list of column swaps to all vector slots. */
    void swaps(int ld,int dim,int slots,int stride,int count,T* vectors,const int* pairs) const;
    /** @brief Fuse the two BiCGSTAB direction or solution updates. */
    void bicg_update(int ld,int dim,int nvec,bool direction,T* out,const T* x,const T* y,const T* a,const T* b) const;
    /** @brief Gather columns using a device destination-to-source map. */
    void gather(int ld, int dim, int nvec, int slots, int stride, const T* in, T* out, const int* map) const;
    void combine(const int n, T* out, const T* x, const T* y, const T a, const T b) const;
    void product(const int n, T* out, const T* x, const T* y) const;
    void swap(const int n, T* x, T* y) const;
    void dot(const int ld, const int dim, const int nvec, const T* x, const T* y, T* out) const;
};

#if defined(__CUDA) || defined(__ROCM)
template <typename T>
class linear_op<T, base_device::DEVICE_GPU>
{
  public:

    /** @brief Apply the same inverse diagonal to every column. */
    void diagonal(int ld, int dim, int nvec, const T* inverse, const T* x, T* y) const;
    /** @brief Fused CGS search direction, alpha step, and final updates. */
    void cgs_direction(int ld, int dim, int nvec, T* p, T* u, const T* r, const T* q, const T* beta) const;
    void cgs_alpha(int ld, int dim, int nvec, T* q, T* direction, const T* u, const T* v, const T* alpha) const;
    void cgs_finish(int ld, int dim, int nvec, T* x, T* r, const T* direction, const T* ad, const T* alpha) const;

    /** @brief Batched updates with optional device coefficients and skip mask. */
    void batch(int ld, int dim, int nvec, T* out, const T* x, const T* y,
               T a, T b, const T* ca, const T* cb, const int* skip) const;
    /** @brief One or two inner products using caller-owned partial sums. */
    void dots(int ld, int dim, int nvec, const T* x, const T* y, const T* z, const T* w,
              T* out, T* partial, int tiles) const;
    /** @brief Apply an ordered list of column swaps to all vector slots. */
    void swaps(int ld,int dim,int slots,int stride,int count,T* vectors,const int* pairs) const;
    /** @brief Fuse the two BiCGSTAB direction or solution updates. */
    void bicg_update(int ld,int dim,int nvec,bool direction,T* out,const T* x,const T* y,const T* a,const T* b) const;
    /** @brief Gather columns using a device destination-to-source map. */
    void gather(int ld, int dim, int nvec, int slots, int stride, const T* in, T* out, const int* map) const;
    void combine(const int n, T* out, const T* x, const T* y, const T a, const T b) const;
    void product(const int n, T* out, const T* x, const T* y) const;
    void swap(const int n, T* x, T* y) const;
    void dot(const int ld, const int dim, const int nvec, const T* x, const T* y, T* out) const;
};
#endif

} // namespace hsolver
#endif
