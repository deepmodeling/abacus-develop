#ifndef HSOLVER_PPCG_KERNEL_OP_H
#define HSOLVER_PPCG_KERNEL_OP_H

#include "source_base/module_device/types.h"
#include <complex>

namespace hsolver
{

// Wavefunction operands reside on Device. Column products return only the
// reduced per-column scalars to host memory for the existing MPI wrappers.
template <typename T, typename Device>
struct PPCGVectorOps
{
    using Real = typename T::value_type;
    static void column_products(
        int dim,
        int ld,
        int count,
        const T *x,
        const T *y,
        double *host_products);
    static void residual(
        int dim,
        int ld,
        int count,
        const T *hx,
        const T *sx,
        const Real *eigenvalues,
        T *output);
    static void precondition(
        int dim,
        int ld,
        int count,
        const Real *diagonal,
        T *values);
    static void scale_columns(
        int dim,
        int ld,
        int count,
        const double *scales,
        T *x,
        T *sx,
        T *hx);
};

#if defined(__CUDA) || defined(__UT_USE_CUDA)
template <typename T>
struct PPCGVectorOps<T, base_device::DEVICE_GPU>
{
    using Real = typename T::value_type;
    static void column_products(
        int dim,
        int ld,
        int count,
        const T *x,
        const T *y,
        double *host_products);
    static void residual(
        int dim,
        int ld,
        int count,
        const T *hx,
        const T *sx,
        const Real *eigenvalues,
        T *output);
    static void precondition(
        int dim,
        int ld,
        int count,
        const Real *diagonal,
        T *values);
    static void scale_columns(
        int dim,
        int ld,
        int count,
        const double *scales,
        T *x,
        T *sx,
        T *hx);
};

#endif

} // namespace hsolver
#endif
