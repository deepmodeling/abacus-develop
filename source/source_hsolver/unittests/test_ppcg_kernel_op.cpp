#include "../kernels/ppcg_kernel_op.h"
#include "../ppcg_workspace.h"

#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <limits>

namespace
{
#ifdef __CUDA
using ResidualDevice = base_device::DEVICE_GPU;
#else
using ResidualDevice = base_device::DEVICE_CPU;
#endif

template <typename Real>
void check_residual_cancellation()
{
    using T = std::complex<Real>;
    using Workspace = typename hsolver::PPCGWorkspaceType<T, ResidualDevice>::type;
    using Values = typename hsolver::PPCGWorkspaceType<Real, ResidualDevice>::type;
    const int exponent = -(std::numeric_limits<Real>::digits / 2 + 1);
    const Real delta = std::ldexp(Real(1), exponent);
    const Real eigenvalue = Real(1) + delta;
    const Real imaginary_h = Real(1) + Real(2) * delta;
    const Real real_s = Real(1) - delta;
    const Real imaginary_s = Real(1) + delta;
    const T h[2] = {T(1, imaginary_h), T(0)};
    const T s[2] = {T(real_s, imaginary_s), T(0)};
    const T padding(123, -45);
    Workspace hx(2, T(0));
    Workspace sx(2, T(0));
    Workspace residual(2, T(0));
    const T initial[2] = {padding, padding};
    Values values(1, Real(0));
    base_device::memory::synchronize_memory_op<T, ResidualDevice, base_device::DEVICE_CPU>()(hx.data(), h, 2);
    base_device::memory::synchronize_memory_op<T, ResidualDevice, base_device::DEVICE_CPU>()(sx.data(), s, 2);
    base_device::memory::synchronize_memory_op<T, ResidualDevice, base_device::DEVICE_CPU>()(
        residual.data(), initial, 2);
    base_device::memory::synchronize_memory_op<Real, ResidualDevice, base_device::DEVICE_CPU>()(
        values.data(), &eigenvalue, 1);
    hsolver::PPCGVectorOps<T, ResidualDevice>::residual(1, 2, 1, hx.data(), sx.data(), values.data(), residual.data());
    T result[2];
    base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, ResidualDevice>()(
        result, residual.data(), 2);
    // The exact residual is (+delta^2, -delta^2). Rounding the product
    // before subtraction would incorrectly report zero in both components.
    const Real square = delta * delta;
    EXPECT_EQ(result[0], T(square, -square));
    EXPECT_EQ(result[1], padding);
}

TEST(PPCGResidualKernel, PreservesDoublePrecisionCancellationRemainder) { check_residual_cancellation<double>(); }

TEST(PPCGResidualKernel, PreservesFloatPrecisionCancellationRemainder) { check_residual_cancellation<float>(); }
} // namespace
