#include "source_hsolver/diago_dav_subspace.h"
#include "source_base/parallel_global.h"
#include "source_base/parallel_reduce.h"

#ifdef ABACUS_DAV_GPU_TEST
#include "source_base/kernels/math_kernel_op.h"
#include <ATen/kernels/blas.h>
#endif

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <type_traits>
#include <vector>

namespace
{
int test_rank = 0;
int test_nproc = 1;
#ifdef ABACUS_DAV_GPU_TEST
using TestDevice = base_device::DEVICE_GPU;
#else
using TestDevice = base_device::DEVICE_CPU;
#endif

template <typename T>
class DiagonalOperator : public hsolver::HSOperator<T, TestDevice>
{
  public:
    using Real = typename GetTypeReal<T>::type;
    explicit DiagonalOperator(const Real overlap) : overlap_(overlap) {}

    void update_k(const int) override {}

    void hpsi(const T* x, T* hx, const int ld, const int nvec) const override
    {
        std::vector<T> values(ld * nvec);
        base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, TestDevice>()(
            values.data(), x, values.size());
        for (int band = 0; band < nvec; ++band)
        {
            for (int i = 0; i < ld; ++i)
            {
                values[band * ld + i] *= overlap(i + test_rank * ld) * eigenvalue(i + test_rank * ld);
            }
        }
        base_device::memory::synchronize_memory_op<T, TestDevice, base_device::DEVICE_CPU>()(
            hx, values.data(), values.size());
    }

    void spsi(const T* x, T* sx, const int ld, const int nvec) const override
    {
        std::vector<T> values(ld * nvec);
        base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, TestDevice>()(
            values.data(), x, values.size());
        for (int i = 0; i < ld * nvec; ++i)
        {
            values[i] *= overlap(i % ld + test_rank * ld);
        }
        base_device::memory::synchronize_memory_op<T, TestDevice, base_device::DEVICE_CPU>()(
            sx, values.data(), values.size());
    }

    Real overlap(const int i) const { return overlap_ == Real(1) ? Real(1) : overlap_ * (Real(1) + Real(0.1) * i); }

    static Real eigenvalue(const int i) { return Real(i - 3); }

  private:
    Real overlap_;
};

template <typename T>
void solve_and_check(const int nbands,
                     const int dim,
                     const int ld,
                     std::vector<T>& psi,
                     const typename GetTypeReal<T>::type overlap,
                     const typename GetTypeReal<T>::type precondition_value,
                     const double tolerance,
                     const double check_tolerance)
{
    using Real = typename GetTypeReal<T>::type;
    std::vector<Real> precondition(dim, precondition_value);
    std::vector<Real> eigenvalues(nbands, 0);
    const std::vector<double> thresholds(nbands, tolerance);
#ifdef __MPI
    const hsolver::diag_comm_info comm(POOL_WORLD, test_rank, test_nproc);
#else
    const hsolver::diag_comm_info comm(test_rank, test_nproc);
#endif
    DiagonalOperator<T> op(overlap);
    hsolver::Diago_DavSubspace<T, TestDevice> solver(precondition, nbands, dim, 2, tolerance, 30, comm, 0, 0);
    T* device_psi = nullptr;
    base_device::memory::resize_memory_op<T, TestDevice>()(device_psi, psi.size());
    base_device::memory::synchronize_memory_op<T, TestDevice, base_device::DEVICE_CPU>()(
        device_psi, psi.data(), psi.size());
    solver.diag(op, device_psi, ld, eigenvalues.data(), thresholds, false, 6);
    base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, TestDevice>()(
        psi.data(), device_psi, psi.size());
    for (int band = 0; band < nbands; ++band)
    {
        EXPECT_NEAR(eigenvalues[band], DiagonalOperator<T>::eigenvalue(band), check_tolerance);
        double norm = 0;
        double residual = 0;
        for (int i = 0; i < dim; ++i)
        {
            EXPECT_TRUE(std::isfinite(std::abs(psi[band * ld + i])));
            const double weight = std::norm(std::complex<double>(psi[band * ld + i]));
            norm += op.overlap(i + test_rank * dim) * weight;
            const double delta = op.overlap(i + test_rank * dim) * (DiagonalOperator<T>::eigenvalue(i + test_rank * dim) - eigenvalues[band]);
            residual += delta * delta * weight;
        }
        Parallel_Reduce::reduce_pool(norm);
        Parallel_Reduce::reduce_pool(residual);
        EXPECT_NEAR(norm, 1.0, check_tolerance);
        EXPECT_LE(std::sqrt(residual), check_tolerance);
        for (int i = dim; i < ld; ++i)
        {
            EXPECT_EQ(psi[band * ld + i], T(0));
        }
    }
    // A second solve must accept the returned eigenvectors without a zero-direction breakdown.
    // Exercise the context-free overload used by LR and Python callers too.
    solver.diag(op, device_psi, ld, eigenvalues.data(), thresholds, false);
    base_device::memory::delete_memory_op<T, TestDevice>()(device_psi);
    for (int band = 0; band < nbands; ++band)
    {
        EXPECT_NEAR(eigenvalues[band], DiagonalOperator<T>::eigenvalue(band), check_tolerance);
    }
}

template <typename T>
class DavSubspaceTest : public testing::Test {};
#if defined(__LCAO) && !defined(ABACUS_DAV_GPU_TEST)
using DavTypes = testing::Types<std::complex<double>, std::complex<float>, double>;
#else
using DavTypes = testing::Types<std::complex<double>, std::complex<float>>;
#endif
TYPED_TEST_SUITE(DavSubspaceTest, DavTypes);

TYPED_TEST(DavSubspaceTest, ExactInitialRitzPairs)
{
    if (test_nproc > 1 && !std::is_same<TypeParam, std::complex<double>>::value) { GTEST_SKIP(); }
    const int dim = 10;
    const int ld = 12;
    std::vector<TypeParam> psi(2 * ld, 0);
    if (test_rank == 0)
    {
        // Rotated basis: retaining the initial Ritz coefficients is essential on early return.
        psi[0] = 0.6;
        psi[1] = 0.8;
        psi[ld] = -0.8;
        psi[ld + 1] = 0.6;
    }
    solve_and_check(2, dim, ld, psi, 1, 1, 1e-6, 2e-5);
}

TYPED_TEST(DavSubspaceTest, MixedConvergedBandsAndOverlap)
{
    if (test_nproc > 1 && !std::is_same<TypeParam, std::complex<double>>::value) { GTEST_SKIP(); }
    const int dim = 10;
    std::vector<TypeParam> psi(3 * dim, 0);
    if (test_rank == 0)
    {
        psi[0] = 0.5;
        psi[dim + 1] = 0.4;
        psi[dim + 3] = 0.3;
        psi[2 * dim + 2] = 0.5;
    }
    solve_and_check(3, dim, dim, psi, 4, 1, 1e-6, 2e-5);
}

TEST(DavSubspace, TinyPositiveResidual)
{
    const int dim = 6;
    std::vector<std::complex<double>> psi(dim, 0);
    if (test_rank == 0)
    {
        psi[0] = 1;
        psi[1] = 1e-14;
    }
    solve_and_check(1, dim, dim, psi, 1, 1, 1e-12, 1e-12);
}

TEST(DavSubspace, TinyResidualAboveRequestedToleranceStillExpands)
{
    const int dim = 6;
    std::vector<std::complex<double>> psi(dim, 0);
    if (test_rank == 0)
    {
        psi[0] = 1;
        psi[1] = 1e-13;
    }
    // The tiny-residual safeguard must not become a floor on the requested tolerance.
    solve_and_check(1, dim, dim, psi, 1, 1, 1e-14, 1e-14);
}

TEST(DavSubspace, SmallPreconditionedDirectionIsNotConvergence)
{
    const int dim = 6;
    std::vector<std::complex<double>> psi(dim, 0);
    if (test_rank == 0)
    {
        psi[0] = 0.6;
        psi[1] = 0.8;
    }
    solve_and_check(1, dim, dim, psi, 1, 1e14, 1e-12, 1e-11);
}

TEST(DavSubspace, LooseEnergyToleranceStillExpandsOrdinaryResidual)
{
    const int dim = 6;
    std::vector<std::complex<double>> psi(dim, 0);
    if (test_rank == 0)
    {
        psi[0] = std::sqrt(1.0 - 0.02 * 0.02);
        psi[1] = 0.02;
    }
    // Smooth band energy tolerances can be loose. Preserve the ordinary first
    // expansion instead of treating a non-negligible residual as a zero direction.
    solve_and_check(1, dim, dim, psi, 1, 1, 0.1, 1e-11);
}

TEST(DavSubspace, RepeatedSubspaceRefresh)
{
    const int dim = 10;
    std::vector<std::complex<double>> psi(dim, 0);
    if (test_rank == 0)
    {
        psi[0] = std::sqrt(0.5);
        psi[1] = std::complex<double>(0, 0.5);
        psi[2] = 0.5;
    }
    // One band and two workspace vectors require refresh to remove both excited components.
    solve_and_check(1, dim, dim, psi, 1, 1, 1e-14, 1e-6);
}

TEST(DavSubspace, NoncontiguousActiveBandsWithSpinorComponents)
{
    const int dim = 12;
    std::vector<std::complex<double>> psi(4 * dim, 0);
    if (test_rank == 0)
    {
        psi[0] = 1;
        psi[dim + 1] = std::sqrt(0.91);
        psi[dim + 6] = std::complex<double>(0, 0.3);
        psi[2 * dim + 2] = 1;
        psi[3 * dim + 3] = std::sqrt(0.91);
        psi[3 * dim + 7] = std::complex<double>(0.18, 0.24);
    }
    // Active bands 1 and 3 have weight in both halves of spinor-shaped storage.
    solve_and_check(4, dim, dim, psi, 1, 1, 1e-12, 1e-11);
}

#ifndef ABACUS_DAV_GPU_TEST
void check_breakdown(const double precondition)
{
    const int dim = 10;
    std::vector<std::complex<double>> psi(2 * dim, 0);
    psi[0] = 1;
    psi[dim + 1] = 0.8;
    psi[dim + 3] = 0.6;
    std::cout.rdbuf(std::cerr.rdbuf());
    solve_and_check(2, dim, dim, psi, 1, precondition, 1e-12, 1e-11);
}

TEST(DavSubspaceDeathTest, ZeroDirectionWithNonzeroResidual)
{
    if (test_nproc > 1) { GTEST_SKIP(); }
    EXPECT_EXIT(check_breakdown(std::numeric_limits<double>::infinity()),
                testing::ExitedWithCode(1),
                "unusable preconditioned direction; band=2, k-point=7.*residual_norm=.*residual_threshold=.*direction_norm=0");
}

TEST(DavSubspaceDeathTest, NonfiniteDirection)
{
    if (test_nproc > 1) { GTEST_SKIP(); }
    EXPECT_EXIT(check_breakdown(std::numeric_limits<double>::quiet_NaN()),
                testing::ExitedWithCode(1),
                "unusable preconditioned direction; band=2, k-point=7.*residual_norm=.*direction_norm=nan");
}
#endif
} // namespace

int main(int argc, char** argv)
{
    int nthread = 1;
    Parallel_Global::read_pal_param(argc, argv, test_nproc, nthread, test_rank);
#ifdef __MPI
    POOL_WORLD = MPI_COMM_WORLD;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
#endif
#ifdef ABACUS_DAV_GPU_TEST
    ModuleBase::createGpuBlasHandle();
    ct::kernels::createGpuBlasHandle();
    ct::kernels::createGpuSolverHandle();
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef ABACUS_DAV_GPU_TEST
    ct::kernels::destroyGpuSolverHandle();
    ct::kernels::destroyGpuBlasHandle();
    ModuleBase::destoryBLAShandle();
#endif
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif
    return result;
}
