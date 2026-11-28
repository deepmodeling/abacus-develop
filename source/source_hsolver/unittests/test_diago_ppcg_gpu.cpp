#include "../diago_ppcg.h"
#include "../hs_operator.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_basis/module_pw/test/test_tool.h"
#include <cmath>
#include <complex>
#include <cuda_runtime_api.h>
#include <gtest/gtest.h>
#include <limits>
#include <random>
#include <sstream>
#include <vector>

namespace
{
using GPU = base_device::DEVICE_GPU;
using CPU = base_device::DEVICE_CPU;

template <typename T>
class DenseDeviceOperator : public hsolver::HSOperator<T, GPU>
{
  public:
    DenseDeviceOperator(const int dimension, const std::vector<T> &h, const std::vector<T> &s) : dimension_(dimension)
    {
        h_.resize(h.size());
        s_.resize(s.size());
        base_device::memory::synchronize_memory_op<T, GPU, CPU>()(h_.data(), h.data(), h.size());
        base_device::memory::synchronize_memory_op<T, GPU, CPU>()(s_.data(), s.data(), s.size());
    }
    void update_k(const int) override {}
    void hpsi(
        const T *x,
        T *hx,
        const int ld,
        const int count) const override
    {
        multiply(h_.data(), x, hx, ld, count);
        ++h_calls;
    }
    void spsi(
        const T *x,
        T *sx,
        const int ld,
        const int count) const override
    {
        multiply(s_.data(), x, sx, ld, count);
        ++s_calls;
    }
    mutable int h_calls = 0;
    mutable int s_calls = 0;

  private:
    void multiply(
        const T *matrix,
        const T *x,
        T *output,
        const int ld,
        const int count) const
    {
        cudaPointerAttributes input_attributes{};
        cudaPointerAttributes output_attributes{};
        ASSERT_EQ(cudaPointerGetAttributes(&input_attributes, x), cudaSuccess);
        ASSERT_EQ(cudaPointerGetAttributes(&output_attributes, output), cudaSuccess);
        ASSERT_EQ(input_attributes.type, cudaMemoryTypeDevice);
        ASSERT_EQ(output_attributes.type, cudaMemoryTypeDevice);
        const int size = ld * count;
        base_device::memory::set_memory_op<T, GPU>()(output, 0, size);
        const T one(1);
        const T zero(0);
        ModuleBase::gemm_op<T, GPU>()(
            'N', 'N', dimension_, count, dimension_, &one, matrix, dimension_, x, ld, &zero, output, ld);
    }
    int dimension_;
    hsolver::PPCGDeviceBuffer<T> h_;
    hsolver::PPCGDeviceBuffer<T> s_;
};

// H = sqrt(S) Q Lambda Q^H sqrt(S), with a complex Fourier Q and
// nonuniform diagonal S. Its generalized eigenvalues are known independently.
template <typename Real>
void check_problem(
    const bool generalized,
    const int rr_step,
    const bool exact_prefix,
    const int maxiter)
{
    using T = std::complex<Real>;
    const int dimension = 24;
    const int bands = 5;
    const int ld = dimension + 7;
    const int matrix_size = dimension * dimension;
    const int vector_size = ld * bands;
    const double pi = 3.14159265358979323846;
    const double threshold = sizeof(Real) == sizeof(float) ? 2e-4 : 1e-8;
    const double check_tolerance = sizeof(Real) == sizeof(float) ? 4e-4 : 3e-8;
    std::vector<T> h(matrix_size, T(0));
    std::vector<T> s(matrix_size, T(0));
    std::vector<Real> overlap(dimension, Real(1));
    std::vector<Real> preconditioner(dimension, Real(1));
    for (int row = 0; row < dimension; ++row)
    {
        if (generalized)
        {
            overlap[row] = Real(0.5 + 0.08 * row);
        }
        s[row + row * dimension] = T(overlap[row]);
    }
    for (int column = 0; column < dimension; ++column)
    {
        for (int row = 0; row < dimension; ++row)
        {
            std::complex<double> sum(0);
            for (int mode = 0; mode < dimension; ++mode)
            {
                const double eigenvalue = 1.0 + 0.2 * mode;
                const double angle = 2.0 * pi * mode * (row - column) / dimension;
                sum += eigenvalue * std::polar(1.0, angle) / double(dimension);
            }
            const double scale = std::sqrt(double(overlap[row]) * double(overlap[column]));
            h[row + column * dimension] = T(sum * scale);
        }
        preconditioner[column] = std::real(h[column + column * dimension]);
    }
    std::vector<T> vectors(vector_size, T(123, -45));
    std::mt19937 generator(1701);
    std::uniform_real_distribution<double> random(-1.0, 1.0);
    for (int band = 0; band < bands; ++band)
    {
        for (int row = 0; row < dimension; ++row)
        {
            const Real real_part = Real(random(generator));
            const Real imaginary_part = Real(random(generator));
            vectors[row + band * ld] = T(real_part, imaginary_part);
        }
    }
    if (exact_prefix)
    {
        for (int row = 0; row < dimension; ++row)
        {
            vectors[row] = Real(1) / std::sqrt(Real(dimension) * overlap[row]);
        }
    }
    DenseDeviceOperator<T> op(dimension, h, s);
    hsolver::PPCGDeviceBuffer<T> device_vectors(vector_size, T(0));
    base_device::memory::synchronize_memory_op<T, GPU, CPU>()(device_vectors.data(), vectors.data(), vector_size);
    std::vector<Real> eigenvalues(bands, Real(0));
    const std::vector<double> thresholds(bands, threshold);
    const Real solver_threshold = Real(threshold);
    hsolver::DiagoPPCG<T, GPU> solver(solver_threshold, maxiter, 2, rr_step);
    const hsolver::PPCGTraceContext context = {1, 1};
    std::ostringstream log;
    solver.diag(
        op, ld, bands, dimension, device_vectors.data(), eigenvalues.data(), thresholds, preconditioner.data(), log,
        context);
    EXPECT_GT(op.h_calls, 0);
    EXPECT_GT(op.s_calls, 0);
    if (maxiter == 1)
    {
        EXPECT_FALSE(solver.converged());
        EXPECT_TRUE(solver.iteration_limit_reached());
        EXPECT_GT(solver.active_band_count(), 0);
        EXPECT_TRUE(std::isfinite(solver.max_residual()));
        return;
    }
    ASSERT_TRUE(solver.converged()) << solver.max_residual();
    base_device::memory::synchronize_memory_op<T, CPU, GPU>()(vectors.data(), device_vectors.data(), vector_size);
    for (int band = 0; band < bands; ++band)
    {
        const double expected_eigenvalue = 1.0 + 0.2 * band;
        EXPECT_NEAR(eigenvalues[band], expected_eigenvalue, check_tolerance);
        double residual_squared = 0.0;
        for (int row = 0; row < dimension; ++row)
        {
            std::complex<double> hx(0);
            for (int column = 0; column < dimension; ++column)
            {
                hx += std::complex<double>(h[row + column * dimension]) *
                      std::complex<double>(vectors[column + band * ld]);
            }
            const auto sx = double(overlap[row]) * std::complex<double>(vectors[row + band * ld]);
            residual_squared += std::norm(hx - double(eigenvalues[band]) * sx);
        }
        EXPECT_LT(std::sqrt(residual_squared), check_tolerance);
        for (int other = 0; other < bands; ++other)
        {
            std::complex<double> product(0);
            for (int row = 0; row < dimension; ++row)
            {
                product += std::conj(std::complex<double>(vectors[row + other * ld])) * double(overlap[row]) *
                           std::complex<double>(vectors[row + band * ld]);
            }
            const double expected_overlap = band == other ? 1.0 : 0.0;
            EXPECT_NEAR(product.real(), expected_overlap, check_tolerance);
            EXPECT_NEAR(product.imag(), 0.0, check_tolerance);
        }
        for (int row = dimension; row < ld; ++row)
        {
            EXPECT_EQ(vectors[row + band * ld], T(0));
        }
    }
}

TEST(DiagoPPCGGPU, StandardComplexDouble) { check_problem<double>(false, 1, false, 200); }
TEST(DiagoPPCGGPU, GeneralizedComplexDouble) { check_problem<double>(true, 1, false, 200); }
TEST(DiagoPPCGGPU, InfrequentRitzAndLockedPrefix) { check_problem<double>(true, 4, true, 200); }
TEST(DiagoPPCGGPU, StandardComplexFloat) { check_problem<float>(false, 1, false, 200); }
TEST(DiagoPPCGGPU, GeneralizedComplexFloat) { check_problem<float>(true, 4, false, 200); }
TEST(DiagoPPCGGPU, ReportsIterationLimit) { check_problem<double>(true, 1, false, 1); }
} // namespace

int main(int argc, char **argv)
{
#ifdef __MPI
    int nproc = 1;
    int rank = 0;
    setupmpi(argc, argv, nproc, rank);
    int nproc_in_pool = 1;
    int kpar = nproc;
    int pool = 0;
    int rank_in_pool = 0;
    divide_pools(nproc, rank, nproc_in_pool, kpar, pool, rank_in_pool);
#endif
    ModuleBase::createGpuBlasHandle();
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    ModuleBase::destoryBLAShandle();
#ifdef __MPI
    finishmpi();
#endif
    return result;
}
