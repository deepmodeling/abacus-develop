#include "source_hsolver/module_linalg/linear_algebra.h"
#include "source_base/parallel_cell.h"
#include "source_base/parallel_device.h"
#include "source_base/kernels/math_kernel_op.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <limits>
#include <string>
#include <vector>

namespace
{
using Wide = std::complex<double>;
#if defined(__CUDA) || defined(__ROCM)
using Device = base_device::DEVICE_GPU;
#else
using Device = base_device::DEVICE_CPU;
#endif

hsolver::diag_comm_info local_comm()
{
#ifdef __MPI
    return hsolver::diag_comm_info(MPI_COMM_SELF, 0, 1);
#else
    return hsolver::diag_comm_info(0, 1);
#endif
}

template <typename T, typename Device>
ct::Tensor upload(const std::vector<T>& values)
{
    const ct::DataType type = ct::DataTypeToEnum<T>::value;
    const ct::DeviceType device = ct::DeviceTypeToEnum<Device>::value;
    const int count = values.size();
    ct::Tensor result(type, device, {count});
    base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>()(
        result.data<T>(), values.data(), count);
    return result;
}

template <typename T, typename Device>
void check_output(const ct::Tensor& output, const std::vector<T>& expected)
{
    std::vector<T> actual(expected.size());
    base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, Device>()(
        actual.data(), output.template data<T>(), actual.size());
    const double tolerance = sizeof(T) == sizeof(std::complex<float>) ? 2e-6 : 1e-12;
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        EXPECT_LT(std::abs(actual[i] - expected[i]), tolerance);
    }
}

template <typename T>
Wide overlap(const int rows, const T* x, const T* y)
{
    Wide result(0);
    for (int k = 0; k < rows; ++k) { result += std::conj(Wide(x[k])) * Wide(y[k]); }
    return result;
}

template <typename T, typename Device>
void check_prepared_coefficients()
{
    const int dim = 5;
    const int ld = 7;
    const int nx = 3;
    const int ny = 2;
    const hsolver::diag_comm_info comm = local_comm();
    hsolver::LinearAlgebra<T, Device> algebra(comm);
    std::vector<Wide> coefficients = {{0.7, -0.2}, {0.1, 0.3}, {-0.4, 0.5},
                                      {0.2, 0.4}, {-0.3, 0.1}, {0.6, -0.2}};
    std::vector<Wide> transient(coefficients);
    const int64_t count = transient.size();
    const T* prepared = algebra.prepare_coefficients(transient.data(), count);
    // Exercise coefficient ownership by releasing the source before rotation.
    std::fill(transient.begin(), transient.end(), Wide(999, -999));
    transient.clear();
    transient.shrink_to_fit();
    for (int image = 0; image < 2; ++image)
    {
        std::vector<T> x(ld * nx, T(77, 9));
        std::vector<T> expected(ld * ny, T(-3, 4));
        for (int j = 0; j < nx; ++j)
        {
            for (int row = 0; row < dim; ++row)
            {
                x[row + j * ld] = T(0.1 * (row + image + 1), 0.2 * (j - row));
            }
        }
        ct::Tensor input = upload<T, Device>(x);
        ct::Tensor output = upload<T, Device>(expected);
        const T beta(0.2, -0.1);
        algebra.expand_prepared(ld, dim, nx, ny, input.template data<T>(), prepared, output.template data<T>(), beta);
        for (int j = 0; j < ny; ++j)
        {
            for (int row = 0; row < dim; ++row)
            {
                T value = beta * expected[row + j * ld];
                for (int k = 0; k < nx; ++k)
                {
                    value += x[row + k * ld] * T(coefficients[k + j * nx]);
                }
                expected[row + j * ld] = value;
            }
        }
        check_output<T, Device>(output, expected);
        // The prepared coefficients remain usable across zero-sized operations.
        algebra.expand_prepared(ld, 0, nx, ny, nullptr, nullptr, nullptr, beta);
    }
    EXPECT_EQ(algebra.prepare_coefficients(nullptr, 0), nullptr);
}

TEST(LinearAlgebra, ReusePreparedCoefficients)
{
    check_prepared_coefficients<std::complex<double>, Device>();
    check_prepared_coefficients<std::complex<float>, Device>();
}

template <typename T, typename Device>
void check_padded_products()
{
    const int dim = 5;
    const int ld = 8;
    const int nx = 3;
    const int ny = 2;
    const hsolver::diag_comm_info comm = local_comm();
    hsolver::LinearAlgebra<T, Device> algebra(comm);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const T padding(nan, nan);
    std::vector<T> x(ld * nx, padding);
    std::vector<T> y(ld * ny, padding);
    for (int j = 0; j < nx; ++j)
    {
        for (int i = 0; i < dim; ++i)
        {
            x[i + j * ld] = T(0.2 * (i + 1), 0.3 * (j - i));
        }
    }
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < dim; ++i)
        {
            y[i + j * ld] = T(0.3 * (j + i), -0.2 * (i + 2));
        }
    }
    ct::Tensor dx = upload<T, Device>(x);
    ct::Tensor dy = upload<T, Device>(y);
    const ct::Tensor gram_tensor = algebra.gram_tensor(ld, dim, nx, dx.template data<T>());
    const ct::Tensor cross_tensor = algebra.cross_tensor(ld, dim, nx, ny, dx.template data<T>(), dy.template data<T>());
    ASSERT_EQ(gram_tensor.device_type(), ct::DeviceType::CpuDevice);
    ASSERT_EQ(cross_tensor.data_type(), ct::DataType::DT_COMPLEX_DOUBLE);
    for (int j = 0; j < nx; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            const Wide expected = overlap(dim, x.data() + i * ld, x.data() + j * ld);
            EXPECT_LT(std::abs(gram_tensor.template data<Wide>()[i + j * nx] - expected), 1e-12);
        }
    }
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            const Wide expected = overlap(dim, x.data() + i * ld, y.data() + j * ld);
            EXPECT_LT(std::abs(cross_tensor.template data<Wide>()[i + j * nx] - expected), 1e-12);
        }
    }
    // Empty local partitions contribute a zero matrix with nx*ny entries.
    const ct::Tensor empty = algebra.cross_tensor(ld, 0, nx, ny, nullptr, nullptr);
    ASSERT_EQ(empty.NumElements(), nx * ny);
    for (int i = 0; i < nx * ny; ++i)
    {
        EXPECT_EQ(empty.template data<Wide>()[i], Wide(0));
    }
}

TEST(LinearAlgebra, PaddedGramAndCross)
{
    check_padded_products<std::complex<double>, Device>();
    check_padded_products<std::complex<float>, Device>();
}

TEST(LinearAlgebra, RootProjectionWithEmptyLocalRows)
{
    const ModuleBase::CommunicationDomain domain = ModuleBase::world_comm_domain();
    const int rank = domain.rank();
    const int ranks = domain.size();
#ifdef __MPI
    const MPI_Comm communicator = domain.communicator();
    const hsolver::diag_comm_info comm(communicator, rank, ranks);
#else
    const hsolver::diag_comm_info comm(0, 1);
#endif
    const int rows = rank == ranks - 1 ? 3 : 0;
    const int ld = 5;
    const int columns = 2;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<Wide> x(ld * columns, Wide(nan, nan));
    for (int j = 0; j < columns; ++j)
    {
        for (int i = 0; i < rows; ++i)
        {
            x[i + j * ld] = Wide(i + 1, j + 1);
        }
    }
    hsolver::LinearAlgebra<Wide, Device> algebra(comm);
    const ct::Tensor input = upload<Wide, Device>(x);
    const Wide* data = input.data<Wide>();
    const ct::Tensor cross = algebra.cross_root_tensor(ld, rows, columns, columns, data, data);
    const ct::Tensor gram = algebra.gram_root_tensor(ld, rows, columns, data);
    if (rank != 0)
    {
        EXPECT_EQ(cross.NumElements(), 0);
        EXPECT_EQ(gram.NumElements(), 0);
        EXPECT_EQ(cross.buffer().GetAllocatedBytes(), 0U);
        EXPECT_EQ(gram.buffer().GetAllocatedBytes(), 0U);
        return;
    }
    ASSERT_EQ(cross.NumElements(), 4);
    ASSERT_EQ(gram.NumElements(), 4);
    for (int j = 0; j < columns; ++j)
    {
        for (int i = 0; i < columns; ++i)
        {
            Wide expected(0);
            for (int k = 0; k < 3; ++k)
            {
                expected += std::conj(Wide(k + 1, i + 1)) * Wide(k + 1, j + 1);
            }
            EXPECT_NEAR(std::abs(cross.data<Wide>()[i + j * columns] - expected), 0.0, 1e-12);
            EXPECT_NEAR(std::abs(gram.data<Wide>()[i + j * columns] - expected), 0.0, 1e-12);
        }
    }
}

#ifdef __MPI
TEST(LinearAlgebra, RootReductionSupportsOtherRootAndSelfCommunicator)
{
    const ModuleBase::CommunicationDomain domain = ModuleBase::world_comm_domain();
    const MPI_Comm communicator = domain.communicator();
    const int rank = domain.rank();
    const int size = domain.size();
    const int root = size - 1;
    Wide value(rank + 1, -rank);
    const Wide before = value;
    Parallel_Common::reduce_data_to_root(&value, 1, communicator, root);
    if (rank == root)
    {
        const Wide expected(size * (size + 1) / 2.0, -size * (size - 1) / 2.0);
        EXPECT_EQ(value, expected);
    }
    else
    {
        EXPECT_EQ(value, before);
    }
    const hsolver::diag_comm_info self = local_comm();
    value = before;
    Parallel_Common::reduce_data_to_root(&value, 1, self.comm, 0);
    EXPECT_EQ(value, before);
    Parallel_Common::reduce_data_to_root(&value, 0, communicator, root);
    EXPECT_EQ(value, before);
}
#endif

#if !defined(__CUDA) && !defined(__ROCM)
TEST(LinearAlgebra, HerkUpperTriangle)
{
    const int rows = 3;
    const int ld = 5;
    const int bands = 2;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<Wide> x(ld * bands, Wide(nan, nan));
    for (int j = 0; j < bands; ++j)
    {
        for (int i = 0; i < rows; ++i) { x[i + j * ld] = Wide(i + 1, j - i); }
    }
    std::vector<Wide> gram(bands * bands, Wide(nan, nan));
    ModuleBase::gram_herk_cpu(ld, rows, bands, x.data(), gram.data());
    EXPECT_TRUE(std::isnan(gram[1].real()));
    EXPECT_EQ(gram[0].imag(), 0.0);
    EXPECT_EQ(gram[3].imag(), 0.0);
    hsolver::linear_complete_hermitian_upper(bands, gram.data());
    for (int j = 0; j < bands; ++j)
    {
        for (int i = 0; i < bands; ++i)
        {
            const Wide expected = overlap(rows, x.data() + i * ld, x.data() + j * ld);
            EXPECT_LT(std::abs(gram[i + j * bands] - expected), 1e-12);
        }
    }
}
#endif

template <typename T, typename Device>
void verify_triangular_apply(const char side, const char uplo, const char trans, const char diag, const T alpha)
{
    SCOPED_TRACE(::testing::Message() << side << uplo << trans << diag << " alpha=" << alpha);
    const int m = 7;
    const int n = 4;
    const int ldb = 10;
    const int order = side == 'L' ? m : n;
    const int lda = order + 2;
    const T poison(73, -41);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<T> triangle(lda * order, T(nan, nan));
    std::vector<Wide> dense(order * order, Wide(0));
    for (int j = 0; j < order; ++j)
    {
        for (int i = 0; i < order; ++i)
        {
            if (diag == 'U' && i == j)
            {
                dense[i + j * order] = 1; // The stored unit diagonal stays NaN.
            }
            else if ((uplo == 'U' && i <= j) || (uplo == 'L' && i >= j))
            {
                const T value(0.1 * (i + 1), 0.05 * (j - i));
                triangle[i + j * lda] = value;
                dense[i + j * order] = Wide(value);
            }
        }
    }
    std::vector<T> input(ldb * n, poison);
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < m; ++i)
        {
            input[i + j * ldb] = T(0.1 * (i + j + 1), 0.07 * (i - j));
        }
    }
    std::vector<T> expected(input);
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < m; ++i)
        {
            Wide sum(0);
            for (int k = 0; k < order; ++k)
            {
                const int row = side == 'L' ? i : k;
                const int column = side == 'L' ? k : j;
                const Wide value = trans == 'N' ? dense[row + column * order] : dense[column + row * order];
                const Wide factor = trans == 'C' ? std::conj(value) : value;
                sum += side == 'L' ? factor * Wide(input[k + j * ldb])
                                   : Wide(input[i + k * ldb]) * factor;
            }
            expected[i + j * ldb] = T(Wide(alpha) * sum);
        }
    }
    ct::Tensor a = upload<T, Device>(triangle);
    ct::Tensor b = upload<T, Device>(input);
    ModuleBase::trmm_op<T, Device>()(side, uplo, trans, diag, m, n, &alpha,
                                   a.template data<T>(), lda, b.template data<T>(), ldb);
    check_output<T, Device>(b, expected);
}


TEST(LinearAlgebra, TriangularProducts)
{
    // side/uplo/trans/diag: exercise both sides and triangles, transpose and unit diagonal.
    for (const std::string mode : {"RUNN", "RUCN", "RUTU", "RLNN", "LUNU", "LLCN"})
    {
        const char side = mode[0];
        const char uplo = mode[1];
        const char trans = mode[2];
        const char diag = mode[3];
        const Wide alpha(.7, -.2);
        const std::complex<float> alpha_float(alpha);
        verify_triangular_apply<Wide, Device>(side, uplo, trans, diag, alpha);
        verify_triangular_apply<std::complex<float>, Device>(side, uplo, trans, diag, alpha_float);
    }
}

TEST(LinearAlgebra, TriangularZeroScaleAndEmptyMatrices)
{
    const Wide zero(0);
    const std::complex<float> zero_float(0);
    verify_triangular_apply<Wide, Device>('R', 'U', 'N', 'N', zero);
    verify_triangular_apply<std::complex<float>, Device>('R', 'U', 'N', 'N', zero_float);
    ModuleBase::trmm_op<Wide, Device>()('R', 'U', 'N', 'N', 0, 4, nullptr, nullptr, 4, nullptr, 10);
    ModuleBase::trmm_op<Wide, Device>()('R', 'U', 'N', 'N', 7, 0, nullptr, nullptr, 1, nullptr, 10);
}
} // namespace
