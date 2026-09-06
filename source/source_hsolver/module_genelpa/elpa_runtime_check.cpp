#include "elpa_runtime_check.h"

#include "elpa_new.h"
#include "source_base/tool_quit.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <utility>
#include <vector>

namespace
{
using Complex = std::complex<double>;

// Keep this case aligned with the standalone ELPA BLOCK2 reproducer and the
// DiagoPxxxgvxElpaTest.ComplexDouble unit test that originally exposed it.
constexpr int kN = 20;
constexpr int kNev = 18;
constexpr int kNblk = 8;
constexpr int kSeed = 0;
constexpr double kPassThreshold = 1.0e-8;

std::size_t index(const int row, const int col, const int leading_dimension)
{
    return static_cast<std::size_t>(row) + static_cast<std::size_t>(col) * leading_dimension;
}

bool set_integer(elpa_t handle, const char* name, const int value, int& status)
{
    elpa_set_integer(handle, name, value, &status);
    return status == ELPA_OK;
}

std::pair<std::vector<Complex>, std::vector<Complex>> generate_hs()
{
    std::mt19937 generator(kSeed);
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);

    std::vector<Complex> h(kN * kN, Complex{0.0, 0.0});
    std::vector<Complex> s(kN * kN, Complex{0.0, 0.0});
    std::vector<Complex> s_tmp(kN * kN, Complex{0.0, 0.0});

    for (int i = 0; i < kN; ++i)
    {
        for (int j = i; j < kN; ++j)
        {
            const double h_real = distribution(generator);
            const double h_imag = distribution(generator);

            h[static_cast<std::size_t>(i) * kN + j] = {h_real, h_imag};
            h[static_cast<std::size_t>(j) * kN + i]
                = i == j ? Complex{h_real, 0.0} : Complex{h_real, -h_imag};

            const double s_real = distribution(generator);
            const double s_imag = distribution(generator);
            s_tmp[static_cast<std::size_t>(i) * kN + j] = {s_real, s_imag};
        }
    }

    for (int i = 0; i < kN; ++i)
    {
        for (int j = 0; j < kN; ++j)
        {
            Complex value{0.0, 0.0};
            for (int k = 0; k < kN; ++k)
            {
                value += s_tmp[static_cast<std::size_t>(i) * kN + k]
                         * std::conj(s_tmp[static_cast<std::size_t>(j) * kN + k]);
            }
            if (i == j)
            {
                value += Complex{2.0, 0.0};
            }
            s[static_cast<std::size_t>(i) * kN + j] = value;
        }
    }

    return {std::move(h), std::move(s)};
}

bool upper_cholesky(const std::vector<Complex>& matrix, std::vector<Complex>& upper)
{
    upper.assign(kN * kN, Complex{0.0, 0.0});

    for (int j = 0; j < kN; ++j)
    {
        double diagonal = matrix[index(j, j, kN)].real();
        for (int k = 0; k < j; ++k)
        {
            diagonal -= std::norm(upper[index(k, j, kN)]);
        }

        if (!(diagonal > 0.0))
        {
            return false;
        }

        upper[index(j, j, kN)] = {std::sqrt(diagonal), 0.0};

        for (int col = j + 1; col < kN; ++col)
        {
            Complex value = matrix[index(j, col, kN)];
            for (int k = 0; k < j; ++k)
            {
                value -= std::conj(upper[index(k, j, kN)]) * upper[index(k, col, kN)];
            }
            upper[index(j, col, kN)] = value / upper[index(j, j, kN)].real();
        }
    }

    return true;
}

std::vector<Complex> invert_upper_triangular(const std::vector<Complex>& upper)
{
    std::vector<Complex> inverse(kN * kN, Complex{0.0, 0.0});

    for (int col = 0; col < kN; ++col)
    {
        for (int row = kN - 1; row >= 0; --row)
        {
            Complex value = row == col ? Complex{1.0, 0.0} : Complex{0.0, 0.0};
            for (int k = row + 1; k < kN; ++k)
            {
                value -= upper[index(row, k, kN)] * inverse[index(k, col, kN)];
            }
            inverse[index(row, col, kN)] = value / upper[index(row, row, kN)];
        }
    }

    return inverse;
}

std::vector<Complex> multiply(const std::vector<Complex>& left, const std::vector<Complex>& right)
{
    std::vector<Complex> result(kN * kN, Complex{0.0, 0.0});
    for (int col = 0; col < kN; ++col)
    {
        for (int row = 0; row < kN; ++row)
        {
            Complex value{0.0, 0.0};
            for (int k = 0; k < kN; ++k)
            {
                value += left[index(row, k, kN)] * right[index(k, col, kN)];
            }
            result[index(row, col, kN)] = value;
        }
    }
    return result;
}

std::vector<Complex> conjugate_transpose(const std::vector<Complex>& matrix)
{
    std::vector<Complex> result(kN * kN);
    for (int col = 0; col < kN; ++col)
    {
        for (int row = 0; row < kN; ++row)
        {
            result[index(row, col, kN)] = std::conj(matrix[index(col, row, kN)]);
        }
    }
    return result;
}

bool generate_a_tilde(std::vector<Complex>& result)
{
    const std::pair<std::vector<Complex>, std::vector<Complex>> hs = generate_hs();
    const std::vector<Complex>& h_row_major = hs.first;
    const std::vector<Complex>& s_row_major = hs.second;

    std::vector<Complex> a(kN * kN);
    std::vector<Complex> s(kN * kN);
    for (int row = 0; row < kN; ++row)
    {
        for (int col = 0; col < kN; ++col)
        {
            a[index(row, col, kN)] = h_row_major[static_cast<std::size_t>(col) * kN + row];
            s[index(row, col, kN)] = s_row_major[static_cast<std::size_t>(col) * kN + row];
        }
    }

    std::vector<Complex> upper;
    if (!upper_cholesky(s, upper))
    {
        return false;
    }

    const auto upper_inverse = invert_upper_triangular(upper);
    const auto right_transformed = multiply(a, upper_inverse);
    result = multiply(conjugate_transpose(upper_inverse), right_transformed);
    return true;
}

double frobenius_norm(const std::vector<Complex>& matrix)
{
    double norm_squared = 0.0;
    for (const Complex value : matrix)
    {
        norm_squared += std::norm(value);
    }
    return std::sqrt(norm_squared);
}

double vector_norm(const Complex* vector, const int size)
{
    double norm_squared = 0.0;
    for (int i = 0; i < size; ++i)
    {
        norm_squared += std::norm(vector[i]);
    }
    return std::sqrt(norm_squared);
}

double relative_residual(const std::vector<Complex>& matrix,
                         const Complex* eigenvector,
                         const double eigenvalue,
                         const double matrix_norm)
{
    std::vector<Complex> residual(kN, Complex{0.0, 0.0});

    for (int row = 0; row < kN; ++row)
    {
        Complex matrix_vector{0.0, 0.0};
        for (int col = 0; col < kN; ++col)
        {
            matrix_vector += matrix[index(row, col, kN)] * eigenvector[col];
        }
        residual[row] = matrix_vector - eigenvalue * eigenvector[row];
    }

    const double eigenvector_norm = vector_norm(eigenvector, kN);
    const double denominator = (matrix_norm + std::abs(eigenvalue)) * eigenvector_norm;
    return denominator > 0.0 ? vector_norm(residual.data(), kN) / denominator
                             : vector_norm(residual.data(), kN);
}

double orthogonality_error(const std::vector<Complex>& eigenvectors)
{
    double error_squared = 0.0;

    for (int i = 0; i < kNev; ++i)
    {
        for (int j = 0; j < kNev; ++j)
        {
            Complex overlap{0.0, 0.0};
            for (int row = 0; row < kN; ++row)
            {
                overlap += std::conj(eigenvectors[index(row, i, kN)]) * eigenvectors[index(row, j, kN)];
            }
            if (i == j)
            {
                overlap -= Complex{1.0, 0.0};
            }
            error_squared += std::norm(overlap);
        }
    }

    return std::sqrt(error_squared);
}

struct ElpaRuntimeCheckResult
{
    bool passed;
    double max_residual;
    double orthogonality;
};

ElpaRuntimeCheckResult failed_result()
{
    return {false, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
}

bool needs_runtime_check(const int kernel)
{
    return kernel == ELPA_2STAGE_COMPLEX_AVX_BLOCK2
           || kernel == ELPA_2STAGE_COMPLEX_AVX2_BLOCK2
           || kernel == ELPA_2STAGE_COMPLEX_AVX512_BLOCK2;
}

ElpaRuntimeCheckResult check_complex_block2_kernel(const int kernel)
{
    std::vector<Complex> original_matrix;
    if (!generate_a_tilde(original_matrix))
    {
        return failed_result();
    }

    std::vector<Complex> matrix = original_matrix;
    std::vector<Complex> eigenvectors(kN * kN, Complex{0.0, 0.0});
    std::vector<double> eigenvalues(kN, 0.0);

    int status = ELPA_OK;
    elpa_t handle = elpa_allocate(&status);
    if (status != ELPA_OK || handle == nullptr)
    {
        return failed_result();
    }

    const auto deallocate = [&handle](int& current_status) {
        int deallocate_status = ELPA_OK;
        elpa_deallocate(handle, &deallocate_status);
        if (current_status == ELPA_OK && deallocate_status != ELPA_OK)
        {
            current_status = deallocate_status;
        }
    };

#ifdef _OPENMP
    if (!set_integer(handle, "omp_threads", 1, status))
    {
        deallocate(status);
        return failed_result();
    }
#endif

    if (!set_integer(handle, "na", kN, status)
        || !set_integer(handle, "nev", kNev, status)
        || !set_integer(handle, "local_nrows", kN, status)
        || !set_integer(handle, "local_ncols", kN, status)
        || !set_integer(handle, "nblk", kNblk, status)
        || !set_integer(handle, "mpi_comm_parent", MPI_Comm_c2f(MPI_COMM_SELF), status)
        || !set_integer(handle, "process_row", 0, status)
        || !set_integer(handle, "process_col", 0, status))
    {
        deallocate(status);
        return failed_result();
    }

    status = elpa_setup(handle);
    if (status != ELPA_OK
        || !set_integer(handle, "solver", ELPA_SOLVER_2STAGE, status)
        || !set_integer(handle, "complex_kernel", kernel, status))
    {
        deallocate(status);
        return failed_result();
    }

    elpa_eigenvectors(handle, matrix.data(), eigenvalues.data(), eigenvectors.data(), &status);
    if (status != ELPA_OK)
    {
        deallocate(status);
        return failed_result();
    }

    const double matrix_norm = frobenius_norm(original_matrix);
    double maximum_residual = 0.0;
    for (int band = 0; band < kNev; ++band)
    {
        maximum_residual
            = std::max(maximum_residual,
                       relative_residual(original_matrix,
                                         &eigenvectors[index(0, band, kN)],
                                         eigenvalues[band],
                                         matrix_norm));
    }

    const double orthogonality = orthogonality_error(eigenvectors);
    const bool passed = maximum_residual < kPassThreshold && orthogonality < kPassThreshold;

    deallocate(status);
    return {passed && status == ELPA_OK, maximum_residual, orthogonality};
}

} // namespace

void validate_elpa_complex_kernel(const int kernel, const int set_status, const MPI_Comm comm)
{
    const int local_check_required = set_status == ELPA_OK && needs_runtime_check(kernel) ? 1 : 0;
    int any_check_required = 0;
    MPI_Allreduce(&local_check_required, &any_check_required, 1, MPI_INT, MPI_MAX, comm);
    if (any_check_required == 0)
    {
        return;
    }

    ElpaRuntimeCheckResult result{true, 0.0, 0.0};
    if (local_check_required != 0)
    {
        static std::map<int, ElpaRuntimeCheckResult> checked_kernels;
        std::map<int, ElpaRuntimeCheckResult>::iterator it = checked_kernels.find(kernel);
        if (it == checked_kernels.end())
        {
            it = checked_kernels.emplace(kernel, check_complex_block2_kernel(kernel)).first;
        }
        result = it->second;
    }

    const int local_failed = result.passed ? 0 : 1;
    int any_failed = 0;
    MPI_Allreduce(&local_failed, &any_failed, 1, MPI_INT, MPI_MAX, comm);
    if (any_failed == 0)
    {
        return;
    }

    double max_residual = result.max_residual;
    double max_orthogonality = result.orthogonality;
    MPI_Allreduce(MPI_IN_PLACE, &max_residual, 1, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(MPI_IN_PLACE, &max_orthogonality, 1, MPI_DOUBLE, MPI_MAX, comm);

    std::ostringstream message;
    message << "The selected ELPA complex 2-stage BLOCK2 kernel failed a runtime correctness check:\n\n "
            << "  Maximum residual = " << std::scientific << std::setprecision(6) << max_residual << "\n "
            << "  Orthogonality error = " << max_orthogonality << "\n\n "
            << "The linked ELPA library may return incorrect eigenvectors.\n\n "
            << "Please use another eigensolver, or rebuild ELPA with \"-fno-tree-slp-vectorize\".\n";
    ModuleBase::WARNING_QUIT("ELPA_Solver::setKernel", message.str());
}
