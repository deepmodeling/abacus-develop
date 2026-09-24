#include "source_lcao/module_ri/sternheimer_delta.h"
#include "source_lcao/module_ri/sternheimer_fd_preconditioner.h"

#include "source_base/module_external/blas_connector.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_lcao/module_ri/sternheimer_fd_solver.h"
#include "source_lcao/module_ri/sternheimer_rpa.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace ModuleRI
{
namespace
{

using Complex = SternheimerFDHamiltonian::Complex;
using Vector = SternheimerFDHamiltonian::Vector;

void check_vector_size(const Vector& vector, const std::size_t expected_size, const std::string& context)
{
    if (vector.size() != expected_size)
    {
        throw std::invalid_argument(context + " vector size does not match the Sternheimer delta grid.");
    }
}

void axpy(const Complex alpha, const Vector& x, Vector& y)
{
    check_vector_size(x, y.size(), "Sternheimer delta axpy");
#pragma omp parallel for schedule(static)
    for (std::size_t ir = 0; ir != y.size(); ++ir)
    {
        y[ir] += alpha * x[ir];
    }
}

double difference_norm(const Vector& lhs, const Vector& rhs, const double volume_element)
{
    check_vector_size(lhs, rhs.size(), "Sternheimer delta difference_norm");
    Vector difference(lhs.size(), Complex(0.0, 0.0));
#pragma omp parallel for schedule(static)
    for (std::size_t ir = 0; ir != lhs.size(); ++ir)
    {
        difference[ir] = lhs[ir] - rhs[ir];
    }
    return sternheimer_fd_grid_norm(difference, volume_element);
}

void scale_vector(Vector& vector, const Complex factor)
{
#pragma omp parallel for schedule(static)
    for (std::size_t ir = 0; ir != vector.size(); ++ir)
    {
        vector[ir] *= factor;
    }
}

void validate_grid_function(const SternheimerDeltaGridFunction& function,
                            const std::size_t grid_size,
                            const std::string& context)
{
    check_vector_size(function.values, grid_size, context + " values");
    for (const Vector& gradient: function.gradients)
    {
        check_vector_size(gradient, grid_size, context + " gradient");
    }
}

void validate_finite_grid_vector(const Vector& values, const std::string& context)
{
    for (const Complex value: values)
    {
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
        {
            throw std::invalid_argument(context + " contains a nonfinite value.");
        }
    }
}

void axpy_grid_function(const Complex alpha,
                        const SternheimerDeltaGridFunction& x,
                        SternheimerDeltaGridFunction& y)
{
    axpy(alpha, x.values, y.values);
    for (int direction = 0; direction != 3; ++direction)
    {
        axpy(alpha, x.gradients[static_cast<std::size_t>(direction)],
             y.gradients[static_cast<std::size_t>(direction)]);
    }
}

void scale_grid_function(SternheimerDeltaGridFunction& function, const Complex factor)
{
    scale_vector(function.values, factor);
    for (Vector& gradient: function.gradients)
    {
        scale_vector(gradient, factor);
    }
}

std::vector<Vector> collect_virtual_orbitals(const std::vector<SternheimerDeltaVirtualState>& virtual_states)
{
    std::vector<Vector> orbitals;
    orbitals.reserve(virtual_states.size());
    for (const SternheimerDeltaVirtualState& state: virtual_states)
    {
        orbitals.push_back(state.orbital);
    }
    return orbitals;
}

std::vector<Vector> collect_fixed_subspace(const std::vector<Vector>& occupied_wavefunctions,
                                           const std::vector<SternheimerDeltaVirtualState>& virtual_states)
{
    std::vector<Vector> subspace = occupied_wavefunctions;
    subspace.reserve(occupied_wavefunctions.size() + virtual_states.size());
    for (const SternheimerDeltaVirtualState& state: virtual_states)
    {
        subspace.push_back(state.orbital);
    }
    return subspace;
}

std::vector<double> diagonalize_complex_hermitian(std::vector<Complex>& matrix, const int size)
{
    char jobz = 'V';
    char uplo = 'U';
    const int lda = size;
    int info = 0;
    std::vector<double> eigenvalues(static_cast<std::size_t>(size), 0.0);
    Complex work_query(0.0, 0.0);
    const int minus_one = -1;
    std::vector<double> rwork(static_cast<std::size_t>(std::max(1, 3 * size - 2)), 0.0);

    zheev_(&jobz,
           &uplo,
           &size,
           matrix.data(),
           &lda,
           eigenvalues.data(),
           &work_query,
           &minus_one,
           rwork.data(),
           &info);
    if (info != 0)
    {
        throw std::runtime_error("Sternheimer complex virtual-subspace diagonalization workspace query failed.");
    }
    const int lwork = std::max(1, static_cast<int>(std::ceil(work_query.real())));
    std::vector<Complex> work(static_cast<std::size_t>(lwork), Complex(0.0, 0.0));
    zheev_(&jobz,
           &uplo,
           &size,
           matrix.data(),
           &lda,
           eigenvalues.data(),
           work.data(),
           &lwork,
           rwork.data(),
           &info);
    if (info != 0)
    {
        throw std::runtime_error("Sternheimer complex virtual-subspace diagonalization failed.");
    }
    return eigenvalues;
}

void validate_postprocess_input(const Vector& standard_delta_wavefunction,
                                const SternheimerDeltaPostprocessInput& input)
{
    if (standard_delta_wavefunction.empty())
    {
        throw std::invalid_argument("Sternheimer delta postprocess requires a non-empty standard solution.");
    }
    if (input.volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer delta postprocess requires a positive grid volume element.");
    }
    if (input.virtual_states.size() != input.perturbation_matrix_elements.size())
    {
        throw std::invalid_argument(
            "Sternheimer delta postprocess requires one perturbation matrix element per virtual state.");
    }

    const std::size_t grid_size = standard_delta_wavefunction.size();
    for (const Vector& occupied: input.occupied_wavefunctions)
    {
        check_vector_size(occupied, grid_size, "Sternheimer delta occupied state");
    }
    for (const SternheimerDeltaVirtualState& state: input.virtual_states)
    {
        check_vector_size(state.orbital, grid_size, "Sternheimer delta virtual orbital");
        check_vector_size(state.residual, grid_size, "Sternheimer delta virtual residual");
    }
}

void validate_virtual_states(const std::vector<SternheimerDeltaVirtualState>& virtual_states,
                             const std::size_t grid_size)
{
    for (const SternheimerDeltaVirtualState& state: virtual_states)
    {
        check_vector_size(state.orbital, grid_size, "Sternheimer delta virtual orbital");
        check_vector_size(state.residual, grid_size, "Sternheimer delta virtual residual");
    }
}

SternheimerDeltaCoefficientComponents compute_delta_coefficient_components(
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const std::vector<Complex>& perturbation_matrix_elements,
    const Vector& out_wavefunction,
    const double occupied_eigenvalue,
    const double omega,
    const double volume_element)
{
    SternheimerDeltaCoefficientComponents components;
    components.sos.assign(virtual_states.size(), Complex(0.0, 0.0));
    components.pulay.assign(virtual_states.size(), Complex(0.0, 0.0));
    components.total.assign(virtual_states.size(), Complex(0.0, 0.0));
#pragma omp parallel for schedule(static)
    for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
    {
        const Complex denominator(occupied_eigenvalue - virtual_states[ia].eigenvalue, -omega);
        if (std::abs(denominator) < 1.0e-30)
        {
            throw std::runtime_error("Sternheimer delta coefficient evaluation found a singular denominator.");
        }
        const Complex residual_overlap
            = sternheimer_fd_grid_dot(virtual_states[ia].residual, out_wavefunction, volume_element);
        components.sos[ia] = perturbation_matrix_elements[ia] / denominator;
        components.pulay[ia] = residual_overlap / denominator;
        components.total[ia] = components.sos[ia] + components.pulay[ia];
    }
    return components;
}

void assemble_delta_wavefunction_components(const std::vector<SternheimerDeltaVirtualState>& virtual_states,
                                            SternheimerDeltaPostprocessResult& result)
{
    result.in_sos_wavefunction.assign(result.out_wavefunction.size(), Complex(0.0, 0.0));
    result.in_pulay_wavefunction.assign(result.out_wavefunction.size(), Complex(0.0, 0.0));
    for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
    {
        axpy(result.sos_coefficients[ia], virtual_states[ia].orbital, result.in_sos_wavefunction);
        axpy(result.pulay_coefficients[ia], virtual_states[ia].orbital, result.in_pulay_wavefunction);
    }

    result.reconstructed_wavefunction = result.out_wavefunction;
    axpy(Complex(1.0, 0.0), result.in_sos_wavefunction, result.reconstructed_wavefunction);
    axpy(Complex(1.0, 0.0), result.in_pulay_wavefunction, result.reconstructed_wavefunction);
}

void evaluate_full_grid_delta_hamiltonian_difference(const SternheimerFDHamiltonian& hamiltonian,
                                                     const double volume_element,
                                                     SternheimerDeltaSubspace& subspace)
{
    double reference_norm_squared = 0.0;
    double difference_norm_squared = 0.0;
    double max_abs_difference = 0.0;
    std::vector<Vector> h_virtual(subspace.virtual_states.size());
    for (std::size_t ket = 0; ket != subspace.virtual_states.size(); ++ket)
    {
        hamiltonian.apply(subspace.virtual_states[ket].orbital, h_virtual[ket]);
    }
    for (std::size_t ket = 0; ket != subspace.virtual_states.size(); ++ket)
    {
        for (std::size_t bra = 0; bra != subspace.virtual_states.size(); ++bra)
        {
            const Complex reference = bra == ket ? Complex(subspace.virtual_states[ket].eigenvalue, 0.0)
                                                 : Complex(0.0, 0.0);
            const Complex full_grid = sternheimer_fd_grid_dot(subspace.virtual_states[bra].orbital,
                                                               h_virtual[ket],
                                                               volume_element);
            const double difference = std::abs(full_grid - reference);
            reference_norm_squared += std::norm(reference);
            difference_norm_squared += difference * difference;
            max_abs_difference = std::max(max_abs_difference, difference);
        }
    }
    const double reference_norm = std::sqrt(reference_norm_squared);
    const double difference_norm = std::sqrt(difference_norm_squared);
    subspace.full_grid_hamiltonian_relative_difference
        = reference_norm > 0.0 ? difference_norm / reference_norm : difference_norm;
    subspace.full_grid_hamiltonian_max_abs_difference = max_abs_difference;
}

constexpr std::size_t delta_grid_block_size = 4096;

const Vector& grid_function_component(const SternheimerDeltaGridFunction& function, const int component)
{
    return component < 0 ? function.values : function.gradients[static_cast<std::size_t>(component)];
}

Vector& grid_function_component(SternheimerDeltaGridFunction& function, const int component)
{
    return component < 0 ? function.values : function.gradients[static_cast<std::size_t>(component)];
}

void pack_grid_function_component(const std::vector<SternheimerDeltaGridFunction>& functions,
                                  const int component,
                                  const std::size_t grid_begin,
                                  const int grid_count,
                                  std::vector<Complex>& packed)
{
    packed.resize(static_cast<std::size_t>(grid_count) * functions.size());
#pragma omp parallel for schedule(static)
    for (std::size_t function = 0; function != functions.size(); ++function)
    {
        const Vector& source = grid_function_component(functions[function], component);
        Complex* destination = packed.data() + static_cast<std::size_t>(grid_count) * function;
        std::copy_n(source.data() + grid_begin, grid_count, destination);
    }
}

void pack_vector_component(const std::vector<const Vector*>& vectors,
                           const std::size_t grid_begin,
                           const int grid_count,
                           std::vector<Complex>& packed)
{
    packed.resize(static_cast<std::size_t>(grid_count) * vectors.size());
#pragma omp parallel for schedule(static)
    for (std::size_t vector = 0; vector != vectors.size(); ++vector)
    {
        Complex* destination = packed.data() + static_cast<std::size_t>(grid_count) * vector;
        std::copy_n(vectors[vector]->data() + grid_begin, grid_count, destination);
    }
}

void accumulate_grid_component_product(
    const std::vector<SternheimerDeltaGridFunction>& left,
    const int left_component,
    const std::vector<SternheimerDeltaGridFunction>& right,
    const int right_component,
    const std::vector<double>* right_weights,
    const Complex alpha,
    std::vector<Complex>& output)
{
    if (left.empty() || right.empty())
    {
        return;
    }
    const std::size_t grid_size = grid_function_component(left.front(), left_component).size();
    const int left_size = static_cast<int>(left.size());
    const int right_size = static_cast<int>(right.size());
    if (output.size() != left.size() * right.size())
    {
        throw std::invalid_argument("Sternheimer blocked grid product output size mismatch.");
    }
    if (right_weights != nullptr && right_weights->size() != grid_size)
    {
        throw std::invalid_argument("Sternheimer blocked grid product weight size mismatch.");
    }

    std::vector<Complex> left_packed;
    std::vector<Complex> right_packed;
    for (std::size_t grid_begin = 0; grid_begin < grid_size; grid_begin += delta_grid_block_size)
    {
        const int grid_count = static_cast<int>(std::min(delta_grid_block_size, grid_size - grid_begin));
        pack_grid_function_component(left, left_component, grid_begin, grid_count, left_packed);
        const bool reuse_left = &left == &right && left_component == right_component
                                && right_weights == nullptr;
        const Complex* right_data = left_packed.data();
        if (!reuse_left)
        {
            pack_grid_function_component(right, right_component, grid_begin, grid_count, right_packed);
            if (right_weights != nullptr)
            {
#pragma omp parallel for schedule(static)
                for (std::size_t function = 0; function != right.size(); ++function)
                {
                    Complex* column
                        = right_packed.data() + static_cast<std::size_t>(grid_count) * function;
                    for (int ir = 0; ir != grid_count; ++ir)
                    {
                        column[ir] *= (*right_weights)[grid_begin + static_cast<std::size_t>(ir)];
                    }
                }
            }
            right_data = right_packed.data();
        }
        BlasConnector::gemm_cm('C',
                               'N',
                               left_size,
                               right_size,
                               grid_count,
                               alpha,
                               left_packed.data(),
                               grid_count,
                               right_data,
                               grid_count,
                               Complex(1.0, 0.0),
                               output.data(),
                               left_size);
    }
}

void project_grid_functions_blocked(
    const std::vector<SternheimerDeltaGridFunction>& occupied_functions,
    std::vector<SternheimerDeltaGridFunction>& candidate_functions,
    const double volume_element)
{
    if (occupied_functions.empty() || candidate_functions.empty())
    {
        return;
    }
    const std::size_t grid_size = candidate_functions.front().values.size();
    const int occupied_count = static_cast<int>(occupied_functions.size());
    const int candidate_count = static_cast<int>(candidate_functions.size());
    std::vector<Complex> coefficients(occupied_functions.size() * candidate_functions.size());
    std::vector<Complex> occupied_packed;
    std::vector<Complex> corrections;

    for (int pass = 0; pass != 2; ++pass)
    {
        std::fill(coefficients.begin(), coefficients.end(), Complex(0.0, 0.0));
        accumulate_grid_component_product(occupied_functions,
                                          -1,
                                          candidate_functions,
                                          -1,
                                          nullptr,
                                          Complex(volume_element, 0.0),
                                          coefficients);
        for (int component = -1; component != 3; ++component)
        {
            for (std::size_t grid_begin = 0; grid_begin < grid_size;
                 grid_begin += delta_grid_block_size)
            {
                const int grid_count
                    = static_cast<int>(std::min(delta_grid_block_size, grid_size - grid_begin));
                pack_grid_function_component(
                    occupied_functions, component, grid_begin, grid_count, occupied_packed);
                corrections.assign(static_cast<std::size_t>(grid_count) * candidate_functions.size(),
                                   Complex(0.0, 0.0));
                BlasConnector::gemm_cm('N',
                                       'N',
                                       grid_count,
                                       candidate_count,
                                       occupied_count,
                                       Complex(-1.0, 0.0),
                                       occupied_packed.data(),
                                       grid_count,
                                       coefficients.data(),
                                       occupied_count,
                                       Complex(0.0, 0.0),
                                       corrections.data(),
                                       grid_count);
#pragma omp parallel for schedule(static)
                for (std::size_t candidate = 0; candidate != candidate_functions.size(); ++candidate)
                {
                    Vector& destination
                        = grid_function_component(candidate_functions[candidate], component);
                    const Complex* correction
                        = corrections.data() + static_cast<std::size_t>(grid_count) * candidate;
                    for (int ir = 0; ir != grid_count; ++ir)
                    {
                        destination[grid_begin + static_cast<std::size_t>(ir)] += correction[ir];
                    }
                }
            }
        }
    }
}

void hermitize_matrix(std::vector<Complex>& matrix, const int size)
{
    for (int column = 0; column != size; ++column)
    {
        matrix[static_cast<std::size_t>(column)
               + static_cast<std::size_t>(size) * static_cast<std::size_t>(column)]
            = Complex(matrix[static_cast<std::size_t>(column)
                             + static_cast<std::size_t>(size) * static_cast<std::size_t>(column)]
                          .real(),
                      0.0);
        for (int row = 0; row != column; ++row)
        {
            const std::size_t upper
                = static_cast<std::size_t>(row)
                  + static_cast<std::size_t>(size) * static_cast<std::size_t>(column);
            const std::size_t lower
                = static_cast<std::size_t>(column)
                  + static_cast<std::size_t>(size) * static_cast<std::size_t>(row);
            const Complex value = 0.5 * (matrix[upper] + std::conj(matrix[lower]));
            matrix[upper] = value;
            matrix[lower] = std::conj(value);
        }
    }
}

SternheimerDeltaGridMatrices assemble_delta_sternheimer_cross_matrices_blocked(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerDeltaGridFunction>& left_functions,
    const std::vector<SternheimerDeltaGridFunction>& right_functions,
    const double volume_element)
{
    const std::size_t max_elements = Vector().max_size();
    const std::size_t max_blas_dimension = static_cast<std::size_t>(std::numeric_limits<int>::max());
    if (left_functions.size() > max_blas_dimension || right_functions.size() > max_blas_dimension
        || (!right_functions.empty() && left_functions.size() > max_elements / right_functions.size()))
    {
        throw std::overflow_error("Sternheimer blocked cross matrix dimensions exceed BLAS or vector limits.");
    }
    const int left_size = static_cast<int>(left_functions.size());
    const int right_size = static_cast<int>(right_functions.size());
    const std::size_t grid_size = left_functions.empty() ? 0 : left_functions.front().values.size();
    const std::size_t packed_rows = std::min(delta_grid_block_size, grid_size);
    if (!right_functions.empty() && packed_rows != 0
        && std::max(left_functions.size(), right_functions.size()) > max_elements / packed_rows)
    {
        throw std::overflow_error("Sternheimer blocked cross matrix packing size exceeds vector limits.");
    }
    const std::size_t matrix_size = left_functions.size() * right_functions.size();
    SternheimerDeltaGridMatrices matrices;
    matrices.row_count = left_functions.size();
    matrices.column_count = right_functions.size();
    matrices.overlap.assign(matrix_size, Complex(0.0, 0.0));
    matrices.kinetic.assign(matrix_size, Complex(0.0, 0.0));
    matrices.local_potential.assign(matrix_size, Complex(0.0, 0.0));
    matrices.nonlocal.assign(matrix_size, Complex(0.0, 0.0));
    matrices.hamiltonian.assign(matrix_size, Complex(0.0, 0.0));
    if (left_functions.empty() || right_functions.empty())
    {
        return matrices;
    }
    accumulate_grid_component_product(left_functions,
                                      -1,
                                      right_functions,
                                      -1,
                                      nullptr,
                                      Complex(volume_element, 0.0),
                                      matrices.overlap);
    for (int direction = 0; direction != 3; ++direction)
    {
        accumulate_grid_component_product(
            left_functions,
            direction,
            right_functions,
            direction,
            nullptr,
            Complex(volume_element * hamiltonian.kinetic_prefactor(), 0.0),
            matrices.kinetic);
    }
    accumulate_grid_component_product(left_functions,
                                      -1,
                                      right_functions,
                                      -1,
                                      &hamiltonian.local_potential(),
                                      Complex(volume_element, 0.0),
                                      matrices.local_potential);

    const SternheimerFDNonlocalProjector* nonlocal_projector = hamiltonian.nonlocal_projector();
    if (nonlocal_projector != nullptr)
    {
        if (std::abs(nonlocal_projector->volume_element() - volume_element)
            > 1.0e-12 * std::max(nonlocal_projector->volume_element(), volume_element))
        {
            throw std::invalid_argument(
                "Sternheimer blocked nonlocal projector uses a different grid volume element.");
        }
        constexpr int nonlocal_block_size = 64;
        std::vector<Complex> basis_packed;
        std::vector<Complex> nonlocal_packed;
        for (int block_begin = 0; block_begin < right_size;)
        {
            const int block_count = std::min(nonlocal_block_size, right_size - block_begin);
            std::vector<Vector> input_vectors;
            input_vectors.reserve(static_cast<std::size_t>(block_count));
            for (int column = 0; column != block_count; ++column)
            {
                input_vectors.push_back(
                    right_functions[static_cast<std::size_t>(block_begin + column)].values);
            }
            std::vector<Vector> nonlocal_vectors;
            nonlocal_projector->apply_batch(input_vectors, nonlocal_vectors);
            std::vector<const Vector*> nonlocal_vector_pointers(static_cast<std::size_t>(block_count));
            for (int column = 0; column != block_count; ++column)
            {
                nonlocal_vector_pointers[static_cast<std::size_t>(column)]
                    = &nonlocal_vectors[static_cast<std::size_t>(column)];
            }
            Complex* output = matrices.nonlocal.data()
                              + static_cast<std::size_t>(left_size)
                                    * static_cast<std::size_t>(block_begin);
            for (std::size_t grid_begin = 0; grid_begin < grid_size;
                 grid_begin += delta_grid_block_size)
            {
                const int grid_count
                    = static_cast<int>(std::min(delta_grid_block_size, grid_size - grid_begin));
                pack_grid_function_component(
                    left_functions, -1, grid_begin, grid_count, basis_packed);
                pack_vector_component(
                    nonlocal_vector_pointers, grid_begin, grid_count, nonlocal_packed);
                BlasConnector::gemm_cm('C',
                                       'N',
                                       left_size,
                                       block_count,
                                       grid_count,
                                       Complex(volume_element, 0.0),
                                       basis_packed.data(),
                                       grid_count,
                                       nonlocal_packed.data(),
                                       grid_count,
                                       Complex(1.0, 0.0),
                                       output,
                                       left_size);
            }
            block_begin += block_count;
        }
    }

#pragma omp parallel for schedule(static)
    for (std::size_t index = 0; index != matrix_size; ++index)
    {
        matrices.hamiltonian[index] = matrices.kinetic[index] + matrices.local_potential[index]
                                      + matrices.nonlocal[index];
    }
    return matrices;
}

SternheimerDeltaGridMatrices assemble_delta_sternheimer_grid_matrices_blocked(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerDeltaGridFunction>& basis_functions,
    const double volume_element)
{
    // Preserve the square path's packing, accumulation order, and S/H cleanup.
    SternheimerDeltaGridMatrices matrices = assemble_delta_sternheimer_cross_matrices_blocked(
        hamiltonian, basis_functions, basis_functions, volume_element);
    hermitize_matrix(matrices.overlap, static_cast<int>(basis_functions.size()));
    hermitize_matrix(matrices.hamiltonian, static_cast<int>(basis_functions.size()));
    return matrices;
}

std::vector<double> diagonalize_complex_generalized_hermitian(std::vector<Complex>& hamiltonian,
                                                              std::vector<Complex>& overlap,
                                                              const int size)
{
    const int itype = 1;
    const char jobz = 'V';
    const char uplo = 'U';
    const int lda = size;
    const int ldb = size;
    int info = 0;
    std::vector<double> eigenvalues(static_cast<std::size_t>(size), 0.0);
    Complex work_query(0.0, 0.0);
    double rwork_query = 0.0;
    int iwork_query = 0;
    const int minus_one = -1;
    zhegvd_(&itype,
            &jobz,
            &uplo,
            &size,
            hamiltonian.data(),
            &lda,
            overlap.data(),
            &ldb,
            eigenvalues.data(),
            &work_query,
            &minus_one,
            &rwork_query,
            &minus_one,
            &iwork_query,
            &minus_one,
            &info);
    if (info != 0)
    {
        throw std::runtime_error(
            "Sternheimer blocked generalized eigensolver workspace query failed.");
    }
    const int lwork = std::max(1, static_cast<int>(std::ceil(work_query.real())));
    const int lrwork = std::max(1, static_cast<int>(std::ceil(rwork_query)));
    const int liwork = std::max(1, iwork_query);
    std::vector<Complex> work(static_cast<std::size_t>(lwork));
    std::vector<double> rwork(static_cast<std::size_t>(lrwork));
    std::vector<int> iwork(static_cast<std::size_t>(liwork));
    zhegvd_(&itype,
            &jobz,
            &uplo,
            &size,
            hamiltonian.data(),
            &lda,
            overlap.data(),
            &ldb,
            eigenvalues.data(),
            work.data(),
            &lwork,
            rwork.data(),
            &lrwork,
            iwork.data(),
            &liwork,
            &info);
    if (info != 0)
    {
        throw std::runtime_error("Sternheimer blocked generalized eigensolver failed with info="
                                 + std::to_string(info) + ".");
    }
    return eigenvalues;
}

void transform_grid_functions_blocked(
    const std::vector<SternheimerDeltaGridFunction>& basis_functions,
    const std::vector<Complex>& coefficients,
    const int output_size,
    const bool retain_grid_functions,
    SternheimerDeltaSubspace& subspace)
{
    const int basis_size = static_cast<int>(basis_functions.size());
    const std::size_t grid_size = basis_functions.front().values.size();
    subspace.virtual_states.resize(static_cast<std::size_t>(output_size));
    for (SternheimerDeltaVirtualState& state: subspace.virtual_states)
    {
        state.orbital.resize(grid_size);
    }
    if (retain_grid_functions)
    {
        subspace.grid_functions.resize(static_cast<std::size_t>(output_size));
        for (SternheimerDeltaGridFunction& function: subspace.grid_functions)
        {
            function.values.resize(grid_size);
            for (Vector& gradient: function.gradients)
            {
                gradient.resize(grid_size);
            }
        }
    }

    std::vector<Complex> basis_packed;
    std::vector<Complex> transformed;
    const int last_component = retain_grid_functions ? 2 : -1;
    for (int component = -1; component <= last_component; ++component)
    {
        for (std::size_t grid_begin = 0; grid_begin < grid_size;
             grid_begin += delta_grid_block_size)
        {
            const int grid_count
                = static_cast<int>(std::min(delta_grid_block_size, grid_size - grid_begin));
            pack_grid_function_component(
                basis_functions, component, grid_begin, grid_count, basis_packed);
            transformed.resize(static_cast<std::size_t>(grid_count)
                               * static_cast<std::size_t>(output_size));
            BlasConnector::gemm_cm('N',
                                   'N',
                                   grid_count,
                                   output_size,
                                   basis_size,
                                   Complex(1.0, 0.0),
                                   basis_packed.data(),
                                   grid_count,
                                   coefficients.data(),
                                   basis_size,
                                   Complex(0.0, 0.0),
                                   transformed.data(),
                                   grid_count);
#pragma omp parallel for schedule(static)
            for (int state_index = 0; state_index != output_size; ++state_index)
            {
                const Complex* source
                    = transformed.data() + static_cast<std::size_t>(grid_count)
                                                   * static_cast<std::size_t>(state_index);
                if (component < 0)
                {
                    std::copy_n(source,
                                grid_count,
                                subspace.virtual_states[static_cast<std::size_t>(state_index)]
                                        .orbital.data()
                                    + grid_begin);
                    if (retain_grid_functions)
                    {
                        std::copy_n(source,
                                    grid_count,
                                    subspace.grid_functions[static_cast<std::size_t>(state_index)]
                                            .values.data()
                                        + grid_begin);
                    }
                }
                else
                {
                    std::copy_n(source,
                                grid_count,
                                subspace.grid_functions[static_cast<std::size_t>(state_index)]
                                        .gradients[static_cast<std::size_t>(component)]
                                        .data()
                                    + grid_begin);
                }
            }
        }
    }
}

void project_vector_batch_blocked(const std::vector<const Vector*>& basis,
                                  const std::vector<Vector*>& vectors,
                                  const double volume_element)
{
    if (basis.empty() || vectors.empty())
    {
        return;
    }
    const std::size_t grid_size = basis.front()->size();
    const int basis_size = static_cast<int>(basis.size());
    const int vector_count = static_cast<int>(vectors.size());
    std::vector<const Vector*> vector_sources(vectors.begin(), vectors.end());
    std::vector<Complex> coefficients(basis.size() * vectors.size(), Complex(0.0, 0.0));
    std::vector<Complex> basis_packed;
    std::vector<Complex> vectors_packed;
    for (std::size_t grid_begin = 0; grid_begin < grid_size; grid_begin += delta_grid_block_size)
    {
        const int grid_count = static_cast<int>(std::min(delta_grid_block_size, grid_size - grid_begin));
        pack_vector_component(basis, grid_begin, grid_count, basis_packed);
        pack_vector_component(vector_sources, grid_begin, grid_count, vectors_packed);
        BlasConnector::gemm_cm('C',
                               'N',
                               basis_size,
                               vector_count,
                               grid_count,
                               Complex(volume_element, 0.0),
                               basis_packed.data(),
                               grid_count,
                               vectors_packed.data(),
                               grid_count,
                               Complex(1.0, 0.0),
                               coefficients.data(),
                               basis_size);
    }
    for (int basis_index = 0; basis_index != basis_size; ++basis_index)
    {
        const Complex norm
            = sternheimer_fd_grid_dot(*basis[static_cast<std::size_t>(basis_index)],
                                      *basis[static_cast<std::size_t>(basis_index)],
                                      volume_element);
        if (std::abs(norm) == 0.0)
        {
            for (int vector = 0; vector != vector_count; ++vector)
            {
                coefficients[static_cast<std::size_t>(basis_index)
                             + basis.size() * static_cast<std::size_t>(vector)]
                    = Complex(0.0, 0.0);
            }
        }
        else
        {
            for (int vector = 0; vector != vector_count; ++vector)
            {
                coefficients[static_cast<std::size_t>(basis_index)
                             + basis.size() * static_cast<std::size_t>(vector)]
                    /= norm;
            }
        }
    }

    std::vector<Complex> corrections;
    for (std::size_t grid_begin = 0; grid_begin < grid_size; grid_begin += delta_grid_block_size)
    {
        const int grid_count = static_cast<int>(std::min(delta_grid_block_size, grid_size - grid_begin));
        pack_vector_component(basis, grid_begin, grid_count, basis_packed);
        corrections.resize(static_cast<std::size_t>(grid_count) * vectors.size());
        BlasConnector::gemm_cm('N',
                               'N',
                               grid_count,
                               vector_count,
                               basis_size,
                               Complex(-1.0, 0.0),
                               basis_packed.data(),
                               grid_count,
                               coefficients.data(),
                               basis_size,
                               Complex(0.0, 0.0),
                               corrections.data(),
                               grid_count);
#pragma omp parallel for schedule(static)
        for (std::size_t vector = 0; vector != vectors.size(); ++vector)
        {
            Complex* destination = vectors[vector]->data() + grid_begin;
            const Complex* correction
                = corrections.data() + static_cast<std::size_t>(grid_count) * vector;
            for (int ir = 0; ir != grid_count; ++ir)
            {
                destination[ir] += correction[ir];
            }
        }
    }
}

SternheimerDeltaSubspace build_complete_reference_delta_sternheimer_subspace_blocked(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerDeltaGridFunction>& occupied_functions,
    std::vector<SternheimerDeltaGridFunction> candidate_functions,
    const double volume_element,
    const SternheimerDeltaSubspaceOptions& options)
{
    project_grid_functions_blocked(occupied_functions, candidate_functions, volume_element);
    SternheimerDeltaGridMatrices matrices = assemble_delta_sternheimer_grid_matrices_blocked(
        hamiltonian, candidate_functions, volume_element);
    const int basis_size = static_cast<int>(candidate_functions.size());
    for (int state = 0; state != basis_size; ++state)
    {
        const double norm = std::sqrt(std::max(
            0.0,
            matrices.overlap[static_cast<std::size_t>(state)
                             + candidate_functions.size() * static_cast<std::size_t>(state)]
                .real()));
        if (norm <= options.norm_tolerance)
        {
            throw std::runtime_error(
                "Sternheimer blocked complete subspace contains a linearly dependent candidate.");
        }
    }
    const std::vector<double> eigenvalues = diagonalize_complex_generalized_hermitian(
        matrices.hamiltonian, matrices.overlap, basis_size);

    SternheimerDeltaSubspace subspace;
    subspace.accepted_candidates = basis_size;
    subspace.discarded_candidates = 0;
    subspace.used_block_generalized_eigensolver = true;
    transform_grid_functions_blocked(candidate_functions,
                                     matrices.hamiltonian,
                                     basis_size,
                                     options.retain_grid_functions,
                                     subspace);
    for (int state = 0; state != basis_size; ++state)
    {
        subspace.virtual_states[static_cast<std::size_t>(state)].eigenvalue
            = eigenvalues[static_cast<std::size_t>(state)];
    }
    std::vector<SternheimerDeltaGridFunction>().swap(candidate_functions);

    if (options.evaluate_full_grid_difference)
    {
        evaluate_full_grid_delta_hamiltonian_difference(hamiltonian, volume_element, subspace);
    }
    else
    {
        subspace.full_grid_hamiltonian_relative_difference = -1.0;
        subspace.full_grid_hamiltonian_max_abs_difference = -1.0;
    }

    std::vector<const Vector*> residual_projector;
    residual_projector.reserve(occupied_functions.size() + subspace.virtual_states.size());
    for (const SternheimerDeltaGridFunction& occupied: occupied_functions)
    {
        residual_projector.push_back(&occupied.values);
    }
    for (const SternheimerDeltaVirtualState& state: subspace.virtual_states)
    {
        residual_projector.push_back(&state.orbital);
    }
    std::vector<Vector*> residuals;
    residuals.reserve(subspace.virtual_states.size());
    for (SternheimerDeltaVirtualState& state: subspace.virtual_states)
    {
        hamiltonian.apply(state.orbital, state.residual);
#pragma omp parallel for schedule(static)
        for (std::size_t ir = 0; ir != state.residual.size(); ++ir)
        {
            state.residual[ir] -= state.eigenvalue * state.orbital[ir];
        }
        residuals.push_back(&state.residual);
    }
    project_vector_batch_blocked(residual_projector, residuals, volume_element);
    return subspace;
}

} // namespace

SternheimerDeltaABlockMode parse_sternheimer_delta_a_block_mode(const std::string& name)
{
    if (name == "reference_value_gradient")
    {
        return SternheimerDeltaABlockMode::ReferenceValueGradient;
    }
    if (name == "grid")
    {
        return SternheimerDeltaABlockMode::FullGrid;
    }
    throw std::invalid_argument("Unknown Sternheimer Delta A-block mode: " + name);
}

const char* sternheimer_delta_a_block_mode_name(const SternheimerDeltaABlockMode mode)
{
    switch (mode)
    {
        case SternheimerDeltaABlockMode::ReferenceValueGradient:
            return "reference_value_gradient";
        case SternheimerDeltaABlockMode::FullGrid:
            return "grid";
    }
    throw std::invalid_argument("Invalid Sternheimer Delta A-block mode.");
}

int sternheimer_delta_virtual_state_limit(const int requested_states,
                                          const int candidate_states,
                                          const int occupied_states)
{
    if (requested_states < 0 || candidate_states < 0 || occupied_states < 0)
    {
        throw std::invalid_argument("Sternheimer delta subspace dimensions must be non-negative.");
    }
    if (occupied_states > candidate_states)
    {
        throw std::invalid_argument("Sternheimer delta occupied dimension exceeds the candidate dimension.");
    }
    const int available_states = candidate_states - occupied_states;
    return requested_states == 0 ? available_states : std::min(requested_states, available_states);
}

SternheimerDeltaGridFunction make_delta_sternheimer_grid_function_with_fd_gradients(
    const SternheimerFDHamiltonian::Vector& values,
    const SternheimerFDHamiltonian::Grid& grid)
{
    if (grid.nx <= 0 || grid.ny <= 0 || grid.nz <= 0 || grid.hx <= 0.0 || grid.hy <= 0.0
        || grid.hz <= 0.0)
    {
        throw std::invalid_argument("Sternheimer delta grid gradient requires a valid uniform grid.");
    }
    if (static_cast<int>(values.size()) != grid.size())
    {
        throw std::invalid_argument("Sternheimer delta grid gradient value size does not match the grid.");
    }

    SternheimerDeltaGridFunction function;
    function.values = values;
    for (Vector& gradient: function.gradients)
    {
        gradient.assign(values.size(), Complex(0.0, 0.0));
    }
    const auto index = [&grid](const int ix, const int iy, const int iz) {
        return (ix * grid.ny + iy) * grid.nz + iz;
    };
    const auto derivative = [&values, &grid, &index](const int ix,
                                                     const int iy,
                                                     const int iz,
                                                     const int direction) {
        const int coordinate = direction == 0 ? ix : (direction == 1 ? iy : iz);
        const int count = direction == 0 ? grid.nx : (direction == 1 ? grid.ny : grid.nz);
        const double spacing = direction == 0 ? grid.hx : (direction == 1 ? grid.hy : grid.hz);
        if (count == 1)
        {
            return Complex(0.0, 0.0);
        }

        int lower = coordinate - 1;
        int upper = coordinate + 1;
        double denominator = 2.0 * spacing;
        if (grid.periodic)
        {
            lower = (lower + count) % count;
            upper %= count;
        }
        else if (coordinate == 0)
        {
            lower = 0;
            upper = 1;
            denominator = spacing;
        }
        else if (coordinate == count - 1)
        {
            lower = count - 2;
            upper = count - 1;
            denominator = spacing;
        }

        const int lower_index = direction == 0 ? index(lower, iy, iz)
                                : direction == 1 ? index(ix, lower, iz)
                                                 : index(ix, iy, lower);
        const int upper_index = direction == 0 ? index(upper, iy, iz)
                                : direction == 1 ? index(ix, upper, iz)
                                                 : index(ix, iy, upper);
        return (values[static_cast<std::size_t>(upper_index)]
                - values[static_cast<std::size_t>(lower_index)])
               / denominator;
    };

    for (int ix = 0; ix != grid.nx; ++ix)
    {
        for (int iy = 0; iy != grid.ny; ++iy)
        {
            for (int iz = 0; iz != grid.nz; ++iz)
            {
                const std::size_t ir = static_cast<std::size_t>(index(ix, iy, iz));
                for (int direction = 0; direction != 3; ++direction)
                {
                    function.gradients[static_cast<std::size_t>(direction)][ir]
                        = derivative(ix, iy, iz, direction);
                }
            }
        }
    }
    return function;
}

void accumulate_delta_sternheimer_bloch_samples(
    const std::vector<double>& sampled_values,
    const std::array<std::vector<double>, 3>& sampled_gradients,
    const int sample_count,
    const int orbital_count,
    const std::size_t grid_begin,
    const std::size_t function_begin,
    const SternheimerReducedKPoint& kpoint,
    const std::array<int, 3>& lattice_translation,
    std::vector<SternheimerDeltaGridFunction>& functions)
{
    if (sample_count < 0 || orbital_count <= 0)
    {
        throw std::invalid_argument("Sternheimer Bloch samples require non-negative samples and positive orbitals.");
    }
    const std::size_t sample_size = static_cast<std::size_t>(sample_count);
    const std::size_t orbital_size = static_cast<std::size_t>(orbital_count);
    const std::size_t buffer_size = sample_size * orbital_size;
    if (sampled_values.size() != buffer_size)
    {
        throw std::invalid_argument("Sternheimer Bloch value buffer size is inconsistent.");
    }
    for (const std::vector<double>& gradient: sampled_gradients)
    {
        if (gradient.size() != buffer_size)
        {
            throw std::invalid_argument("Sternheimer Bloch gradient buffer size is inconsistent.");
        }
    }
    if (function_begin > functions.size() || orbital_size > functions.size() - function_begin)
    {
        throw std::invalid_argument("Sternheimer Bloch sample function range is out of bounds.");
    }
    for (std::size_t orbital = 0; orbital != orbital_size; ++orbital)
    {
        const SternheimerDeltaGridFunction& function = functions[function_begin + orbital];
        if (grid_begin > function.values.size() || sample_size > function.values.size() - grid_begin)
        {
            throw std::invalid_argument("Sternheimer Bloch sample grid range is out of bounds.");
        }
        validate_grid_function(function, function.values.size(), "Sternheimer Bloch sample");
    }

    const Complex phase = sternheimer_bloch_phase(kpoint, lattice_translation);
    for (std::size_t sample = 0; sample != sample_size; ++sample)
    {
        const std::size_t grid_index = grid_begin + sample;
        for (std::size_t orbital = 0; orbital != orbital_size; ++orbital)
        {
            const std::size_t buffer_index = sample * orbital_size + orbital;
            SternheimerDeltaGridFunction& function = functions[function_begin + orbital];
            function.values[grid_index] += phase * sampled_values[buffer_index];
            for (std::size_t direction = 0; direction != function.gradients.size(); ++direction)
            {
                function.gradients[direction][grid_index]
                    += phase * sampled_gradients[direction][buffer_index];
            }
        }
    }
}

void accumulate_delta_sternheimer_lcao_state_samples(
    const std::vector<double>& sampled_values,
    const std::array<std::vector<double>, 3>& sampled_gradients,
    const int sample_count,
    const int orbital_count,
    const std::size_t grid_begin,
    const std::size_t ao_begin,
    const std::vector<std::vector<SternheimerFDHamiltonian::Complex>>& coefficients,
    const SternheimerReducedKPoint& kpoint,
    const std::array<int, 3>& lattice_translation,
    std::vector<SternheimerDeltaGridFunction>& functions)
{
    if (sample_count < 0 || orbital_count <= 0 || coefficients.size() != functions.size())
    {
        throw std::invalid_argument(
            "Sternheimer direct LCAO sampling requires valid sample, orbital, coefficient, and state counts.");
    }
    const std::size_t sample_size = static_cast<std::size_t>(sample_count);
    const std::size_t orbital_size = static_cast<std::size_t>(orbital_count);
    const std::size_t buffer_size = sample_size * orbital_size;
    if (sampled_values.size() != buffer_size)
    {
        throw std::invalid_argument("Sternheimer direct LCAO value buffer size is inconsistent.");
    }
    for (const std::vector<double>& gradient: sampled_gradients)
    {
        if (gradient.size() != buffer_size)
        {
            throw std::invalid_argument("Sternheimer direct LCAO gradient buffer size is inconsistent.");
        }
    }
    for (std::size_t state_index = 0; state_index != functions.size(); ++state_index)
    {
        if (ao_begin > coefficients[state_index].size()
            || orbital_size > coefficients[state_index].size() - ao_begin)
        {
            throw std::invalid_argument("Sternheimer direct LCAO coefficient AO range is out of bounds.");
        }
        const SternheimerDeltaGridFunction& function = functions[state_index];
        if (grid_begin > function.values.size() || sample_size > function.values.size() - grid_begin)
        {
            throw std::invalid_argument("Sternheimer direct LCAO sample grid range is out of bounds.");
        }
        validate_grid_function(function, function.values.size(), "Sternheimer direct LCAO sample");
    }

    const Complex image_phase = sternheimer_bloch_phase(kpoint, lattice_translation);
#pragma omp parallel for schedule(static)
    for (std::size_t sample = 0; sample != sample_size; ++sample)
    {
        const std::size_t grid_index = grid_begin + sample;
        for (std::size_t state_index = 0; state_index != functions.size(); ++state_index)
        {
            SternheimerDeltaGridFunction& state = functions[state_index];
            Complex value(0.0, 0.0);
            std::array<Complex, 3> gradient{{Complex(0.0, 0.0),
                                             Complex(0.0, 0.0),
                                             Complex(0.0, 0.0)}};
            for (std::size_t orbital = 0; orbital != orbital_size; ++orbital)
            {
                const std::size_t buffer_index = sample * orbital_size + orbital;
                const Complex coefficient = coefficients[state_index][ao_begin + orbital];
                value += coefficient * sampled_values[buffer_index];
                for (std::size_t direction = 0; direction != gradient.size(); ++direction)
                {
                    gradient[direction] += coefficient * sampled_gradients[direction][buffer_index];
                }
            }
            state.values[grid_index] += image_phase * value;
            for (std::size_t direction = 0; direction != gradient.size(); ++direction)
            {
                state.gradients[direction][grid_index] += image_phase * gradient[direction];
            }
        }
    }
}

SternheimerDeltaGridFunction linear_combination_delta_sternheimer_grid_functions(
    const std::vector<SternheimerDeltaGridFunction>& basis_functions,
    const std::vector<SternheimerFDHamiltonian::Complex>& coefficients)
{
    if (basis_functions.size() != coefficients.size())
    {
        throw std::invalid_argument("Sternheimer LCAO grid combination requires one coefficient per AO.");
    }

    SternheimerDeltaGridFunction state;
    if (basis_functions.empty())
    {
        return state;
    }

    const std::size_t grid_size = basis_functions.front().values.size();
    state.values.assign(grid_size, Complex(0.0, 0.0));
    for (Vector& gradient: state.gradients)
    {
        gradient.assign(grid_size, Complex(0.0, 0.0));
    }
    for (std::size_t ia = 0; ia != basis_functions.size(); ++ia)
    {
        validate_grid_function(basis_functions[ia], grid_size, "Sternheimer LCAO AO grid function");
        axpy_grid_function(coefficients[ia], basis_functions[ia], state);
    }
    return state;
}

std::vector<SternheimerDeltaGridFunction> orthonormalize_delta_sternheimer_grid_functions(
    const std::vector<SternheimerDeltaGridFunction>& functions,
    const double volume_element,
    const double norm_tolerance)
{
    if (volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer grid-function orthonormalization requires a positive volume element.");
    }
    if (norm_tolerance < 0.0)
    {
        throw std::invalid_argument("Sternheimer grid-function orthonormalization requires a non-negative tolerance.");
    }
    if (functions.empty())
    {
        return {};
    }

    const std::size_t grid_size = functions.front().values.size();
    auto dot = [volume_element](const Vector& lhs, const Vector& rhs) {
        return sternheimer_fd_grid_dot(lhs, rhs, volume_element);
    };

    std::vector<SternheimerDeltaGridFunction> orthonormal;
    orthonormal.reserve(functions.size());
    for (SternheimerDeltaGridFunction function: functions)
    {
        validate_grid_function(function, grid_size, "Sternheimer grid-function orthonormalization");
        for (int pass = 0; pass != 2; ++pass)
        {
            for (const SternheimerDeltaGridFunction& accepted: orthonormal)
            {
                axpy_grid_function(-dot(accepted.values, function.values), accepted, function);
            }
        }
        const double norm = sternheimer_fd_grid_norm(function.values, volume_element);
        if (norm <= norm_tolerance)
        {
            throw std::runtime_error(
                "Sternheimer grid-function orthonormalization found a linearly dependent state.");
        }
        scale_grid_function(function, Complex(1.0 / norm, 0.0));
        orthonormal.push_back(std::move(function));
    }
    return orthonormal;
}

std::vector<std::array<int, 3>> enumerate_delta_sternheimer_periodic_images(
    const SternheimerFDHamiltonian::Grid& grid,
    const std::array<double, 3>& atom_position,
    const double cutoff_radius)
{
    if (grid.nx <= 0 || grid.ny <= 0 || grid.nz <= 0 || grid.hx <= 0.0 || grid.hy <= 0.0
        || grid.hz <= 0.0)
    {
        throw std::invalid_argument("Sternheimer delta periodic images require a valid uniform grid.");
    }
    if (!std::isfinite(cutoff_radius) || cutoff_radius < 0.0)
    {
        throw std::invalid_argument("Sternheimer delta periodic images require a non-negative finite cutoff.");
    }
    for (const double coordinate: atom_position)
    {
        if (!std::isfinite(coordinate))
        {
            throw std::invalid_argument("Sternheimer delta periodic images require finite atom coordinates.");
        }
    }
    if (!grid.periodic)
    {
        return {{0, 0, 0}};
    }

    const SternheimerFDLatticeVectors dual = sternheimer_fd_grid_dual_vectors(grid);
    std::array<int, 3> image_min{};
    std::array<int, 3> image_max{};
    for (int direction = 0; direction != 3; ++direction)
    {
        double reduced_position = 0.0;
        double dual_norm_squared = 0.0;
        for (int component = 0; component != 3; ++component)
        {
            reduced_position += dual[direction][component] * atom_position[component];
            dual_norm_squared += dual[direction][component] * dual[direction][component];
        }
        reduced_position -= std::floor(reduced_position);
        const double reduced_cutoff = std::sqrt(dual_norm_squared) * cutoff_radius;
        image_min[static_cast<std::size_t>(direction)]
            = static_cast<int>(std::ceil(-reduced_position - reduced_cutoff));
        image_max[static_cast<std::size_t>(direction)]
            = static_cast<int>(std::floor(1.0 - reduced_position + reduced_cutoff));
    }

    std::vector<std::array<int, 3>> images;
    for (int ix = image_min[0]; ix <= image_max[0]; ++ix)
    {
        for (int iy = image_min[1]; iy <= image_max[1]; ++iy)
        {
            for (int iz = image_min[2]; iz <= image_max[2]; ++iz)
            {
                images.push_back({ix, iy, iz});
            }
        }
    }
    return images;
}

std::vector<std::vector<Complex>> recover_delta_grid_ao_coefficients(
    const std::vector<SternheimerDeltaGridFunction>& basis,
    const std::vector<const Vector*>& states,
    const double volume_element,
    const double reconstruction_tolerance,
    double& maximum_relative_error,
    double* gram_reciprocal_condition)
{
    if (basis.empty() || states.empty() || !std::isfinite(volume_element) || volume_element <= 0.
        || !std::isfinite(reconstruction_tolerance) || reconstruction_tolerance <= 0.)
    {
        throw std::invalid_argument("Invalid AO coordinate reconstruction request.");
    }
    const int nao = static_cast<int>(basis.size());
    const int nstates = static_cast<int>(states.size());
    const std::size_t grid_size = basis.front().values.size();
    if (grid_size == 0)
    {
        throw std::invalid_argument("Empty AO reconstruction grid.");
    }
    for (const auto& function: basis)
    {
        if (function.values.size() != grid_size)
        {
            throw std::invalid_argument("AO reconstruction grid dimensions differ.");
        }
    }
    for (const auto* state: states)
    {
        if (state == nullptr || state->size() != grid_size)
        {
            throw std::invalid_argument("AO reconstruction state dimensions differ.");
        }
    }
    std::vector<Complex> gram(basis.size() * basis.size(), 0.);
    accumulate_grid_component_product(basis, -1, basis, -1, nullptr, Complex(volume_element, 0.), gram);
    for (const auto value: gram)
    {
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
        {
            throw std::runtime_error("Non-finite sampled AO Gram matrix.");
        }
    }
    // A small wavefunction residual cannot bound coordinates in a nearly null AO direction.
    auto gram_eigenvectors = gram;
    const auto gram_eigenvalues = diagonalize_complex_hermitian(gram_eigenvectors, nao);
    const double reciprocal_condition = gram_eigenvalues.front() / gram_eigenvalues.back();
    if (!std::isfinite(reciprocal_condition) || gram_eigenvalues.back() <= 0. || reciprocal_condition < 1.e-8)
    {
        throw std::runtime_error("Sampled AO Gram matrix is ill-conditioned (rcond < 1e-8); "
                                 "refine the grid or remove near-dependent orbitals.");
    }
    if (gram_reciprocal_condition != nullptr)
    {
        *gram_reciprocal_condition = reciprocal_condition;
    }
    std::vector<Complex> coordinates(basis.size() * states.size(), 0.);
    std::vector<Complex> packed_basis, packed_states;
    for (std::size_t begin = 0; begin < grid_size; begin += delta_grid_block_size)
    {
        const int count = static_cast<int>(std::min(delta_grid_block_size, grid_size - begin));
        pack_grid_function_component(basis, -1, begin, count, packed_basis);
        pack_vector_component(states, begin, count, packed_states);
        BlasConnector::gemm_cm('C',
                               'N',
                               nao,
                               nstates,
                               count,
                               Complex(volume_element, 0.),
                               packed_basis.data(),
                               count,
                               packed_states.data(),
                               count,
                               Complex(1., 0.),
                               coordinates.data(),
                               nao);
    }
    std::vector<int> pivots(basis.size());
    int info = 0;
    zgesv_(&nao, &nstates, gram.data(), &nao, pivots.data(), coordinates.data(), &nao, &info);
    if (info != 0)
    {
        throw std::runtime_error("AO coordinate reconstruction Gram solve failed.");
    }
    std::vector<double> errors(states.size(), 0.), norms(states.size(), 0.);
    for (std::size_t begin = 0; begin < grid_size; begin += delta_grid_block_size)
    {
        const int count = static_cast<int>(std::min(delta_grid_block_size, grid_size - begin));
        pack_grid_function_component(basis, -1, begin, count, packed_basis);
        pack_vector_component(states, begin, count, packed_states);
        std::vector<Complex> reconstructed(static_cast<std::size_t>(count) * states.size());
        BlasConnector::gemm_cm('N',
                               'N',
                               count,
                               nstates,
                               nao,
                               Complex(1., 0.),
                               packed_basis.data(),
                               count,
                               coordinates.data(),
                               nao,
                               Complex(0., 0.),
                               reconstructed.data(),
                               count);
        for (int state = 0; state < nstates; ++state)
        {
            for (int ir = 0; ir < count; ++ir)
            {
                const std::size_t index = static_cast<std::size_t>(state) * count + ir;
                errors[state] += std::norm(reconstructed[index] - packed_states[index]);
                norms[state] += std::norm(packed_states[index]);
            }
        }
    }
    maximum_relative_error = 0.;
    std::vector<std::vector<Complex>> result(states.size());
    for (int state = 0; state < nstates; ++state)
    {
        const double error = std::sqrt(errors[state] / norms[state]);
        if (!std::isfinite(error) || error > reconstruction_tolerance)
        {
            throw std::runtime_error("Delta state cannot be reconstructed in the supplied AO span.");
        }
        maximum_relative_error = std::max(maximum_relative_error, error);
        result[state].assign(coordinates.begin() + state * nao, coordinates.begin() + (state + 1) * nao);
    }
    return result;
}

SternheimerDeltaCoefficientComponents solve_delta_sternheimer_subspace_coefficients(
    const std::vector<SternheimerFDHamiltonian::Complex>& hamiltonian_matrix,
    const std::vector<SternheimerFDHamiltonian::Complex>& overlap_matrix,
    const std::vector<SternheimerFDHamiltonian::Complex>& rhs,
    const std::vector<SternheimerFDHamiltonian::Complex>& hamiltonian_out_coupling,
    const SternheimerFDHamiltonian::Complex shift)
{
    const int size = static_cast<int>(rhs.size());
    if (hamiltonian_out_coupling.size() != rhs.size())
    {
        throw std::invalid_argument("Sternheimer delta subspace coupling size does not match the right-hand side.");
    }
    const std::size_t matrix_size = static_cast<std::size_t>(size) * static_cast<std::size_t>(size);
    if (hamiltonian_matrix.size() != matrix_size || overlap_matrix.size() != matrix_size)
    {
        throw std::invalid_argument("Sternheimer delta subspace matrices do not match the right-hand side size.");
    }

    SternheimerDeltaCoefficientComponents result;
    result.sos.assign(rhs.size(), Complex(0.0, 0.0));
    result.pulay.assign(rhs.size(), Complex(0.0, 0.0));
    result.total.assign(rhs.size(), Complex(0.0, 0.0));
    if (size == 0)
    {
        return result;
    }

    std::vector<Complex> shifted_operator(matrix_size, Complex(0.0, 0.0));
    for (std::size_t index = 0; index != matrix_size; ++index)
    {
        shifted_operator[index] = hamiltonian_matrix[index] - shift * overlap_matrix[index];
    }
    constexpr int right_hand_side_count = 2;
    std::vector<Complex> solutions(static_cast<std::size_t>(right_hand_side_count) * rhs.size(), Complex(0.0, 0.0));
    for (int index = 0; index != size; ++index)
    {
        solutions[static_cast<std::size_t>(index)] = rhs[static_cast<std::size_t>(index)];
        solutions[static_cast<std::size_t>(size + index)] = -hamiltonian_out_coupling[static_cast<std::size_t>(index)];
    }

    std::vector<int> pivots(static_cast<std::size_t>(size), 0);
    int info = 0;
    zgesv_(&size,
           &right_hand_side_count,
           shifted_operator.data(),
           &size,
           pivots.data(),
           solutions.data(),
           &size,
           &info);
    if (info < 0)
    {
        throw std::runtime_error("Sternheimer delta subspace solve received an invalid LAPACK argument.");
    }
    if (info > 0)
    {
        throw std::runtime_error("Sternheimer delta subspace shifted matrix is singular.");
    }

    for (int index = 0; index != size; ++index)
    {
        const std::size_t coefficient_index = static_cast<std::size_t>(index);
        result.sos[coefficient_index] = solutions[coefficient_index];
        result.pulay[coefficient_index] = solutions[static_cast<std::size_t>(size + index)];
        result.total[coefficient_index] = result.sos[coefficient_index] + result.pulay[coefficient_index];
    }
    return result;
}

SternheimerDeltaPostprocessResult postprocess_delta_sternheimer_solution(
    const SternheimerFDHamiltonian::Vector& standard_delta_wavefunction,
    const SternheimerDeltaPostprocessInput& input)
{
    validate_postprocess_input(standard_delta_wavefunction, input);

    auto dot = [&input](const Vector& lhs, const Vector& rhs) {
        return sternheimer_fd_grid_dot(lhs, rhs, input.volume_element);
    };

    SternheimerDeltaPostprocessResult result;
    result.out_wavefunction = standard_delta_wavefunction;
    SternheimerRPA::project_out_subspace(input.occupied_wavefunctions, dot, result.out_wavefunction);
    const std::vector<Vector> virtual_orbitals = collect_virtual_orbitals(input.virtual_states);
    SternheimerRPA::project_out_subspace(virtual_orbitals, dot, result.out_wavefunction);

    const SternheimerDeltaCoefficientComponents components
        = compute_delta_coefficient_components(input.virtual_states,
                                               input.perturbation_matrix_elements,
                                               result.out_wavefunction,
                                               input.occupied_eigenvalue,
                                               input.omega,
                                               input.volume_element);
    result.sos_coefficients = components.sos;
    result.pulay_coefficients = components.pulay;
    result.coefficients = components.total;
    assemble_delta_wavefunction_components(input.virtual_states, result);
    result.out_norm = sternheimer_fd_grid_norm(result.out_wavefunction, input.volume_element);
    result.reconstruction_error
        = difference_norm(result.reconstructed_wavefunction, standard_delta_wavefunction, input.volume_element);
    return result;
}

SternheimerDeltaPulayOperatorComponents decompose_delta_sternheimer_pulay_operator_terms(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<double>& fixed_local_potential,
    const std::vector<Vector>& occupied_wavefunctions,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const Vector& out_wavefunction,
    const double occupied_eigenvalue,
    const double omega,
    const double volume_element)
{
    const std::size_t grid_size = static_cast<std::size_t>(hamiltonian.grid().size());
    if (fixed_local_potential.size() != grid_size || out_wavefunction.size() != grid_size)
    {
        throw std::invalid_argument(
            "Sternheimer Pulay operator decomposition requires full-grid potentials and wavefunctions.");
    }
    if (volume_element <= 0.0)
    {
        throw std::invalid_argument(
            "Sternheimer Pulay operator decomposition requires a positive grid volume element.");
    }
    validate_virtual_states(virtual_states, grid_size);
    for (const Vector& occupied: occupied_wavefunctions)
    {
        check_vector_size(occupied, grid_size, "Sternheimer Pulay operator occupied state");
    }

    SternheimerDeltaPulayOperatorComponents result;
    Vector* outputs[] = {&result.kinetic,
                         &result.fixed_local,
                         &result.hxc_local,
                         &result.nonlocal,
                         &result.eigenvalue,
                         &result.total};
    for (Vector* output: outputs)
    {
        output->assign(grid_size, Complex(0.0, 0.0));
    }

    const auto dot = [volume_element](const Vector& lhs, const Vector& rhs) {
        return sternheimer_fd_grid_dot(lhs, rhs, volume_element);
    };
    std::vector<Vector> fixed_subspace = occupied_wavefunctions;
    const std::vector<Vector> virtual_orbitals = collect_virtual_orbitals(virtual_states);
    fixed_subspace.insert(fixed_subspace.end(), virtual_orbitals.begin(), virtual_orbitals.end());
    const std::vector<double>& full_local_potential = hamiltonian.local_potential();

    for (const SternheimerDeltaVirtualState& state: virtual_states)
    {
        const Complex denominator(occupied_eigenvalue - state.eigenvalue, -omega);
        if (std::abs(denominator) < 1.0e-30)
        {
            throw std::runtime_error(
                "Sternheimer Pulay operator decomposition found a singular denominator.");
        }

        std::array<Vector, 5> residual_terms;
        for (Vector& term: residual_terms)
        {
            term.assign(grid_size, Complex(0.0, 0.0));
        }
        hamiltonian.apply_kinetic(state.orbital, residual_terms[0]);
        hamiltonian.apply_nonlocal(state.orbital, residual_terms[3]);
#pragma omp parallel for schedule(static)
        for (std::size_t ir = 0; ir != grid_size; ++ir)
        {
            residual_terms[1][ir] = fixed_local_potential[ir] * state.orbital[ir];
            residual_terms[2][ir]
                = (full_local_potential[ir] - fixed_local_potential[ir]) * state.orbital[ir];
            residual_terms[4][ir] = -state.eigenvalue * state.orbital[ir];
        }

        Vector* component_outputs[] = {&result.kinetic,
                                       &result.fixed_local,
                                       &result.hxc_local,
                                       &result.nonlocal,
                                       &result.eigenvalue};
        for (std::size_t iterm = 0; iterm != residual_terms.size(); ++iterm)
        {
            SternheimerRPA::project_out_subspace(fixed_subspace, dot, residual_terms[iterm]);
            const Complex coefficient = dot(residual_terms[iterm], out_wavefunction) / denominator;
            axpy(coefficient, state.orbital, *component_outputs[iterm]);
        }

        const Complex total_coefficient = dot(state.residual, out_wavefunction) / denominator;
        axpy(total_coefficient, state.orbital, result.total);
    }
    return result;
}

SternheimerDeltaGridMatrices assemble_delta_sternheimer_grid_matrices_fast(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerDeltaGridFunction>& basis_functions,
    const double volume_element)
{
    if (!std::isfinite(volume_element) || volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer reference delta matrices require a finite positive grid volume element.");
    }
    for (const auto& function: basis_functions)
    {
        validate_grid_function(function, hamiltonian.grid().size(), "Sternheimer reference delta basis");
    }
    return assemble_delta_sternheimer_grid_matrices_blocked(hamiltonian, basis_functions, volume_element);
}

SternheimerDeltaGridMatrices assemble_delta_sternheimer_cross_matrices_fast(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerDeltaGridFunction>& left_functions,
    const std::vector<SternheimerDeltaGridFunction>& right_functions,
    const double volume_element)
{
    if (!std::isfinite(volume_element) || volume_element <= 0.0
        || !std::isfinite(hamiltonian.kinetic_prefactor()))
    {
        throw std::invalid_argument(
            "Sternheimer cross matrices require a finite positive grid volume element and finite kinetic prefactor.");
    }
    if (!std::isfinite(volume_element * hamiltonian.kinetic_prefactor()))
    {
        throw std::overflow_error("Sternheimer cross matrix kinetic quadrature weight overflow.");
    }
    const std::size_t grid_size = static_cast<std::size_t>(hamiltonian.grid().size());
    for (const auto* functions: {&left_functions, &right_functions})
    {
        const std::string context = functions == &left_functions
                                        ? "Sternheimer cross matrix left basis"
                                        : "Sternheimer cross matrix right basis";
        for (const auto& function: *functions)
        {
            validate_grid_function(function, grid_size, context);
            for (int component = -1; component != 3; ++component)
            {
                validate_finite_grid_vector(grid_function_component(function, component), context);
            }
        }
    }
    for (const double value: hamiltonian.local_potential())
    {
        if (!std::isfinite(value))
        {
            throw std::invalid_argument("Sternheimer cross matrix local potential contains a nonfinite value.");
        }
    }
    const SternheimerFDNonlocalProjector* projector = hamiltonian.nonlocal_projector();
    if (projector != nullptr)
    {
        if (!std::isfinite(projector->volume_element()) || projector->volume_element() <= 0.0
            || std::abs(projector->volume_element() - volume_element)
                   > 1.0e-12 * std::max(projector->volume_element(), volume_element))
        {
            throw std::invalid_argument(
                "Sternheimer cross matrix nonlocal projector uses an invalid or different grid volume element.");
        }
        for (const auto& block: projector->blocks())
        {
            for (const Vector& values: block.projectors)
            {
                validate_finite_grid_vector(values, "Sternheimer cross matrix nonlocal projector");
            }
            for (const Vector& row: block.d_matrix)
            {
                validate_finite_grid_vector(row, "Sternheimer cross matrix nonlocal D");
            }
        }
    }
    SternheimerDeltaGridMatrices matrices = assemble_delta_sternheimer_cross_matrices_blocked(
        hamiltonian, left_functions, right_functions, volume_element);
    for (const Vector* component: {&matrices.overlap, &matrices.kinetic, &matrices.local_potential,
                                   &matrices.nonlocal, &matrices.hamiltonian})
    {
        for (const Complex value: *component)
        {
            if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            {
                throw std::overflow_error("Sternheimer cross matrix quadrature produced a nonfinite result.");
            }
        }
    }
    return matrices;
}

SternheimerDeltaGridMatrices assemble_delta_sternheimer_grid_matrices(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerDeltaGridFunction>& basis_functions,
    const double volume_element)
{
    if (volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer reference delta matrices require a positive grid volume element.");
    }
    const std::size_t grid_size = static_cast<std::size_t>(hamiltonian.grid().size());
    for (const SternheimerDeltaGridFunction& function: basis_functions)
    {
        validate_grid_function(function, grid_size, "Sternheimer reference delta basis");
    }

    const SternheimerFDNonlocalProjector* nonlocal_projector = hamiltonian.nonlocal_projector();
    if (nonlocal_projector != nullptr
        && std::abs(nonlocal_projector->volume_element() - volume_element)
               > 1.0e-12 * std::max(nonlocal_projector->volume_element(), volume_element))
    {
        throw std::invalid_argument(
            "Sternheimer reference delta nonlocal projector uses a different grid volume element.");
    }

    const int basis_size = static_cast<int>(basis_functions.size());
    const std::size_t matrix_size
        = static_cast<std::size_t>(basis_size) * static_cast<std::size_t>(basis_size);
    SternheimerDeltaGridMatrices matrices;
    matrices.row_count = basis_functions.size();
    matrices.column_count = basis_functions.size();
    matrices.overlap.assign(matrix_size, Complex(0.0, 0.0));
    matrices.kinetic.assign(matrix_size, Complex(0.0, 0.0));
    matrices.local_potential.assign(matrix_size, Complex(0.0, 0.0));
    matrices.nonlocal.assign(matrix_size, Complex(0.0, 0.0));
    matrices.hamiltonian.assign(matrix_size, Complex(0.0, 0.0));
    if (basis_size == 0)
    {
        return matrices;
    }

    std::vector<Vector> nonlocal_basis(static_cast<std::size_t>(basis_size));
    if (nonlocal_projector != nullptr)
    {
        for (int ib = 0; ib != basis_size; ++ib)
        {
            nonlocal_projector->apply(basis_functions[static_cast<std::size_t>(ib)].values,
                                      nonlocal_basis[static_cast<std::size_t>(ib)]);
        }
    }

    const std::vector<double>& local_potential = hamiltonian.local_potential();
    const double kinetic_prefactor = hamiltonian.kinetic_prefactor();
    for (int ib = 0; ib != basis_size; ++ib)
    {
        const SternheimerDeltaGridFunction& ket = basis_functions[static_cast<std::size_t>(ib)];
        for (int ia = 0; ia != basis_size; ++ia)
        {
            const SternheimerDeltaGridFunction& bra = basis_functions[static_cast<std::size_t>(ia)];
            Complex overlap(0.0, 0.0);
            Complex kinetic_element(0.0, 0.0);
            Complex local_potential_element(0.0, 0.0);
            Complex nonlocal_element(0.0, 0.0);
            for (std::size_t ir = 0; ir != grid_size; ++ir)
            {
                overlap += volume_element * std::conj(bra.values[ir]) * ket.values[ir];
                Complex gradient_dot(0.0, 0.0);
                for (int direction = 0; direction != 3; ++direction)
                {
                    gradient_dot += std::conj(bra.gradients[static_cast<std::size_t>(direction)][ir])
                                    * ket.gradients[static_cast<std::size_t>(direction)][ir];
                }
                kinetic_element += volume_element * kinetic_prefactor * gradient_dot;
                local_potential_element += volume_element * std::conj(bra.values[ir])
                                           * local_potential[ir] * ket.values[ir];
            }
            if (nonlocal_projector != nullptr)
            {
                nonlocal_element = sternheimer_fd_grid_dot(
                    bra.values, nonlocal_basis[static_cast<std::size_t>(ib)], volume_element);
            }
            const std::size_t index
                = static_cast<std::size_t>(ia) + static_cast<std::size_t>(basis_size) * static_cast<std::size_t>(ib);
            matrices.overlap[index] = overlap;
            matrices.kinetic[index] = kinetic_element;
            matrices.local_potential[index] = local_potential_element;
            matrices.nonlocal[index] = nonlocal_element;
            matrices.hamiltonian[index] = kinetic_element + local_potential_element + nonlocal_element;
        }
    }
    return matrices;
}

SternheimerDeltaSubspace build_reference_delta_sternheimer_subspace(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerDeltaGridFunction>& occupied_functions,
    std::vector<SternheimerDeltaGridFunction> candidate_functions,
    const double volume_element,
    const SternheimerDeltaSubspaceOptions& options)
{
    if (volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer reference delta subspace requires a positive grid volume element.");
    }
    if (options.max_virtual_states < 0)
    {
        throw std::invalid_argument("Sternheimer reference delta max_virtual_states must be non-negative.");
    }
    if (options.norm_tolerance < 0.0)
    {
        throw std::invalid_argument("Sternheimer reference delta norm_tolerance must be non-negative.");
    }

    const std::size_t grid_size = static_cast<std::size_t>(hamiltonian.grid().size());
    for (const SternheimerDeltaGridFunction& occupied: occupied_functions)
    {
        validate_grid_function(occupied, grid_size, "Sternheimer reference delta occupied function");
    }
    for (const SternheimerDeltaGridFunction& candidate: candidate_functions)
    {
        validate_grid_function(candidate, grid_size, "Sternheimer reference delta candidate function");
    }

    const int input_candidate_count = static_cast<int>(candidate_functions.size());
    const bool complete_requested_space
        = options.max_virtual_states == 0 || options.max_virtual_states >= input_candidate_count;
    if (options.use_block_generalized_eigensolver && complete_requested_space
        && input_candidate_count > 0)
    {
        return build_complete_reference_delta_sternheimer_subspace_blocked(
            hamiltonian,
            occupied_functions,
            std::move(candidate_functions),
            volume_element,
            options);
    }

    auto dot = [volume_element](const Vector& lhs, const Vector& rhs) {
        return sternheimer_fd_grid_dot(lhs, rhs, volume_element);
    };
    std::vector<SternheimerDeltaGridFunction> residual_candidates = std::move(candidate_functions);
    for (SternheimerDeltaGridFunction& candidate: residual_candidates)
    {
        for (int pass = 0; pass != 2; ++pass)
        {
            for (const SternheimerDeltaGridFunction& occupied: occupied_functions)
            {
                axpy_grid_function(-dot(occupied.values, candidate.values), occupied, candidate);
            }
        }
    }

    const int candidate_limit = options.max_virtual_states > 0
                                    ? std::min(options.max_virtual_states,
                                               static_cast<int>(residual_candidates.size()))
                                    : static_cast<int>(residual_candidates.size());
    std::vector<bool> selected(residual_candidates.size(), false);
    std::vector<SternheimerDeltaGridFunction> orthonormal_candidates;
    orthonormal_candidates.reserve(static_cast<std::size_t>(candidate_limit));
    while (static_cast<int>(orthonormal_candidates.size()) < candidate_limit)
    {
        int pivot = -1;
        double pivot_norm = options.norm_tolerance;
        for (std::size_t candidate_index = 0; candidate_index != residual_candidates.size(); ++candidate_index)
        {
            if (selected[candidate_index])
            {
                continue;
            }
            const double norm
                = sternheimer_fd_grid_norm(residual_candidates[candidate_index].values, volume_element);
            if (norm > pivot_norm)
            {
                pivot = static_cast<int>(candidate_index);
                pivot_norm = norm;
            }
        }
        if (pivot < 0)
        {
            break;
        }

        const std::size_t pivot_index = static_cast<std::size_t>(pivot);
        selected[pivot_index] = true;
        SternheimerDeltaGridFunction accepted = std::move(residual_candidates[pivot_index]);
        scale_grid_function(accepted, Complex(1.0 / pivot_norm, 0.0));
        orthonormal_candidates.push_back(std::move(accepted));
        const SternheimerDeltaGridFunction& newest = orthonormal_candidates.back();
        for (std::size_t candidate_index = 0; candidate_index != residual_candidates.size(); ++candidate_index)
        {
            if (selected[candidate_index])
            {
                continue;
            }
            for (int pass = 0; pass != 2; ++pass)
            {
                SternheimerDeltaGridFunction& candidate = residual_candidates[candidate_index];
                axpy_grid_function(-dot(newest.values, candidate.values), newest, candidate);
            }
        }
    }

    SternheimerDeltaSubspace subspace;
    subspace.accepted_candidates = static_cast<int>(orthonormal_candidates.size());
    subspace.discarded_candidates
        = input_candidate_count - subspace.accepted_candidates;
    if (orthonormal_candidates.empty())
    {
        return subspace;
    }

    const SternheimerDeltaGridMatrices matrices
        = assemble_delta_sternheimer_grid_matrices(hamiltonian, orthonormal_candidates, volume_element);
    const int basis_size = static_cast<int>(orthonormal_candidates.size());
    std::vector<Complex> eigenvectors = matrices.hamiltonian;
    const std::vector<double> eigenvalues = diagonalize_complex_hermitian(eigenvectors, basis_size);

    subspace.grid_functions.reserve(static_cast<std::size_t>(basis_size));
    subspace.virtual_states.reserve(static_cast<std::size_t>(basis_size));
    for (int ia = 0; ia != basis_size; ++ia)
    {
        SternheimerDeltaGridFunction eigenfunction;
        eigenfunction.values.assign(grid_size, Complex(0.0, 0.0));
        for (Vector& gradient: eigenfunction.gradients)
        {
            gradient.assign(grid_size, Complex(0.0, 0.0));
        }
        for (int ib = 0; ib != basis_size; ++ib)
        {
            const Complex coefficient
                = eigenvectors[static_cast<std::size_t>(ib)
                               + static_cast<std::size_t>(basis_size) * static_cast<std::size_t>(ia)];
            axpy_grid_function(coefficient,
                               orthonormal_candidates[static_cast<std::size_t>(ib)],
                               eigenfunction);
        }

        SternheimerDeltaVirtualState state;
        state.orbital = options.retain_grid_functions ? eigenfunction.values
                                                      : std::move(eigenfunction.values);
        state.eigenvalue = eigenvalues[static_cast<std::size_t>(ia)];
        if (options.retain_grid_functions)
        {
            subspace.grid_functions.push_back(std::move(eigenfunction));
        }
        subspace.virtual_states.push_back(std::move(state));
    }

    if (options.evaluate_full_grid_difference)
    {
        evaluate_full_grid_delta_hamiltonian_difference(hamiltonian, volume_element, subspace);
    }
    else
    {
        subspace.full_grid_hamiltonian_relative_difference = -1.0;
        subspace.full_grid_hamiltonian_max_abs_difference = -1.0;
    }

    std::vector<Vector> residual_projector;
    residual_projector.reserve(occupied_functions.size() + subspace.virtual_states.size());
    for (const SternheimerDeltaGridFunction& occupied: occupied_functions)
    {
        residual_projector.push_back(occupied.values);
    }
    for (const SternheimerDeltaVirtualState& state: subspace.virtual_states)
    {
        residual_projector.push_back(state.orbital);
    }
    for (SternheimerDeltaVirtualState& state: subspace.virtual_states)
    {
        hamiltonian.apply(state.orbital, state.residual);
        for (std::size_t ir = 0; ir != state.residual.size(); ++ir)
        {
            state.residual[ir] -= state.eigenvalue * state.orbital[ir];
        }
        SternheimerRPA::project_out_subspace(residual_projector, dot, state.residual);
    }
    return subspace;
}

SternheimerDeltaSubspace build_delta_sternheimer_subspace(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerFDHamiltonian::Vector>& occupied_wavefunctions,
    std::vector<SternheimerFDHamiltonian::Vector> candidate_orbitals,
    const double volume_element,
    const SternheimerDeltaSubspaceOptions& options)
{
    const int grid_size = hamiltonian.grid().size();
    if (volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer delta subspace requires a positive grid volume element.");
    }
    if (options.max_virtual_states < 0)
    {
        throw std::invalid_argument("Sternheimer delta subspace max_virtual_states must be non-negative.");
    }
    if (options.norm_tolerance < 0.0)
    {
        throw std::invalid_argument("Sternheimer delta subspace norm_tolerance must be non-negative.");
    }
    for (const Vector& occupied: occupied_wavefunctions)
    {
        check_vector_size(occupied, static_cast<std::size_t>(grid_size), "Sternheimer delta occupied state");
    }

    auto dot = [volume_element](const Vector& lhs, const Vector& rhs) {
        return sternheimer_fd_grid_dot(lhs, rhs, volume_element);
    };

    std::vector<Vector> orthonormal_candidates;
    const int input_candidate_count = static_cast<int>(candidate_orbitals.size());
    for (Vector& candidate: candidate_orbitals)
    {
        check_vector_size(candidate, static_cast<std::size_t>(grid_size), "Sternheimer delta candidate orbital");
        SternheimerRPA::project_out_subspace(occupied_wavefunctions, dot, candidate);
        SternheimerRPA::project_out_subspace(orthonormal_candidates, dot, candidate);
        const double norm = sternheimer_fd_grid_norm(candidate, volume_element);
        if (norm <= options.norm_tolerance)
        {
            continue;
        }
        scale_vector(candidate, Complex(1.0 / norm, 0.0));
        orthonormal_candidates.push_back(std::move(candidate));
        if (options.max_virtual_states > 0
            && static_cast<int>(orthonormal_candidates.size()) >= options.max_virtual_states)
        {
            break;
        }
    }

    SternheimerDeltaSubspace subspace;
    subspace.accepted_candidates = static_cast<int>(orthonormal_candidates.size());
    subspace.discarded_candidates
        = input_candidate_count - subspace.accepted_candidates;
    if (orthonormal_candidates.empty())
    {
        return subspace;
    }

    const int nvirtual = static_cast<int>(orthonormal_candidates.size());
    std::vector<Vector> h_candidates(nvirtual);
    for (int j = 0; j != nvirtual; ++j)
    {
        hamiltonian.apply(orthonormal_candidates[j], h_candidates[j]);
    }

    std::vector<Complex> h_matrix(static_cast<std::size_t>(nvirtual) * static_cast<std::size_t>(nvirtual),
                                  Complex(0.0, 0.0));
    for (int j = 0; j != nvirtual; ++j)
    {
        for (int i = 0; i != nvirtual; ++i)
        {
            h_matrix[static_cast<std::size_t>(i) + static_cast<std::size_t>(nvirtual) * static_cast<std::size_t>(j)]
                = dot(orthonormal_candidates[i], h_candidates[j]);
        }
    }
    const std::vector<double> eigenvalues = diagonalize_complex_hermitian(h_matrix, nvirtual);

    subspace.virtual_states.reserve(static_cast<std::size_t>(nvirtual));
    for (int ia = 0; ia != nvirtual; ++ia)
    {
        Vector eta(static_cast<std::size_t>(grid_size), Complex(0.0, 0.0));
        for (int j = 0; j != nvirtual; ++j)
        {
            const Complex coefficient
                = h_matrix[static_cast<std::size_t>(j) + static_cast<std::size_t>(nvirtual) * static_cast<std::size_t>(ia)];
            axpy(coefficient, orthonormal_candidates[j], eta);
        }

        SternheimerRPA::project_out_subspace(occupied_wavefunctions, dot, eta);
        const double eta_norm = sternheimer_fd_grid_norm(eta, volume_element);
        if (eta_norm <= options.norm_tolerance)
        {
            ++subspace.discarded_candidates;
            continue;
        }
        scale_vector(eta, Complex(1.0 / eta_norm, 0.0));

        Vector residual;
        hamiltonian.apply(eta, residual);
        for (std::size_t ir = 0; ir != residual.size(); ++ir)
        {
            residual[ir] -= eigenvalues[static_cast<std::size_t>(ia)] * eta[ir];
        }

        SternheimerDeltaVirtualState state;
        state.orbital = std::move(eta);
        state.eigenvalue = eigenvalues[static_cast<std::size_t>(ia)];
        subspace.virtual_states.push_back(std::move(state));
    }

    if (options.evaluate_full_grid_difference)
    {
        evaluate_full_grid_delta_hamiltonian_difference(hamiltonian, volume_element, subspace);
    }
    else
    {
        subspace.full_grid_hamiltonian_relative_difference = -1.0;
        subspace.full_grid_hamiltonian_max_abs_difference = -1.0;
    }

    const std::vector<Vector> virtual_orbitals = collect_virtual_orbitals(subspace.virtual_states);
    std::vector<Vector> residual_projector = occupied_wavefunctions;
    residual_projector.insert(residual_projector.end(), virtual_orbitals.begin(), virtual_orbitals.end());
    for (SternheimerDeltaVirtualState& state: subspace.virtual_states)
    {
        hamiltonian.apply(state.orbital, state.residual);
        for (std::size_t ir = 0; ir != state.residual.size(); ++ir)
        {
            state.residual[ir] -= state.eigenvalue * state.orbital[ir];
        }
        SternheimerRPA::project_out_subspace(residual_projector, dot, state.residual);
    }
    return subspace;
}

SternheimerDeltaSubspace build_delta_sternheimer_subspace_by_mode(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerDeltaGridFunction>& occupied_functions,
    std::vector<SternheimerDeltaGridFunction> candidate_functions,
    const double volume_element,
    const SternheimerDeltaSubspaceOptions& options,
    const SternheimerDeltaABlockMode mode)
{
    if (mode == SternheimerDeltaABlockMode::ReferenceValueGradient)
    {
        return build_reference_delta_sternheimer_subspace(
            hamiltonian, occupied_functions, std::move(candidate_functions), volume_element, options);
    }
    if (mode == SternheimerDeltaABlockMode::FullGrid)
    {
        std::vector<Vector> occupied_values;
        occupied_values.reserve(occupied_functions.size());
        for (const SternheimerDeltaGridFunction& function: occupied_functions)
        {
            occupied_values.push_back(function.values);
        }
        std::vector<Vector> candidate_values;
        candidate_values.reserve(candidate_functions.size());
        for (SternheimerDeltaGridFunction& function: candidate_functions)
        {
            candidate_values.push_back(std::move(function.values));
        }
        return build_delta_sternheimer_subspace(
            hamiltonian, occupied_values, std::move(candidate_values), volume_element, options);
    }
    throw std::invalid_argument("Invalid Sternheimer Delta A-block mode.");
}

std::vector<SternheimerFDHamiltonian::Complex> delta_sternheimer_perturbation_matrix_elements(
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const std::vector<double>& perturbation_potential,
    const SternheimerFDHamiltonian::Vector& occupied_wavefunction,
    const double volume_element)
{
    if (volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer delta perturbation elements require a positive grid volume element.");
    }
    if (perturbation_potential.size() != occupied_wavefunction.size())
    {
        throw std::invalid_argument("Sternheimer delta perturbation potential size does not match wavefunction.");
    }
    validate_virtual_states(virtual_states, occupied_wavefunction.size());

    std::vector<Complex> elements(virtual_states.size(), Complex(0.0, 0.0));
#pragma omp parallel for schedule(static)
    for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
    {
        Complex value(0.0, 0.0);
        for (std::size_t ir = 0; ir != occupied_wavefunction.size(); ++ir)
        {
            value += std::conj(virtual_states[ia].orbital[ir]) * perturbation_potential[ir]
                * occupied_wavefunction[ir];
        }
        elements[ia] = volume_element * value;
    }
    return elements;
}

std::vector<SternheimerFDHamiltonian::Complex> delta_sternheimer_perturbation_matrix_elements(
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const SternheimerFDHamiltonian::Vector& perturbation_potential,
    const SternheimerFDHamiltonian::Vector& occupied_wavefunction,
    const double volume_element)
{
    if (volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer delta perturbation elements require a positive grid volume element.");
    }
    if (perturbation_potential.size() != occupied_wavefunction.size())
    {
        throw std::invalid_argument("Sternheimer delta perturbation potential size does not match wavefunction.");
    }
    validate_virtual_states(virtual_states, occupied_wavefunction.size());

    std::vector<Complex> elements(virtual_states.size(), Complex(0.0, 0.0));
#pragma omp parallel for schedule(static)
    for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
    {
        Complex value(0.0, 0.0);
        for (std::size_t ir = 0; ir != occupied_wavefunction.size(); ++ir)
        {
            value += std::conj(virtual_states[ia].orbital[ir]) * perturbation_potential[ir]
                * occupied_wavefunction[ir];
        }
        elements[ia] = volume_element * value;
    }
    return elements;
}

SternheimerFDHamiltonian::Vector build_delta_sternheimer_sos_wavefunction(
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const std::vector<SternheimerFDHamiltonian::Complex>& perturbation_matrix_elements,
    const double occupied_eigenvalue,
    const double omega)
{
    if (virtual_states.empty())
    {
        throw std::invalid_argument("Sternheimer direct SOS requires at least one virtual state.");
    }
    if (virtual_states.size() != perturbation_matrix_elements.size())
    {
        throw std::invalid_argument(
            "Sternheimer direct SOS requires one perturbation matrix element per virtual state.");
    }

    const std::size_t grid_size = virtual_states.front().orbital.size();
    if (grid_size == 0)
    {
        throw std::invalid_argument("Sternheimer direct SOS virtual states must be non-empty.");
    }
    Vector response(grid_size, Complex(0.0, 0.0));
    for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
    {
        check_vector_size(virtual_states[ia].orbital, grid_size, "Sternheimer direct SOS virtual state");
        const Complex denominator(occupied_eigenvalue - virtual_states[ia].eigenvalue, -omega);
        if (std::abs(denominator) < 1.0e-30)
        {
            throw std::runtime_error("Sternheimer direct SOS found a singular virtual-state denominator.");
        }
        axpy(perturbation_matrix_elements[ia] / denominator, virtual_states[ia].orbital, response);
    }
    return response;
}

SternheimerDeltaFixedSubspace build_delta_sternheimer_fixed_subspace(
    const std::vector<SternheimerFDHamiltonian::Vector>& occupied_wavefunctions,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states)
{
    SternheimerDeltaFixedSubspace fixed_subspace;
    fixed_subspace.functions = collect_fixed_subspace(occupied_wavefunctions, virtual_states);
    return fixed_subspace;
}

SternheimerDeltaLinearResponse solve_delta_sternheimer_linear_response(
    const SternheimerFDHamiltonian& hamiltonian,
    const SternheimerDeltaFixedSubspace& fixed_subspace,
    const double reference_eigenvalue,
    const SternheimerFDHamiltonian::Vector& rhs,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const std::vector<SternheimerFDHamiltonian::Complex>& perturbation_matrix_elements,
    const double omega,
    const double volume_element,
    const SternheimerRPA::SolverOptions& options)
{
    for (const Vector& function: fixed_subspace.functions)
    {
        check_vector_size(function, static_cast<std::size_t>(hamiltonian.grid().size()),
                          "Sternheimer delta fixed-subspace state");
    }
    const SternheimerSubspaceProjector projector(fixed_subspace.functions, volume_element);
    return solve_delta_sternheimer_linear_response(hamiltonian, projector, reference_eigenvalue,
        rhs, virtual_states, perturbation_matrix_elements, omega, volume_element, options);
}

SternheimerDeltaLinearResponse solve_delta_sternheimer_linear_response(
    const SternheimerFDHamiltonian& hamiltonian,
    const SternheimerSubspaceProjector& fixed_projector,
    const double reference_eigenvalue,
    const SternheimerFDHamiltonian::Vector& rhs,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const std::vector<SternheimerFDHamiltonian::Complex>& perturbation_matrix_elements,
    const double omega,
    const double volume_element,
    const SternheimerRPA::SolverOptions& options)
{
    const int grid_size = hamiltonian.grid().size();
    check_vector_size(rhs, static_cast<std::size_t>(grid_size), "Sternheimer delta rhs");
    if (virtual_states.size() != perturbation_matrix_elements.size())
    {
        throw std::invalid_argument(
            "Sternheimer delta linear response requires one perturbation matrix element per virtual state.");
    }
    if (volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer delta linear response requires a positive grid volume element.");
    }
    validate_virtual_states(virtual_states, static_cast<std::size_t>(grid_size));

    auto dot = [volume_element](const Vector& lhs, const Vector& rhs_vec) {
        return sternheimer_fd_grid_dot(lhs, rhs_vec, volume_element);
    };

    std::vector<Complex> denominators(virtual_states.size(), Complex(0.0, 0.0));
    for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
    {
        denominators[ia] = Complex(reference_eigenvalue - virtual_states[ia].eigenvalue, -omega);
        if (std::abs(denominators[ia]) < 1.0e-30)
        {
            throw std::runtime_error("Sternheimer delta linear response found a singular virtual-state denominator.");
        }
    }

    Vector projected_rhs = rhs;
    fixed_projector.project(projected_rhs);
    for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
    {
        axpy(-perturbation_matrix_elements[ia] / denominators[ia], virtual_states[ia].residual, projected_rhs);
    }

    SternheimerRPA::LinearProblem problem;
    problem.dot = dot;
    problem.apply = [&hamiltonian, &fixed_projector, &virtual_states, &denominators, reference_eigenvalue, omega, dot](
                        const Vector& input,
                        Vector& output) {
        Vector q_input = input;
        fixed_projector.project(q_input);

        hamiltonian.apply(q_input, output);
        const Complex shift(-reference_eigenvalue, omega);
#pragma omp parallel for schedule(static)
        for (std::size_t ir = 0; ir != output.size(); ++ir)
        {
            output[ir] += shift * q_input[ir];
        }
        fixed_projector.project(output);
        for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
        {
            const Complex coupling = dot(virtual_states[ia].residual, q_input) / denominators[ia];
            axpy(coupling, virtual_states[ia].residual, output);
        }
    };
    std::shared_ptr<SternheimerFDSpectralPreconditioner> spectral_preconditioner;
    if (options.use_fd_spectral_preconditioner)
    {
        spectral_preconditioner = get_thread_local_sternheimer_fd_spectral_preconditioner(
            hamiltonian,
            reference_eigenvalue,
            omega,
            options.fd_spectral_preconditioner_regularization);
        problem.precondition = [spectral_preconditioner, &fixed_projector](
                                   const Vector& input,
                                   Vector& output) {
            spectral_preconditioner->apply(input, output);
            fixed_projector.project(output);
        };
    }

    SternheimerDeltaLinearResponse result;
    result.response.out_wavefunction.assign(static_cast<std::size_t>(grid_size), Complex(0.0, 0.0));
    if (fixed_projector.basis_size() >= static_cast<std::size_t>(grid_size))
    {
        const double projected_rhs_norm = sternheimer_fd_grid_norm(projected_rhs, volume_element);
        const double rhs_norm = sternheimer_fd_grid_norm(rhs, volume_element);
        const double zero_tolerance
            = std::max(options.breakdown_tol, options.residual_tol) * std::max(1.0, rhs_norm);
        if (projected_rhs_norm > zero_tolerance)
        {
            throw std::runtime_error(
                "Sternheimer delta fixed subspace fills the grid but leaves a nonzero Q-space rhs.");
        }
        result.solver.converged = true;
        result.solver.absolute_residual = projected_rhs_norm;
        result.solver.relative_residual = projected_rhs_norm / std::max(1.0, rhs_norm);
    }
    else
    {
        result.solver
            = SternheimerRPA::solve_gmres(problem, projected_rhs, result.response.out_wavefunction, options);
    }
    fixed_projector.project(result.response.out_wavefunction);

    const SternheimerDeltaCoefficientComponents components
        = compute_delta_coefficient_components(virtual_states,
                                               perturbation_matrix_elements,
                                               result.response.out_wavefunction,
                                               reference_eigenvalue,
                                               omega,
                                               volume_element);
    result.response.sos_coefficients = components.sos;
    result.response.pulay_coefficients = components.pulay;
    result.response.coefficients = components.total;
    assemble_delta_wavefunction_components(virtual_states, result.response);
    result.response.out_norm = sternheimer_fd_grid_norm(result.response.out_wavefunction, volume_element);

    Vector q_residual;
    hamiltonian.apply(result.response.out_wavefunction, q_residual);
    const Complex shift(-reference_eigenvalue, omega);
#pragma omp parallel for schedule(static)
    for (std::size_t ir = 0; ir != q_residual.size(); ++ir)
    {
        q_residual[ir] += shift * result.response.out_wavefunction[ir];
    }
    fixed_projector.project(q_residual);
    for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
    {
        axpy(result.response.coefficients[ia], virtual_states[ia].residual, q_residual);
    }

    Vector q_rhs = rhs;
    fixed_projector.project(q_rhs);
#pragma omp parallel for schedule(static)
    for (std::size_t ir = 0; ir != q_residual.size(); ++ir)
    {
        q_residual[ir] -= q_rhs[ir];
    }
    const double q_residual_norm = sternheimer_fd_grid_norm(q_residual, volume_element);
    double residual_norm_squared = q_residual_norm * q_residual_norm;
    for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
    {
        const Complex delta_block(virtual_states[ia].eigenvalue - reference_eigenvalue, omega);
        const Complex eta_residual = delta_block * result.response.coefficients[ia]
                                     + dot(virtual_states[ia].residual, result.response.out_wavefunction)
                                     + perturbation_matrix_elements[ia];
        residual_norm_squared += std::norm(eta_residual);
    }
    result.residual_norm = std::sqrt(residual_norm_squared);
    result.response.reconstruction_error = result.residual_norm;
    return result;
}

std::vector<SternheimerDeltaLinearResponse> solve_delta_sternheimer_linear_response_batch(
    const SternheimerFDHamiltonian& hamiltonian,
    const SternheimerDeltaFixedSubspace& fixed_subspace,
    const double reference_eigenvalue,
    const SternheimerFDHamiltonian::Matrix& rhs,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const std::vector<std::vector<SternheimerFDHamiltonian::Complex>>& perturbation_matrix_elements,
    const double omega,
    const double volume_element,
    const SternheimerRPA::SolverOptions& options)
{
    if (rhs.empty())
    {
        return {};
    }
    for (const Vector& function: fixed_subspace.functions)
    {
        check_vector_size(function, static_cast<std::size_t>(hamiltonian.grid().size()),
                          "Sternheimer delta batch fixed-subspace state");
    }
    const SternheimerSubspaceProjector projector(fixed_subspace.functions, volume_element, true);
    return solve_delta_sternheimer_linear_response_batch(hamiltonian, projector, reference_eigenvalue,
        rhs, virtual_states, perturbation_matrix_elements, omega, volume_element, options);
}

std::vector<SternheimerDeltaLinearResponse> solve_delta_sternheimer_linear_response_batch(
    const SternheimerFDHamiltonian& hamiltonian,
    const SternheimerSubspaceProjector& fixed_projector,
    const double reference_eigenvalue,
    const SternheimerFDHamiltonian::Matrix& rhs,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const std::vector<std::vector<SternheimerFDHamiltonian::Complex>>& perturbation_matrix_elements,
    const double omega,
    const double volume_element,
    const SternheimerRPA::SolverOptions& options)
{
    using Matrix = SternheimerFDHamiltonian::Matrix;
    if (rhs.empty())
    {
        return {};
    }
    if (rhs.size() != perturbation_matrix_elements.size())
    {
        throw std::invalid_argument(
            "Sternheimer delta batch response requires one perturbation-element vector per rhs.");
    }
    const int grid_size = hamiltonian.grid().size();
    for (std::size_t column = 0; column != rhs.size(); ++column)
    {
        check_vector_size(rhs[column], static_cast<std::size_t>(grid_size), "Sternheimer delta batch rhs");
        if (virtual_states.size() != perturbation_matrix_elements[column].size())
        {
            throw std::invalid_argument(
                "Sternheimer delta batch response requires one perturbation matrix element per virtual state.");
        }
    }
    if (volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer delta batch response requires a positive grid volume element.");
    }
    validate_virtual_states(virtual_states, static_cast<std::size_t>(grid_size));

    const auto dot = [volume_element](const Vector& lhs, const Vector& rhs_vec) {
        return sternheimer_fd_grid_dot(lhs, rhs_vec, volume_element);
    };

    std::vector<Complex> denominators(virtual_states.size(), Complex(0.0, 0.0));
    for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
    {
        denominators[ia] = Complex(reference_eigenvalue - virtual_states[ia].eigenvalue, -omega);
        if (std::abs(denominators[ia]) < 1.0e-30)
        {
            throw std::runtime_error(
                "Sternheimer delta batch response found a singular virtual-state denominator.");
        }
    }

    Matrix projected_rhs = rhs;
    fixed_projector.project_batch(projected_rhs);
    for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
    {
        for (std::size_t column = 0; column != rhs.size(); ++column)
        {
            axpy(-perturbation_matrix_elements[column][ia] / denominators[ia],
                 virtual_states[ia].residual,
                 projected_rhs[column]);
        }
    }

    SternheimerRPA::BatchLinearProblem problem;
    problem.dot = dot;
    problem.apply = [&hamiltonian,
                     &fixed_projector,
                     &virtual_states,
                     &denominators,
                     reference_eigenvalue,
                     omega,
                     dot](const Matrix& input, Matrix& output) {
        Matrix q_input = input;
        fixed_projector.project_batch(q_input);
        hamiltonian.apply_batch(q_input, output);
        const Complex shift(-reference_eigenvalue, omega);
        for (std::size_t column = 0; column != output.size(); ++column)
        {
#pragma omp parallel for schedule(static)
            for (std::size_t ir = 0; ir != output[column].size(); ++ir)
            {
                output[column][ir] += shift * q_input[column][ir];
            }
        }
        fixed_projector.project_batch(output);
        for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
        {
            for (std::size_t column = 0; column != output.size(); ++column)
            {
                const Complex coupling = dot(virtual_states[ia].residual, q_input[column]) / denominators[ia];
                axpy(coupling, virtual_states[ia].residual, output[column]);
            }
        }
    };
    std::shared_ptr<SternheimerFDSpectralPreconditioner> spectral_preconditioner;
    if (options.use_fd_spectral_preconditioner)
    {
        spectral_preconditioner = get_thread_local_sternheimer_fd_spectral_preconditioner(
            hamiltonian,
            reference_eigenvalue,
            omega,
            options.fd_spectral_preconditioner_regularization);
        problem.precondition = [spectral_preconditioner, &fixed_projector](const Matrix& input, Matrix& output) {
            spectral_preconditioner->apply_batch(input, output);
            fixed_projector.project_batch(output);
        };
    }

    std::vector<SternheimerDeltaLinearResponse> results(rhs.size());
    Matrix out_wavefunctions(rhs.size(), Vector(static_cast<std::size_t>(grid_size), Complex(0.0, 0.0)));
    if (fixed_projector.basis_size() >= static_cast<std::size_t>(grid_size))
    {
        for (std::size_t column = 0; column != rhs.size(); ++column)
        {
            const double projected_rhs_norm = sternheimer_fd_grid_norm(projected_rhs[column], volume_element);
            const double rhs_norm = sternheimer_fd_grid_norm(rhs[column], volume_element);
            const double zero_tolerance
                = std::max(options.breakdown_tol, options.residual_tol) * std::max(1.0, rhs_norm);
            if (projected_rhs_norm > zero_tolerance)
            {
                throw std::runtime_error(
                    "Sternheimer delta batch fixed subspace fills the grid but leaves a nonzero Q-space rhs.");
            }
            results[column].solver.converged = true;
            results[column].solver.absolute_residual = projected_rhs_norm;
            results[column].solver.relative_residual = projected_rhs_norm / std::max(1.0, rhs_norm);
        }
    }
    else
    {
        const std::vector<SternheimerRPA::SolverResult> solver_results
            = SternheimerRPA::solve_gmres_batch(problem, projected_rhs, out_wavefunctions, options);
        for (std::size_t column = 0; column != rhs.size(); ++column)
        {
            results[column].solver = solver_results[column];
        }
    }
    fixed_projector.project_batch(out_wavefunctions);
    for (std::size_t column = 0; column != rhs.size(); ++column)
    {
        results[column].response.out_wavefunction = out_wavefunctions[column];
        const SternheimerDeltaCoefficientComponents components
            = compute_delta_coefficient_components(virtual_states,
                                                   perturbation_matrix_elements[column],
                                                   results[column].response.out_wavefunction,
                                                   reference_eigenvalue,
                                                   omega,
                                                   volume_element);
        results[column].response.sos_coefficients = components.sos;
        results[column].response.pulay_coefficients = components.pulay;
        results[column].response.coefficients = components.total;
        assemble_delta_wavefunction_components(virtual_states, results[column].response);
        results[column].response.out_norm
            = sternheimer_fd_grid_norm(results[column].response.out_wavefunction, volume_element);
    }

    Matrix q_residual;
    hamiltonian.apply_batch(out_wavefunctions, q_residual);
    const Complex shift(-reference_eigenvalue, omega);
    for (std::size_t column = 0; column != q_residual.size(); ++column)
    {
#pragma omp parallel for schedule(static)
        for (std::size_t ir = 0; ir != q_residual[column].size(); ++ir)
        {
            q_residual[column][ir] += shift * out_wavefunctions[column][ir];
        }
    }
    fixed_projector.project_batch(q_residual);
    for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
    {
        for (std::size_t column = 0; column != q_residual.size(); ++column)
        {
            axpy(results[column].response.coefficients[ia], virtual_states[ia].residual, q_residual[column]);
        }
    }

    Matrix q_rhs = rhs;
    fixed_projector.project_batch(q_rhs);
    for (std::size_t column = 0; column != rhs.size(); ++column)
    {
#pragma omp parallel for schedule(static)
        for (std::size_t ir = 0; ir != q_residual[column].size(); ++ir)
        {
            q_residual[column][ir] -= q_rhs[column][ir];
        }
        const double q_residual_norm = sternheimer_fd_grid_norm(q_residual[column], volume_element);
        double residual_norm_squared = q_residual_norm * q_residual_norm;
        for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
        {
            const Complex delta_block(virtual_states[ia].eigenvalue - reference_eigenvalue, omega);
            const Complex eta_residual = delta_block * results[column].response.coefficients[ia]
                                         + dot(virtual_states[ia].residual, results[column].response.out_wavefunction)
                                         + perturbation_matrix_elements[column][ia];
            residual_norm_squared += std::norm(eta_residual);
        }
        results[column].residual_norm = std::sqrt(residual_norm_squared);
        results[column].response.reconstruction_error = results[column].residual_norm;
    }
    return results;
}

SternheimerDeltaLinearResponse solve_delta_sternheimer_linear_response(
    const SternheimerFDHamiltonian& hamiltonian,
    const std::vector<SternheimerFDHamiltonian::Vector>& occupied_wavefunctions,
    const double reference_eigenvalue,
    const SternheimerFDHamiltonian::Vector& rhs,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const std::vector<SternheimerFDHamiltonian::Complex>& perturbation_matrix_elements,
    const double omega,
    const double volume_element,
    const SternheimerRPA::SolverOptions& options)
{
    const SternheimerDeltaFixedSubspace fixed_subspace
        = build_delta_sternheimer_fixed_subspace(occupied_wavefunctions, virtual_states);
    return solve_delta_sternheimer_linear_response(hamiltonian,
                                                    fixed_subspace,
                                                    reference_eigenvalue,
                                                    rhs,
                                                    virtual_states,
                                                    perturbation_matrix_elements,
                                                    omega,
                                                    volume_element,
                                                    options);
}

SternheimerFDHamiltonian::Complex accumulate_delta_sternheimer_response(
    const SternheimerFDHamiltonian::Vector& probe_wavefunction,
    const SternheimerDeltaPostprocessResult& response,
    const std::vector<SternheimerDeltaVirtualState>& virtual_states,
    const double volume_element)
{
    if (volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer delta response requires a positive grid volume element.");
    }
    if (response.coefficients.size() != virtual_states.size())
    {
        throw std::invalid_argument("Sternheimer delta response coefficient count does not match virtual states.");
    }
    check_vector_size(response.out_wavefunction, probe_wavefunction.size(), "Sternheimer delta response out state");

    Complex value = sternheimer_fd_grid_dot(probe_wavefunction, response.out_wavefunction, volume_element);
    for (std::size_t ia = 0; ia != virtual_states.size(); ++ia)
    {
        check_vector_size(virtual_states[ia].orbital, probe_wavefunction.size(), "Sternheimer delta response virtual");
        value += response.coefficients[ia]
            * sternheimer_fd_grid_dot(probe_wavefunction, virtual_states[ia].orbital, volume_element);
    }
    return value;
}

} // namespace ModuleRI
