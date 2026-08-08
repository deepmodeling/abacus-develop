#include "source_lcao/module_ri/sternheimer_grid_diagnostics.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace ModuleRI
{
namespace
{

using Complex = std::complex<double>;

std::size_t checked_product(const int first, const int second, const char* label)
{
    if (first < 0 || second < 0)
    {
        throw std::invalid_argument(std::string(label) + " dimensions must be non-negative.");
    }
    const std::size_t lhs = static_cast<std::size_t>(first);
    const std::size_t rhs = static_cast<std::size_t>(second);
    if (rhs != 0 && lhs > std::numeric_limits<std::size_t>::max() / rhs)
    {
        throw std::overflow_error(std::string(label) + " dimensions overflow size_t.");
    }
    return lhs * rhs;
}

std::size_t checked_product(const int first, const int second, const int third, const char* label)
{
    const std::size_t first_second = checked_product(first, second, label);
    if (third < 0)
    {
        throw std::invalid_argument(std::string(label) + " dimensions must be non-negative.");
    }
    const std::size_t rhs = static_cast<std::size_t>(third);
    if (rhs != 0 && first_second > std::numeric_limits<std::size_t>::max() / rhs)
    {
        throw std::overflow_error(std::string(label) + " dimensions overflow size_t.");
    }
    return first_second * rhs;
}

void validate_metadata(const SternheimerGridDiagnosticMetadata& metadata)
{
    if (metadata.nx <= 0 || metadata.ny <= 0 || metadata.nz <= 0)
    {
        throw std::invalid_argument("Sternheimer grid diagnostic dimensions must be positive.");
    }
    if (metadata.spin < 0 || metadata.occupied < 0 || metadata.virtuals < 0 || metadata.auxiliaries < 0)
    {
        throw std::invalid_argument("Sternheimer grid diagnostic state counts must be non-negative.");
    }
    if (metadata.volume_element <= 0.0)
    {
        throw std::invalid_argument("Sternheimer grid diagnostic volume element must be positive.");
    }
}

void validate_operator_sizes(const SternheimerDeltaGridMatrices& matrices, const std::size_t expected_size)
{
    if (matrices.overlap.size() != expected_size || matrices.kinetic.size() != expected_size
        || matrices.local_potential.size() != expected_size || matrices.nonlocal.size() != expected_size
        || matrices.hamiltonian.size() != expected_size)
    {
        throw std::invalid_argument("Sternheimer grid diagnostic operator matrix dimensions do not match.");
    }
}

double relative_norm(const double difference_squared, const double reference_squared)
{
    const double difference = std::sqrt(difference_squared);
    const double reference = std::sqrt(reference_squared);
    return reference > 0.0 ? difference / reference : difference;
}

void write_metadata(std::ofstream& output, const SternheimerGridDiagnosticMetadata& metadata, const char* kind)
{
    output << "ABACUS_STERNHEIMER_GRID_DIAGNOSTICS 1\n";
    output << "kind " << kind << '\n';
    output << "grid " << metadata.nx << ' ' << metadata.ny << ' ' << metadata.nz << '\n';
    output << "spin " << metadata.spin << '\n';
    output << "occupied " << metadata.occupied << '\n';
    output << "virtuals " << metadata.virtuals << '\n';
    output << "auxiliaries " << metadata.auxiliaries << '\n';
    output << "volume_element " << std::setprecision(17) << std::scientific << metadata.volume_element << '\n';
}

void write_column_major_matrix(std::ofstream& output,
                               const char* label,
                               const std::vector<Complex>& matrix,
                               const int dimension)
{
    output << "matrix " << label << '\n';
    for (int row = 0; row != dimension; ++row)
    {
        for (int column = 0; column != dimension; ++column)
        {
            const Complex value = matrix[static_cast<std::size_t>(row)
                                         + static_cast<std::size_t>(dimension) * static_cast<std::size_t>(column)];
            output << row << ' ' << column << ' ' << value.real() << ' ' << value.imag() << '\n';
        }
    }
}

} // namespace

SternheimerPerturbationTensor::SternheimerPerturbationTensor(const int occupied_count,
                                                             const int virtual_count,
                                                             const int auxiliary_count)
    : occupied(occupied_count), virtuals(virtual_count), auxiliaries(auxiliary_count),
      values(checked_product(occupied_count, virtual_count, auxiliary_count, "Sternheimer perturbation tensor"),
             Complex(0.0, 0.0))
{
}

SternheimerPerturbationTensor::Complex& SternheimerPerturbationTensor::at(const int occupied_index,
                                                                          const int virtual_index,
                                                                          const int auxiliary_index)
{
    return const_cast<Complex&>(
        static_cast<const SternheimerPerturbationTensor&>(*this).at(occupied_index, virtual_index, auxiliary_index));
}

const SternheimerPerturbationTensor::Complex& SternheimerPerturbationTensor::at(const int occupied_index,
                                                                                const int virtual_index,
                                                                                const int auxiliary_index) const
{
    if (occupied_index < 0 || occupied_index >= occupied || virtual_index < 0 || virtual_index >= virtuals
        || auxiliary_index < 0 || auxiliary_index >= auxiliaries)
    {
        throw std::out_of_range("Sternheimer perturbation tensor index is out of range.");
    }
    const std::size_t index = (static_cast<std::size_t>(occupied_index) * static_cast<std::size_t>(virtuals)
                               + static_cast<std::size_t>(virtual_index))
                                  * static_cast<std::size_t>(auxiliaries)
                              + static_cast<std::size_t>(auxiliary_index);
    return values[index];
}

double relative_component_reconstruction_error(const std::vector<Complex>& total,
                                               const std::vector<Complex>& sos,
                                               const std::vector<Complex>& pulay,
                                               const std::vector<Complex>& qspace)
{
    if (total.size() != sos.size() || total.size() != pulay.size() || total.size() != qspace.size())
    {
        throw std::invalid_argument("Sternheimer response component dimensions do not match.");
    }
    double difference_squared = 0.0;
    double reference_squared = 0.0;
    for (std::size_t index = 0; index != total.size(); ++index)
    {
        difference_squared += std::norm(total[index] - sos[index] - pulay[index] - qspace[index]);
        reference_squared += std::norm(total[index]);
    }
    return relative_norm(difference_squared, reference_squared);
}

double relative_operator_reconstruction_error(const SternheimerDeltaGridMatrices& matrices)
{
    const std::size_t size = matrices.hamiltonian.size();
    validate_operator_sizes(matrices, size);
    double difference_squared = 0.0;
    double reference_squared = 0.0;
    for (std::size_t index = 0; index != size; ++index)
    {
        difference_squared += std::norm(matrices.hamiltonian[index] - matrices.kinetic[index]
                                        - matrices.local_potential[index] - matrices.nonlocal[index]);
        reference_squared += std::norm(matrices.hamiltonian[index]);
    }
    return relative_norm(difference_squared, reference_squared);
}

void write_delta_grid_matrices(const std::string& path,
                               const SternheimerGridDiagnosticMetadata& metadata,
                               const SternheimerDeltaGridMatrices& matrices,
                               const std::vector<Complex>& occupied_virtual_overlap,
                               const double reconstruction_tolerance)
{
    validate_metadata(metadata);
    if (reconstruction_tolerance < 0.0)
    {
        throw std::invalid_argument("Sternheimer operator reconstruction tolerance must be non-negative.");
    }
    const std::size_t matrix_size
        = checked_product(metadata.virtuals, metadata.virtuals, "Sternheimer grid diagnostic matrix");
    validate_operator_sizes(matrices, matrix_size);
    if (occupied_virtual_overlap.size()
        != checked_product(metadata.occupied, metadata.virtuals, "Sternheimer occupied-virtual overlap"))
    {
        throw std::invalid_argument("Sternheimer occupied-virtual overlap dimensions do not match.");
    }
    if (relative_operator_reconstruction_error(matrices) > reconstruction_tolerance)
    {
        throw std::runtime_error("Sternheimer grid operator components do not reconstruct the Hamiltonian.");
    }

    std::ofstream output(path.c_str(), std::ios::out | std::ios::trunc);
    if (!output)
    {
        throw std::runtime_error("Failed to open Sternheimer grid diagnostic output: " + path);
    }
    write_metadata(output, metadata, "delta_grid_matrices");
    output << "storage column_major\n";
    write_column_major_matrix(output, "overlap", matrices.overlap, metadata.virtuals);
    write_column_major_matrix(output, "kinetic", matrices.kinetic, metadata.virtuals);
    write_column_major_matrix(output, "local_potential", matrices.local_potential, metadata.virtuals);
    write_column_major_matrix(output, "nonlocal", matrices.nonlocal, metadata.virtuals);
    write_column_major_matrix(output, "hamiltonian", matrices.hamiltonian, metadata.virtuals);
    output << "matrix occupied_virtual_overlap\n";
    for (int occupied = 0; occupied != metadata.occupied; ++occupied)
    {
        for (int virtual_index = 0; virtual_index != metadata.virtuals; ++virtual_index)
        {
            const Complex value
                = occupied_virtual_overlap[static_cast<std::size_t>(occupied) * metadata.virtuals + virtual_index];
            output << occupied << ' ' << virtual_index << ' ' << value.real() << ' ' << value.imag() << '\n';
        }
    }
    if (!output)
    {
        throw std::runtime_error("Failed to write Sternheimer grid diagnostic output: " + path);
    }
}

void write_delta_perturbation_tensor(const std::string& path,
                                     const SternheimerGridDiagnosticMetadata& metadata,
                                     const SternheimerPerturbationTensor& tensor)
{
    validate_metadata(metadata);
    if (tensor.occupied != metadata.occupied || tensor.virtuals != metadata.virtuals
        || tensor.auxiliaries != metadata.auxiliaries
        || tensor.values.size()
               != checked_product(tensor.occupied,
                                  tensor.virtuals,
                                  tensor.auxiliaries,
                                  "Sternheimer perturbation tensor"))
    {
        throw std::invalid_argument("Sternheimer perturbation tensor dimensions do not match metadata.");
    }

    std::ofstream output(path.c_str(), std::ios::out | std::ios::trunc);
    if (!output)
    {
        throw std::runtime_error("Failed to open Sternheimer perturbation diagnostic output: " + path);
    }
    write_metadata(output, metadata, "delta_perturbation_tensor");
    output << "storage row_major\n";
    output << "tensor perturbation occupied virtual auxiliary\n";
    for (int occupied = 0; occupied != tensor.occupied; ++occupied)
    {
        for (int virtual_index = 0; virtual_index != tensor.virtuals; ++virtual_index)
        {
            for (int auxiliary = 0; auxiliary != tensor.auxiliaries; ++auxiliary)
            {
                const Complex value = tensor.at(occupied, virtual_index, auxiliary);
                output << occupied << ' ' << virtual_index << ' ' << auxiliary << ' ' << value.real() << ' '
                       << value.imag() << '\n';
            }
        }
    }
    if (!output)
    {
        throw std::runtime_error("Failed to write Sternheimer perturbation diagnostic output: " + path);
    }
}

} // namespace ModuleRI
