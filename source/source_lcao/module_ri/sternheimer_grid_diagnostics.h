#ifndef STERNHEIMER_GRID_DIAGNOSTICS_H
#define STERNHEIMER_GRID_DIAGNOSTICS_H

#include "source_lcao/module_ri/sternheimer_delta.h"

#include <complex>
#include <string>
#include <vector>

namespace ModuleRI
{

struct SternheimerGridDiagnosticMetadata
{
    int nx = 0;
    int ny = 0;
    int nz = 0;
    int spin = 0;
    int occupied = 0;
    int virtuals = 0;
    int auxiliaries = 0;
    double volume_element = 0.0;
};

struct SternheimerPerturbationTensor
{
    using Complex = std::complex<double>;

    SternheimerPerturbationTensor() = default;
    SternheimerPerturbationTensor(int occupied_count, int virtual_count, int auxiliary_count);

    Complex& at(int occupied_index, int virtual_index, int auxiliary_index);
    const Complex& at(int occupied_index, int virtual_index, int auxiliary_index) const;

    int occupied = 0;
    int virtuals = 0;
    int auxiliaries = 0;
    std::vector<Complex> values;
};

double relative_component_reconstruction_error(const std::vector<std::complex<double>>& total,
                                               const std::vector<std::complex<double>>& sos,
                                               const std::vector<std::complex<double>>& pulay,
                                               const std::vector<std::complex<double>>& qspace);

double relative_operator_reconstruction_error(const SternheimerDeltaGridMatrices& matrices);

void write_delta_grid_matrices(const std::string& path,
                               const SternheimerGridDiagnosticMetadata& metadata,
                               const SternheimerDeltaGridMatrices& matrices,
                               const std::vector<std::complex<double>>& occupied_virtual_overlap,
                               double reconstruction_tolerance);

void write_delta_perturbation_tensor(const std::string& path,
                                     const SternheimerGridDiagnosticMetadata& metadata,
                                     const SternheimerPerturbationTensor& tensor);

} // namespace ModuleRI

#endif
