#ifndef STERNHEIMER_ABFS_PERTURBATION_H
#define STERNHEIMER_ABFS_PERTURBATION_H

#include "source_base/vector3.h"
#include "source_basis/module_ao/ORB_atomic_lm.h"
#include "source_lcao/module_ri/sternheimer_fd_hamiltonian.h"

#include <complex>
#include <string>
#include <vector>

namespace ModuleRI
{

constexpr int sternheimer_abfs_transform_grid_chunk = 1024;

struct SternheimerRadialPerturbation
{
    int type_index = -1;
    int angular_momentum = 0;
    int radial_index = 0;
    std::string label;
    std::vector<double> radial_grid;
    std::vector<double> radial_values;
};

struct SternheimerABFGridChannel
{
    int channel_index = -1;
    int atom_index = -1;
    int atom_local_index = -1;
    int type_index = -1;
    int angular_momentum = 0;
    int radial_index = 0;
    int magnetic_index = 0;
    std::string label;
    std::vector<double> potential_r;
    double max_abs = 0.0;
};

struct SternheimerABFBlochGridChannel
{
    int channel_index = -1;
    int atom_index = -1;
    int atom_local_index = -1;
    int type_index = -1;
    int angular_momentum = 0;
    int radial_index = 0;
    int magnetic_index = 0;
    std::string label;
    std::vector<std::complex<double>> potential_r;
    double max_abs = 0.0;
};

struct SternheimerCoulombProjectionDiagnostic
{
    double relative_error = 0.0;
};

std::vector<std::vector<SternheimerRadialPerturbation>> make_sternheimer_radial_perturbations_from_orbitals(
    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orbitals);

std::vector<SternheimerABFGridChannel> describe_sternheimer_abf_grid_channels(
    const std::vector<std::vector<SternheimerRadialPerturbation>>& radials_by_type,
    const std::vector<int>& atom_types,
    const std::vector<ModuleBase::Vector3<double>>& atom_positions,
    int max_channels = -1);

std::vector<SternheimerABFGridChannel> sample_sternheimer_abf_grid_channels(
    const std::vector<std::vector<SternheimerRadialPerturbation>>& radials_by_type,
    const std::vector<int>& atom_types,
    const std::vector<ModuleBase::Vector3<double>>& atom_positions,
    const SternheimerFDHamiltonian::Grid& grid,
    int max_channels = -1);

std::vector<std::vector<double>> sample_sternheimer_abf_grid_channel_transform(
    const std::vector<std::vector<SternheimerRadialPerturbation>>& radials_by_type,
    const std::vector<int>& atom_types,
    const std::vector<ModuleBase::Vector3<double>>& atom_positions,
    const SternheimerFDHamiltonian::Grid& grid,
    std::vector<SternheimerABFGridChannel>& raw_channels,
    const std::vector<double>& raw_to_output,
    int output_channels);

} // namespace ModuleRI

#endif
