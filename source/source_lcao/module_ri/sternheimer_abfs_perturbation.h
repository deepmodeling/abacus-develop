#ifndef STERNHEIMER_ABFS_PERTURBATION_H
#define STERNHEIMER_ABFS_PERTURBATION_H

#include "source_base/vector3.h"
#include "source_basis/module_ao/ORB_atomic_lm.h"
#include "source_lcao/module_ri/sternheimer_fd_hamiltonian.h"

#include <complex>
#include <limits>
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
    bool continue_coulomb_tail = false;
    double coulomb_tail_coefficient = std::numeric_limits<double>::quiet_NaN();
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

using SternheimerOrbitalSet = std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>;

struct SternheimerCoulombProjectionDiagnostic
{
    double relative_error = 0.0;
};

std::vector<std::vector<SternheimerRadialPerturbation>> make_sternheimer_radial_perturbations_from_orbitals(
    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orbitals);

// Set the exterior r^(-l-1) coefficient from the finalized compact auxiliary charge.
void set_sternheimer_coulomb_tail_from_charge(SternheimerRadialPerturbation& potential,
                                              const Numerical_Orbital_Lm& charge);

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

// max_channels limits the total returned count, not the count per atom.
// start_channel skips global descriptors without sampling their fine-grid values.
std::vector<SternheimerABFBlochGridChannel> sample_sternheimer_abf_bloch_grid_channels(
    const std::vector<std::vector<SternheimerRadialPerturbation>>& radials_by_type,
    const std::vector<int>& atom_types,
    const std::vector<ModuleBase::Vector3<double>>& atom_positions,
    const SternheimerFDHamiltonian::Grid& grid,
    const SternheimerReducedKPoint& qpoint,
    int max_channels = -1,
    int start_channel = 0);

// Solve the periodic Poisson equation for Bloch auxiliary densities.
// Input channel values are densities; output channel values are Hartree potentials in Ha.
// gamma_inverse_k2 replaces 1/|G+q|^2 only for the Gamma zero mode and must be
// zero for non-Gamma q points.
std::vector<SternheimerABFBlochGridChannel> solve_sternheimer_abf_periodic_full_coulomb(
    const std::vector<SternheimerABFBlochGridChannel>& density_channels,
    const SternheimerFDHamiltonian::Grid& grid,
    const SternheimerReducedKPoint& qpoint,
    double gamma_inverse_k2);

// Replace each input density by its periodic Hartree potential. This is the
// production interface for large auxiliary spaces because it does not retain
// a second full-grid channel set.
void solve_sternheimer_abf_periodic_full_coulomb_in_place(std::vector<SternheimerABFBlochGridChannel>& density_channels,
                                                          const SternheimerFDHamiltonian::Grid& grid,
                                                          const SternheimerReducedKPoint& qpoint,
                                                          double gamma_inverse_k2);

std::vector<std::complex<double>> sternheimer_grid_projected_matrix(
    const std::vector<SternheimerABFBlochGridChannel>& densities,
    const std::vector<SternheimerABFBlochGridChannel>& potentials,
    double volume_element);

SternheimerCoulombProjectionDiagnostic compare_sternheimer_periodic_coulomb_projection(
    const std::vector<SternheimerABFBlochGridChannel>& densities,
    const std::vector<SternheimerABFBlochGridChannel>& potentials,
    const std::vector<std::complex<double>>& target_coulomb,
    double volume_element);

} // namespace ModuleRI

#endif
