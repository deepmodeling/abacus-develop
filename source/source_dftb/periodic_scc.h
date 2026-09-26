#ifndef SOURCE_DFTB_PERIODIC_SCC_H
#define SOURCE_DFTB_PERIODIC_SCC_H

#include "source_dftb/eigensolver.h"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace ModuleDFTB
{

struct DftbWeightedKPoint
{
    std::array<double, 3> fractional{};
    double weight = 0.0;
};

struct DftbBandKPoint
{
    std::array<double, 3> fractional{};
    std::string label;
};

struct DftbBandPoint
{
    std::array<double, 3> fractional{};
    std::string label;
    double distance_inverse_bohr = 0.0;
    std::vector<double> eigenvalues_hartree;
};

struct DftbPeriodicInput
{
    std::vector<DftbSpAtom> atoms;
    std::vector<DftbPairParameters> pair_parameters;
    // Lattice vectors are Cartesian Bohr vectors, one vector per column.
    std::array<std::array<double, 3>, 3> lattice_bohr{};
    std::vector<DftbWeightedKPoint> kpoints;
    // Optional non-self-consistent band path evaluated at the converged SCC potential.
    std::vector<DftbBandKPoint> band_kpoints;
    // dU/d(delta-q) for each species; delta-q is electron population minus
    // neutral reference population, matching DFTB+ internal convention.
    std::vector<double> hubbard_derivative;
    double total_electrons = 0.0;
    double thermal_energy_hartree = 0.0;
    double scc_tolerance = 1.0e-6;
    int maximum_scc_iterations = 200;
    double mixing_parameter = 0.2;
    bool third_order = true;
};

struct DftbPeriodicResult
{
    std::vector<double> electron_excess_charges;
    std::vector<double> scc_residual_history;
    std::vector<DftbBandPoint> band_structure;
    double fermi_energy_hartree = 0.0;
    double band_energy_hartree = 0.0;
    double band_free_energy_hartree = 0.0;
    double h0_energy_hartree = 0.0;
    double scc_energy_hartree = 0.0;
    double third_order_energy_hartree = 0.0;
    double repulsive_energy_hartree = 0.0;
    double total_free_energy_hartree = 0.0;
    double maximum_charge_residual = 0.0;
    int scc_iterations = 0;
    bool converged = false;
    std::vector<double> atomic_electron_populations;
};

/**
 * Native periodic DFTB2/DFTB3 SCC engine for s+p valence bases.
 *
 * Implements Bloch H/S assembly, 3D-periodic Ewald electrostatics with the
 * SKF short-range gamma correction, Mulliken populations, atomic SCC shifts,
 * full third-order atomic-charge terms, finite-temperature filling, and SKF
 * repulsive energy. Forces, stress, spin, and shell-resolved Hubbard data are
 * not yet implemented.
 */
DftbPeriodicResult solve_periodic_dftb(const DftbPeriodicInput& input);

} // namespace ModuleDFTB

#endif // SOURCE_DFTB_PERIODIC_SCC_H
