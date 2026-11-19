#ifndef RELAX_GEOMETRY_H
#define RELAX_GEOMETRY_H

#include <iosfwd>
#include <string>

class UnitCell;
namespace ModuleBase
{
class matrix;
class Matrix3;
} // namespace ModuleBase

namespace relax_geometry
{
// Apply the atomic move flags to forces in Ry/Bohr. The output and its
// largest absolute component (the return value) are both in eV/Angstrom.
double project_ionic_forces(
    const UnitCell &ucell, const ModuleBase::matrix &force, ModuleBase::matrix &projected_force);

// Apply the cell constraints to stress in Ry/Bohr^3. Fill the projected
// volume times stress in eV and return the largest absolute stress in kbar.
double project_cell_stress(
    const UnitCell &ucell,
    const ModuleBase::matrix &stress,
    const std::string &fixed_axes,
    ModuleBase::matrix &projected_virial);

// Add a displacement in latvec units, enforcing locked rows, volume and Bravais constraints.
void apply_lattice_step(
    UnitCell &ucell,
    ModuleBase::Matrix3 displacement,
    bool fixed_volume,
    bool fixed_ibrav);

// Convert a fraction of a Cartesian displacement (Angstrom) with the pre-step GT,
// apply move flags and symmetry, and update the atoms using the current lattice.
void apply_ionic_step(
    UnitCell &ucell,
    const ModuleBase::Matrix3 &old_gt,
    ModuleBase::matrix displacement,
    double fraction);

// Refresh the reciprocal cell, volume and the cell data used by the next SCF.
void refresh_cell(UnitCell &ucell, int nspin, std::ofstream &ofs_running);
} // namespace relax_geometry

#endif
