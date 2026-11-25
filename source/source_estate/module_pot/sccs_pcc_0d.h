#ifndef SCCS_PCC_0D_ADAPTER_H
#define SCCS_PCC_0D_ADAPTER_H

#include <memory>
#include <vector>

class UnitCell;
namespace ModuleBase
{
template <typename T> class Vector3;
}
namespace ModulePW
{
class PW_Basis;
}
namespace ModuleSccs
{
class CoulombOperator;
}
namespace unitcell
{
struct AtomData;
struct OrthogonalCell;
}
namespace elecstate
{
struct Pcc0dParameters;

// Cubic PCC cell with the mass-weighted ionic origin, as PotPcc; a cell that
// is not an equal-edge cube stops the run.
void pcc_0d_geometry(const UnitCell& cell,
                     const std::vector<unitcell::AtomData>& atoms,
                     unitcell::OrthogonalCell& geometry,
                     Pcc0dParameters& parameters);

// Minimum-image displacements of Cartesian positions from the origin, in place.
void pcc_0d_relative(const unitcell::OrthogonalCell& geometry,
                     std::vector<ModuleBase::Vector3<double>>& positions);

// Assemble an SCCS grid-charge operator using the existing cell geometry and
// mass-weighted ionic origin, on the Cartesian grid coordinates grid_positions
// (e.g. from SccsFixedSources). An unsupported cell stops the run.
void pcc_0d_operator(const UnitCell& cell,
                     const ModulePW::PW_Basis& basis,
                     const std::vector<unitcell::AtomData>& atoms,
                     const std::vector<ModuleBase::Vector3<double>>& grid_positions,
                     std::unique_ptr<ModuleSccs::CoulombOperator>& coulomb);
} // namespace elecstate

#endif
