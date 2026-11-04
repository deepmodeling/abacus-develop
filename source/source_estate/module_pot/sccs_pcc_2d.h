#ifndef SCCS_PCC_2D_ADAPTER_H
#define SCCS_PCC_2D_ADAPTER_H

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
struct SlabCell;
}
namespace elecstate
{
struct Pcc2dParameters;

// Slab PCC cell with the mass-weighted ionic origin along the normal, with
// the PotPcc alignment tolerance; a tilted open axis stops the run.
void pcc_2d_geometry(const UnitCell& cell,
                     const std::vector<unitcell::AtomData>& atoms,
                     int open_axis,
                     unitcell::SlabCell& geometry,
                     Pcc2dParameters& parameters);

// Replace Cartesian positions, in place, by the normal times their
// minimum-image coordinate from the slab origin.
void pcc_2d_project(const unitcell::SlabCell& geometry,
                    std::vector<ModuleBase::Vector3<double>>& positions);

// Assemble an SCCS grid-charge operator using the existing cell geometry and
// mass-weighted ionic origin, on the Cartesian grid coordinates grid_positions
// (e.g. from SccsFixedSources). An unsupported cell stops the run.
void pcc_2d_operator(const UnitCell& cell,
                     const ModulePW::PW_Basis& basis,
                     const std::vector<unitcell::AtomData>& atoms,
                     int open_axis,
                     const std::vector<ModuleBase::Vector3<double>>& grid_positions,
                     std::unique_ptr<ModuleSccs::CoulombOperator>& coulomb);
} // namespace elecstate

#endif
