#include "sccs_pcc_2d.h"

#include "source_base/tool_quit.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/cell_geometry.h"
#include "source_cell/cell_tools.h"
#include "source_cell/unitcell.h"
#include "source_hamilt/module_sccs/sccs_coul_2d.h"

namespace elecstate
{
void pcc_2d_geometry(const UnitCell& cell,
                     const std::vector<unitcell::AtomData>& atoms,
                     int open_axis,
                     unitcell::SlabCell& geometry,
                     Pcc2dParameters& parameters)
{
    // Same alignment tolerance as PotPcc, loose enough for lattice vectors
    // typed with rounded digits.
    const double alignment_tolerance = 1e-6;
    const bool perpendicular = unitcell::make_slab_cell(cell.latvec, cell.lat0, open_axis, alignment_tolerance,
                                                        geometry);
    if (!perpendicular)
    {
        ModuleBase::WARNING_QUIT("pcc_2d_geometry",
                                 "PCC 2D requires the open lattice vector to be perpendicular to the periodic plane");
    }
    parameters = make_pcc_2d_parameters(geometry);
    std::vector<ModuleBase::Vector3<double>> positions;
    std::vector<double> masses;
    for (const unitcell::AtomData& atom : atoms)
    {
        positions.push_back(atom.position);
        masses.push_back(atom.mass);
    }
    geometry.origin = unitcell::weighted_center(positions, masses, geometry);
}

void pcc_2d_project(const unitcell::SlabCell& geometry,
                    std::vector<ModuleBase::Vector3<double>>& positions)
{
    for (ModuleBase::Vector3<double>& position : positions)
    {
        const double coordinate = unitcell::relative_coordinate(position, geometry);
        position = geometry.normal * coordinate;
    }
}

void pcc_2d_operator(const UnitCell& cell,
                     const ModulePW::PW_Basis& basis,
                     const std::vector<unitcell::AtomData>& atoms,
                     int open_axis,
                     const std::vector<ModuleBase::Vector3<double>>& grid_positions,
                     std::unique_ptr<ModuleSccs::CoulombOperator>& coulomb)
{
    unitcell::SlabCell geometry;
    Pcc2dParameters parameters;
    pcc_2d_geometry(cell, atoms, open_axis, geometry, parameters);
    std::vector<ModuleBase::Vector3<double>> projected_positions = grid_positions;
    pcc_2d_project(geometry, projected_positions);
    coulomb.reset(new ModuleSccs::Pcc2dCoulombOperator(basis, cell.tpiba, projected_positions, geometry.normal,
                                                       parameters));
}
} // namespace elecstate
