#include "sccs_pcc_0d.h"

#include "source_base/tool_quit.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/cell_geometry.h"
#include "source_cell/cell_tools.h"
#include "source_cell/unitcell.h"
#include "source_hamilt/module_sccs/sccs_coul_0d.h"

namespace elecstate
{
void pcc_0d_geometry(const UnitCell& cell,
                     const std::vector<unitcell::AtomData>& atoms,
                     unitcell::OrthogonalCell& geometry,
                     Pcc0dParameters& parameters)
{
    // Same orthogonality and equal-edge tolerance as PotPcc.
    const double cell_tolerance = 1e-10;
    const bool orthogonal = unitcell::make_orthogonal_cell(cell.latvec, cell.lat0, cell_tolerance, geometry);
    const bool cubic = orthogonal && make_pcc_0d_parameters(geometry, cell_tolerance, parameters);
    if (!cubic)
    {
        ModuleBase::WARNING_QUIT("pcc_0d_geometry", "PCC 0D requires an equal-edge cubic cell");
    }
    std::vector<ModuleBase::Vector3<double>> positions;
    std::vector<double> masses;
    for (const unitcell::AtomData& atom : atoms)
    {
        positions.push_back(atom.position);
        masses.push_back(atom.mass);
    }
    geometry.origin = unitcell::weighted_center(positions, masses, geometry);
}

void pcc_0d_relative(const unitcell::OrthogonalCell& geometry,
                     std::vector<ModuleBase::Vector3<double>>& positions)
{
    for (ModuleBase::Vector3<double>& position : positions)
    {
        position = unitcell::relative_position(position, geometry);
    }
}

void pcc_0d_operator(const UnitCell& cell,
                     const ModulePW::PW_Basis& basis,
                     const std::vector<unitcell::AtomData>& atoms,
                     const std::vector<ModuleBase::Vector3<double>>& grid_positions,
                     std::unique_ptr<ModuleSccs::CoulombOperator>& coulomb)
{
    unitcell::OrthogonalCell geometry;
    Pcc0dParameters parameters;
    pcc_0d_geometry(cell, atoms, geometry, parameters);
    std::vector<ModuleBase::Vector3<double>> relative_positions = grid_positions;
    pcc_0d_relative(geometry, relative_positions);
    coulomb.reset(new ModuleSccs::Pcc0dCoulombOperator(basis, cell.tpiba, relative_positions, parameters));
}
} // namespace elecstate
