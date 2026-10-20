#include "sccs_output.h"
#include "sccs_pcc_0d.h"
#include "sccs_pcc_2d.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/cell_geometry.h"
#include "source_cell/cell_tools.h"
#include "source_cell/unitcell.h"
#include "source_estate/pcc_0d.h"
#include "source_estate/pcc_2d.h"
#include "source_hamilt/module_sccs/sccs_diag.h"
#include "source_hamilt/module_sccs/sccs_parameters.h"
#include "source_hamilt/module_sccs/sccs_moments.h"
#include "source_hamilt/module_sccs/sccs_response.h"

namespace elecstate
{
void collect_output(const UnitCell& cell,
                     const ModulePW::PW_Basis& basis,
                     const std::vector<unitcell::AtomData>& atoms,
                     const std::vector<ModuleBase::Vector3<double>>& grid_positions,
                     const std::vector<double>& ions,
                     const std::vector<double>& charge,
                     const ModuleSccs::SccsConfig& config,
                     const ModuleSccs::SccsResponse& response,
                     SccsOutput& result)
{
    std::vector<double> polarization;
    ModuleSccs::cont_pol_charge(charge, response, basis, cell.tpiba, polarization);
    PccOutput candidate;
    candidate.slab = config.boundary == ModuleSccs::Boundary::Pcc2d;
    candidate.axis = config.pcc_2d_axis;
    unitcell::OrthogonalCell geometry;
    unitcell::SlabCell slab;
    Pcc0dParameters parameters;
    Pcc2dParameters slab_parameters;
    std::vector<ModuleBase::Vector3<double>> positions;
    std::vector<double> ionic_charges;
    for (const unitcell::AtomData& atom : atoms)
    {
        positions.push_back(atom.position);
        ionic_charges.push_back(atom.valence_charge);
    }
    std::vector<ModuleBase::Vector3<double>> grid = grid_positions;
    // The same geometry as the SCCS operator of this evaluation.
    if (candidate.slab)
    {
        pcc_2d_geometry(cell, atoms, candidate.axis, slab, slab_parameters);
        candidate.coordinate = slab.origin;
        candidate.normal = slab.normal;
        pcc_2d_project(slab, grid);
        pcc_2d_project(slab, positions);
    }
    else
    {
        pcc_0d_geometry(cell, atoms, geometry, parameters);
        candidate.origin = geometry.origin;
        pcc_0d_relative(geometry, grid);
        pcc_0d_relative(geometry, positions);
    }
    std::vector<double> electrons(charge.size());
    for (std::size_t i = 0; i < charge.size(); ++i)
    {
        electrons[i] = charge[i] - ions[i];
    }
    candidate.moments[0] = ModuleSccs::pool_moments(charge, grid, basis);
    candidate.moments[1] = ModuleSccs::pool_moments(electrons, grid, basis);
    candidate.moments[2] = ModuleSccs::pool_moments(polarization, grid, basis);
    const int atom_count = atoms.size();
    const double* ionic_charge_data = ionic_charges.data();
    const ModuleBase::Vector3<double>* ionic_position_data = positions.data();
    const ChargeMoments point_ions = charge_moments(ionic_charge_data, ionic_position_data, atom_count, 1.0);
    candidate.moments[1] = add_charge_moments(candidate.moments[1], point_ions);
    candidate.moments[3] = add_charge_moments(candidate.moments[0], candidate.moments[2]);
    if (candidate.slab)
    {
        candidate.smooth_energy = pcc_2d_energy(candidate.moments[0], slab_parameters);
    }
    else
    {
        candidate.smooth_energy = pcc_0d_energy(candidate.moments[0], parameters);
    }
    result.pcc = candidate;
    result.expected_charge = -(1.0 - 1.0 / config.cavity.epsilon_bulk) * candidate.moments[0].charge;
}
} // namespace elecstate
