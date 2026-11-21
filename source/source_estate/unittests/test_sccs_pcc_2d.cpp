#include "source_estate/module_pot/pot_sccs.h"
#include "source_estate/module_pot/pot_pcc.h"
#include "source_estate/module_pot/sccs_pcc_2d.h"
#include "source_basis/module_pw/pw_grid_geometry.h"
#include "source_cell/cell_tools.h"
#include "source_hamilt/module_sccs/sccs_coulomb.h"
#include "source_hamilt/module_sccs/sccs_functional.h"
#include "source_hamilt/module_sccs/sccs_ion_charge.h"
#include "source_hamilt/module_sccs/sccs_ion_force.h"
#include "source_hamilt/module_sccs/sccs_response.h"
#include "sccs_cell_test.h"
#include "source_io/module_parameter/input_parameter.h"

namespace
{
const elecstate::SccsResume no_resume = elecstate::SccsResume();
}


class SccsPcc2dTest : public SccsTest::CellTest
{
protected:
    void SetUp() override
    {
        cell.latvec.e33 = 1.5;
        set_up_basis(cell.latvec);
        set_up_cell();
    }
};

// With epsilon 1 the solvent adds nothing, so the point-ion PCC is applied
// exactly once; with a dielectric it adds only the reaction and cavity terms.
TEST_F(SccsPcc2dTest, SolventAddsOnlyReactionAndCavityTermsToThePointIonPcc)
{
    const double epsilons[] = {1.0, 5.0};
    for (double epsilon : epsilons)
    {
        Input_para input;
        input.assume_isolated = "pcc_2d";
        input.sccs_epsilon = epsilon;
        input.sccs_rho_min = 1e-3;
        ModuleSccs::SccsConfig config;
        ModuleSccs::PolarizationSolverParameters solver;
        elecstate::sccs_cfg_input(input, config, solver);
        const std::vector<unitcell::AtomData> atoms = unitcell::get_atom_data(cell.atoms, cell.ntype, cell.lat0);
        std::vector<ModuleBase::Vector3<double>> grid;
        ModulePW::grid_positions(basis, cell.latvec, cell.lat0, grid);
        std::unique_ptr<ModuleSccs::CoulombOperator> coulomb;
        elecstate::pcc_2d_operator(cell, basis, atoms, input.pcc_2d_axis, grid, coulomb);
        std::vector<double> ions;
        ModuleSccs::gauss_ion_rho(atoms, basis, tpiba,
                                  ModuleSccs::gaussian_ion_spread, ions);
        std::vector<double> solute_charge(basis.nrxx);
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            solute_charge[ir] = ions[ir] - density[ir];
        }
        const std::vector<double> cold;
        ModuleSccs::SccsResponse response;
        ModuleSccs::solve_response(density, solute_charge, config.cavity, SccsTest::no_probe, solver, cold,
                                   basis, tpiba, *coulomb, response);
        ModuleSccs::FunctionalResult functional;
        ModuleSccs::eval_functional(solute_charge, response, config, SccsTest::no_probe, basis, tpiba,
                                    *coulomb, functional);
        elecstate::PotSccs solvent(&basis, config, solver, 0.6, no_resume);
        elecstate::PotPcc correction(&basis, elecstate::PotPcc::Dimension::slab, input.pcc_2d_axis, 0.6);
        ModuleBase::matrix potential(1, basis.nrxx);
        correction.cal_v_eff(&charge, &cell, potential);
        ModuleBase::matrix point_potential = potential;
        const double point_energy = correction.get_energy();
        solvent.cal_v_eff(&charge, &cell, potential);
        const double reaction_energy_rydberg = 2.0 * functional.reaction_energy;
        EXPECT_NEAR(solvent.get_energy(), reaction_energy_rydberg, 1e-12);
        const double combined_energy = solvent.get_energy() + correction.get_energy();
        const double expected_energy = point_energy + 2.0 * functional.reaction_energy;
        EXPECT_NEAR(combined_energy, expected_energy, 1e-12);
        std::vector<ModuleBase::Vector3<double>> reaction_force;
        ModuleSccs::gauss_ion_force(atoms, functional.reaction_potential, basis,
                                    tpiba, ModuleSccs::gaussian_ion_spread, reaction_force);
        ModuleBase::matrix point_force(cell.nat, 3);
        correction.add_force(cell, point_force);
        ModuleBase::matrix combined_force = point_force;
        solvent.add_sol_force(cell, combined_force);
        for (int axis = 0; axis < 3; ++axis)
        {
            const double expected = point_force(0, axis) + 2.0 * reaction_force[0][axis];
            EXPECT_NEAR(combined_force(0, axis), expected, 1e-12);
        }
        const std::vector<double>* electrostatic = solvent.sol_elec_pot();
        ASSERT_NE(electrostatic, nullptr);
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            const double expected = point_potential(0, ir) + 2.0 * functional.electron_potential[ir];
            EXPECT_NEAR(potential(0, ir), expected, 1e-12);
            const double expected_potential = -2.0 * functional.reaction_potential[ir];
            EXPECT_NEAR((*electrostatic)[ir], expected_potential, 1e-12);
        }
    }
}

TEST_F(SccsPcc2dTest, AcceptsThePccSlabAlignmentTolerance)
{
    // A tilt of 1e-8 relative to the open vector, as from rounded lattice
    // digits, passes both PCC and the SCCS operator; a real tilt stops both.
    cell.latvec.e31 = 1.5e-8;
    const std::vector<unitcell::AtomData> atoms = unitcell::get_atom_data(cell.atoms, cell.ntype, cell.lat0);
    std::vector<ModuleBase::Vector3<double>> grid;
    ModulePW::grid_positions(basis, cell.latvec, cell.lat0, grid);
    std::unique_ptr<ModuleSccs::CoulombOperator> coulomb;
    elecstate::pcc_2d_operator(cell, basis, atoms, 2, grid, coulomb);
    EXPECT_NE(coulomb.get(), nullptr);
    elecstate::PotPcc correction(&basis, elecstate::PotPcc::Dimension::slab, 2, 0.6);
    ModuleBase::matrix potential(1, basis.nrxx);
    correction.cal_v_eff(&charge, &cell, potential);
    cell.latvec.e31 = 1.5e-4;
    ModulePW::grid_positions(basis, cell.latvec, cell.lat0, grid);
    std::unique_ptr<ModuleSccs::CoulombOperator> tilted;
    EXPECT_EXIT(elecstate::pcc_2d_operator(cell, basis, atoms, 2, grid, tilted),
                testing::ExitedWithCode(1),
                "");
}
