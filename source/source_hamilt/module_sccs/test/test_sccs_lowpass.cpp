#include "sccs_test.h"
#include "../sccs_lowpass.h"
#include "../sccs_cav_deriv.h"
#include "../sccs_functional.h"
#include "../sccs_parameters.h"
#include "../sccs_coul_0d.h"
#include "../sccs_coul_2d.h"
#include "../sccs_response.h"

#include "source_base/parallel_reduce.h"
#include "source_basis/module_pw/pw_grid_geometry.h"


class SccsLowpassTest : public SccsTest::PwTest
{
protected:
    void check_electron_derivative(ModuleSccs::CoulombOperator& coulomb, ModuleSccs::Boundary boundary)
    {
        ModuleSccs::SccsConfig config = ModuleSccs::preset_config(ModuleSccs::Preset::WaterNeutral);
        config.boundary = boundary;
        config.cavity.epsilon_bulk = 5.0;
        // Strong filtering on the first shell makes omitting either transpose
        // detectable with a small FFT grid.
        config.cavity.lowpass_p1 = basis.ggecut;
        config.cavity.lowpass_p2 = 0.5;
        config.surface_tension = 0.0;
        config.pressure = 0.0;
        ModuleSccs::PolarizationSolverParameters solver;
        solver.tolerance_rms = 1e-13;
        solver.tolerance_max = 1e-13;
        const std::vector<double> mode = cosine_mode(0);
        std::vector<double> density(basis.nrxx);
        std::vector<double> charge(basis.nrxx);
        std::vector<double> direction(basis.nrxx);
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            density[ir] = 1e-3 + 1e-4 * mode[ir];
            charge[ir] = 5e-4 + 2e-4 * mode[ir];
            direction[ir] = 1e-4 * mode[ir];
        }
        const std::vector<double> cold;
        ModuleSccs::SccsResponse center_response;
        ModuleSccs::solve_response(density, charge, config.cavity, solver, cold,
                                   basis, tpiba, coulomb, center_response);
        EXPECT_TRUE(center_response.polarization.gradient.empty());
        ModuleSccs::FunctionalResult center;
        ModuleSccs::eval_functional(charge, center_response, config, basis,
                                    tpiba, coulomb, center);
        double predicted = 0.0;
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            predicted += center.electron_potential[ir] * direction[ir];
        }
        Parallel_Reduce::reduce_pool(predicted);
        predicted *= basis.omega / basis.nxyz;
        const double steps[] = {1e-3, 5e-4};
        for (double step : steps)
        {
            double energies[2];
            for (int side = 0; side < 2; ++side)
            {
                const double sign = side == 0 ? 1.0 : -1.0;
                std::vector<double> displaced_density(basis.nrxx);
                std::vector<double> displaced_charge(basis.nrxx);
                for (int ir = 0; ir < basis.nrxx; ++ir)
                {
                    const double change = sign * step * direction[ir];
                    displaced_density[ir] = density[ir] + change;
                    displaced_charge[ir] = charge[ir] - change;
                }
                ModuleSccs::SccsResponse response;
                ModuleSccs::solve_response(displaced_density, displaced_charge,
                                           config.cavity, solver, cold, basis, tpiba, coulomb, response);
                ModuleSccs::FunctionalResult functional;
                ModuleSccs::eval_functional(displaced_charge, response, config,
                                            basis, tpiba, coulomb, functional);
                energies[side] = functional.reaction_energy;
            }
            const double finite_difference = (energies[0] - energies[1]) / (2.0 * step);
            EXPECT_GT(std::abs(predicted), 1e-5);
            EXPECT_NEAR(predicted, finite_difference, 1e-8);
        }
    }

    std::vector<ModuleBase::Vector3<double>> centered_positions()
    {
        const ModuleBase::Matrix3 lattice;
        std::vector<ModuleBase::Vector3<double>> positions;
        ModulePW::grid_positions(basis, lattice, length, positions);
        const double half_length = 0.5 * length;
        const ModuleBase::Vector3<double> center(half_length, half_length, half_length);
        for (ModuleBase::Vector3<double>& position : positions)
        {
            position -= center;
        }
        return positions;
    }
};

TEST_F(SccsLowpassTest, PccElectronPotentialDifferentiatesDiscreteEnergy)
{
    elecstate::Pcc0dParameters molecule;
    molecule.length = length;
    const std::vector<ModuleBase::Vector3<double>> positions = centered_positions();
    ModuleSccs::Pcc0dCoulombOperator coulomb(basis, tpiba, positions, molecule);
    check_electron_derivative(coulomb, ModuleSccs::Boundary::Pcc0d);
    elecstate::Pcc2dParameters slab;
    slab.length = length;
    slab.area = length * length;
    std::vector<ModuleBase::Vector3<double>> projected = centered_positions();
    for (ModuleBase::Vector3<double>& position : projected)
    {
        position.x = 0.0;
        position.y = 0.0;
    }
    const ModuleBase::Vector3<double> normal(0.0, 0.0, 1.0);
    ModuleSccs::Pcc2dCoulombOperator slab_coulomb(basis, tpiba, projected, normal, slab);
    check_electron_derivative(slab_coulomb, ModuleSccs::Boundary::Pcc2d);
}

TEST_F(SccsLowpassTest, FilterUsesDensityCutoff)
{
    ModuleSccs::SccsConfig config = ModuleSccs::preset_config(ModuleSccs::Preset::WaterNeutral);
    std::vector<double> filter = ModuleSccs::switch_filter(config.cavity, basis);
    EXPECT_TRUE(filter.empty());
    config.cavity.lowpass_p1 = 10.0;
    config.cavity.lowpass_p2 = 5.0;
    filter = ModuleSccs::switch_filter(config.cavity, basis);
    ASSERT_EQ(filter.size(), static_cast<std::size_t>(basis.npw));
    const double half_cutoff = basis.ggecut / 2.0;
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        if (basis.gg[ig] == 0.0)
        {
            EXPECT_NEAR(filter[ig], 1.0, 1e-12);
        }
        if (basis.gg[ig] > half_cutoff)
        {
            EXPECT_LT(filter[ig], 0.5);
        }
        else
        {
            EXPECT_GE(filter[ig], 0.5);
        }
    }
}
