#include "sccs_test.h"
#include "../sccs_lowpass.h"
#include "../sccs_cav_deriv.h"
#include "../sccs_functional.h"
#include "../sccs_parameters.h"
#include "../sccs_coul_0d.h"
#include "../sccs_coul_2d.h"
#include "../sccs_response.h"
#include "../sccs_sa.h"

#include "source_base/parallel_reduce.h"
#include "source_basis/module_pw/pw_grid_geometry.h"


class SccsLowpassTest : public SccsTest::PwTest
{
protected:
    void check_electron_derivative(ModuleSccs::CoulombOperator& coulomb, ModuleSccs::Boundary boundary)
    {
        const bool filled = false;
        check_derivative(coulomb, boundary, filled);
    }

    // A filled cavity puts the solvent-aware adjoint between dE/ds_sa and dE/dn.
    void check_derivative(ModuleSccs::CoulombOperator& coulomb, ModuleSccs::Boundary boundary, bool filled)
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
        std::vector<double> probe_kernel;
        if (filled)
        {
            // The probe solute fraction is about 0.67 everywhere, inside the filling step.
            config.cavity.solvent_aware.solvent_radius = 1.5;
            config.cavity.solvent_aware.filling_threshold = 0.67;
            config.cavity.solvent_aware.filling_spread = 0.05;
            const ModuleBase::Matrix3 lattice;
            probe_kernel = ModuleSccs::sa_probe_kernel(basis, lattice, length, config.cavity.solvent_aware);
        }
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
        ModuleSccs::solve_response(density, charge, config.cavity, probe_kernel, solver, cold,
                                   basis, tpiba, coulomb, center_response);
        EXPECT_TRUE(center_response.polarization.gradient.empty());
        if (filled)
        {
            EXPECT_GT(center_response.solvent_aware.volume, 50.0);
        }
        ModuleSccs::FunctionalResult center;
        ModuleSccs::eval_functional(charge, center_response, config, probe_kernel, basis,
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
                                           config.cavity, probe_kernel, solver, cold, basis, tpiba, coulomb,
                                           response);
                ModuleSccs::FunctionalResult functional;
                ModuleSccs::eval_functional(displaced_charge, response, config, probe_kernel,
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

// PCC 0D with the solvent-aware filling.
TEST_F(SccsLowpassTest, FilledPccElectronPotentialDifferentiatesDiscreteEnergy)
{
    elecstate::Pcc0dParameters molecule;
    molecule.length = length;
    const std::vector<ModuleBase::Vector3<double>> positions = centered_positions();
    ModuleSccs::Pcc0dCoulombOperator coulomb(basis, tpiba, positions, molecule);
    const bool filled = true;
    check_derivative(coulomb, ModuleSccs::Boundary::Pcc0d, filled);
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

// The periodic chain path takes no FFT derivative of s, so lowpass parameters
// set in code (the INPUT reader rejects them without PCC) leave the periodic
// response unchanged.
TEST_F(SccsLowpassTest, PeriodicResponseIgnoresTheLowpass)
{
    ModuleSccs::CavityParameters cavity;
    cavity.density_min = 1e-4;
    cavity.density_max = 5e-3;
    cavity.epsilon_bulk = 5.0;
    ModuleSccs::PolarizationSolverParameters solver;
    const std::vector<double> mode = cosine_mode(0);
    std::vector<double> density(basis.nrxx);
    std::vector<double> charge(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        density[ir] = 1e-3 + 1e-4 * mode[ir];
        charge[ir] = 2e-4 * mode[ir];
    }
    const std::vector<double> cold;
    ModuleSccs::SccsResponse plain;
    ModuleSccs::solve_response(density, charge, cavity, solver, cold, basis, tpiba, plain);
    cavity.lowpass_p1 = basis.ggecut;
    cavity.lowpass_p2 = 0.5;
    ModuleSccs::SccsResponse filtered;
    ModuleSccs::solve_response(density, charge, cavity, solver, cold, basis, tpiba, filtered);
    EXPECT_FALSE(filtered.polarization.gradient.empty());
    ASSERT_EQ(filtered.cavity_potential.size(), plain.cavity_potential.size());
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_EQ(filtered.cavity_potential[ir], plain.cavity_potential[ir]);
    }
}
