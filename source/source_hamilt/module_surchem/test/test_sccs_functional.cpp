#ifdef __MPI
#include "source_base/parallel_global.h"
#include <mpi.h>
#endif

#include "../sccs/sccs_functional.h"
#include "../sccs/sccs_pw_coulomb.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <vector>

namespace
{

ModuleSccs::ElectrostaticFunctionalResult fixed_dielectric_functional(
    const std::vector<double>& solute_charge,
    const double epsilon,
    const double volume_element,
    const ModuleSccs::PeriodicCoulombOperator& coulomb)
{
    ModuleSccs::ElectrostaticField vacuum_field;
    coulomb.apply(solute_charge, vacuum_field);
    std::vector<double> screened_charge(solute_charge.size());
    for (std::size_t index = 0; index < solute_charge.size(); ++index)
    {
        screened_charge[index] = solute_charge[index] / epsilon;
    }
    ModuleSccs::ElectrostaticField dielectric_field;
    coulomb.apply(screened_charge, dielectric_field);
    const ModuleSccs::SerialChargeReduction reduction;
    return ModuleSccs::evaluate_electrostatic_functional(solute_charge,
                                                          dielectric_field,
                                                          vacuum_field,
                                                          std::vector<double>(solute_charge.size(), 0.0),
                                                          volume_element,
                                                          reduction);
}

TEST(SccsFunctional, FixedDielectricEnergyMatchesElectronPotentialDerivative)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    const double cell_length = 10.0;
    basis.initgrids(cell_length, lattice, 20.0);
    basis.initparameters(false, 20.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();

    std::vector<std::complex<double>> source_g(basis.npw);
    std::vector<std::complex<double>> direction_g(basis.npw);
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        if (std::abs(std::abs(basis.gdirect[ig].x) - 1.0) < 1.0e-12
            && std::abs(basis.gdirect[ig].y) < 1.0e-12
            && std::abs(basis.gdirect[ig].z) < 1.0e-12)
        {
            source_g[ig] = 2.0e-3;
            direction_g[ig] = 7.0e-4;
        }
    }
    std::vector<double> solute_charge(basis.nrxx);
    std::vector<double> electron_direction(basis.nrxx);
    basis.recip2real(source_g.data(), solute_charge.data());
    basis.recip2real(direction_g.data(), electron_direction.data());

    const double volume = cell_length * cell_length * cell_length;
    const double volume_element = volume / static_cast<double>(basis.nxyz);
    const double epsilon = 5.0;
    const ModuleSccs::PeriodicCoulombOperator coulomb(basis,
                                                      ModuleBase::TWO_PI / cell_length);
    const ModuleSccs::ElectrostaticFunctionalResult center
        = fixed_dielectric_functional(solute_charge, epsilon, volume_element, coulomb);

    const double step = 1.0e-4;
    std::vector<double> plus(solute_charge.size());
    std::vector<double> minus(solute_charge.size());
    double potential_derivative = 0.0;
    for (std::size_t index = 0; index < solute_charge.size(); ++index)
    {
        plus[index] = solute_charge[index] - step * electron_direction[index];
        minus[index] = solute_charge[index] + step * electron_direction[index];
        potential_derivative
            += center.electron_potential[index] * electron_direction[index] * volume_element;
    }
    const double finite_difference
        = (fixed_dielectric_functional(plus, epsilon, volume_element, coulomb).reaction_energy
           - fixed_dielectric_functional(minus, epsilon, volume_element, coulomb).reaction_energy)
          / (2.0 * step);

    EXPECT_LT(center.reaction_energy, 0.0);
    EXPECT_NEAR(finite_difference, potential_derivative, 1.0e-11);
}

TEST(SccsFunctional, IncludesDielectricCavityDerivativePotential)
{
    const std::vector<double> charge(1, 0.0);
    ModuleSccs::ElectrostaticField dielectric;
    dielectric.potential.assign(1, 2.0);
    dielectric.gradient.assign(1, ModuleBase::Vector3<double>(1.0, 2.0, 2.0));
    ModuleSccs::ElectrostaticField vacuum = dielectric;
    const ModuleSccs::SerialChargeReduction reduction;
    const ModuleSccs::ElectrostaticFunctionalResult result
        = ModuleSccs::evaluate_electrostatic_functional(charge,
                                                        dielectric,
                                                        vacuum,
                                                        std::vector<double>(1, -4.0),
                                                        1.0,
                                                        reduction);
    EXPECT_NEAR(result.electron_potential[0], 36.0 / (8.0 * ModuleBase::PI), 1.0e-14);
}

} // namespace

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_split(MPI_COMM_WORLD, 0, 1, &POOL_WORLD);
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Comm_free(&POOL_WORLD);
    MPI_Finalize();
#endif
    return result;
}
