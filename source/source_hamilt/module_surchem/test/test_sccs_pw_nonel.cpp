#ifdef __MPI
#include "source_base/parallel_global.h"
#include <mpi.h>
#endif

#include "../sccs/sccs_pw_nonel.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <vector>

namespace
{

TEST(SccsPwNonel, EnergyDerivativeMatchesPotential)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    const double length = 10.0;
    basis.initgrids(length, lattice, 30.0);
    basis.initparameters(false, 30.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();

    std::vector<std::complex<double>> solute_g(basis.npw);
    std::vector<std::complex<double>> direction_g(basis.npw);
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        if (std::abs(std::abs(basis.gdirect[ig].x) - 1.0) < 1.0e-12
            && std::abs(basis.gdirect[ig].y) < 1.0e-12
            && std::abs(basis.gdirect[ig].z) < 1.0e-12)
        {
            solute_g[ig] = 0.12;
            direction_g[ig] = 0.03;
        }
    }
    if (basis.ig_gge0 >= 0)
    {
        solute_g[basis.ig_gge0] = 0.5;
        direction_g[basis.ig_gge0] = 0.02;
    }
    std::vector<double> solute(basis.nrxx);
    std::vector<double> direction(basis.nrxx);
    basis.recip2real(solute_g.data(), solute.data());
    basis.recip2real(direction_g.data(), direction.data());

    const double volume_element = length * length * length / static_cast<double>(basis.nxyz);
    ModuleSccs::NonElectrostaticParameters parameters;
    parameters.surface_tension = 0.7;
    parameters.pressure = -0.04;
    parameters.surface_regularization = 1.0e-3;
    const ModuleSccs::SerialChargeReduction reduction;
    const std::vector<double> derivative(solute.size(), 1.0);
    const ModuleSccs::NonElectrostaticResult center
        = ModuleSccs::evaluate_pw_non_electrostatic(basis,
                                                    ModuleBase::TWO_PI / length,
                                                    volume_element,
                                                    parameters,
                                                    solute,
                                                    derivative,
                                                    reduction);

    const double step = 1.0e-5;
    std::vector<double> plus(solute.size());
    std::vector<double> minus(solute.size());
    double analytic = 0.0;
    for (std::size_t index = 0; index < solute.size(); ++index)
    {
        plus[index] = solute[index] + step * direction[index];
        minus[index] = solute[index] - step * direction[index];
        analytic += center.density_potential[index] * direction[index] * volume_element;
    }
    const ModuleSccs::NonElectrostaticResult plus_result
        = ModuleSccs::evaluate_pw_non_electrostatic(basis,
                                                    ModuleBase::TWO_PI / length,
                                                    volume_element,
                                                    parameters,
                                                    plus,
                                                    derivative,
                                                    reduction);
    const ModuleSccs::NonElectrostaticResult minus_result
        = ModuleSccs::evaluate_pw_non_electrostatic(basis,
                                                    ModuleBase::TWO_PI / length,
                                                    volume_element,
                                                    parameters,
                                                    minus,
                                                    derivative,
                                                    reduction);
    const double finite_difference
        = ((plus_result.surface_energy + plus_result.volume_energy)
           - (minus_result.surface_energy + minus_result.volume_energy))
          / (2.0 * step);
    EXPECT_NEAR(finite_difference, analytic, 1.0e-7);
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
