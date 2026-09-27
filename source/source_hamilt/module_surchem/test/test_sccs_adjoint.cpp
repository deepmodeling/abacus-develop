#include "../sccs/sccs_adjoint.h"
#include "../sccs/sccs_functional.h"
#include "../sccs/sccs_periodic.h"
#include "../sccs/sccs_pw_charge.h"
#include "../pcc/sccs_pcc_coulomb.h"
#include "../pcc/sccs_pcc_2d_coulomb.h"
#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_base/parallel_global.h"
#include "source_basis/module_pw/pw_basis.h"

#include <gtest/gtest.h>
#include <cmath>
#include <memory>

namespace
{

class SccsAdjointTest : public testing::Test
{
  protected:
    void SetUp() override
    {
        basis.reset(new ModulePW::PW_Basis("cpu", "double"));
#ifdef __MPI
        basis->initmpi(1, 0, POOL_WORLD);
#endif
        const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        const double length = 12.0;
        basis->initgrids(length, lattice, 30.0);
        basis->initparameters(false, 30.0, 1, false);
        basis->setuptransform();
        basis->collect_local_pw();
        tpiba = ModuleBase::TWO_PI / length;
        dv = length * length * length / basis->nxyz;
        positions = ModuleSccs::pw_grid_positions(*basis, lattice, length);
        geometry = ModuleSccs::pcc_geometry(lattice, length, 1.0e-10);
        geometry_2d = ModuleSccs::pcc_2d_geometry(lattice, length, 1.0e-10);
        density.resize(basis->nrxx);
        charge.resize(basis->nrxx);
        direction.resize(basis->nrxx);
        for (int i = 0; i < basis->nrxx; ++i)
        {
            const double x = positions[i].x - 6.0;
            const double y = positions[i].y - 6.0;
            const double z = positions[i].z - 6.0;
            const double r2 = x * x + y * y + z * z;
            density[i] = 0.2 * std::exp(-r2 / 3.0);
            charge[i] = -0.02 * std::exp(-((x - 0.7) * (x - 0.7) + y * y + z * z) / 2.0);
            direction[i] = 0.01 * (0.5 + x + y) * std::exp(-r2 / 3.0);
        }
        cavity.density_min = 0.0024;
        cavity.density_max = 0.0155;
        cavity.epsilon_bulk = 78.3;
        solver.max_iterations = 400;
        solver.mixing_method = "pulay";
        solver.mixing = 0.3;
        solver.tolerance_rms = 1.0e-13;
        solver.tolerance_max = 1.0e-11;
    }

    std::unique_ptr<ModuleSccs::CoulombOperator> make_operator(const int boundary)
    {
        if (boundary == 0)
        {
            return std::unique_ptr<ModuleSccs::CoulombOperator>(
                new ModuleSccs::PeriodicCoulombOperator(*basis, tpiba));
        }
        if (boundary == 1)
        {
            return std::unique_ptr<ModuleSccs::CoulombOperator>(
                new ModuleSccs::PccCoulombOperator(*basis, tpiba, positions, dv, geometry, charge_reduction));
        }
        return std::unique_ptr<ModuleSccs::CoulombOperator>(
            new ModuleSccs::Pcc2dCoulombOperator(*basis, tpiba, positions, dv, geometry_2d, charge_reduction));
    }

    double energy(const std::vector<double>& rho,
                  const std::vector<double>& source,
                  const ModuleSccs::CoulombOperator& op,
                  const double coefficient)
    {
        const std::vector<double> initial;
        const ModuleSccs::PeriodicSccsResult response = ModuleSccs::solve_sccs_response(
            rho, source, cavity, solver, initial, *basis, tpiba, op, polarization_reduction);
        EXPECT_EQ(response.polarization.status, ModuleSccs::PolarizationStatus::Converged);
        ModuleSccs::ElectrostaticField vacuum;
        op.apply(source, vacuum);
        const ModuleSccs::ElectrostaticFunctionalResult value = ModuleSccs::evaluate_electrostatic_functional(
            source, response.polarization.field, vacuum, response.depsilon_drho, dv, charge_reduction);
        double total = value.reaction_energy;
        for (const double polarization : response.polarization.polarization_charge)
        {
            total += coefficient * polarization * dv;
        }
        return total;
    }

    std::unique_ptr<ModulePW::PW_Basis> basis;
    std::vector<ModuleBase::Vector3<double>> positions;
    std::vector<double> density;
    std::vector<double> charge;
    std::vector<double> direction;
    ModuleSccs::PccGeometry geometry;
    ModuleSccs::Pcc2dGeometry geometry_2d;
    ModuleSccs::CavityParameters cavity;
    ModuleSccs::PolarizationSolverParameters solver;
    ModuleSccs::SerialChargeReduction charge_reduction;
    ModuleSccs::SerialPolarizationReduction polarization_reduction;
    double tpiba = 0.0;
    double dv = 0.0;
};

TEST_F(SccsAdjointTest, CoulombGradientAdjointsPreserveInnerProducts)
{
    for (int boundary = 0; boundary < 3; ++boundary)
    {
        SCOPED_TRACE(boundary);
        const auto op = make_operator(boundary);
        ModuleSccs::ElectrostaticField field;
        op->apply(charge, field);
        std::vector<ModuleBase::Vector3<double>> probe(basis->nrxx);
        for (int i = 0; i < basis->nrxx; ++i)
        {
            probe[i].x = 0.3 + std::sin(positions[i].x);
            probe[i].y = -0.5 + std::cos(positions[i].y);
            probe[i].z = positions[i].z - 6.0;
        }
        std::vector<double> transpose;
        op->apply_gradient_adjoint(probe, transpose);
        double left = 0.0;
        double right = 0.0;
        for (int i = 0; i < basis->nrxx; ++i)
        {
            left += (field.gradient[i] * probe[i]) * dv;
            right += charge[i] * transpose[i] * dv;
        }
        EXPECT_NEAR(left, right, 2.0e-12);
    }
}

TEST_F(SccsAdjointTest, GradientOnlyMatchesFullFieldForAllBoundaries)
{
    for (int boundary = 0; boundary < 3; ++boundary)
    {
        SCOPED_TRACE(boundary);
        const auto op = make_operator(boundary);
        ModuleSccs::ElectrostaticField field;
        std::vector<ModuleBase::Vector3<double>> gradient;
        for (const double scale : {1.0, 0.0, -0.7})
        {
            auto source = charge;
            for (double& value : source)
            {
                value *= scale;
            }
            op->apply_gradient(source, gradient);
            op->apply(source, field);
            ASSERT_EQ(gradient.size(), field.gradient.size());
            for (std::size_t i = 0; i < gradient.size(); ++i)
            {
                for (int d = 0; d < 3; ++d)
                {
                    EXPECT_DOUBLE_EQ(gradient[i][d], field.gradient[i][d]);
                }
            }
        }
    }
}

TEST_F(SccsAdjointTest, SharpWaterCavityAndSourceDerivativesMatchEnergy)
{
    for (int boundary = 0; boundary < 3; ++boundary)
    {
        SCOPED_TRACE(boundary);
        const auto op = make_operator(boundary);
        const double coefficient = boundary == 2 ? 0.017 : 0.0;
        const std::vector<double> initial;
        const auto response = ModuleSccs::solve_sccs_response(
            density, charge, cavity, solver, initial, *basis, tpiba, *op, polarization_reduction);
        ASSERT_EQ(response.polarization.status, ModuleSccs::PolarizationStatus::Converged);
        ModuleSccs::ElectrostaticField vacuum;
        op->apply(charge, vacuum);
        auto functional = ModuleSccs::evaluate_electrostatic_functional(
            charge, response.polarization.field, vacuum, response.depsilon_drho, dv, charge_reduction);
        const auto adjoint = ModuleSccs::evaluate_discrete_electrostatic_derivative(
            charge, response, vacuum, *basis, tpiba, coefficient, solver, initial,
            *op, polarization_reduction, functional);
        EXPECT_LE(adjoint.residual_rms, solver.tolerance_rms);
        EXPECT_LE(adjoint.residual_max, solver.tolerance_max);
        for (int mode = 0; mode < 3; ++mode)
        {
            SCOPED_TRACE(mode);
            double analytic = 0.0;
            for (int i = 0; i < basis->nrxx; ++i)
            {
                double potential = functional.electron_potential[i];
                if (mode == 0)
                {
                    potential = functional.charge_potential[i];
                }
                else if (mode == 1)
                {
                    potential += functional.charge_potential[i];
                }
                analytic += potential * direction[i] * dv;
            }
            double differences[2];
            for (int resolution = 0; resolution < 2; ++resolution)
            {
                const double h = resolution == 0 ? 0.002 : 0.001;
                double energies[2];
                for (int side = 0; side < 2; ++side)
                {
                    auto rho = density;
                    auto source = charge;
                    const double shift = (2 * side - 1) * h;
                    for (int i = 0; i < basis->nrxx; ++i)
                    {
                        if (mode == 0)
                        {
                            source[i] += shift * direction[i];
                        }
                        else
                        {
                            rho[i] += shift * direction[i];
                            if (mode == 2)
                            {
                                source[i] -= shift * direction[i];
                            }
                        }
                    }
                    energies[side] = energy(rho, source, *op, coefficient);
                }
                differences[resolution] = (energies[1] - energies[0]) / (2.0 * h);
            }
            const double extrapolated = (4.0 * differences[1] - differences[0]) / 3.0;
            EXPECT_NEAR(analytic, extrapolated, 2.0e-8);
        }
    }
}

} // namespace

int main(int argc, char** argv)
{
#ifdef __MPI
    int process_count = 1;
    int thread_count = 1;
    int rank = 0;
    Parallel_Global::read_pal_param(argc, argv, process_count, thread_count, rank);
    POOL_WORLD = MPI_COMM_WORLD;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif
    return result;
}
