#ifdef __MPI
#include "source_base/parallel_global.h"
#endif

#include "../sccs_pcc_2d_coulomb.h"
#include "../sccs_functional.h"
#include "../sccs_periodic.h"
#include "../sccs_pw_charge.h"
#include "../sccs_pw_reduction.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace
{

int test_process_count = 1;
int test_rank = 0;

class PresetArrayReduction : public ModuleSccs::ChargeReduction
{
  public:
    explicit PresetArrayReduction(const std::vector<double>& remote) : remote_(remote)
    {
    }

    void reduce_sum(double&) const override
    {
    }

    void reduce_sum(double* values, const int count) const override
    {
        if (count != static_cast<int>(remote_.size()))
        {
            throw std::invalid_argument("preset reduction size mismatch");
        }
        for (int index = 0; index < count; ++index)
        {
            values[index] += remote_[index];
        }
    }

  private:
    std::vector<double> remote_;
};

class SccsPcc2dCoulombTest : public testing::Test
{
  protected:
    SccsPcc2dCoulombTest()
        : basis_("cpu", "double"), polarization_reduction_(test_process_count)
    {
    }

    void SetUp() override
    {
#ifdef __MPI
        basis_.initmpi(test_process_count, test_rank, POOL_WORLD);
#endif
        const ModuleBase::Matrix3 lattice(0.8, 0.0, 0.0,
                                          0.0, 1.2, 0.0,
                                          0.0, 0.0, 1.0);
        lattice_scale_ = 10.0;
        basis_.initgrids(lattice_scale_, lattice, 80.0);
        basis_.initparameters(false, 80.0, 1, false);
        basis_.setuptransform();
        basis_.collect_local_pw();
        geometry_ = ModuleSccs::pcc_2d_geometry(lattice, lattice_scale_, 1.0e-10);
        volume_ = geometry_.parameters.periodic_area
                  * geometry_.parameters.cell_length_y;
        volume_element_ = volume_ / static_cast<double>(basis_.nxyz);
        positions_ = ModuleSccs::pw_grid_positions(basis_, lattice, lattice_scale_);
    }

    ModuleSccs::ElectrostaticFunctionalResult evaluate_reaction(
        const std::vector<double>& electron_density,
        const std::vector<double>& ionic_density,
        const ModuleSccs::CavityParameters& cavity,
        const ModuleSccs::Pcc2dCoulombOperator& coulomb) const
    {
        std::vector<double> solute_charge(electron_density.size());
        for (std::size_t index = 0; index < electron_density.size(); ++index)
        {
            solute_charge[index] = ionic_density[index] - electron_density[index];
        }
        ModuleSccs::PolarizationSolverParameters solver;
        solver.max_iterations = 1000;
        solver.mixing = 0.2;
        solver.tolerance_rms = 1.0e-12;
        solver.tolerance_max = 1.0e-11;
        const ModuleSccs::PeriodicSccsResult response
            = ModuleSccs::solve_sccs_response(electron_density,
                                               solute_charge,
                                               cavity,
                                               solver,
                                               std::vector<double>(),
                                               basis_,
                                               ModuleBase::TWO_PI / lattice_scale_,
                                               coulomb,
                                               polarization_reduction_);
        if (response.polarization.status != ModuleSccs::PolarizationStatus::Converged)
        {
            throw std::runtime_error("two-dimensional SCCS test response did not converge");
        }
        ModuleSccs::ElectrostaticField vacuum_field;
        coulomb.apply(solute_charge, vacuum_field);
        return ModuleSccs::evaluate_electrostatic_functional(solute_charge,
                                                              response.polarization.field,
                                                              vacuum_field,
                                                              response.depsilon_drho,
                                                              volume_element_,
                                                              reduction_);
    }

    ModulePW::PW_Basis basis_;
    double lattice_scale_ = 0.0;
    double volume_ = 0.0;
    double volume_element_ = 0.0;
    ModuleSccs::Pcc2dGeometry geometry_;
    std::vector<ModuleBase::Vector3<double>> positions_;
    ModuleSccs::PoolChargeReduction reduction_;
    ModuleSccs::PoolPolarizationReduction polarization_reduction_;
};

TEST_F(SccsPcc2dCoulombTest, ReducesYMomentsAndBuildsPlaneAverages)
{
    const std::vector<double> uniform(basis_.nrxx, 1.0 / volume_);
    const ModuleSccs::Pcc2dMoments moments
        = ModuleSccs::reduced_pcc_2d_density_moments(uniform,
                                                     positions_,
                                                     volume_element_,
                                                     geometry_,
                                                     reduction_);
    EXPECT_NEAR(moments.charge, 1.0, 1.0e-12);
    EXPECT_NEAR(moments.dipole_y, 0.0, 1.0e-14);

    std::vector<double> y_values(positions_.size());
    for (std::size_t index = 0; index < positions_.size(); ++index)
    {
        y_values[index] = positions_[index].y;
    }
    const std::vector<double> average
        = ModuleSccs::pcc_2d_plane_average(y_values, basis_, reduction_);
    ASSERT_EQ(average.size(), static_cast<std::size_t>(basis_.ny));
    for (int iy = 0; iy < basis_.ny; ++iy)
    {
        const double expected = geometry_.parameters.cell_length_y
                                * (static_cast<double>(iy) + 0.5)
                                / static_cast<double>(basis_.ny);
        EXPECT_NEAR(average[iy], expected, 2.0e-13);
    }
}

TEST_F(SccsPcc2dCoulombTest, CombinesDistributedZFragmentsWithoutChangingYProfiles)
{
    const auto fragment_geometry = [](const double origin_y) {
        ModuleSccs::Pcc2dGeometry value;
        value.parameters.periodic_area = 10.0;
        value.parameters.cell_length_y = 20.0;
        value.origin_y = origin_y;
        return value;
    };
    std::vector<double> local_density;
    std::vector<double> remote_density;
    std::vector<ModuleBase::Vector3<double>> local_positions;
    std::vector<ModuleBase::Vector3<double>> remote_positions;
    const int nx = 2;
    const int ny = 3;
    const int nz = 5;
    const int split_z = 2;
    const double fragment_volume_element = 0.4;
    for (int ix = 0; ix < nx; ++ix)
    {
        for (int iy = 0; iy < ny; ++iy)
        {
            for (int iz = 0; iz < nz; ++iz)
            {
                const ModuleBase::Vector3<double> position(
                    static_cast<double>(ix) + 0.5,
                    2.0 * (static_cast<double>(iy) + 0.5),
                    static_cast<double>(iz) + 0.5);
                const double value = 1.0e-3 * (1.0 + position.y);
                if (iz < split_z)
                {
                    local_density.push_back(value);
                    local_positions.push_back(position);
                }
                else
                {
                    remote_density.push_back(value);
                    remote_positions.push_back(position);
                }
            }
        }
    }
    const ModuleSccs::Pcc2dMoments remote_moments
        = ModuleSccs::pcc_2d_density_moments(remote_density,
                                             remote_positions,
                                             fragment_volume_element,
                                             fragment_geometry(3.0));
    const std::vector<double> remote_values{remote_moments.charge,
                                            remote_moments.dipole_y,
                                            remote_moments.quadrupole_yy};
    const PresetArrayReduction moment_reduction(remote_values);
    const ModuleSccs::Pcc2dMoments reduced
        = ModuleSccs::reduced_pcc_2d_density_moments(local_density,
                                                     local_positions,
                                                     fragment_volume_element,
                                                     fragment_geometry(3.0),
                                                     moment_reduction);
    std::vector<double> full_density = local_density;
    std::vector<ModuleBase::Vector3<double>> full_positions = local_positions;
    full_density.insert(full_density.end(), remote_density.begin(), remote_density.end());
    full_positions.insert(full_positions.end(), remote_positions.begin(), remote_positions.end());
    const ModuleSccs::Pcc2dMoments expected
        = ModuleSccs::pcc_2d_density_moments(full_density,
                                             full_positions,
                                             fragment_volume_element,
                                             fragment_geometry(3.0));
    EXPECT_NEAR(reduced.charge, expected.charge, 1.0e-12);
    EXPECT_NEAR(reduced.dipole_y, expected.dipole_y, 1.0e-12);
    EXPECT_NEAR(reduced.quadrupole_yy, expected.quadrupole_yy, 1.0e-11);

    ModulePW::PW_Basis fragment("cpu", "double");
    fragment.nx = 3;
    fragment.ny = 5;
    fragment.nz = 7;
    fragment.nplane = 3;
    fragment.startz_current = 2;
    fragment.nrxx = fragment.nx * fragment.ny * fragment.nplane;
    std::vector<double> fragment_values(fragment.nrxx);
    std::vector<double> remote_plane_sums(fragment.ny);
    for (int ix = 0; ix < fragment.nx; ++ix)
    {
        for (int iy = 0; iy < fragment.ny; ++iy)
        {
            const double profile = 0.25 + 0.5 * static_cast<double>(iy);
            remote_plane_sums[iy]
                = static_cast<double>(fragment.nx * (fragment.nz - fragment.nplane))
                  * profile;
            for (int iz_local = 0; iz_local < fragment.nplane; ++iz_local)
            {
                const int index
                    = (ix * fragment.ny + iy) * fragment.nplane + iz_local;
                fragment_values[index] = profile;
            }
        }
    }
    const PresetArrayReduction plane_reduction(remote_plane_sums);
    const std::vector<double> profile
        = ModuleSccs::pcc_2d_plane_average(fragment_values, fragment, plane_reduction);
    for (int iy = 0; iy < fragment.ny; ++iy)
    {
        EXPECT_DOUBLE_EQ(profile[iy], 0.25 + 0.5 * static_cast<double>(iy));
    }
}

TEST_F(SccsPcc2dCoulombTest, AddsCorrectionFromCurrentChargeOnEveryApplication)
{
    std::vector<double> charge(basis_.nrxx);
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        const double relative_y = ModuleSccs::pcc_2d_relative_y(positions_[ir].y, geometry_);
        charge[ir] = 2.0e-3 * std::exp(-relative_y * relative_y / 3.0)
                     - 7.0e-4 * relative_y * std::exp(-relative_y * relative_y / 2.0);
    }
    const ModuleSccs::Pcc2dCoulombOperator coulomb(basis_,
                                                   ModuleBase::TWO_PI / lattice_scale_,
                                                   positions_,
                                                   volume_element_,
                                                   geometry_,
                                                   reduction_);
    const ModuleSccs::PeriodicCoulombOperator periodic(
        basis_, ModuleBase::TWO_PI / lattice_scale_);
    ModuleSccs::ElectrostaticField field;
    ModuleSccs::ElectrostaticField periodic_field;
    coulomb.apply(charge, field);
    periodic.apply(charge, periodic_field);
    const ModuleSccs::Pcc2dMoments moments
        = ModuleSccs::reduced_pcc_2d_density_moments(charge,
                                                     positions_,
                                                     volume_element_,
                                                     geometry_,
                                                     reduction_);
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        const double relative_y = ModuleSccs::pcc_2d_relative_y(positions_[ir].y, geometry_);
        EXPECT_NEAR(field.potential[ir] - periodic_field.potential[ir],
                    ModuleSccs::pcc_2d_potential(moments,
                                                 relative_y,
                                                 geometry_.parameters),
                    2.0e-12);
        EXPECT_NEAR(field.gradient[ir].x - periodic_field.gradient[ir].x, 0.0, 1.0e-14);
        EXPECT_NEAR(field.gradient[ir].y - periodic_field.gradient[ir].y,
                    ModuleSccs::pcc_2d_potential_gradient(moments,
                                                          relative_y,
                                                          geometry_.parameters)
                        .y,
                    2.0e-12);
        EXPECT_NEAR(field.gradient[ir].z - periodic_field.gradient[ir].z, 0.0, 1.0e-14);
    }

    std::vector<double> opposite(charge.size());
    for (std::size_t index = 0; index < charge.size(); ++index)
    {
        opposite[index] = -charge[index];
    }
    ModuleSccs::ElectrostaticField opposite_field;
    coulomb.apply(opposite, opposite_field);
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        EXPECT_NEAR(opposite_field.potential[ir], -field.potential[ir], 2.0e-12);
        EXPECT_NEAR(opposite_field.gradient[ir].y, -field.gradient[ir].y, 2.0e-12);
    }
}

TEST_F(SccsPcc2dCoulombTest, UniformDielectricScreensChargeAndField)
{
    const ModuleSccs::Pcc2dCoulombOperator coulomb(basis_,
                                                   ModuleBase::TWO_PI / lattice_scale_,
                                                   positions_,
                                                   volume_element_,
                                                   geometry_,
                                                   reduction_);
    const std::vector<double> solute_charge(basis_.nrxx, 1.0 / volume_);
    const double epsilon_value = 5.0;
    const std::vector<double> epsilon(basis_.nrxx, epsilon_value);
    const std::vector<ModuleBase::Vector3<double>> grad_log_epsilon(basis_.nrxx);
    ModuleSccs::PolarizationSolverParameters solver;
    solver.max_iterations = 100;
    solver.mixing = 0.7;
    solver.tolerance_rms = 1.0e-14;
    solver.tolerance_max = 1.0e-14;
    const ModuleSccs::PolarizationResult result
        = ModuleSccs::solve_polarization(solute_charge,
                                         epsilon,
                                         grad_log_epsilon,
                                         std::vector<double>(),
                                         solver,
                                         coulomb,
                                         polarization_reduction_);
    ASSERT_EQ(result.status, ModuleSccs::PolarizationStatus::Converged);

    std::vector<double> screened_charge(solute_charge.size());
    for (std::size_t index = 0; index < solute_charge.size(); ++index)
    {
        screened_charge[index] = solute_charge[index] / epsilon_value;
        EXPECT_NEAR(result.polarization_charge[index],
                    -(1.0 - 1.0 / epsilon_value) * solute_charge[index],
                    1.0e-14);
    }
    ModuleSccs::ElectrostaticField expected;
    coulomb.apply(screened_charge, expected);
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        EXPECT_NEAR(result.field.potential[ir], expected.potential[ir], 2.0e-12);
        EXPECT_NEAR(result.field.gradient[ir].y, expected.gradient[ir].y, 2.0e-12);
    }
}

TEST_F(SccsPcc2dCoulombTest, SmoothLayeredDielectricMatchesOpenOneDimensionalField)
{
    const ModuleSccs::Pcc2dCoulombOperator coulomb(basis_,
                                                   ModuleBase::TWO_PI / lattice_scale_,
                                                   positions_,
                                                   volume_element_,
                                                   geometry_,
                                                   reduction_);
    const double width = 2.0;
    const double source_width = 1.3;
    const double amplitude = 3.0e-3;
    const double half_length = 0.5 * geometry_.parameters.cell_length_y;
    std::vector<double> solute_charge(basis_.nrxx);
    std::vector<double> epsilon(basis_.nrxx);
    std::vector<ModuleBase::Vector3<double>> grad_log_epsilon(basis_.nrxx);
    std::vector<double> reference_gradient(basis_.nrxx);
    std::vector<double> reference_polarization(basis_.nrxx);
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        const double y = ModuleSccs::pcc_2d_relative_y(positions_[ir].y, geometry_);
        const double source_exponential
            = std::exp(-y * y / (source_width * source_width));
        solute_charge[ir] = amplitude * y * source_exponential;
        const double dielectric_exponential = std::exp(-y * y / (width * width));
        epsilon[ir] = 1.0 + 3.0 * (1.0 - dielectric_exponential);
        const double grad_epsilon
            = 6.0 * y * dielectric_exponential / (width * width);
        grad_log_epsilon[ir].y = grad_epsilon / epsilon[ir];
        const double boundary_exponential
            = std::exp(-half_length * half_length / (source_width * source_width));
        const double cumulative
            = -0.5 * amplitude * source_width * source_width
              * (source_exponential - boundary_exponential);
        reference_gradient[ir] = -ModuleBase::FOUR_PI * cumulative / epsilon[ir];
        reference_polarization[ir]
            = (1.0 / epsilon[ir] - 1.0) * solute_charge[ir]
              + grad_log_epsilon[ir].y * reference_gradient[ir] / ModuleBase::FOUR_PI;
    }

    ModuleSccs::PolarizationSolverParameters solver;
    solver.max_iterations = 1000;
    solver.mixing = 0.2;
    solver.tolerance_rms = 1.0e-12;
    solver.tolerance_max = 1.0e-11;
    const ModuleSccs::PolarizationResult result
        = ModuleSccs::solve_polarization(solute_charge,
                                         epsilon,
                                         grad_log_epsilon,
                                         std::vector<double>(),
                                         solver,
                                         coulomb,
                                         polarization_reduction_);
    ASSERT_EQ(result.status, ModuleSccs::PolarizationStatus::Converged);

    double maximum_gradient_error = 0.0;
    double maximum_polarization_error = 0.0;
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        maximum_gradient_error
            = std::max(maximum_gradient_error,
                       std::abs(result.field.gradient[ir].y - reference_gradient[ir]));
        maximum_polarization_error
            = std::max(maximum_polarization_error,
                       std::abs(result.polarization_charge[ir]
                                - reference_polarization[ir]));
    }
    EXPECT_LT(maximum_gradient_error, 2.0e-4);
    EXPECT_LT(maximum_polarization_error, 3.0e-5);
}

TEST_F(SccsPcc2dCoulombTest, NonuniformCavityEnergyMatchesNeutralDensityDerivative)
{
    const ModuleSccs::Pcc2dCoulombOperator coulomb(basis_,
                                                   ModuleBase::TWO_PI / lattice_scale_,
                                                   positions_,
                                                   volume_element_,
                                                   geometry_,
                                                   reduction_);
    std::vector<double> electron_density(basis_.nrxx);
    std::vector<double> ionic_density(basis_.nrxx);
    std::vector<double> direction(basis_.nrxx);
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        const double y = ModuleSccs::pcc_2d_relative_y(positions_[ir].y, geometry_);
        electron_density[ir] = 4.0e-3 + 2.5e-2 * std::exp(-y * y / 4.0);
        const double solute_charge = 2.0e-3 * y * std::exp(-y * y / 2.25);
        ionic_density[ir] = electron_density[ir] + solute_charge;
        direction[ir]
            = 1.0e-3
              * (std::exp(-(y - 0.6) * (y - 0.6) / 2.0)
                 - std::exp(-(y + 0.6) * (y + 0.6) / 2.0));
    }
    ModuleSccs::CavityParameters cavity;
    cavity.density_min = 5.0e-3;
    cavity.density_max = 2.0e-2;
    cavity.epsilon_bulk = 4.0;
    const ModuleSccs::ElectrostaticFunctionalResult center
        = evaluate_reaction(electron_density, ionic_density, cavity, coulomb);

    double analytic = 0.0;
    double direction_charge = 0.0;
    for (std::size_t index = 0; index < electron_density.size(); ++index)
    {
        analytic += center.electron_potential[index] * direction[index] * volume_element_;
        direction_charge += direction[index] * volume_element_;
    }
    reduction_.reduce_sum(analytic);
    reduction_.reduce_sum(direction_charge);
    EXPECT_NEAR(direction_charge, 0.0, 1.0e-12);
    const double steps[4] = {4.0e-2, 2.0e-2, 1.0e-2, 5.0e-3};
    for (int step_index = 0; step_index < 4; ++step_index)
    {
        const double step = steps[step_index];
        std::vector<double> plus(electron_density.size());
        std::vector<double> minus(electron_density.size());
        for (std::size_t index = 0; index < electron_density.size(); ++index)
        {
            plus[index] = electron_density[index] + step * direction[index];
            minus[index] = electron_density[index] - step * direction[index];
        }
        const double finite_difference
            = (evaluate_reaction(plus, ionic_density, cavity, coulomb).reaction_energy
               - evaluate_reaction(minus, ionic_density, cavity, coulomb).reaction_energy)
              / (2.0 * step);
        EXPECT_NEAR(finite_difference, analytic, 1.0e-7)
            << "step " << step;
    }
}

TEST_F(SccsPcc2dCoulombTest, NonuniformCavityEnergyMatchesFixedChargeDensityDerivative)
{
    const ModuleSccs::Pcc2dCoulombOperator coulomb(basis_,
                                                   ModuleBase::TWO_PI / lattice_scale_,
                                                   positions_,
                                                   volume_element_,
                                                   geometry_,
                                                   reduction_);
    std::vector<double> electron_density(basis_.nrxx);
    std::vector<double> ionic_density(basis_.nrxx);
    std::vector<double> direction(basis_.nrxx);
    double net_charge = 0.0;
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        const double y = ModuleSccs::pcc_2d_relative_y(positions_[ir].y, geometry_);
        electron_density[ir] = 4.0e-3 + 2.5e-2 * std::exp(-y * y / 4.0);
        const double solute_charge
            = 2.0e-3 * (0.4 + y) * std::exp(-y * y / 2.25);
        ionic_density[ir] = electron_density[ir] + solute_charge;
        direction[ir]
            = 1.0e-3
              * (std::exp(-(y - 0.6) * (y - 0.6) / 2.0)
                 - std::exp(-(y + 0.6) * (y + 0.6) / 2.0));
        net_charge += solute_charge * volume_element_;
    }
    reduction_.reduce_sum(net_charge);
    ASSERT_GT(std::abs(net_charge), 1.0e-3);

    ModuleSccs::CavityParameters cavity;
    cavity.density_min = 5.0e-3;
    cavity.density_max = 2.0e-2;
    cavity.epsilon_bulk = 4.0;
    const ModuleSccs::ElectrostaticFunctionalResult center
        = evaluate_reaction(electron_density, ionic_density, cavity, coulomb);

    double analytic = 0.0;
    double direction_charge = 0.0;
    for (std::size_t index = 0; index < electron_density.size(); ++index)
    {
        analytic += center.electron_potential[index] * direction[index] * volume_element_;
        direction_charge += direction[index] * volume_element_;
    }
    reduction_.reduce_sum(analytic);
    reduction_.reduce_sum(direction_charge);
    EXPECT_NEAR(direction_charge, 0.0, 1.0e-12);
    const double steps[4] = {4.0e-2, 2.0e-2, 1.0e-2, 5.0e-3};
    for (int step_index = 0; step_index < 4; ++step_index)
    {
        const double step = steps[step_index];
        std::vector<double> plus(electron_density.size());
        std::vector<double> minus(electron_density.size());
        for (std::size_t index = 0; index < electron_density.size(); ++index)
        {
            plus[index] = electron_density[index] + step * direction[index];
            minus[index] = electron_density[index] - step * direction[index];
        }
        const double finite_difference
            = (evaluate_reaction(plus, ionic_density, cavity, coulomb).reaction_energy
               - evaluate_reaction(minus, ionic_density, cavity, coulomb).reaction_energy)
              / (2.0 * step);
        EXPECT_NEAR(finite_difference, analytic, 1.0e-7)
            << "step " << step;
    }
}

} // namespace

int main(int argc, char** argv)
{
#ifdef __MPI
    int thread_count = 1;
    Parallel_Global::read_pal_param(argc,
                                    argv,
                                    test_process_count,
                                    thread_count,
                                    test_rank);
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
