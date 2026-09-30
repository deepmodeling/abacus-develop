#ifdef __MPI
#include "source_base/parallel_global.h"
#endif

#include "../pcc/pcc_moments.h"
#include "../sccs/sccs_pcc_2d_coulomb.h"
#include "../sccs/sccs_response.h"
#include "../sccs/sccs_pw_coulomb.h"
#include "../common/pw_grid.h"
#include "../common/charge_reduction.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{

int test_process_count = 1;
int test_rank = 0;

class PresetArrayReduction : public ModuleSurchem::ChargeReduction
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

    void reduce_max(double&) const override
    {
    }

  private:
    std::vector<double> remote_;
};

class SccsPcc2dCoulombTest : public testing::Test
{
  protected:
    SccsPcc2dCoulombTest()
        : basis_("cpu", "double"), reduction_(test_process_count)
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
        geometry_ = ModulePcc::pcc_2d_geometry(lattice, lattice_scale_, 1, 1.0e-10);
        volume_ = geometry_.parameters.periodic_area
                  * geometry_.parameters.cell_length;
        volume_element_ = volume_ / static_cast<double>(basis_.nxyz);
        positions_ = ModuleSurchem::pw_grid_positions(basis_, lattice, lattice_scale_);
    }

    ModulePW::PW_Basis basis_;
    double lattice_scale_ = 0.0;
    double volume_ = 0.0;
    double volume_element_ = 0.0;
    ModulePcc::Pcc2dGeometry geometry_;
    std::vector<ModuleBase::Vector3<double>> positions_;
    ModuleSurchem::PoolChargeReduction reduction_;
};

TEST_F(SccsPcc2dCoulombTest, ReducesYMoments)
{
    const std::vector<double> uniform(basis_.nrxx, 1.0 / volume_);
    const ModulePcc::Pcc2dMoments moments
        = ModulePcc::reduced_pcc_2d_density_moments(uniform,
                                                     positions_,
                                                     volume_element_,
                                                     geometry_,
                                                     reduction_);
    EXPECT_NEAR(moments.charge, 1.0, 1.0e-12);
    // Integer FFT nodes sample [0, Ly), so their mean is half a step below Ly/2.
    const double uniform_dipole = -geometry_.parameters.cell_length / (2.0 * basis_.ny);
    EXPECT_NEAR(moments.dipole, uniform_dipole, 1.0e-14);
}

TEST_F(SccsPcc2dCoulombTest, CombinesDistributedZFragmentMoments)
{
    const auto fragment_geometry = [](const double origin_y) {
        ModulePcc::Pcc2dGeometry value;
        value.parameters.periodic_area = 10.0;
        value.parameters.cell_length = 20.0;
        value.origin = origin_y;
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
    const ModulePcc::Pcc2dMoments remote_moments
        = ModulePcc::pcc_2d_density_moments(remote_density,
                                             remote_positions,
                                             fragment_volume_element,
                                             fragment_geometry(3.0));
    const std::vector<double> remote_values{remote_moments.charge,
                                            remote_moments.dipole,
                                            remote_moments.quadrupole};
    const PresetArrayReduction moment_reduction(remote_values);
    const ModulePcc::Pcc2dMoments reduced
        = ModulePcc::reduced_pcc_2d_density_moments(local_density,
                                                     local_positions,
                                                     fragment_volume_element,
                                                     fragment_geometry(3.0),
                                                     moment_reduction);
    std::vector<double> full_density = local_density;
    std::vector<ModuleBase::Vector3<double>> full_positions = local_positions;
    full_density.insert(full_density.end(), remote_density.begin(), remote_density.end());
    full_positions.insert(full_positions.end(), remote_positions.begin(), remote_positions.end());
    const ModulePcc::Pcc2dMoments expected
        = ModulePcc::pcc_2d_density_moments(full_density,
                                             full_positions,
                                             fragment_volume_element,
                                             fragment_geometry(3.0));
    EXPECT_NEAR(reduced.charge, expected.charge, 1.0e-12);
    EXPECT_NEAR(reduced.dipole, expected.dipole, 1.0e-12);
    EXPECT_NEAR(reduced.quadrupole, expected.quadrupole, 1.0e-11);
}

TEST_F(SccsPcc2dCoulombTest, PotentialUsesOnePeriodicTransformPair)
{
    const double tpiba = ModuleBase::TWO_PI / lattice_scale_;
    const ModuleSccs::Pcc2dCoulombOperator coulomb(basis_, tpiba, positions_,
                                                 volume_element_, geometry_, reduction_);
    const std::vector<double> charge(basis_.nrxx, 1.0 / volume_);
    std::vector<double> potential;
    coulomb.apply_potential(charge, potential);
    EXPECT_EQ(coulomb.transform_counts().forward_calls, 1);
    EXPECT_EQ(coulomb.transform_counts().inverse_calls, 1);
    const std::vector<double> zero(basis_.nrxx, 0.0);
    coulomb.apply_potential(zero, potential);
    for (double value : potential)
    {
        EXPECT_DOUBLE_EQ(value, 0.0);
    }
}

TEST_F(SccsPcc2dCoulombTest, AddsCorrectionFromCurrentChargeOnEveryApplication)
{
    std::vector<double> charge(basis_.nrxx);
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        const double relative_y = ModulePcc::pcc_2d_relative_coordinate(positions_[ir], geometry_);
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
    std::vector<double> potential;
    std::vector<double> periodic_potential;
    coulomb.apply_potential(charge, potential);
    periodic.apply_potential(charge, periodic_potential);
    const ModulePcc::Pcc2dMoments moments
        = ModulePcc::reduced_pcc_2d_density_moments(charge,
                                                     positions_,
                                                     volume_element_,
                                                     geometry_,
                                                     reduction_);
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        const double relative_y = ModulePcc::pcc_2d_relative_coordinate(positions_[ir], geometry_);
        const double correction = ModulePcc::pcc_2d_potential(moments, relative_y, geometry_.parameters);
        EXPECT_NEAR(potential[ir] - periodic_potential[ir], correction, 2.0e-12);
    }

    std::vector<double> opposite(charge.size());
    for (std::size_t index = 0; index < charge.size(); ++index)
    {
        opposite[index] = -charge[index];
    }
    std::vector<double> opposite_potential;
    coulomb.apply_potential(opposite, opposite_potential);
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        EXPECT_NEAR(opposite_potential[ir], -potential[ir], 2.0e-12);
    }
}

// The production sqrt-CG solves eps^-1/2 C_PCC eps^-1/2 exactly in one step for
// a uniform dielectric and keeps the PCC2D gauge instead of the periodic zero mean.
TEST_F(SccsPcc2dCoulombTest, SqrtCgKeepsChargedUniformDielectricPccGauge)
{
    const double tpiba = ModuleBase::TWO_PI / lattice_scale_;
    const ModuleSccs::Pcc2dCoulombOperator coulomb(basis_,
                                                   tpiba,
                                                   positions_,
                                                   volume_element_,
                                                   geometry_,
                                                   reduction_);
    // A zero electron density puts the whole cell in bulk solvent.
    const std::vector<double> cavity_density(basis_.nrxx, 0.0);
    const std::vector<double> solute_charge(basis_.nrxx, 1.0 / volume_);
    ModuleSccs::CavityParameters cavity;
    cavity.density_min = 1.0e-4;
    cavity.density_max = 5.0e-3;
    cavity.epsilon_bulk = 5.0;
    ModuleSccs::PolarizationSolverParameters solver;
    solver.max_iterations = 10;
    solver.tolerance_rms = 1.0e-14;
    solver.tolerance_max = 1.0e-14;
    const std::vector<double> cold_start;
    const ModuleSccs::SccsResponse result
        = ModuleSccs::solve_sccs_response(cavity_density,
                                          solute_charge,
                                          cavity,
                                          solver,
                                          cold_start,
                                          basis_,
                                          tpiba,
                                          coulomb,
                                          reduction_);
    EXPECT_EQ(result.polarization.iterations, 1);
    EXPECT_NEAR(result.far_field_polarization_charge,
                -(1.0 - 1.0 / cavity.epsilon_bulk),
                1.0e-12);

    std::vector<double> screened_charge(solute_charge.size());
    for (std::size_t index = 0; index < solute_charge.size(); ++index)
    {
        screened_charge[index] = solute_charge[index] / cavity.epsilon_bulk;
    }
    std::vector<double> expected;
    coulomb.apply_potential(screened_charge, expected);
    double local_mean = 0.0;
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        local_mean += expected[ir];
        EXPECT_NEAR(result.polarization.field.potential[ir], expected[ir], 2.0e-12);
        EXPECT_DOUBLE_EQ(result.restart_potential[ir], result.polarization.field.potential[ir]);
        EXPECT_NEAR(result.polarization.polarization_charge[ir],
                    -(1.0 - 1.0 / cavity.epsilon_bulk) * solute_charge[ir],
                    1.0e-14);
    }
    // The ENVIRON monopole constant keeps a nonzero cell average, so a
    // periodic zero-mean shift would fail the pointwise comparison above.
    reduction_.reduce_sum(local_mean);
    const double mean = local_mean / static_cast<double>(basis_.nxyz);
    EXPECT_GT(std::abs(mean), 1.0e-3);
}

// Open 1D reference for the production sqrt-CG: a neutral dipolar layer in a
// cavity uniform in x and z. D = 4 pi int q dy vanishes outside the layer, so
// dv/dy = -4 pi cumulative(y) / eps(y). Central differences of the solved
// potential check the solution itself. With the chain-rule factsqrt the
// sampled f v source (f drops steeply to zero at density_min, v carries the
// dipole plateau) left a 1.4% field error here that did not converge with the
// grid; the FFT derivatives of the switching function used with PCC converge.
// The FFT gradient of the corrected potential (Environ de_dboundary, used by
// the default continuum cavity potential) rings from the dipole step at the
// open boundary; this cavity reaches within 1.2 bohr of it, so the ringing
// also enters the transition region. It is printed, not asserted: the switching
// lowpass mode does not use this gradient.
TEST(SccsPcc2dSqrtCg, LayeredCavityMatchesOpenOneDimensionalField)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(0.4, 0.0, 0.0,
                                      0.0, 1.2, 0.0,
                                      0.0, 0.0, 0.4);
    const double scale = 10.0;
    const double tpiba = ModuleBase::TWO_PI / scale;
    basis.initgrids(scale, lattice, 320.0);
    basis.initparameters(false, 320.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();
    const ModulePcc::Pcc2dGeometry geometry = ModulePcc::pcc_2d_geometry(lattice, scale, 1, 1.0e-10);
    const double volume = geometry.parameters.periodic_area * geometry.parameters.cell_length;
    const double volume_element = volume / static_cast<double>(basis.nxyz);
    const std::vector<ModuleBase::Vector3<double>> positions
        = ModuleSurchem::pw_grid_positions(basis, lattice, scale);
    const ModuleSurchem::SerialChargeReduction charge_reduction;
    const ModuleSccs::Pcc2dCoulombOperator coulomb(basis,
                                                   tpiba,
                                                   positions,
                                                   volume_element,
                                                   geometry,
                                                   charge_reduction);
    ModuleSccs::CavityParameters cavity;
    cavity.density_min = 1.0e-4;
    cavity.density_max = 5.0e-2;
    cavity.epsilon_bulk = 5.0;
    const double density_peak = 5.0e-2;
    const double density_width = 2.0;
    const double source_width = 1.3;
    const double amplitude = 3.0e-3;
    const double half_length = 0.5 * geometry.parameters.cell_length;
    const double boundary_exponential
        = std::exp(-half_length * half_length / (source_width * source_width));
    std::vector<double> cavity_density(basis.nrxx);
    std::vector<double> solute_charge(basis.nrxx);
    std::vector<double> reference_gradient(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double y = ModulePcc::pcc_2d_relative_coordinate(positions[ir], geometry);
        cavity_density[ir] = density_peak * std::exp(-y * y / (density_width * density_width));
        const double source_exponential = std::exp(-y * y / (source_width * source_width));
        solute_charge[ir] = amplitude * y * source_exponential;
        const double cumulative = -0.5 * amplitude * source_width * source_width
                                  * (source_exponential - boundary_exponential);
        const ModuleSccs::CavityPoint point = ModuleSccs::evaluate_cavity(cavity_density[ir], cavity);
        reference_gradient[ir] = -ModuleBase::FOUR_PI * cumulative / point.epsilon;
    }
    ModuleSccs::PolarizationSolverParameters solver;
    solver.max_iterations = 200;
    solver.tolerance_rms = 1.0e-13;
    solver.tolerance_max = 1.0e-12;
    const std::vector<double> cold_start;
    const ModuleSccs::SccsResponse result
        = ModuleSccs::solve_sccs_response(cavity_density,
                                          solute_charge,
                                          cavity,
                                          solver,
                                          cold_start,
                                          basis,
                                          tpiba,
                                          coulomb,
                                          charge_reduction);
    EXPECT_GT(result.polarization.iterations, 1);

    // The potential depends on y only; grid index = (ix ny + iy) nplane + iz.
    const double step_y = geometry.parameters.cell_length / static_cast<double>(basis.ny);
    const double interior = half_length - 2.5 * step_y;
    double maximum_difference_error = 0.0;
    double maximum_gradient = 0.0;
    double maximum_transverse_gradient = 0.0;
    double maximum_transition_fft_error = 0.0;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double y = ModulePcc::pcc_2d_relative_coordinate(positions[ir], geometry);
        maximum_gradient = std::max(maximum_gradient, std::abs(reference_gradient[ir]));
        const double transverse = std::max(std::abs(result.polarization.field.gradient[ir].x),
                                           std::abs(result.polarization.field.gradient[ir].z));
        maximum_transverse_gradient = std::max(maximum_transverse_gradient, transverse);
        if (cavity_density[ir] > cavity.density_min && cavity_density[ir] < cavity.density_max)
        {
            const double fft_error
                = std::abs(result.polarization.field.gradient[ir].y - reference_gradient[ir]);
            maximum_transition_fft_error = std::max(maximum_transition_fft_error, fft_error);
        }
        const int iy = (ir / basis.nplane) % basis.ny;
        if (std::abs(y) >= interior || iy == 0 || iy == basis.ny - 1)
        {
            continue;
        }
        const double upper = result.polarization.field.potential[ir + basis.nplane];
        const double lower = result.polarization.field.potential[ir - basis.nplane];
        const double difference = (upper - lower) / (2.0 * step_y);
        const double error = std::abs(difference - reference_gradient[ir]);
        maximum_difference_error = std::max(maximum_difference_error, error);
    }
    std::cout << "PCC2D_SQRT_CG_LAYERED max_difference_error " << maximum_difference_error
              << " max_gradient " << maximum_gradient << " transition_fft_error "
              << maximum_transition_fft_error
              << " far_field_polarization_charge " << result.far_field_polarization_charge
              << " iterations " << result.polarization.iterations << std::endl;
    // Measured at this grid: central-difference field error 1.7e-4 (0.55% of
    // the peak field, mostly the O(step^2) difference error); the FFT gradient
    // misses by 0.023 in the transition region.
    const double difference_tolerance = 1.0e-2 * maximum_gradient;
    EXPECT_LT(maximum_difference_error, difference_tolerance);
    EXPECT_LT(maximum_transverse_gradient, 1.0e-12);
    // The open boundary only weakly pins the constant potential mode; rounding
    // in the source leaves a gauge offset whose far-field trace is 3.4e-5 here.
    EXPECT_NEAR(result.far_field_polarization_charge, 0.0, 1.0e-4);
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
