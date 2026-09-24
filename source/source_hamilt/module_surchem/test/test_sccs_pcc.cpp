#include "../sccs_pcc.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "gtest/gtest.h"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace
{

TEST(SccsPcc, AccumulatesPointChargeMomentsAboutTheRequestedOrigin)
{
    std::vector<ModuleSccs::PointCharge> charges(2);
    charges[0].charge = 2.0;
    charges[0].position = ModuleBase::Vector3<double>(2.0, 0.0, 0.0);
    charges[1].charge = -1.0;
    charges[1].position = ModuleBase::Vector3<double>(0.0, 2.0, 0.0);

    const ModuleSccs::MultipoleMoments moments
        = ModuleSccs::point_charge_moments(charges, ModuleBase::Vector3<double>(1.0, 1.0, 0.0));
    EXPECT_DOUBLE_EQ(moments.charge, 1.0);
    EXPECT_DOUBLE_EQ(moments.dipole.x, 3.0);
    EXPECT_DOUBLE_EQ(moments.dipole.y, -3.0);
    EXPECT_DOUBLE_EQ(moments.dipole.z, 0.0);
    EXPECT_DOUBLE_EQ(moments.quadrupole_trace, 2.0);
}

TEST(SccsPcc, BuildsGeometryForRotatedCubicCell)
{
    const ModuleBase::Matrix3 lattice(0.0, 1.0, 0.0,
                                      -1.0, 0.0, 0.0,
                                      0.0, 0.0, 1.0);
    const ModuleSccs::PccGeometry geometry
        = ModuleSccs::pcc_geometry(lattice, 10.0, 1.0e-10);
    EXPECT_DOUBLE_EQ(geometry.parameters.cube_length, 10.0);
    EXPECT_DOUBLE_EQ(geometry.origin.x, -5.0);
    EXPECT_DOUBLE_EQ(geometry.origin.y, 5.0);
    EXPECT_DOUBLE_EQ(geometry.origin.z, 5.0);

    const ModuleBase::Vector3<double> position(
        geometry.origin.x + 6.0 * geometry.axis_a.x,
        geometry.origin.y + 6.0 * geometry.axis_a.y,
        geometry.origin.z + 6.0 * geometry.axis_a.z);
    const ModuleBase::Vector3<double> relative
        = ModuleSccs::pcc_relative_position(position, geometry);
    EXPECT_NEAR(relative.x, -4.0 * geometry.axis_a.x, 1.0e-14);
    EXPECT_NEAR(relative.y, -4.0 * geometry.axis_a.y, 1.0e-14);
    EXPECT_NEAR(relative.z, -4.0 * geometry.axis_a.z, 1.0e-14);
}

TEST(SccsPcc, SystemCenterUnwrapsAtomsAcrossPeriodicBoundaries)
{
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    const ModuleSccs::PccGeometry geometry
        = ModuleSccs::pcc_geometry(lattice, 10.0, 1.0e-10);
    const std::vector<ModuleBase::Vector3<double>> positions{
        ModuleBase::Vector3<double>(9.8, 9.7, 9.6),
        ModuleBase::Vector3<double>(0.2, 0.3, 0.4)};
    const std::vector<double> masses{1.0, 1.0};
    const ModuleBase::Vector3<double> center
        = ModuleSccs::pcc_system_center(positions, masses, geometry);
    EXPECT_NEAR(center.x, 0.0, 1.0e-14);
    EXPECT_NEAR(center.y, 0.0, 1.0e-14);
    EXPECT_NEAR(center.z, 0.0, 1.0e-14);
}

TEST(SccsPcc, SystemCenterMakesWrappedRigidTranslationInvariant)
{
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    ModuleSccs::PccGeometry first_geometry
        = ModuleSccs::pcc_geometry(lattice, 10.0, 1.0e-10);
    std::vector<ModuleSccs::PointCharge> first(2);
    first[0].charge = 1.2;
    first[0].position = ModuleBase::Vector3<double>(9.8, 9.6, 0.1);
    first[1].charge = -0.3;
    first[1].position = ModuleBase::Vector3<double>(0.3, 0.4, 9.7);
    const std::vector<double> masses{2.0, 1.0};
    const std::vector<ModuleBase::Vector3<double>> first_positions{
        first[0].position, first[1].position};
    first_geometry.origin
        = ModuleSccs::pcc_system_center(first_positions, masses, first_geometry);

    std::vector<ModuleSccs::PointCharge> shifted = first;
    for (std::size_t index = 0; index < shifted.size(); ++index)
    {
        shifted[index].position.x = std::fmod(shifted[index].position.x + 1.4, 10.0);
        shifted[index].position.y = std::fmod(shifted[index].position.y + 1.4, 10.0);
        shifted[index].position.z = std::fmod(shifted[index].position.z + 1.4, 10.0);
    }
    ModuleSccs::PccGeometry shifted_geometry
        = ModuleSccs::pcc_geometry(lattice, 10.0, 1.0e-10);
    const std::vector<ModuleBase::Vector3<double>> shifted_positions{
        shifted[0].position, shifted[1].position};
    shifted_geometry.origin
        = ModuleSccs::pcc_system_center(shifted_positions, masses, shifted_geometry);

    const ModuleSccs::MultipoleMoments first_moments
        = ModuleSccs::point_charge_moments(first, first_geometry);
    const ModuleSccs::MultipoleMoments shifted_moments
        = ModuleSccs::point_charge_moments(shifted, shifted_geometry);
    EXPECT_NEAR(first_moments.charge, shifted_moments.charge, 1.0e-14);
    EXPECT_NEAR(first_moments.dipole.x, shifted_moments.dipole.x, 1.0e-14);
    EXPECT_NEAR(first_moments.dipole.y, shifted_moments.dipole.y, 1.0e-14);
    EXPECT_NEAR(first_moments.dipole.z, shifted_moments.dipole.z, 1.0e-14);
    EXPECT_NEAR(first_moments.quadrupole_trace,
                shifted_moments.quadrupole_trace,
                1.0e-14);
    EXPECT_NEAR(ModuleSccs::pcc_self_energy(first_moments,
                                            first_geometry.parameters),
                ModuleSccs::pcc_self_energy(shifted_moments,
                                            shifted_geometry.parameters),
                1.0e-14);
}

TEST(SccsPcc, PotentialGradientMatchesCentralDifferences)
{
    ModuleSccs::MultipoleMoments moments;
    moments.charge = -1.2;
    moments.dipole = ModuleBase::Vector3<double>(0.3, -0.7, 0.2);
    moments.quadrupole_trace = 4.1;
    ModuleSccs::PccParameters parameters;
    parameters.cube_length = 18.0;
    const ModuleBase::Vector3<double> position(0.9, -1.3, 0.4);
    const ModuleBase::Vector3<double> gradient
        = ModuleSccs::pcc_potential_gradient(moments, position, parameters);
    const double step = 1.0e-5;

    for (int direction = 0; direction < 3; ++direction)
    {
        ModuleBase::Vector3<double> plus = position;
        ModuleBase::Vector3<double> minus = position;
        if (direction == 0)
        {
            plus.x += step;
            minus.x -= step;
        }
        else if (direction == 1)
        {
            plus.y += step;
            minus.y -= step;
        }
        else
        {
            plus.z += step;
            minus.z -= step;
        }
        const double finite = (ModuleSccs::pcc_potential(moments, plus, parameters)
                               - ModuleSccs::pcc_potential(moments, minus, parameters))
                              / (2.0 * step);
        const double analytic = direction == 0 ? gradient.x : (direction == 1 ? gradient.y : gradient.z);
        EXPECT_NEAR(analytic, finite, 1.0e-11);
    }
}

TEST(SccsPcc, IntegratesDensityMomentsWithTheVolumeElement)
{
    const std::vector<double> density{0.5, -0.25};
    const std::vector<ModuleBase::Vector3<double>> positions{
        ModuleBase::Vector3<double>(1.0, 0.0, 0.0),
        ModuleBase::Vector3<double>(0.0, 2.0, 0.0)};
    const ModuleSccs::MultipoleMoments moments
        = ModuleSccs::density_moments(density,
                                      positions,
                                      2.0,
                                      ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    EXPECT_DOUBLE_EQ(moments.charge, 0.5);
    EXPECT_DOUBLE_EQ(moments.dipole.x, 1.0);
    EXPECT_DOUBLE_EQ(moments.dipole.y, -1.0);
    EXPECT_DOUBLE_EQ(moments.quadrupole_trace, -1.0);
}

TEST(SccsPcc, BilinearKernelIsSymmetric)
{
    ModuleSccs::MultipoleMoments left;
    left.charge = 1.0;
    left.dipole = ModuleBase::Vector3<double>(0.1, 0.2, 0.3);
    left.quadrupole_trace = 2.0;
    ModuleSccs::MultipoleMoments right;
    right.charge = -2.0;
    right.dipole = ModuleBase::Vector3<double>(-0.4, 0.7, 0.2);
    right.quadrupole_trace = -1.0;
    ModuleSccs::PccParameters parameters;
    parameters.cube_length = 20.0;

    EXPECT_DOUBLE_EQ(ModuleSccs::pcc_bilinear_energy(left, right, parameters),
                     ModuleSccs::pcc_bilinear_energy(right, left, parameters));
}

TEST(SccsPcc, SelfEnergyHasExpectedMakovPayneForm)
{
    ModuleSccs::MultipoleMoments moments;
    moments.charge = 1.5;
    moments.dipole = ModuleBase::Vector3<double>(0.2, -0.4, 0.7);
    moments.quadrupole_trace = 3.2;
    ModuleSccs::PccParameters parameters;
    parameters.cube_length = 14.0;
    const double volume = std::pow(parameters.cube_length, 3);
    const double dipole_squared = 0.2 * 0.2 + 0.4 * 0.4 + 0.7 * 0.7;
    const double expected = parameters.madelung * moments.charge * moments.charge
                                / (2.0 * parameters.cube_length)
                            - 2.0 * ModuleBase::PI
                                  * (moments.charge * moments.quadrupole_trace - dipole_squared)
                                  / (3.0 * volume);
    EXPECT_NEAR(ModuleSccs::pcc_self_energy(moments, parameters), expected, 1.0e-15);
}

TEST(SccsPcc, SelfEnergyIsIndependentOfTheMultipoleOrigin)
{
    std::vector<ModuleSccs::PointCharge> charges(2);
    charges[0].charge = 1.7;
    charges[0].position = ModuleBase::Vector3<double>(1.2, -0.5, 0.8);
    charges[1].charge = -0.4;
    charges[1].position = ModuleBase::Vector3<double>(-0.3, 0.6, 1.1);
    const ModuleSccs::MultipoleMoments first
        = ModuleSccs::point_charge_moments(charges, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    const ModuleSccs::MultipoleMoments second
        = ModuleSccs::point_charge_moments(charges, ModuleBase::Vector3<double>(0.7, -1.0, 0.2));
    ModuleSccs::PccParameters parameters;
    parameters.cube_length = 24.0;

    EXPECT_NEAR(ModuleSccs::pcc_self_energy(first, parameters),
                ModuleSccs::pcc_self_energy(second, parameters),
                1.0e-15);
}

TEST(SccsPcc, PointIonVacuumEnergyDerivativeMatchesElectronPotential)
{
    const ModuleBase::Vector3<double> origin(0.0, 0.0, 0.0);
    const ModuleBase::Vector3<double> electron_position(1.1, -0.8, 0.4);
    const double volume_element = 0.3;
    const double electron_density = 0.7;
    const double step = 1.0e-6;
    ModuleSccs::PccParameters parameters;
    parameters.cube_length = 19.0;

    const auto energy = [&](const double density) {
        std::vector<ModuleSccs::PointCharge> charges(2);
        charges[0].charge = 1.0;
        charges[0].position = ModuleBase::Vector3<double>(-0.2, 0.5, -0.1);
        charges[1].charge = -density * volume_element;
        charges[1].position = electron_position;
        return ModuleSccs::pcc_self_energy(ModuleSccs::point_charge_moments(charges, origin),
                                            parameters);
    };

    std::vector<ModuleSccs::PointCharge> center_charges(2);
    center_charges[0].charge = 1.0;
    center_charges[0].position = ModuleBase::Vector3<double>(-0.2, 0.5, -0.1);
    center_charges[1].charge = -electron_density * volume_element;
    center_charges[1].position = electron_position;
    const ModuleSccs::MultipoleMoments center
        = ModuleSccs::point_charge_moments(center_charges, origin);
    const double analytic = -volume_element
                            * ModuleSccs::pcc_potential(center,
                                                       electron_position,
                                                       parameters);
    const double finite = (energy(electron_density + step)
                           - energy(electron_density - step))
                          / (2.0 * step);
    EXPECT_NEAR(analytic, finite, 1.0e-11);
}

TEST(SccsPcc, PointChargeForceMatchesSelfEnergyFiniteDifference)
{
    std::vector<ModuleSccs::PointCharge> charges(2);
    charges[0].charge = 1.2;
    charges[0].position = ModuleBase::Vector3<double>(0.7, -0.4, 0.2);
    charges[1].charge = -0.3;
    charges[1].position = ModuleBase::Vector3<double>(-0.6, 0.5, 0.9);
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    ModuleSccs::PccGeometry geometry
        = ModuleSccs::pcc_geometry(lattice, 17.0, 1.0e-10);
    geometry.origin = ModuleBase::Vector3<double>();
    const ModuleSccs::MultipoleMoments moments
        = ModuleSccs::point_charge_moments(charges, geometry);
    const ModuleBase::Vector3<double> analytic
        = ModuleSccs::pcc_point_charge_force(moments, charges[0], geometry);
    const double step = 1.0e-5;

    for (int direction = 0; direction < 3; ++direction)
    {
        std::vector<ModuleSccs::PointCharge> plus = charges;
        std::vector<ModuleSccs::PointCharge> minus = charges;
        plus[0].position[direction] += step;
        minus[0].position[direction] -= step;
        const double energy_plus
            = ModuleSccs::pcc_self_energy(
                ModuleSccs::point_charge_moments(plus, geometry),
                geometry.parameters);
        const double energy_minus
            = ModuleSccs::pcc_self_energy(
                ModuleSccs::point_charge_moments(minus, geometry),
                geometry.parameters);
        const double finite_force = -(energy_plus - energy_minus) / (2.0 * step);
        EXPECT_NEAR(analytic[direction], finite_force, 1.0e-11);
    }
}

TEST(SccsPcc, RejectsInvalidCubeLength)
{
    ModuleSccs::PccParameters parameters;
    EXPECT_THROW(ModuleSccs::pcc_self_energy(ModuleSccs::MultipoleMoments(), parameters),
                 std::invalid_argument);
}

} // namespace
