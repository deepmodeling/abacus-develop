#include "../pcc/pcc_2d.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "gtest/gtest.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{

ModulePcc::Pcc2dParameters parameters()
{
    ModulePcc::Pcc2dParameters value;
    value.periodic_area = 31.0;
    value.cell_length = 17.0;
    return value;
}

ModulePcc::Pcc2dGeometry geometry(const double origin_y)
{
    ModulePcc::Pcc2dGeometry value;
    value.parameters = parameters();
    value.origin = origin_y;
    return value;
}

TEST(SccsPcc2d, BuildsGeometryForCubicCell)
{
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    const ModulePcc::Pcc2dGeometry geometry
        = ModulePcc::pcc_2d_geometry(lattice, 10.0, 1, 1.0e-10);

    EXPECT_DOUBLE_EQ(geometry.parameters.periodic_area, 100.0);
    EXPECT_DOUBLE_EQ(geometry.parameters.cell_length, 10.0);
    EXPECT_DOUBLE_EQ(geometry.origin, 5.0);
}

TEST(SccsPcc2d, BuildsGeometryForNonOrthogonalPeriodicPlane)
{
    const ModuleBase::Matrix3 lattice(2.0, 0.0, 0.0,
                                      0.0, 5.0, 0.0,
                                      1.0, 0.0, 3.0);
    const ModulePcc::Pcc2dGeometry geometry
        = ModulePcc::pcc_2d_geometry(lattice, 4.0, 1, 1.0e-10);

    EXPECT_DOUBLE_EQ(geometry.parameters.periodic_area, 96.0);
    EXPECT_DOUBLE_EQ(geometry.parameters.cell_length, 20.0);
    EXPECT_DOUBLE_EQ(geometry.origin, 10.0);
}

TEST(SccsPcc2d, GeometryAcceptsNumericalAlignmentNoise)
{
    const ModuleBase::Matrix3 lattice(1.0, 1.0e-12, 0.0,
                                      -1.0e-12, 2.0, 1.0e-12,
                                      0.3, -1.0e-12, 1.0);
    EXPECT_NO_THROW(ModulePcc::pcc_2d_geometry(lattice, 7.0, 1, 1.0e-10));
}

// The normal follows the open lattice vector, so a reversed vector only
// flips the sign convention of the coordinate along it.
TEST(SccsPcc2d, NormalFollowsReversedOpenVector)
{
    const ModuleBase::Matrix3 reversed_normal(1.0, 0.0, 0.0,
                                              0.0, -2.0, 0.0,
                                              0.0, 0.0, 1.0);
    const ModulePcc::Pcc2dGeometry geometry
        = ModulePcc::pcc_2d_geometry(reversed_normal, 1.0, 1, 1.0e-10);
    EXPECT_DOUBLE_EQ(geometry.normal.y, -1.0);
    EXPECT_DOUBLE_EQ(geometry.parameters.cell_length, 2.0);
    EXPECT_DOUBLE_EQ(geometry.parameters.periodic_area, 1.0);
}

TEST(SccsPcc2d, BuildsGeometryForEveryOpenAxis)
{
    const ModuleBase::Matrix3 lattice(2.0, 0.0, 0.0,
                                      0.0, 5.0, 0.0,
                                      0.0, 0.0, 3.0);
    const ModulePcc::Pcc2dGeometry along_x = ModulePcc::pcc_2d_geometry(lattice, 4.0, 0, 1.0e-10);
    EXPECT_DOUBLE_EQ(along_x.normal.x, 1.0);
    EXPECT_DOUBLE_EQ(along_x.parameters.cell_length, 8.0);
    EXPECT_DOUBLE_EQ(along_x.parameters.periodic_area, 240.0);
    EXPECT_DOUBLE_EQ(along_x.origin, 4.0);
    const ModulePcc::Pcc2dGeometry along_z = ModulePcc::pcc_2d_geometry(lattice, 4.0, 2, 1.0e-10);
    EXPECT_DOUBLE_EQ(along_z.normal.z, 1.0);
    EXPECT_DOUBLE_EQ(along_z.parameters.cell_length, 12.0);
    EXPECT_DOUBLE_EQ(along_z.parameters.periodic_area, 160.0);
    EXPECT_DOUBLE_EQ(along_z.origin, 6.0);
    // A hexagonal periodic plane with the open vector along z.
    const ModuleBase::Matrix3 hexagonal(1.0, 0.0, 0.0,
                                        -0.5, 0.5 * std::sqrt(3.0), 0.0,
                                        0.0, 0.0, 4.0);
    const ModulePcc::Pcc2dGeometry slab = ModulePcc::pcc_2d_geometry(hexagonal, 2.0, 2, 1.0e-10);
    EXPECT_DOUBLE_EQ(slab.normal.z, 1.0);
    EXPECT_NEAR(slab.parameters.periodic_area, 2.0 * std::sqrt(3.0), 1.0e-14);
    EXPECT_DOUBLE_EQ(slab.parameters.cell_length, 8.0);
    const ModuleBase::Vector3<double> position(0.3, -0.2, 7.0);
    EXPECT_DOUBLE_EQ(ModulePcc::pcc_2d_coordinate(position, slab), 7.0);
    EXPECT_THROW(ModulePcc::pcc_2d_geometry(lattice, 4.0, 3, 1.0e-10), std::invalid_argument);
}

TEST(SccsPcc2d, GeometryRejectsUnsupportedSlabOrientations)
{
    const ModuleBase::Matrix3 tilted_normal(1.0, 0.0, 0.0,
                                             0.2, 2.0, 0.0,
                                             0.0, 0.0, 1.0);
    const ModuleBase::Matrix3 tilted_plane(1.0, 0.1, 0.0,
                                           0.0, 2.0, 0.0,
                                           0.0, 0.0, 1.0);
    const ModuleBase::Matrix3 reversed_normal(1.0, 0.0, 0.0,
                                              0.0, -2.0, 0.0,
                                              0.0, 0.0, 1.0);
    const ModuleBase::Matrix3 degenerate_plane(1.0, 0.0, 0.0,
                                               0.0, 2.0, 0.0,
                                               2.0, 0.0, 0.0);

    EXPECT_THROW(ModulePcc::pcc_2d_geometry(tilted_normal, 1.0, 1, 1.0e-10),
                 std::invalid_argument);
    EXPECT_THROW(ModulePcc::pcc_2d_geometry(tilted_plane, 1.0, 1, 1.0e-10),
                 std::invalid_argument);
    EXPECT_THROW(ModulePcc::pcc_2d_geometry(degenerate_plane, 1.0, 1, 1.0e-10),
                 std::invalid_argument);
}

TEST(SccsPcc2d, SystemCenterUnwrapsAcrossYBoundary)
{
    const std::vector<double> positions_y{9.8, 0.2};
    const std::vector<double> weights{1.0, 3.0};
    EXPECT_NEAR(ModulePcc::pcc_2d_system_center(positions_y, weights, 10.0),
                0.1,
                1.0e-14);

    ModulePcc::Pcc2dGeometry value = geometry(0.1);
    value.parameters.cell_length = 10.0;
    EXPECT_NEAR(ModulePcc::pcc_2d_relative_coordinate(ModuleBase::Vector3<double>(0.0, 9.8, 0.0), value), -0.3, 1.0e-14);
    EXPECT_NEAR(ModulePcc::pcc_2d_relative_coordinate(ModuleBase::Vector3<double>(0.0, 0.2, 0.0), value), 0.1, 1.0e-14);
}

TEST(SccsPcc2d, MomentsRemainInvariantWhenChargesCrossYBoundary)
{
    std::vector<ModulePcc::PointCharge> original(2);
    original[0].charge = 1.2;
    original[0].position.y = 9.8;
    original[1].charge = -0.4;
    original[1].position.y = 0.2;
    const std::vector<double> weights{1.0, 3.0};
    const std::vector<double> original_y{9.8, 0.2};
    ModulePcc::Pcc2dGeometry first = geometry(
        ModulePcc::pcc_2d_system_center(original_y, weights, 10.0));
    first.parameters.cell_length = 10.0;

    std::vector<ModulePcc::PointCharge> shifted = original;
    shifted[0].position.y = 3.5;
    shifted[1].position.y = 3.9;
    const std::vector<double> shifted_y{3.5, 3.9};
    ModulePcc::Pcc2dGeometry second = geometry(
        ModulePcc::pcc_2d_system_center(shifted_y, weights, 10.0));
    second.parameters.cell_length = 10.0;

    const ModulePcc::Pcc2dMoments original_moments
        = ModulePcc::pcc_2d_point_charge_moments(original, first);
    const ModulePcc::Pcc2dMoments shifted_moments
        = ModulePcc::pcc_2d_point_charge_moments(shifted, second);
    EXPECT_NEAR(shifted_moments.charge, original_moments.charge, 1.0e-15);
    EXPECT_NEAR(shifted_moments.dipole, original_moments.dipole, 1.0e-14);
    EXPECT_NEAR(shifted_moments.quadrupole,
                original_moments.quadrupole,
                1.0e-14);
}

TEST(SccsPcc2d, AccumulatesOnlyYMoments)
{
    std::vector<ModulePcc::PointCharge> charges(2);
    charges[0].charge = 2.0;
    charges[0].position = ModuleBase::Vector3<double>(4.0, 3.0, -7.0);
    charges[1].charge = -1.0;
    charges[1].position = ModuleBase::Vector3<double>(-8.0, -1.0, 9.0);

    const ModulePcc::Pcc2dMoments moments
        = ModulePcc::pcc_2d_point_charge_moments(charges, geometry(1.0));
    EXPECT_DOUBLE_EQ(moments.charge, 1.0);
    EXPECT_DOUBLE_EQ(moments.dipole, 6.0);
    EXPECT_DOUBLE_EQ(moments.quadrupole, 4.0);
}

TEST(SccsPcc2d, SumMomentsEqualTheMomentsOfTheCombinedCharges)
{
    std::vector<ModulePcc::PointCharge> ions(2);
    ions[0].charge = 2.0;
    ions[0].position = ModuleBase::Vector3<double>(4.0, 3.0, -7.0);
    ions[1].charge = 1.0;
    ions[1].position = ModuleBase::Vector3<double>(0.5, -2.0, 1.0);
    std::vector<ModulePcc::PointCharge> electrons(1);
    electrons[0].charge = -2.5;
    electrons[0].position = ModuleBase::Vector3<double>(-8.0, 0.5, 9.0);
    const ModulePcc::Pcc2dGeometry slab = geometry(1.0);

    const ModulePcc::Pcc2dMoments ionic = ModulePcc::pcc_2d_point_charge_moments(ions, slab);
    const ModulePcc::Pcc2dMoments electronic = ModulePcc::pcc_2d_point_charge_moments(electrons, slab);
    const ModulePcc::Pcc2dMoments sum = ModulePcc::pcc_2d_sum_moments(ionic, electronic);
    std::vector<ModulePcc::PointCharge> combined = ions;
    combined.push_back(electrons[0]);
    const ModulePcc::Pcc2dMoments expected = ModulePcc::pcc_2d_point_charge_moments(combined, slab);
    EXPECT_DOUBLE_EQ(sum.charge, expected.charge);
    EXPECT_DOUBLE_EQ(sum.dipole, expected.dipole);
    EXPECT_DOUBLE_EQ(sum.quadrupole, expected.quadrupole);
}

TEST(SccsPcc2d, DensityMomentsIncludeTheVolumeElement)
{
    const std::vector<double> density{0.5, -0.25};
    const std::vector<ModuleBase::Vector3<double>> positions{
        ModuleBase::Vector3<double>(7.0, 1.0, -4.0),
        ModuleBase::Vector3<double>(-2.0, 3.0, 6.0)};
    const ModulePcc::Pcc2dMoments moments
        = ModulePcc::pcc_2d_density_moments(density, positions, 2.0, geometry(0.0));
    EXPECT_DOUBLE_EQ(moments.charge, 0.5);
    EXPECT_DOUBLE_EQ(moments.dipole, -0.5);
    EXPECT_DOUBLE_EQ(moments.quadrupole, -3.5);
}

// Planar open kernel -2*pi*|u|/A plus a constant gauge shift, minus the
// zero-mean periodic planar kernel, summed over sources closer than L_y/2.
double direct_planar_correction(const std::vector<double>& locations,
                                const std::vector<double>& charges,
                                const double y,
                                const double gauge_shift,
                                const ModulePcc::Pcc2dParameters& value)
{
    double direct = 0.0;
    for (std::size_t source = 0; source < charges.size(); ++source)
    {
        const double separation = y - locations[source];
        const double open_kernel
            = -2.0 * ModuleBase::PI * std::abs(separation) / value.periodic_area
              + gauge_shift;
        const double periodic_kernel
            = 2.0 * ModuleBase::PI / value.periodic_area
              * (separation * separation / value.cell_length
                 - std::abs(separation) + value.cell_length / 6.0);
        direct += charges[source] * (open_kernel - periodic_kernel);
    }
    return direct;
}

void expect_planar_correction(const ModulePcc::Pcc2dParameters& value,
                              const double gauge_shift)
{
    const std::vector<double> locations{-3.4, -1.1, 0.5, 2.8};
    const std::vector<double> charges{0.7, -0.2, 0.3, -0.1};
    std::vector<ModulePcc::PointCharge> points(charges.size());
    for (std::size_t index = 0; index < charges.size(); ++index)
    {
        points[index].charge = charges[index];
        points[index].position.y = locations[index];
    }
    ModulePcc::Pcc2dGeometry origin_geometry;
    origin_geometry.parameters = value;
    origin_geometry.origin = 0.0;
    const ModulePcc::Pcc2dMoments moments
        = ModulePcc::pcc_2d_point_charge_moments(points, origin_geometry);
    const std::vector<double> evaluation_points{-5.2, -0.8, 1.4, 5.1};
    for (std::size_t evaluation = 0; evaluation < evaluation_points.size(); ++evaluation)
    {
        const double y = evaluation_points[evaluation];
        for (std::size_t source = 0; source < locations.size(); ++source)
        {
            ASSERT_LE(std::abs(y - locations[source]), 0.5 * value.cell_length);
        }
        const double direct
            = direct_planar_correction(locations, charges, y, gauge_shift, value);
        EXPECT_NEAR(ModulePcc::pcc_2d_potential(moments, y, value), direct, 1.0e-14);
    }
}

TEST(SccsPcc2d, MomentPotentialMatchesIndependentPlanarGreenFunctions)
{
    const ModulePcc::Pcc2dParameters value = parameters();
    // ENVIRON's monopole constant -pi*q/(3 L_y) equals this shift of the open kernel.
    const double gauge_shift
        = ModuleBase::PI * value.cell_length / (3.0 * value.periodic_area)
          - ModuleBase::PI / (3.0 * value.cell_length);
    expect_planar_correction(value, gauge_shift);
}

TEST(SccsPcc2d, SquareCellCorrectionUsesUnshiftedOpenPlanarKernel)
{
    ModulePcc::Pcc2dParameters value;
    value.cell_length = 17.0;
    value.periodic_area = value.cell_length * value.cell_length;
    expect_planar_correction(value, 0.0);
}

TEST(SccsPcc2d, ChargedMonopoleUsesTheEnvironGauge)
{
    ModulePcc::Pcc2dMoments moments;
    moments.charge = 1.7;
    const ModulePcc::Pcc2dParameters value = parameters();
    const double expected_potential
        = -ModuleBase::PI * moments.charge / (3.0 * value.cell_length);
    const double expected_energy
        = -ModuleBase::PI * moments.charge * moments.charge
          / (6.0 * value.cell_length);

    EXPECT_NEAR(ModulePcc::pcc_2d_potential(moments, 0.0, value),
                expected_potential,
                1.0e-15);
    EXPECT_NEAR(ModulePcc::pcc_2d_self_energy(moments, value),
                expected_energy,
                1.0e-15);
}

TEST(SccsPcc2d, PositiveAndNegativeMonopolesHaveOddPotentialAndEqualEnergy)
{
    ModulePcc::Pcc2dMoments positive;
    positive.charge = 0.8;
    positive.dipole = -0.3;
    positive.quadrupole = 1.4;
    ModulePcc::Pcc2dMoments negative;
    negative.charge = -positive.charge;
    negative.dipole = -positive.dipole;
    negative.quadrupole = -positive.quadrupole;
    const ModulePcc::Pcc2dParameters value = parameters();
    const double y = 2.6;

    EXPECT_DOUBLE_EQ(ModulePcc::pcc_2d_potential(negative, y, value),
                     -ModulePcc::pcc_2d_potential(positive, y, value));
    EXPECT_DOUBLE_EQ(ModulePcc::pcc_2d_self_energy(negative, value),
                     ModulePcc::pcc_2d_self_energy(positive, value));
}

TEST(SccsPcc2d, GradientFollowsNormalAndMatchesCentralDifference)
{
    ModulePcc::Pcc2dMoments moments;
    moments.charge = -1.2;
    moments.dipole = 0.7;
    moments.quadrupole = 4.1;
    const ModulePcc::Pcc2dParameters value = parameters();
    ModulePcc::Pcc2dGeometry geometry;
    geometry.parameters = value;
    const double y = -1.3;
    const double step = 1.0e-5;
    const ModuleBase::Vector3<double> gradient
        = ModulePcc::pcc_2d_potential_gradient(moments, y, geometry);
    const double finite
        = (ModulePcc::pcc_2d_potential(moments, y + step, value)
           - ModulePcc::pcc_2d_potential(moments, y - step, value))
          / (2.0 * step);

    EXPECT_DOUBLE_EQ(gradient.x, 0.0);
    EXPECT_NEAR(gradient.y, finite, 1.0e-11);
    EXPECT_DOUBLE_EQ(gradient.z, 0.0);

    geometry.axis = 2;
    geometry.normal = ModuleBase::Vector3<double>(0.0, 0.0, 1.0);
    const ModuleBase::Vector3<double> along_z
        = ModulePcc::pcc_2d_potential_gradient(moments, y, geometry);
    EXPECT_DOUBLE_EQ(along_z.x, 0.0);
    EXPECT_DOUBLE_EQ(along_z.y, 0.0);
    EXPECT_NEAR(along_z.z, finite, 1.0e-11);
}

TEST(SccsPcc2d, BilinearKernelIsSymmetric)
{
    ModulePcc::Pcc2dMoments left;
    left.charge = 1.0;
    left.dipole = 0.2;
    left.quadrupole = 2.0;
    ModulePcc::Pcc2dMoments right;
    right.charge = -2.0;
    right.dipole = 0.7;
    right.quadrupole = -1.0;
    const ModulePcc::Pcc2dParameters value = parameters();

    EXPECT_DOUBLE_EQ(ModulePcc::pcc_2d_bilinear_energy(left, right, value),
                     ModulePcc::pcc_2d_bilinear_energy(right, left, value));
}

TEST(SccsPcc2d, NeutralSelfEnergyHasExpectedDipoleForm)
{
    ModulePcc::Pcc2dMoments moments;
    moments.dipole = 0.7;
    moments.quadrupole = 3.2;
    const ModulePcc::Pcc2dParameters value = parameters();
    const double expected = 2.0 * ModuleBase::PI * moments.dipole * moments.dipole
                            / (value.periodic_area * value.cell_length);

    EXPECT_NEAR(ModulePcc::pcc_2d_self_energy(moments, value), expected, 1.0e-15);
}

TEST(SccsPcc2d, PotentialAndEnergyAreIndependentOfMomentOrigin)
{
    std::vector<ModulePcc::PointCharge> charges(2);
    charges[0].charge = 1.7;
    charges[0].position = ModuleBase::Vector3<double>(1.2, -0.5, 0.8);
    charges[1].charge = -0.4;
    charges[1].position = ModuleBase::Vector3<double>(-0.3, 0.6, 1.1);
    const double shifted_origin = 0.7;
    const double evaluation_y = 2.1;
    const ModulePcc::Pcc2dMoments first
        = ModulePcc::pcc_2d_point_charge_moments(charges, geometry(0.0));
    const ModulePcc::Pcc2dMoments second
        = ModulePcc::pcc_2d_point_charge_moments(charges, geometry(shifted_origin));
    const ModulePcc::Pcc2dParameters value = parameters();

    EXPECT_NEAR(ModulePcc::pcc_2d_potential(first, evaluation_y, value),
                ModulePcc::pcc_2d_potential(second,
                                             evaluation_y - shifted_origin,
                                             value),
                1.0e-15);
    EXPECT_NEAR(ModulePcc::pcc_2d_self_energy(first, value),
                ModulePcc::pcc_2d_self_energy(second, value),
                1.0e-15);
}

TEST(SccsPcc2d, DensityDerivativeMatchesPotential)
{
    const double origin_y = 0.0;
    const double electron_y = -0.8;
    const double volume_element = 0.3;
    const double electron_density = 0.7;
    const double step = 1.0e-6;
    const ModulePcc::Pcc2dParameters value = parameters();

    const auto energy = [&](const double density) {
        std::vector<ModulePcc::PointCharge> charges(2);
        charges[0].charge = 1.0;
        charges[0].position.y = 0.5;
        charges[1].charge = -density * volume_element;
        charges[1].position.y = electron_y;
        const ModulePcc::Pcc2dMoments moments
            = ModulePcc::pcc_2d_point_charge_moments(charges, geometry(origin_y));
        return ModulePcc::pcc_2d_self_energy(moments, value);
    };

    std::vector<ModulePcc::PointCharge> center_charges(2);
    center_charges[0].charge = 1.0;
    center_charges[0].position.y = 0.5;
    center_charges[1].charge = -electron_density * volume_element;
    center_charges[1].position.y = electron_y;
    const ModulePcc::Pcc2dMoments center
        = ModulePcc::pcc_2d_point_charge_moments(center_charges, geometry(origin_y));
    const double analytic = -volume_element
                            * ModulePcc::pcc_2d_potential(center, electron_y, value);
    const double finite = (energy(electron_density + step)
                           - energy(electron_density - step))
                          / (2.0 * step);
    EXPECT_NEAR(analytic, finite, 1.0e-11);
}

TEST(SccsPcc2d, PointChargeForceMatchesSelfEnergyFiniteDifference)
{
    const double origin_y = 0.0;
    std::vector<ModulePcc::PointCharge> charges(2);
    charges[0].charge = 1.2;
    charges[0].position = ModuleBase::Vector3<double>(0.7, -0.4, 0.2);
    charges[1].charge = -0.3;
    charges[1].position = ModuleBase::Vector3<double>(-0.6, 0.5, 0.9);
    const ModulePcc::Pcc2dParameters value = parameters();
    const ModulePcc::Pcc2dMoments moments
        = ModulePcc::pcc_2d_point_charge_moments(charges, geometry(origin_y));
    const ModuleBase::Vector3<double> analytic
        = ModulePcc::pcc_2d_point_charge_force(moments,
                                                charges[0],
                                                geometry(origin_y));
    const double step = 1.0e-5;

    std::vector<ModulePcc::PointCharge> plus = charges;
    std::vector<ModulePcc::PointCharge> minus = charges;
    plus[0].position.y += step;
    minus[0].position.y -= step;
    const double energy_plus
        = ModulePcc::pcc_2d_self_energy(
            ModulePcc::pcc_2d_point_charge_moments(plus, geometry(origin_y)),
            value);
    const double energy_minus
        = ModulePcc::pcc_2d_self_energy(
            ModulePcc::pcc_2d_point_charge_moments(minus, geometry(origin_y)),
            value);
    const double finite_force = -(energy_plus - energy_minus) / (2.0 * step);
    EXPECT_DOUBLE_EQ(analytic.x, 0.0);
    EXPECT_NEAR(analytic.y, finite_force, 1.0e-11);
    EXPECT_DOUBLE_EQ(analytic.z, 0.0);
}

TEST(SccsPcc2d, RejectsInvalidAndNonFiniteInputs)
{
    ModulePcc::Pcc2dParameters value;
    EXPECT_THROW(ModulePcc::pcc_2d_self_energy(ModulePcc::Pcc2dMoments(), value),
                 std::invalid_argument);
    value = parameters();
    ModulePcc::Pcc2dMoments moments;
    moments.charge = std::numeric_limits<double>::quiet_NaN();
    ModulePcc::Pcc2dGeometry geometry;
    geometry.parameters = value;
    ModulePcc::PointCharge ion;
    ion.charge = 1.0;
    EXPECT_THROW(ModulePcc::pcc_2d_point_charge_force(moments, ion, geometry), std::domain_error);
    geometry.origin = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(ModulePcc::validate_pcc_2d_geometry(geometry), std::invalid_argument);
}

} // namespace
