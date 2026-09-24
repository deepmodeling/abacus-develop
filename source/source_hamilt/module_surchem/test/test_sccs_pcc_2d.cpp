#include "../sccs_pcc_2d.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "gtest/gtest.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{

ModuleSccs::Pcc2dParameters parameters()
{
    ModuleSccs::Pcc2dParameters value;
    value.periodic_area = 31.0;
    value.cell_length_y = 17.0;
    return value;
}

ModuleSccs::Pcc2dGeometry geometry(const double origin_y)
{
    ModuleSccs::Pcc2dGeometry value;
    value.parameters = parameters();
    value.origin_y = origin_y;
    return value;
}

TEST(SccsPcc2d, BuildsGeometryForCubicCell)
{
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    const ModuleSccs::Pcc2dGeometry geometry
        = ModuleSccs::pcc_2d_geometry(lattice, 10.0, 1.0e-10);

    EXPECT_DOUBLE_EQ(geometry.parameters.periodic_area, 100.0);
    EXPECT_DOUBLE_EQ(geometry.parameters.cell_length_y, 10.0);
    EXPECT_DOUBLE_EQ(geometry.origin_y, 5.0);
}

TEST(SccsPcc2d, BuildsGeometryForNonOrthogonalPeriodicPlane)
{
    const ModuleBase::Matrix3 lattice(2.0, 0.0, 0.0,
                                      0.0, 5.0, 0.0,
                                      1.0, 0.0, 3.0);
    const ModuleSccs::Pcc2dGeometry geometry
        = ModuleSccs::pcc_2d_geometry(lattice, 4.0, 1.0e-10);

    EXPECT_DOUBLE_EQ(geometry.parameters.periodic_area, 96.0);
    EXPECT_DOUBLE_EQ(geometry.parameters.cell_length_y, 20.0);
    EXPECT_DOUBLE_EQ(geometry.origin_y, 10.0);
}

TEST(SccsPcc2d, GeometryAcceptsNumericalAlignmentNoise)
{
    const ModuleBase::Matrix3 lattice(1.0, 1.0e-12, 0.0,
                                      -1.0e-12, 2.0, 1.0e-12,
                                      0.3, -1.0e-12, 1.0);
    EXPECT_NO_THROW(ModuleSccs::pcc_2d_geometry(lattice, 7.0, 1.0e-10));
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

    EXPECT_THROW(ModuleSccs::pcc_2d_geometry(tilted_normal, 1.0, 1.0e-10),
                 std::invalid_argument);
    EXPECT_THROW(ModuleSccs::pcc_2d_geometry(tilted_plane, 1.0, 1.0e-10),
                 std::invalid_argument);
    EXPECT_THROW(ModuleSccs::pcc_2d_geometry(reversed_normal, 1.0, 1.0e-10),
                 std::invalid_argument);
    EXPECT_THROW(ModuleSccs::pcc_2d_geometry(degenerate_plane, 1.0, 1.0e-10),
                 std::invalid_argument);
}

TEST(SccsPcc2d, SystemCenterUnwrapsAcrossYBoundary)
{
    const std::vector<double> positions_y{9.8, 0.2};
    const std::vector<double> weights{1.0, 3.0};
    EXPECT_NEAR(ModuleSccs::pcc_2d_system_center_y(positions_y, weights, 10.0),
                0.1,
                1.0e-14);

    ModuleSccs::Pcc2dGeometry value = geometry(0.1);
    value.parameters.cell_length_y = 10.0;
    EXPECT_NEAR(ModuleSccs::pcc_2d_relative_y(9.8, value), -0.3, 1.0e-14);
    EXPECT_NEAR(ModuleSccs::pcc_2d_relative_y(0.2, value), 0.1, 1.0e-14);
}

TEST(SccsPcc2d, MomentsRemainInvariantWhenChargesCrossYBoundary)
{
    std::vector<ModuleSccs::PointCharge> original(2);
    original[0].charge = 1.2;
    original[0].position.y = 9.8;
    original[1].charge = -0.4;
    original[1].position.y = 0.2;
    const std::vector<double> weights{1.0, 3.0};
    const std::vector<double> original_y{9.8, 0.2};
    ModuleSccs::Pcc2dGeometry first = geometry(
        ModuleSccs::pcc_2d_system_center_y(original_y, weights, 10.0));
    first.parameters.cell_length_y = 10.0;

    std::vector<ModuleSccs::PointCharge> shifted = original;
    shifted[0].position.y = 3.5;
    shifted[1].position.y = 3.9;
    const std::vector<double> shifted_y{3.5, 3.9};
    ModuleSccs::Pcc2dGeometry second = geometry(
        ModuleSccs::pcc_2d_system_center_y(shifted_y, weights, 10.0));
    second.parameters.cell_length_y = 10.0;

    const ModuleSccs::Pcc2dMoments original_moments
        = ModuleSccs::pcc_2d_point_charge_moments(original, first);
    const ModuleSccs::Pcc2dMoments shifted_moments
        = ModuleSccs::pcc_2d_point_charge_moments(shifted, second);
    EXPECT_NEAR(shifted_moments.charge, original_moments.charge, 1.0e-15);
    EXPECT_NEAR(shifted_moments.dipole_y, original_moments.dipole_y, 1.0e-14);
    EXPECT_NEAR(shifted_moments.quadrupole_yy,
                original_moments.quadrupole_yy,
                1.0e-14);
}

TEST(SccsPcc2d, AccumulatesOnlyYMoments)
{
    std::vector<ModuleSccs::PointCharge> charges(2);
    charges[0].charge = 2.0;
    charges[0].position = ModuleBase::Vector3<double>(4.0, 3.0, -7.0);
    charges[1].charge = -1.0;
    charges[1].position = ModuleBase::Vector3<double>(-8.0, -1.0, 9.0);

    const ModuleSccs::Pcc2dMoments moments
        = ModuleSccs::pcc_2d_point_charge_moments(charges, geometry(1.0));
    EXPECT_DOUBLE_EQ(moments.charge, 1.0);
    EXPECT_DOUBLE_EQ(moments.dipole_y, 6.0);
    EXPECT_DOUBLE_EQ(moments.quadrupole_yy, 4.0);
}

TEST(SccsPcc2d, DensityMomentsIncludeTheVolumeElement)
{
    const std::vector<double> density{0.5, -0.25};
    const std::vector<ModuleBase::Vector3<double>> positions{
        ModuleBase::Vector3<double>(7.0, 1.0, -4.0),
        ModuleBase::Vector3<double>(-2.0, 3.0, 6.0)};
    const ModuleSccs::Pcc2dMoments moments
        = ModuleSccs::pcc_2d_density_moments(density, positions, 2.0, geometry(0.0));
    EXPECT_DOUBLE_EQ(moments.charge, 0.5);
    EXPECT_DOUBLE_EQ(moments.dipole_y, -0.5);
    EXPECT_DOUBLE_EQ(moments.quadrupole_yy, -3.5);
}

TEST(SccsPcc2d, MomentPotentialMatchesIndependentPlanarGreenFunctions)
{
    const ModuleSccs::Pcc2dParameters value = parameters();
    const std::vector<double> locations{-3.4, -1.1, 0.5, 2.8};
    const std::vector<double> charges{0.7, -0.2, 0.3, -0.1};
    std::vector<ModuleSccs::PointCharge> points(charges.size());
    for (std::size_t index = 0; index < charges.size(); ++index)
    {
        points[index].charge = charges[index];
        points[index].position.y = locations[index];
    }
    const ModuleSccs::Pcc2dMoments moments
        = ModuleSccs::pcc_2d_point_charge_moments(points, geometry(0.0));
    const std::vector<double> evaluation_points{-5.2, -0.8, 1.4, 5.1};

    for (std::size_t evaluation = 0; evaluation < evaluation_points.size(); ++evaluation)
    {
        const double y = evaluation_points[evaluation];
        double direct = 0.0;
        for (std::size_t source = 0; source < charges.size(); ++source)
        {
            const double separation = y - locations[source];
            ASSERT_LE(std::abs(separation), 0.5 * value.cell_length_y);
            const double gauge_shift
                = ModuleBase::PI / (3.0 * value.cell_length_y)
                  + ModuleBase::PI * value.cell_length_y
                        / (3.0 * value.periodic_area);
            const double open_kernel
                = -2.0 * ModuleBase::PI * std::abs(separation) / value.periodic_area
                  + gauge_shift;
            const double periodic_kernel
                = 2.0 * ModuleBase::PI / value.periodic_area
                  * (separation * separation / value.cell_length_y
                     - std::abs(separation) + value.cell_length_y / 6.0);
            direct += charges[source] * (open_kernel - periodic_kernel);
        }
        EXPECT_NEAR(ModuleSccs::pcc_2d_potential(moments, y, value), direct, 1.0e-14);
    }
}

TEST(SccsPcc2d, ChargedMonopoleUsesTheLiteratureGauge)
{
    ModuleSccs::Pcc2dMoments moments;
    moments.charge = 1.7;
    const ModuleSccs::Pcc2dParameters value = parameters();
    const double expected_potential
        = ModuleBase::PI * moments.charge / (3.0 * value.cell_length_y);
    const double expected_energy
        = ModuleBase::PI * moments.charge * moments.charge
          / (6.0 * value.cell_length_y);

    EXPECT_NEAR(ModuleSccs::pcc_2d_potential(moments, 0.0, value),
                expected_potential,
                1.0e-15);
    EXPECT_NEAR(ModuleSccs::pcc_2d_self_energy(moments, value),
                expected_energy,
                1.0e-15);
}

TEST(SccsPcc2d, PositiveAndNegativeMonopolesHaveOddPotentialAndEqualEnergy)
{
    ModuleSccs::Pcc2dMoments positive;
    positive.charge = 0.8;
    positive.dipole_y = -0.3;
    positive.quadrupole_yy = 1.4;
    ModuleSccs::Pcc2dMoments negative;
    negative.charge = -positive.charge;
    negative.dipole_y = -positive.dipole_y;
    negative.quadrupole_yy = -positive.quadrupole_yy;
    const ModuleSccs::Pcc2dParameters value = parameters();
    const double y = 2.6;

    EXPECT_DOUBLE_EQ(ModuleSccs::pcc_2d_potential(negative, y, value),
                     -ModuleSccs::pcc_2d_potential(positive, y, value));
    EXPECT_DOUBLE_EQ(ModuleSccs::pcc_2d_self_energy(negative, value),
                     ModuleSccs::pcc_2d_self_energy(positive, value));
}

TEST(SccsPcc2d, GradientIsYOnlyAndMatchesCentralDifference)
{
    ModuleSccs::Pcc2dMoments moments;
    moments.charge = -1.2;
    moments.dipole_y = 0.7;
    moments.quadrupole_yy = 4.1;
    const ModuleSccs::Pcc2dParameters value = parameters();
    const double y = -1.3;
    const double step = 1.0e-5;
    const ModuleBase::Vector3<double> gradient
        = ModuleSccs::pcc_2d_potential_gradient(moments, y, value);
    const double finite
        = (ModuleSccs::pcc_2d_potential(moments, y + step, value)
           - ModuleSccs::pcc_2d_potential(moments, y - step, value))
          / (2.0 * step);

    EXPECT_DOUBLE_EQ(gradient.x, 0.0);
    EXPECT_NEAR(gradient.y, finite, 1.0e-11);
    EXPECT_DOUBLE_EQ(gradient.z, 0.0);
}

TEST(SccsPcc2d, BilinearKernelIsSymmetric)
{
    ModuleSccs::Pcc2dMoments left;
    left.charge = 1.0;
    left.dipole_y = 0.2;
    left.quadrupole_yy = 2.0;
    ModuleSccs::Pcc2dMoments right;
    right.charge = -2.0;
    right.dipole_y = 0.7;
    right.quadrupole_yy = -1.0;
    const ModuleSccs::Pcc2dParameters value = parameters();

    EXPECT_DOUBLE_EQ(ModuleSccs::pcc_2d_bilinear_energy(left, right, value),
                     ModuleSccs::pcc_2d_bilinear_energy(right, left, value));
}

TEST(SccsPcc2d, NeutralSelfEnergyHasExpectedDipoleForm)
{
    ModuleSccs::Pcc2dMoments moments;
    moments.dipole_y = 0.7;
    moments.quadrupole_yy = 3.2;
    const ModuleSccs::Pcc2dParameters value = parameters();
    const double expected = 2.0 * ModuleBase::PI * moments.dipole_y * moments.dipole_y
                            / (value.periodic_area * value.cell_length_y);

    EXPECT_NEAR(ModuleSccs::pcc_2d_self_energy(moments, value), expected, 1.0e-15);
}

TEST(SccsPcc2d, PotentialAndEnergyAreIndependentOfMomentOrigin)
{
    std::vector<ModuleSccs::PointCharge> charges(2);
    charges[0].charge = 1.7;
    charges[0].position = ModuleBase::Vector3<double>(1.2, -0.5, 0.8);
    charges[1].charge = -0.4;
    charges[1].position = ModuleBase::Vector3<double>(-0.3, 0.6, 1.1);
    const double shifted_origin = 0.7;
    const double evaluation_y = 2.1;
    const ModuleSccs::Pcc2dMoments first
        = ModuleSccs::pcc_2d_point_charge_moments(charges, geometry(0.0));
    const ModuleSccs::Pcc2dMoments second
        = ModuleSccs::pcc_2d_point_charge_moments(charges, geometry(shifted_origin));
    const ModuleSccs::Pcc2dParameters value = parameters();

    EXPECT_NEAR(ModuleSccs::pcc_2d_potential(first, evaluation_y, value),
                ModuleSccs::pcc_2d_potential(second,
                                             evaluation_y - shifted_origin,
                                             value),
                1.0e-15);
    EXPECT_NEAR(ModuleSccs::pcc_2d_self_energy(first, value),
                ModuleSccs::pcc_2d_self_energy(second, value),
                1.0e-15);
}

TEST(SccsPcc2d, DensityDerivativeMatchesPotential)
{
    const double origin_y = 0.0;
    const double electron_y = -0.8;
    const double volume_element = 0.3;
    const double electron_density = 0.7;
    const double step = 1.0e-6;
    const ModuleSccs::Pcc2dParameters value = parameters();

    const auto energy = [&](const double density) {
        std::vector<ModuleSccs::PointCharge> charges(2);
        charges[0].charge = 1.0;
        charges[0].position.y = 0.5;
        charges[1].charge = -density * volume_element;
        charges[1].position.y = electron_y;
        const ModuleSccs::Pcc2dMoments moments
            = ModuleSccs::pcc_2d_point_charge_moments(charges, geometry(origin_y));
        return ModuleSccs::pcc_2d_self_energy(moments, value);
    };

    std::vector<ModuleSccs::PointCharge> center_charges(2);
    center_charges[0].charge = 1.0;
    center_charges[0].position.y = 0.5;
    center_charges[1].charge = -electron_density * volume_element;
    center_charges[1].position.y = electron_y;
    const ModuleSccs::Pcc2dMoments center
        = ModuleSccs::pcc_2d_point_charge_moments(center_charges, geometry(origin_y));
    const double analytic = -volume_element
                            * ModuleSccs::pcc_2d_potential(center, electron_y, value);
    const double finite = (energy(electron_density + step)
                           - energy(electron_density - step))
                          / (2.0 * step);
    EXPECT_NEAR(analytic, finite, 1.0e-11);
}

TEST(SccsPcc2d, IonicShapeEnergyMatchesAppendixA2)
{
    const ModuleSccs::Pcc2dParameters value = parameters();
    const double polarization_charge = -0.8;
    ModuleSccs::Pcc2dMoments smooth_ionic;
    smooth_ionic.charge = 2.0;
    smooth_ionic.dipole_y = 6.0;
    smooth_ionic.quadrupole_yy = 6.5;
    ModuleSccs::Pcc2dMoments point_ionic;
    point_ionic.charge = 2.0;
    point_ionic.dipole_y = 6.0;
    point_ionic.quadrupole_yy = 5.2;
    const double expected = ModuleBase::PI * polarization_charge
                            * (smooth_ionic.quadrupole_yy
                               - point_ionic.quadrupole_yy)
                            / (value.periodic_area * value.cell_length_y);
    EXPECT_NEAR(ModuleSccs::pcc_2d_ionic_shape_energy(polarization_charge,
                                                       smooth_ionic,
                                                       point_ionic,
                                                       value),
                expected,
                1.0e-15);
    EXPECT_THROW(ModuleSccs::pcc_2d_ionic_shape_energy(
                     std::numeric_limits<double>::infinity(),
                     smooth_ionic,
                     point_ionic,
                     value),
                 std::domain_error);
    point_ionic.charge = 0.0;
    EXPECT_THROW(ModuleSccs::pcc_2d_ionic_shape_energy(polarization_charge,
                                                       smooth_ionic,
                                                       point_ionic,
                                                       value),
                 std::domain_error);
}

TEST(SccsPcc2d, IonicShapeEnergyIsInvariantUnderRigidTranslation)
{
    const ModuleSccs::Pcc2dParameters value = parameters();
    ModuleSccs::Pcc2dMoments smooth_ionic;
    smooth_ionic.charge = 4.0;
    smooth_ionic.dipole_y = -2.0;
    smooth_ionic.quadrupole_yy = 7.5;
    ModuleSccs::Pcc2dMoments point_ionic;
    point_ionic.charge = 4.0;
    point_ionic.dipole_y = -2.4;
    point_ionic.quadrupole_yy = 5.1;
    const double reference
        = ModuleSccs::pcc_2d_ionic_shape_energy(-0.7,
                                                smooth_ionic,
                                                point_ionic,
                                                value);
    const double translation = 3.2;
    const auto translate = [translation](const ModuleSccs::Pcc2dMoments& moments) {
        ModuleSccs::Pcc2dMoments shifted;
        shifted.charge = moments.charge;
        shifted.dipole_y = moments.dipole_y + translation * moments.charge;
        shifted.quadrupole_yy
            = moments.quadrupole_yy + 2.0 * translation * moments.dipole_y
              + translation * translation * moments.charge;
        return shifted;
    };
    EXPECT_NEAR(ModuleSccs::pcc_2d_ionic_shape_energy(-0.7,
                                                       translate(smooth_ionic),
                                                       translate(point_ionic),
                                                       value),
                reference,
                1.0e-15);
}

TEST(SccsPcc2d, PointChargeForceMatchesSelfEnergyFiniteDifference)
{
    const double origin_y = 0.0;
    std::vector<ModuleSccs::PointCharge> charges(2);
    charges[0].charge = 1.2;
    charges[0].position = ModuleBase::Vector3<double>(0.7, -0.4, 0.2);
    charges[1].charge = -0.3;
    charges[1].position = ModuleBase::Vector3<double>(-0.6, 0.5, 0.9);
    const ModuleSccs::Pcc2dParameters value = parameters();
    const ModuleSccs::Pcc2dMoments moments
        = ModuleSccs::pcc_2d_point_charge_moments(charges, geometry(origin_y));
    const ModuleBase::Vector3<double> analytic
        = ModuleSccs::pcc_2d_point_charge_force(moments,
                                                charges[0],
                                                geometry(origin_y));
    const double step = 1.0e-5;

    std::vector<ModuleSccs::PointCharge> plus = charges;
    std::vector<ModuleSccs::PointCharge> minus = charges;
    plus[0].position.y += step;
    minus[0].position.y -= step;
    const double energy_plus
        = ModuleSccs::pcc_2d_self_energy(
            ModuleSccs::pcc_2d_point_charge_moments(plus, geometry(origin_y)),
            value);
    const double energy_minus
        = ModuleSccs::pcc_2d_self_energy(
            ModuleSccs::pcc_2d_point_charge_moments(minus, geometry(origin_y)),
            value);
    const double finite_force = -(energy_plus - energy_minus) / (2.0 * step);
    EXPECT_DOUBLE_EQ(analytic.x, 0.0);
    EXPECT_NEAR(analytic.y, finite_force, 1.0e-11);
    EXPECT_DOUBLE_EQ(analytic.z, 0.0);
}

TEST(SccsPcc2d, RejectsInvalidAndNonFiniteInputs)
{
    ModuleSccs::Pcc2dParameters value;
    EXPECT_THROW(ModuleSccs::pcc_2d_self_energy(ModuleSccs::Pcc2dMoments(), value),
                 std::invalid_argument);
    value = parameters();
    ModuleSccs::Pcc2dMoments moments;
    moments.charge = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(ModuleSccs::pcc_2d_potential(moments, 0.0, value), std::domain_error);
}

} // namespace
