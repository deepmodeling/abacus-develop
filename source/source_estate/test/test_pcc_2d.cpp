#include "source_estate/pcc_2d.h"
#include "source_base/constants.h"
#include "source_cell/cell_geometry.h"

#include <gtest/gtest.h>

#include <cmath>

namespace
{
/// Hartree energy per cell of periodic Gaussian sheets (surface charge
/// `density`, normal width sigma) with the G = 0 term excluded. One sheet at
/// the origin, or a +density/-density pair separated by `separation` along the
/// normal when `dipole_pair` is true. Only G parallel to the normal contributes.
double periodic_sheet_energy(const double area,
                             const double length,
                             const double sigma,
                             const double density,
                             const double separation,
                             const bool dipole_pair)
{
    const double step = 2.0 * ModuleBase::PI / length;
    const int max_index = static_cast<int>(8.0 / (sigma * step)) + 1;
    double sum = 0.0;
    for (int k = 1; k <= max_index; ++k)
    {
        const double g = step * k;
        double structure = 1.0;
        if (dipole_pair)
        {
            const double phase = std::sin(0.5 * g * separation);
            structure = 4.0 * phase * phase;
        }
        sum += density * density * structure * std::exp(-g * g * sigma * sigma) / (g * g);
    }
    return 4.0 * ModuleBase::PI * area * sum / length;
}
} // namespace

TEST(Pcc2d, MatchesChargedSheetGaugeAndNeutralDipoleEnergy)
{
    unitcell::SlabCell cell;
    cell.area = 6.0;
    cell.length = 8.0;
    const elecstate::Pcc2dParameters parameters = elecstate::make_pcc_2d_parameters(cell);
    elecstate::ChargeMoments moments;
    moments.charge = 2.0;
    const double charged = elecstate::pcc_2d_energy(moments, parameters);
    const double expected_charged = -ModuleBase::PI * 8.0 * 4.0 / 36.0;
    EXPECT_NEAR(charged, expected_charged, 1.0e-14);
    moments.charge = 0.0;
    moments.dipole.x = 2.0;
    const double neutral = elecstate::pcc_2d_energy(moments, parameters);
    const double expected_neutral = 2.0 * ModuleBase::PI * 4.0 / 48.0;
    EXPECT_NEAR(neutral, expected_neutral, 1.0e-14);
}

TEST(Pcc2d, NormalForceMatchesEnergyDerivativeAndHasNoTangentialComponent)
{
    elecstate::Pcc2dParameters parameters;
    parameters.area = 6.0;
    parameters.length = 8.0;
    const ModuleBase::Vector3<double> normal(0.0, 1.0, 0.0);
    const double charges[2] = {1.3, -0.4};
    ModuleBase::Vector3<double> positions[2] = {
        ModuleBase::Vector3<double>(0.0, 0.2, 0.0),
        ModuleBase::Vector3<double>(0.0, -0.5, 0.0)};
    elecstate::ChargeMoments moments;
    moments = elecstate::charge_moments(charges, positions, 2, 1.0);
    const double coordinate = positions[0].y;
    const ModuleBase::Vector3<double> force = elecstate::pcc_2d_force(moments, charges[0], coordinate, normal, parameters);
    EXPECT_DOUBLE_EQ(force.x, 0.0);
    EXPECT_DOUBLE_EQ(force.z, 0.0);
    const double step = 1.0e-5;
    positions[0].y = coordinate + step;
    moments = elecstate::charge_moments(charges, positions, 2, 1.0);
    const double plus = elecstate::pcc_2d_energy(moments, parameters);
    positions[0].y = coordinate - step;
    moments = elecstate::charge_moments(charges, positions, 2, 1.0);
    const double minus = elecstate::pcc_2d_energy(moments, parameters);
    const double numerical_force = -(plus - minus) / (2.0 * step);
    EXPECT_NEAR(force.y, numerical_force, 1.0e-10);
}

TEST(Pcc2d, ElectronicPotentialMatchesEnergyDerivative)
{
    elecstate::Pcc2dParameters parameters;
    parameters.area = 6.0;
    parameters.length = 8.0;
    const ModuleBase::Vector3<double> normal(1.0, 0.0, 0.0);
    const double coordinate = 0.3;
    elecstate::ChargeMoments moments;
    moments.charge = 1.2;
    moments.dipole.x = 0.2;
    moments.second_moment = 0.8;
    const double positive_potential = elecstate::pcc_2d_potential(moments, coordinate, normal, parameters);
    const double expected = -2.0 * positive_potential;
    const double step = 1.0e-5;
    elecstate::ChargeMoments plus = moments;
    plus.charge -= step;
    plus.dipole.x -= coordinate * step;
    plus.second_moment -= coordinate * coordinate * step;
    elecstate::ChargeMoments minus = moments;
    minus.charge += step;
    minus.dipole.x += coordinate * step;
    minus.second_moment += coordinate * coordinate * step;
    const double plus_energy = 2.0 * elecstate::pcc_2d_energy(plus, parameters);
    const double minus_energy = 2.0 * elecstate::pcc_2d_energy(minus, parameters);
    const double derivative = (plus_energy - minus_energy) / (2.0 * step);
    EXPECT_NEAR(derivative, expected, 1.0e-10);
}

// Physical check of the charged-slab gauge: a charged Gaussian layer in cells
// of growing vacuum. periodic + PCC must give the open-boundary energy with the
// potential referenced to zero on the plane of the charge, -2 sqrt(pi) A s^2 sigma,
// independently of the vacuum size.
TEST(Pcc2d, RecoversOpenBoundaryEnergyOfChargedLayerForAnyVacuum)
{
    const double area = 50.0;
    const double sigma = 1.0;
    const double density = 0.1;
    const double charge = area * density;
    const double isolated = -2.0 * std::sqrt(ModuleBase::PI) * area * density * density * sigma;
    const double lengths[3] = {20.0, 30.0, 40.0};
    for (int index = 0; index < 3; ++index)
    {
        elecstate::Pcc2dParameters parameters;
        parameters.area = area;
        parameters.length = lengths[index];
        elecstate::ChargeMoments layer;
        layer.charge = charge;
        layer.second_moment = charge * sigma * sigma;
        const double periodic = periodic_sheet_energy(area, lengths[index], sigma, density, 0.0, false);
        const double correction = elecstate::pcc_2d_energy(layer, parameters);
        EXPECT_NEAR(periodic + correction, isolated, 1.0e-10);
    }
}

// Polar neutral slab modelled as two opposite Gaussian layers. Periodic images
// impose a spurious zero-average field that changes the energy by ~1/L; PCC
// restores the isolated dipole-layer energy for every vacuum size.
TEST(Pcc2d, RecoversIsolatedEnergyOfPolarNeutralSlabForAnyVacuum)
{
    const double area = 50.0;
    const double sigma = 1.0;
    const double density = 0.1;
    const double separation = 3.0;
    const double same_layer = 2.0 * sigma / std::sqrt(ModuleBase::PI);
    const double gap = 0.5 * separation / sigma;
    const double cross_layer = separation * std::erf(gap) + same_layer * std::exp(-gap * gap);
    const double isolated = -2.0 * ModuleBase::PI * area * density * density * (same_layer - cross_layer);
    const double lengths[3] = {20.0, 30.0, 40.0};
    for (int index = 0; index < 3; ++index)
    {
        elecstate::Pcc2dParameters parameters;
        parameters.area = area;
        parameters.length = lengths[index];
        elecstate::ChargeMoments slab;
        slab.dipole.z = area * density * separation;
        const double periodic = periodic_sheet_energy(area, lengths[index], sigma, density, separation, true);
        const double correction = elecstate::pcc_2d_energy(slab, parameters);
        const double periodic_error = std::abs(periodic - isolated);
        EXPECT_GT(periodic_error, 0.1);
        EXPECT_NEAR(periodic + correction, isolated, 1.0e-10);
    }
}
