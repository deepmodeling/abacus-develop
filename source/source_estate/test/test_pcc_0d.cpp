#include "source_estate/pcc_0d.h"
#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_cell/cell_geometry.h"

#include <gtest/gtest.h>

#include <cmath>

namespace
{
/// Hartree energy per cell of periodic Gaussian charges (width sigma) under the
/// usual G = 0 exclusion. One charge at the origin, or a +q/-q pair separated
/// by `separation` along z when `dipole_pair` is true.
double periodic_gaussian_energy(const double length,
                                const double sigma,
                                const double charge,
                                const double separation,
                                const bool dipole_pair)
{
    const double step = 2.0 * ModuleBase::PI / length;
    const int max_index = static_cast<int>(8.0 / (sigma * step)) + 1;
    const double volume = length * length * length;
    double sum = 0.0;
    for (int i = -max_index; i <= max_index; ++i)
    {
        for (int j = -max_index; j <= max_index; ++j)
        {
            for (int k = -max_index; k <= max_index; ++k)
            {
                if (i == 0 && j == 0 && k == 0)
                {
                    continue;
                }
                const double gx = step * i;
                const double gy = step * j;
                const double gz = step * k;
                const double g2 = gx * gx + gy * gy + gz * gz;
                double structure = 1.0;
                if (dipole_pair)
                {
                    const double phase = std::sin(0.5 * gz * separation);
                    structure = 4.0 * phase * phase;
                }
                sum += charge * charge * structure * std::exp(-g2 * sigma * sigma) / g2;
            }
        }
    }
    return 2.0 * ModuleBase::PI * sum / volume;
}

class Pcc0dTest : public testing::Test
{
  protected:
    void SetUp() override
    {
        const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                         0.0, 1.0, 0.0,
                                         0.0, 0.0, 1.0);
        ASSERT_TRUE(unitcell::make_orthogonal_cell(lattice, 10.0, 1.0e-10, cell));
        ASSERT_TRUE(elecstate::make_pcc_0d_parameters(cell, 1.0e-10, parameters));
        moments.charge = 1.2;
        moments.dipole = ModuleBase::Vector3<double>(0.1, -0.2, 0.3);
        moments.second_moment = 0.8;
    }

    unitcell::OrthogonalCell cell;
    elecstate::Pcc0dParameters parameters;
    elecstate::ChargeMoments moments;
};
} // namespace

TEST_F(Pcc0dTest, SeparatesCubicRestrictionFromOrthogonalGeometry)
{
    EXPECT_DOUBLE_EQ(parameters.length, 10.0);
    const ModuleBase::Matrix3 rectangular(1.0, 0.0, 0.0,
                                         0.0, 2.0, 0.0,
                                         0.0, 0.0, 1.0);
    ASSERT_TRUE(unitcell::make_orthogonal_cell(rectangular, 10.0, 1.0e-10, cell));
    EXPECT_FALSE(elecstate::make_pcc_0d_parameters(cell, 1.0e-10, parameters));
}

// The gradient, the ionic force and the electron potential are derivatives of
// the potential and of the energy with respect to a point charge.
TEST_F(Pcc0dTest, DerivativesMatchFiniteDifferences)
{
    const double step = 1.0e-5;
    const ModuleBase::Vector3<double> position(0.4, -0.5, 0.6);
    const ModuleBase::Vector3<double> gradient = elecstate::pcc_0d_gradient(moments, position, parameters);
    const double charges[2] = {1.3, -0.4};
    ModuleBase::Vector3<double> positions[2] = {position, ModuleBase::Vector3<double>(-0.5, 0.6, -0.7)};
    const elecstate::ChargeMoments ions = elecstate::charge_moments(charges, positions, 2, 1.0);
    const ModuleBase::Vector3<double> force = elecstate::pcc_0d_force(ions, charges[0], position, parameters);
    for (int axis = 0; axis < 3; ++axis)
    {
        ModuleBase::Vector3<double> plus = position;
        ModuleBase::Vector3<double> minus = position;
        plus[axis] += step;
        minus[axis] -= step;
        const double plus_potential = elecstate::pcc_0d_potential(moments, plus, parameters);
        const double minus_potential = elecstate::pcc_0d_potential(moments, minus, parameters);
        const double potential_slope = (plus_potential - minus_potential) / (2.0 * step);
        EXPECT_NEAR(potential_slope, gradient[axis], 1.0e-10);
        positions[0] = plus;
        const elecstate::ChargeMoments plus_ions = elecstate::charge_moments(charges, positions, 2, 1.0);
        positions[0] = minus;
        const elecstate::ChargeMoments minus_ions = elecstate::charge_moments(charges, positions, 2, 1.0);
        const double plus_energy = elecstate::pcc_0d_energy(plus_ions, parameters);
        const double minus_energy = elecstate::pcc_0d_energy(minus_ions, parameters);
        const double numerical_force = -(plus_energy - minus_energy) / (2.0 * step);
        EXPECT_NEAR(numerical_force, force[axis], 1.0e-10);
    }
    // Adding electrons -step at position changes the energy by the electron potential (Ry).
    const double positive_potential = elecstate::pcc_0d_potential(moments, position, parameters);
    const double electron_potential_ry = -2.0 * positive_potential;
    elecstate::ChargeMoments plus = moments;
    plus.charge -= step;
    plus.dipole -= position * step;
    plus.second_moment -= step * position.norm2();
    elecstate::ChargeMoments minus = moments;
    minus.charge += step;
    minus.dipole += position * step;
    minus.second_moment += step * position.norm2();
    const double plus_energy = 2.0 * elecstate::pcc_0d_energy(plus, parameters);
    const double minus_energy = 2.0 * elecstate::pcc_0d_energy(minus, parameters);
    const double energy_slope = (plus_energy - minus_energy) / (2.0 * step);
    EXPECT_NEAR(energy_slope, electron_potential_ry, 1.0e-10);
}

TEST_F(Pcc0dTest, EnergyIsIndependentOfMultipoleOrigin)
{
    const double charges[2] = {1.2, -0.3};
    ModuleBase::Vector3<double> positions[2] = {
        ModuleBase::Vector3<double>(0.3, -0.2, 0.1),
        ModuleBase::Vector3<double>(-0.5, 0.6, -0.7)};
    moments = elecstate::charge_moments(charges, positions, 2, 1.0);
    const double original_energy = elecstate::pcc_0d_energy(moments, parameters);
    const ModuleBase::Vector3<double> shift(0.4, -0.3, 0.2);
    positions[0] -= shift;
    positions[1] -= shift;
    moments = elecstate::charge_moments(charges, positions, 2, 1.0);
    const double shifted_energy = elecstate::pcc_0d_energy(moments, parameters);
    EXPECT_NEAR(original_energy, shifted_energy, 1.0e-14);
}

// Physical check of the model: a charged Gaussian ion (H3O+-like) in cubic
// cells of growing size. The periodic Hartree energy drifts as 1/L, while
// periodic + PCC reproduces the isolated self-energy q^2 / (2 sigma sqrt(pi)).
TEST_F(Pcc0dTest, RecoversIsolatedEnergyOfChargedGaussianForAnyCellSize)
{
    const double sigma = 1.0;
    const double charge = 1.0;
    const double isolated = charge * charge / (2.0 * sigma * std::sqrt(ModuleBase::PI));
    const double lengths[3] = {16.0, 20.0, 24.0};
    for (int index = 0; index < 3; ++index)
    {
        const double length = lengths[index];
        parameters.length = length;
        elecstate::ChargeMoments gaussian;
        gaussian.charge = charge;
        gaussian.second_moment = 3.0 * sigma * sigma * charge;
        const double periodic = periodic_gaussian_energy(length, sigma, charge, 0.0, false);
        const double correction = elecstate::pcc_0d_energy(gaussian, parameters);
        const double periodic_error = std::abs(periodic - isolated);
        EXPECT_GT(periodic_error, 1.0e-2);
        EXPECT_NEAR(periodic + correction, isolated, 1.0e-10);
    }
}

// Neutral polar molecule modelled as a +q/-q Gaussian pair. PCC removes the
// leading 1/L^3 dipole interaction; the residual error comes from higher
// multipoles and must be far smaller and decay faster than the uncorrected one.
TEST_F(Pcc0dTest, RemovesDipoleImageInteractionOfNeutralGaussianPair)
{
    const double sigma = 1.0;
    const double charge = 1.0;
    const double separation = 2.0;
    const double self_energy = charge * charge / (2.0 * sigma * std::sqrt(ModuleBase::PI));
    const double attraction = charge * charge * std::erf(0.5 * separation / sigma) / separation;
    const double isolated = 2.0 * self_energy - attraction;
    elecstate::ChargeMoments pair;
    pair.dipole.z = charge * separation;
    const double lengths[2] = {16.0, 24.0};
    double corrected_errors[2] = {0.0, 0.0};
    for (int index = 0; index < 2; ++index)
    {
        const double length = lengths[index];
        parameters.length = length;
        const double periodic = periodic_gaussian_energy(length, sigma, charge, separation, true);
        const double correction = elecstate::pcc_0d_energy(pair, parameters);
        const double periodic_error = std::abs(periodic - isolated);
        const double corrected_error = std::abs(periodic + correction - isolated);
        EXPECT_LT(corrected_error, 0.05 * periodic_error);
        corrected_errors[index] = corrected_error;
    }
    const double ratio = lengths[1] / lengths[0];
    const double dipole_decay = std::pow(ratio, 3.0);
    EXPECT_GT(corrected_errors[0] / corrected_errors[1], dipole_decay);
}
