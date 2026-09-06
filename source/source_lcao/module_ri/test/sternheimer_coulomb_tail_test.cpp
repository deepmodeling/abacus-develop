#include "source_basis/module_ao/ORB_atomic_lm.h"
#include "source_lcao/module_ri/conv_coulomb_pot_k.h"
#include "source_lcao/module_ri/sternheimer_abfs_perturbation.h"

#include <cmath>
#include <gtest/gtest.h>
#include <limits>

namespace
{
auto sample_tail(const int l, const bool enabled, const bool periodic = false)
{
    ModuleRI::SternheimerRadialPerturbation radial;
    radial.angular_momentum = l;
    radial.radial_grid = {0.0, 0.5, 1.0};
    radial.radial_values = {0.0, 2.0, 3.0};
    radial.continue_coulomb_tail = enabled;
    radial.coulomb_tail_coefficient = 3.0;
    ModuleRI::SternheimerFDHamiltonian::Grid grid;
    grid.nx = 3;
    grid.ny = grid.nz = 1;
    grid.hx = grid.hy = grid.hz = 1.0;
    grid.periodic = periodic;
    return ModuleRI::sample_sternheimer_abf_grid_channels({{radial}},
                                                          {0},
                                                          {ModuleBase::Vector3<double>(0.0, 0.0, 0.0)},
                                                          grid);
}
} // namespace

TEST(SternheimerCoulombTail, PreservesInnerSamplesAndExtendsMultipoles)
{
    for (const int l: {0, 1, 2})
    {
        const auto old = sample_tail(l, false);
        const auto extended = sample_tail(l, true);
        ASSERT_EQ(old.size(), extended.size());
        double endpoint_norm = 0.0;
        for (std::size_t m = 0; m < old.size(); ++m)
        {
            EXPECT_EQ(extended[m].potential_r[0], old[m].potential_r[0]);
            EXPECT_EQ(extended[m].potential_r[1], old[m].potential_r[1]);
            EXPECT_EQ(old[m].potential_r[2], 0.0);
            EXPECT_NEAR(extended[m].potential_r[2], old[m].potential_r[1] / std::pow(2.0, l + 1), 1.e-14);
            endpoint_norm += std::abs(old[m].potential_r[1]);
        }
        EXPECT_GT(endpoint_norm, 0.1);
    }
}

TEST(SternheimerCoulombTail, RejectsPeriodicImageSum)
{
    EXPECT_THROW(sample_tail(0, true, true), std::invalid_argument);
    EXPECT_NO_THROW(sample_tail(0, false, true));
}

TEST(SternheimerCoulombTail, ZeroMomentChargeDoesNotAcquireAnExteriorPotential)
{
    // P_l=r^l (1-r^2)^2 (1-a*r^2) has exactly zero l-th moment on [0,1].
    for (const int l: {0, 1, 2})
    {
        const auto integral = [l](const int extra) {
            return 1.0 / (2 * l + 3 + extra) - 2.0 / (2 * l + 5 + extra) + 1.0 / (2 * l + 7 + extra);
        };
        constexpr int nr = 1001;
        std::vector<double> r(nr), rab(nr, 1.0 / (nr - 1)), p(nr);
        const double a = integral(0) / integral(2);
        for (int i = 0; i < nr; ++i)
        {
            r[i] = i * rab[i];
            p[i] = std::pow(r[i], l) * std::pow(1 - r[i] * r[i], 2) * (1 - a * r[i] * r[i]);
        }
        Numerical_Orbital_Lm charge;
        charge.set_orbital_info("X",
                                0,
                                l,
                                0,
                                nr,
                                rab.data(),
                                r.data(),
                                Numerical_Orbital_Lm::Psi_Type::Psi,
                                p.data(),
                                2001,
                                0.05,
                                0.001,
                                false,
                                true,
                                false);
        const auto potential = Conv_Coulomb_Pot_K::cal_orbs_ccp(
            charge,
            {{Conv_Coulomb_Pot_K::Coulomb_Type::Fock, {{{"alpha", "1"}, {"singularity_correction", "limits"}}}}},
            1.0);
        std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> orbitals(1);
        orbitals[0].resize(l + 1);
        orbitals[0][l].push_back(potential);
        auto radials = ModuleRI::make_sternheimer_radial_perturbations_from_orbitals(orbitals);
        ModuleRI::set_sternheimer_coulomb_tail_from_charge(radials[0][0], charge);
        ModuleRI::SternheimerFDHamiltonian::Grid grid;
        grid.nx = 2;
        grid.ny = grid.nz = 1;
        grid.hx = 2.0;
        grid.hy = grid.hz = 1.0;
        grid.periodic = false;
        const auto channels = ModuleRI::sample_sternheimer_abf_grid_channels(radials,
                                                                             {0},
                                                                             {ModuleBase::Vector3<double>(0, 0, 0)},
                                                                             grid);
        for (const auto& channel: channels)
        {
            EXPECT_NEAR(channel.potential_r[1], 0.0, 1.e-11) << "l=" << l;
        }
    }
}

TEST(SternheimerCoulombTail, ChargeMomentHasAnalyticNormalizationAndRejectsTruncation)
{
    for (const int l: {0, 1, 2, 3})
    {
        constexpr int nr = 1001;
        std::vector<double> r(nr), rab(nr, 0.001), p(nr);
        for (int i = 0; i < nr; ++i)
        {
            r[i] = i * 0.001;
            p[i] = std::pow(r[i], l) * std::pow(1 - r[i] * r[i], 2);
        }
        Numerical_Orbital_Lm charge;
        charge.set_orbital_info("X",
                                0,
                                l,
                                0,
                                nr,
                                rab.data(),
                                r.data(),
                                Numerical_Orbital_Lm::Psi_Type::Psi,
                                p.data(),
                                101,
                                0.05,
                                0.001,
                                false,
                                true,
                                false);
        ModuleRI::SternheimerRadialPerturbation potential;
        potential.type_index = 0;
        potential.angular_momentum = l;
        potential.radial_grid = {0.0, 0.5, 1.0};
        potential.radial_values = {7.0, 8.0, 9.0};
        const auto original_values = potential.radial_values;
        ModuleRI::set_sternheimer_coulomb_tail_from_charge(potential, charge);
        const double moment = 1.0 / (2 * l + 3) - 2.0 / (2 * l + 5) + 1.0 / (2 * l + 7);
        EXPECT_NEAR(potential.coulomb_tail_coefficient, 4 * std::acos(-1.0) * moment / (2 * l + 1), 1.e-11);
        EXPECT_EQ(potential.radial_values, original_values);
        // The charge, not an arbitrary last potential value, sets the tail.
        potential.radial_values.back() = 1.e4;
        ModuleRI::set_sternheimer_coulomb_tail_from_charge(potential, charge);
        EXPECT_NEAR(potential.coulomb_tail_coefficient, 4 * std::acos(-1.0) * moment / (2 * l + 1), 1.e-11);
        potential.radial_grid.back() = 0.9;
        EXPECT_THROW(ModuleRI::set_sternheimer_coulomb_tail_from_charge(potential, charge), std::invalid_argument);
    }
}

TEST(SternheimerCoulombTail, RejectsUnspecifiedMoment)
{
    ModuleRI::SternheimerRadialPerturbation radial;
    radial.radial_grid = {0, 0.5, 1};
    radial.radial_values = {0, 1, 2};
    radial.continue_coulomb_tail = true;
    ModuleRI::SternheimerFDHamiltonian::Grid grid;
    grid.nx = grid.ny = grid.nz = 1;
    grid.hx = grid.hy = grid.hz = 1;
    grid.periodic = false;
    EXPECT_THROW(
        ModuleRI::sample_sternheimer_abf_grid_channels({{radial}}, {0}, {ModuleBase::Vector3<double>(0, 0, 0)}, grid),
        std::invalid_argument);
}
