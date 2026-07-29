#include "source_lcao/module_ri/sternheimer_supercell_perturbation.h"

#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <vector>

namespace
{

using Complex = std::complex<double>;
using Channel = ModuleRI::SternheimerABFBlochGridChannel;

Channel make_channel(const int atom_index, const Complex value)
{
    Channel channel;
    channel.channel_index = atom_index;
    channel.atom_index = atom_index;
    channel.atom_local_index = atom_index;
    channel.type_index = 0;
    channel.angular_momentum = 0;
    channel.radial_index = 0;
    channel.magnetic_index = 0;
    channel.label = "H_s";
    channel.potential_r = {value, 2.0 * value};
    channel.max_abs = 2.0 * std::abs(value);
    return channel;
}

} // namespace

TEST(SternheimerSupercellPerturbation, ParsesCommensurateTranslationSum)
{
    const auto config = ModuleRI::parse_sternheimer_supercell_translation_sum(
        "repeats=2x2x2,primitive_q=0.5:0:0,atoms_per_primitive=2,basis_atom=0,channel_within_atom=0");
    EXPECT_EQ(config.repeats, (std::array<int, 3>{2, 2, 2}));
    EXPECT_EQ(config.primitive_qpoint, (std::array<double, 3>{0.5, 0.0, 0.0}));
    EXPECT_EQ(config.atoms_per_primitive, 2);
    EXPECT_EQ(config.basis_atom, 0);
    EXPECT_EQ(config.channel_within_atom, 0);

    EXPECT_THROW(ModuleRI::parse_sternheimer_supercell_translation_sum(
                     "repeats=2x2x2,primitive_q=0.25:0:0,atoms_per_primitive=2,basis_atom=0,channel_within_atom=0"),
                 std::invalid_argument);
}

TEST(SternheimerSupercellPerturbation, CombinesEquivalentAtomsWithBlochPhase)
{
    ModuleRI::SternheimerSupercellTranslationSum config;
    config.repeats = {2, 1, 1};
    config.primitive_qpoint = {0.5, 0.0, 0.0};
    config.atoms_per_primitive = 2;
    config.basis_atom = 0;
    config.channel_within_atom = 0;

    const std::vector<Channel> channels{make_channel(0, Complex(1.0, 0.0)),
                                        make_channel(1, Complex(20.0, 0.0)),
                                        make_channel(2, Complex(3.0, 0.0)),
                                        make_channel(3, Complex(40.0, 0.0))};
    const Channel combined = ModuleRI::combine_sternheimer_supercell_translation_channel(channels, config);
    const double inverse_sqrt_two = 1.0 / std::sqrt(2.0);

    EXPECT_EQ(combined.channel_index, 0);
    EXPECT_EQ(combined.atom_index, 0);
    ASSERT_EQ(combined.potential_r.size(), 2U);
    EXPECT_NEAR(combined.potential_r[0].real(), -2.0 * inverse_sqrt_two, 1.0e-14);
    EXPECT_NEAR(combined.potential_r[0].imag(), 0.0, 1.0e-14);
    EXPECT_NEAR(combined.potential_r[1].real(), -4.0 * inverse_sqrt_two, 1.0e-14);
    EXPECT_NEAR(combined.max_abs, 4.0 * inverse_sqrt_two, 1.0e-14);
}

TEST(SternheimerSupercellPerturbation, RejectsMissingEquivalentAtomChannel)
{
    ModuleRI::SternheimerSupercellTranslationSum config;
    config.repeats = {2, 1, 1};
    config.primitive_qpoint = {0.5, 0.0, 0.0};
    config.atoms_per_primitive = 2;
    config.basis_atom = 0;
    config.channel_within_atom = 0;

    const std::vector<Channel> channels{make_channel(0, Complex(1.0, 0.0)),
                                        make_channel(1, Complex(20.0, 0.0)),
                                        make_channel(3, Complex(40.0, 0.0))};
    EXPECT_THROW(ModuleRI::combine_sternheimer_supercell_translation_channel(channels, config),
                 std::invalid_argument);
}
