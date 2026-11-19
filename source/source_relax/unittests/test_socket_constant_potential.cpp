#include "source_relax/socket_constant_potential.h"

#include "gtest/gtest.h"

#include <limits>
#include <locale>
#include <stdexcept>
#include <string>

namespace
{
const std::string metadata = "NELECT:8.25,neutral_electrons:8";
}

TEST(SocketConstantPotentialTest, OrdinaryInitDoesNotEnableElectronControl)
{
    EXPECT_FALSE(SocketCP::parse_init("").enabled);
    EXPECT_FALSE(SocketCP::parse_init("label:ordinary,unused:123").enabled);
    EXPECT_THROW(SocketCP::make_control(SocketCP::Configuration(), 8.0, 1), std::runtime_error);
}

TEST(SocketConstantPotentialTest, FractionalElectronRequestIsPreserved)
{
    const auto config = SocketCP::parse_init(metadata);
    const auto control = SocketCP::make_control(config, 8.3125, 0);
    EXPECT_TRUE(config.enabled);
    EXPECT_DOUBLE_EQ(8.3125, control.nelect);
}

TEST(SocketConstantPotentialTest, ConstantPotentialInitRequiresCompleteMetadata)
{
    EXPECT_THROW(SocketCP::parse_init("NELECT:8"), std::runtime_error);
    EXPECT_THROW(SocketCP::parse_init("neutral_electrons:8"), std::runtime_error);
    EXPECT_TRUE(SocketCP::parse_init("NELECT:8,neutral_electrons:8").enabled);
    EXPECT_THROW(SocketCP::parse_init(metadata + ",NELECT:9"), std::runtime_error);
    EXPECT_THROW(SocketCP::parse_init(metadata + ",neutral_electrons:8"), std::runtime_error);
}

TEST(SocketConstantPotentialTest, InitRejectsInvalidNumericalValues)
{
    EXPECT_THROW(SocketCP::parse_init("NELECT:8junk,neutral_electrons:8"), std::runtime_error);
    EXPECT_THROW(SocketCP::parse_init("NELECT:0,neutral_electrons:8"), std::runtime_error);
    EXPECT_THROW(SocketCP::parse_init("NELECT:8,neutral_electrons:-1"), std::runtime_error);
    EXPECT_THROW(SocketCP::parse_init("NELECT:8,neutral_electrons:8,solvation_update_stride:0"), std::runtime_error);
    EXPECT_THROW(SocketCP::parse_init("NELECT:8,neutral_electrons:8,solvation_update_stride:1.5"), std::runtime_error);
    EXPECT_THROW(SocketCP::parse_init("NELECT:8,neutral_electrons:8,solvation_update_stride:2147483648"), std::runtime_error);
}

TEST(SocketConstantPotentialTest, PotentialWindowsSupportEachCartesianAxis)
{
    for (int axis = 0; axis < 3; ++axis)
    {
        const auto config = SocketCP::parse_init(metadata + ",potential_average_axis:"
                                                + std::string(1, 'x' + axis)
                                                + ",potential_average_min_A:12.5,potential_average_max_A:15");
        EXPECT_EQ(axis, config.window.axis);
        EXPECT_DOUBLE_EQ(12.5, config.window.minimum_A);
        EXPECT_DOUBLE_EQ(15.0, config.window.maximum_A);
    }
}

TEST(SocketConstantPotentialTest, PotentialWindowRejectsInvalidOrMissingBounds)
{
    EXPECT_THROW(SocketCP::parse_init(metadata + ",potential_average_axis:xy"), std::runtime_error);
    EXPECT_THROW(SocketCP::parse_init(metadata + ",potential_average_axis:z,potential_average_min_A:1"), std::runtime_error);
    EXPECT_THROW(SocketCP::parse_init(metadata + ",potential_average_axis:z,potential_average_min_A:2,potential_average_max_A:1"), std::runtime_error);
    EXPECT_THROW(SocketCP::parse_init(metadata + ",potential_average_axis:z,potential_average_min_A:-1,potential_average_max_A:2"), std::runtime_error);
}

TEST(SocketConstantPotentialTest, RejectsRemovedStrideEvenForOrdinaryInit)
{
    EXPECT_THROW(SocketCP::parse_init("solvation_update_stride:1"), std::runtime_error);
    EXPECT_THROW(SocketCP::parse_init(metadata + ",solvation_update_stride:5"), std::runtime_error);
}

TEST(SocketConstantPotentialTest, RejectsNeElectrodeAndUnknownFlags)
{
    const auto config = SocketCP::parse_init(metadata);
    EXPECT_THROW(SocketCP::make_control(config, 8.0, 1), std::runtime_error);
    EXPECT_THROW(SocketCP::make_control(config, 8.0, 2), std::runtime_error);
    EXPECT_THROW(SocketCP::make_control(config, 8.0, 3), std::runtime_error);
    EXPECT_THROW(SocketCP::make_control(config, 8.0, 4), std::runtime_error);
    EXPECT_THROW(SocketCP::make_control(config, 8.0, -1), std::runtime_error);
}

TEST(SocketConstantPotentialTest, RejectsInvalidElectronRequests)
{
    const auto config = SocketCP::parse_init(metadata);
    EXPECT_THROW(SocketCP::make_control(config, 0.0, 0), std::runtime_error);
    EXPECT_THROW(SocketCP::make_control(config, -1.0, 0), std::runtime_error);
    EXPECT_THROW(SocketCP::make_control(config, std::numeric_limits<double>::quiet_NaN(), 0), std::runtime_error);
    EXPECT_THROW(SocketCP::make_control(config, std::numeric_limits<double>::infinity(), 0), std::runtime_error);
}

TEST(SocketConstantPotentialTest, SerializesTheExistingElectronicExtrasKeys)
{
    SocketCP::Observables values;
    values.nelect = 8.125;
    values.fermi_eV = -0.25;
    EXPECT_EQ("{\"fermi_level_eV\":-0.25,\"nelect\":8.125}",
              SocketCP::serialize_observables(values));
    values.has_workfunction = true;
    values.reference_eV = 4.0;
    values.workfunction_eV = 4.25;
    EXPECT_EQ("{\"fermi_level_eV\":-0.25,\"nelect\":8.125,\"workfunction_eV\":4.25,"
              "\"electrostatic_reference_eV\":4}", SocketCP::serialize_observables(values));
}

namespace
{
class DecimalComma : public std::numpunct<char>
{
    char do_decimal_point() const override { return ','; }
};
}

TEST(SocketConstantPotentialTest, JsonNumbersIgnoreTheGlobalDecimalSeparator)
{
    SocketCP::Observables values;
    values.nelect = 8.125;
    values.fermi_eV = -0.25;
    const std::locale previous = std::locale();
    std::locale::global(std::locale(previous, new DecimalComma));
    const std::string json = SocketCP::serialize_observables(values);
    std::locale::global(previous);
    EXPECT_EQ("{\"fermi_level_eV\":-0.25,\"nelect\":8.125}", json);
}
