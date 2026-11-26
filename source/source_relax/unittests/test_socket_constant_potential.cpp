#include "source_base/constants.h"
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
    const SocketCP::Configuration config = SocketCP::parse_init(metadata);
    const elecstate::ElectronicControl control = SocketCP::make_control(config, 8.3125, 0);
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
}

TEST(SocketConstantPotentialTest, PotentialWindowsSupportEachCartesianAxis)
{
    for (int axis = 0; axis < 3; ++axis)
    {
        const SocketCP::Configuration config = SocketCP::parse_init(
            metadata + ",potential_average_axis:" + std::string(1, 'x' + axis) +
            ",potential_average_min_A:12.5,potential_average_max_A:15");
        EXPECT_EQ(axis, config.window.axis);
        EXPECT_DOUBLE_EQ(12.5, config.window.minimum_A);
        EXPECT_DOUBLE_EQ(15.0, config.window.maximum_A);
    }
}

TEST(SocketConstantPotentialTest, PotentialWindowRejectsInvalidOrMissingBounds)
{
    EXPECT_THROW(SocketCP::parse_init(metadata + ",potential_average_axis:xy"), std::runtime_error);
    EXPECT_THROW(
        SocketCP::parse_init(metadata + ",potential_average_axis:z,potential_average_min_A:1"), std::runtime_error);
    EXPECT_THROW(
        SocketCP::parse_init(
            metadata + ",potential_average_axis:z,potential_average_min_A:2,potential_average_max_A:1"),
        std::runtime_error);
    EXPECT_THROW(
        SocketCP::parse_init(
            metadata + ",potential_average_axis:z,potential_average_min_A:-1,potential_average_max_A:2"),
        std::runtime_error);
}

TEST(SocketConstantPotentialTest, RejectsRemovedStrideEvenForOrdinaryInit)
{
    EXPECT_THROW(SocketCP::parse_init("solvation_update_stride:1"), std::runtime_error);
    EXPECT_THROW(SocketCP::parse_init(metadata + ",solvation_update_stride:5"), std::runtime_error);
}

TEST(SocketConstantPotentialTest, RejectsNeElectrodeAndUnknownFlags)
{
    const SocketCP::Configuration config = SocketCP::parse_init(metadata);
    EXPECT_THROW(SocketCP::make_control(config, 8.0, 1), std::runtime_error);
    EXPECT_THROW(SocketCP::make_control(config, 8.0, 2), std::runtime_error);
    EXPECT_THROW(SocketCP::make_control(config, 8.0, 3), std::runtime_error);
    EXPECT_THROW(SocketCP::make_control(config, 8.0, 4), std::runtime_error);
    EXPECT_THROW(SocketCP::make_control(config, 8.0, -1), std::runtime_error);
}

TEST(SocketConstantPotentialTest, RejectsInvalidElectronRequests)
{
    const SocketCP::Configuration config = SocketCP::parse_init(metadata);
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
    EXPECT_EQ("{\"fermi_level_eV\":-0.25,\"nelect\":8.125}", SocketCP::serialize_observables(values));
    values.has_workfunction = true;
    values.reference_eV = 4.0;
    values.workfunction_eV = 4.25;
    EXPECT_EQ(
        "{\"fermi_level_eV\":-0.25,\"nelect\":8.125,\"workfunction_eV\":4.25,"
        "\"electrostatic_reference_eV\":4}",
        SocketCP::serialize_observables(values));
}

namespace
{
class DecimalComma : public std::numpunct<char>
{
    char do_decimal_point() const override { return ','; }
};
} // namespace

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

TEST(SocketElectronicFrameTest, BindingIsScopedAndRejectsNestedRequests)
{
    elecstate::ElectronicStep *slot = nullptr;
    elecstate::ElectronicStep first;
    elecstate::ElectronicStep second;
    EXPECT_EQ(nullptr, slot);
    {
        SocketCP::ScopedElectronicStep binding(&slot, &first);
        EXPECT_EQ(&first, slot);
        EXPECT_THROW(SocketCP::ScopedElectronicStep nested(&slot, &second), std::runtime_error);
        EXPECT_EQ(&first, slot);
    }
    EXPECT_EQ(nullptr, slot);
    try
    {
        SocketCP::ScopedElectronicStep binding(&slot, &second);
        throw std::runtime_error("evaluation failed");
    }
    catch (const std::runtime_error &)
    {
        EXPECT_EQ(nullptr, slot);
    }
}

TEST(SocketElectronicFrameTest, MissingSlotOnlyAcceptsOrdinaryEvaluation)
{
    elecstate::ElectronicStep step;
    EXPECT_NO_THROW(SocketCP::ScopedElectronicStep ordinary(nullptr, nullptr));
    EXPECT_THROW(SocketCP::ScopedElectronicStep missing(nullptr, &step), std::runtime_error);
}

TEST(SocketElectronicFrameTest, SeparateSolversDoNotShareRequestState)
{
    elecstate::ElectronicStep *first_slot = nullptr;
    elecstate::ElectronicStep *second_slot = nullptr;
    elecstate::ElectronicStep first;
    elecstate::ElectronicStep second;
    first.control.nelect = 8.25;
    second.control.nelect = 7.75;
    SocketCP::ScopedElectronicStep first_binding(&first_slot, &first);
    {
        SocketCP::ScopedElectronicStep second_binding(&second_slot, &second);
        first_slot->request_applied = true;
        EXPECT_FALSE(second_slot->request_applied);
        EXPECT_DOUBLE_EQ(8.25, first_slot->control.nelect);
        EXPECT_DOUBLE_EQ(7.75, second_slot->control.nelect);
    }
    EXPECT_EQ(nullptr, second_slot);
    EXPECT_EQ(&first, first_slot);
}

TEST(SocketElectronicFrameTest, CollectionConvertsPhysicalUnitsAndKeepsZeroFermi)
{
    elecstate::ElectronicStep step;
    step.request_applied = true;
    step.result_ready = true;
    step.result.nelect = 8.125;
    step.result.fermi_ry = 0.0;
    step.result.reference_ry = 0.25;
    const SocketCP::Configuration config = SocketCP::parse_init(
        "NELECT:8,neutral_electrons:8,potential_average_axis:y,potential_average_min_A:2,potential_average_max_A:3");
    step.window = config.window;
    const SocketCP::Observables values = SocketCP::collect_observables(config, step);
    EXPECT_DOUBLE_EQ(0.0, values.fermi_eV);
    EXPECT_DOUBLE_EQ(8.125, values.nelect);
    EXPECT_TRUE(values.has_workfunction);
    EXPECT_DOUBLE_EQ(0.25 * ModuleBase::Ry_to_eV, values.reference_eV);
    EXPECT_DOUBLE_EQ(values.reference_eV, values.workfunction_eV);
}

TEST(SocketElectronicFrameTest, CollectionRejectsMissingOrInvalidPhysicalObservables)
{
    elecstate::ElectronicStep step;
    const SocketCP::Configuration config = SocketCP::parse_init("NELECT:8,neutral_electrons:8");
    const SocketCP::Configuration ordinary;
    EXPECT_THROW(SocketCP::collect_observables(ordinary, step), std::runtime_error);
    EXPECT_THROW(SocketCP::collect_observables(config, step), std::runtime_error);
    step.result_ready = true;
    EXPECT_THROW(SocketCP::collect_observables(config, step), std::runtime_error);
    step.request_applied = true;
    step.result.nelect = 0.0;
    EXPECT_THROW(SocketCP::collect_observables(config, step), std::runtime_error);
    step.result.nelect = 8.0;
    step.result.fermi_ry = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(SocketCP::collect_observables(config, step), std::runtime_error);
}
