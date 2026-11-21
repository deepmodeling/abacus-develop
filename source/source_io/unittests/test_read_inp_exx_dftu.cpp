#include "source_io/module_parameter/read_input.h"
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace
{
ModuleIO::Input_Item item_named(const ModuleIO::ReadInput& reader, const std::string& name)
{
    const auto& items = reader.get_input_lists();
    const auto found = std::find_if(items.begin(), items.end(),
                                  [&name](const std::pair<std::string, ModuleIO::Input_Item>& entry) {
                                      return entry.first == name;
                                  });
    if (found == items.end())
    {
        throw std::runtime_error("Missing INPUT registration: " + name);
    }
    return found->second;
}

void read_value(const ModuleIO::ReadInput& reader, Parameter& parameter,
                const std::string& name, const std::string& value)
{
    auto item = item_named(reader, name);
    item.str_values = {value};
    item.read_value(item, parameter);
}
} // namespace

TEST(ReadInpExxDftu, RegisteredDefaultsAndScalarReader)
{
    ModuleIO::ReadInput reader(0);
    Parameter parameter;
    auto dc = item_named(reader, "dftu_dc");
    auto alpha = item_named(reader, "dftu_alpha");
    EXPECT_EQ(parameter.inp.dftu_dc, "fll");
    ASSERT_TRUE(alpha.check_value);
    alpha.check_value(alpha, parameter); // Default FLL is not restricted.
    dc.str_values = {"mixed"};
    dc.read_value(dc, parameter);
    alpha.str_values = {"0.25"};
    alpha.read_value(alpha, parameter);
    EXPECT_EQ(parameter.inp.dftu_dc, "mixed");
    EXPECT_DOUBLE_EQ(parameter.inp.dftu_alpha, 0.25);
    alpha.get_final_value(alpha, parameter);
    EXPECT_EQ(alpha.final_value.str(), "0.25");
    alpha.str_values = {"0.25garbage"};
    EXPECT_THROW(alpha.read_value(alpha, parameter), std::exception);
}

TEST(ReadInpExxDftu, RegisteredCrossParameterChecks)
{
    ModuleIO::ReadInput reader(0);
    Parameter parameter;
    read_value(reader, parameter, "dft_plus_u", "1");
    read_value(reader, parameter, "basis_type", "lcao");
    read_value(reader, parameter, "ks_solver", "genelpa");
    read_value(reader, parameter, "dft_functional", "PBE");
    read_value(reader, parameter, "symmetry", "0");
    read_value(reader, parameter, "dftu_dc", "mixed");
    for (const std::string name : {"exx_fock_alpha", "exx_erfc_alpha"})
    {
        auto item = item_named(reader, name);
        item.reset_value(item, parameter);
    }
    auto alpha = item_named(reader, "dftu_alpha");
    EXPECT_EXIT(alpha.check_value(alpha, parameter), ::testing::ExitedWithCode(1), "");
    alpha.str_values = {"0"};
    alpha.read_value(alpha, parameter);
    alpha.check_value(alpha, parameter);
    read_value(reader, parameter, "cal_force", "1");
    EXPECT_EXIT(alpha.check_value(alpha, parameter), ::testing::ExitedWithCode(1), "");
    // An explicit FLL request does not inherit the AMF/mixed force restriction.
    read_value(reader, parameter, "dftu_dc", "fll");
    alpha.str_values.clear();
    alpha.check_value(alpha, parameter);
}

TEST(ReadInpExxDftu, AlphaWriterPreservesDoubleAfterCopyAndRepeatedUse)
{
    ModuleIO::ReadInput reader(0);
    Parameter parameter;
    auto registered = item_named(reader, "dftu_alpha");
    ModuleIO::Input_Item copied(registered);
    const double below_one = std::nextafter(1.0, 0.0);
    const double values[] = {0.0, 1e-8, 0.12345678901234566, 0.5, 0.99999999, below_one, 1.0};
    for (const double value : values)
    {
        SCOPED_TRACE(value);
        for (int repeat = 0; repeat < 2; ++repeat)
        {
            std::ostringstream input_text;
            input_text << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
            copied.str_values = {input_text.str()};
            copied.read_value(copied, parameter);
            ASSERT_EQ(parameter.inp.dftu_alpha, value);
            copied.final_value.str("");
            copied.final_value.clear();
            copied.final_value.precision(6);
            copied.get_final_value(copied, parameter);
            const std::string emitted = copied.final_value.str();
            const double restored = std::stod(emitted);
            EXPECT_EQ(restored, value);
            if (value > 0.0 && value < 1.0)
            {
                EXPECT_GT(restored, 0.0);
                EXPECT_LT(restored, 1.0);
            }
        }
    }
}
