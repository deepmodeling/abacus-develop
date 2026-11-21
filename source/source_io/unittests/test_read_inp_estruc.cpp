#include <algorithm>
#include <string>
#include <vector>
#include <utility>

#include "gtest/gtest.h"
#include "source_io/module_parameter/read_input.h"

class TestParameters
{
  public:
    static Input_para& input(Parameter& param) { return param.input; }
};

class LobpcgInputTest : public testing::Test
{
  protected:
    using Items = std::vector<std::pair<std::string, ModuleIO::Input_Item>>;
    Items::iterator find_label(const std::string& label, Items& items)
    {
        return std::find_if(items.begin(), items.end(),
            [&label](const Items::value_type& item) { return item.first == label; });
    }
    template <typename ReadValue>
    void check_option(const std::string& label, const std::string& default_value, const std::string& type,
                      const std::vector<std::string>& valid, const std::vector<std::string>& invalid, ReadValue value)
    {
        ModuleIO::ReadInput reader(0);
        Parameter param;
        auto items = reader.get_input_lists();
        auto item = find_label(label, items);
        ASSERT_NE(item, items.end());
        auto& option = item->second;
        EXPECT_EQ(value(param), default_value);
        EXPECT_EQ(option.default_value, default_value);
        EXPECT_EQ(option.type, type);
        EXPECT_EQ(option.get_availability(), "basis_type==pw and ks_solver==lobpcg");
        for (const auto& token : valid)
        {
            SCOPED_TRACE(token);
            option.str_values = {token};
            EXPECT_TRUE(option.domain.validate_explicit(label, option.str_values).empty());
            option.read_value(option, param);
            EXPECT_EQ(value(param), token);
            option.final_value.str("");
            option.get_final_value(option, param);
            EXPECT_EQ(option.final_value.str(), token);
        }
        for (const auto& token : invalid)
        {
            const std::vector<std::string> tokens = {token};
            EXPECT_FALSE(option.domain.validate_explicit(label, tokens).empty());
        }
    }
    void check_item(const ModuleIO::Input_Item& item, const Parameter& param)
    {
        ModuleIO::ReadInput schema(-1);
        ModuleIO::InputValues values;
        for (const auto& entry : schema.get_input_lists())
        {
            if (entry.second.domain.bound())
            {
                values[entry.first] = entry.second.domain.value(param);
            }
        }
        EXPECT_TRUE(item.domain.validate_resolved(item.label, values).empty());
        EXPECT_TRUE(item.domain.validate_requirements(item.label, values).empty());
        if (item.check_value)
        {
            item.check_value(item, param);
        }
    }
};

TEST_F(LobpcgInputTest, LobpcgExtraBands)
{
    check_option("lobpcg_extra_bands", "0", "Integer", {"0", "5", "10"}, {"-1", "1.5", "invalid", "2147483648"},
                 [](const Parameter& p) { return std::to_string(p.inp.lobpcg_extra_bands); });
}

TEST_F(LobpcgInputTest, LobpcgLocking)
{
    check_option("lobpcg_locking", "hard", "String", {"hard", "soft"}, {"auto", "none", "1"},
                 [](const Parameter& p) { return p.inp.lobpcg_locking; });
}

TEST_F(LobpcgInputTest, LobpcgDeviceSupport)
{
    ModuleIO::ReadInput reader(0);
    Parameter param;
    auto items = reader.get_input_lists();
    auto item = find_label("ks_solver", items);
    ASSERT_NE(item, items.end());
    TestParameters::input(param).basis_type = "pw";
    TestParameters::input(param).ks_solver = "lobpcg";
    for (const std::string device : {"cpu", "gpu"})
    {
        TestParameters::input(param).device = device;
        check_item(item->second, param);
    }
}
