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
    ModuleIO::ReadInput reader{0};
    Parameter param;
    Items items;
    void SetUp() override { items = Items(reader.get_input_lists()); }

};

TEST_F(LobpcgInputTest, PwDiagNmaxMustBePositive)
{
    Items::iterator item = find_label("pw_diag_nmax", items);
    ASSERT_NE(item, items.end());
    ModuleIO::Input_Item& option = item->second;
    EXPECT_EQ(option.default_value, "50");
    EXPECT_EQ(option.type, "Integer");
    EXPECT_TRUE(option.domain.validate_explicit(option.label, {"37"}).empty());
    EXPECT_FALSE(option.domain.validate_explicit(option.label, {"0"}).empty());
    EXPECT_FALSE(option.domain.validate_explicit(option.label, {"-1"}).empty());
}

TEST_F(LobpcgInputTest, LobpcgLocking)
{
    Items::iterator found = find_label("lobpcg_locking", items);
    ASSERT_NE(found, items.end());
    ModuleIO::Input_Item& item = found->second;
    EXPECT_EQ(param.inp.lobpcg_locking, "hard");
    EXPECT_EQ(item.default_value, "hard");
    EXPECT_TRUE(item.domain.validate_explicit(item.label, {"hard"}).empty());
    EXPECT_TRUE(item.domain.validate_explicit(item.label, {"soft"}).empty());
    EXPECT_FALSE(item.domain.validate_explicit(item.label, {"invalid"}).empty());
    item.str_values = {"soft"};
    item.read_value(item, param);
    EXPECT_EQ(param.inp.lobpcg_locking, "soft");
    item.get_final_value(item, param);
    EXPECT_EQ(item.final_value.str(), "soft");
}

TEST_F(LobpcgInputTest, LobpcgDeviceSupport)
{
    Items::iterator item = find_label("ks_solver", items);
    ASSERT_NE(item, items.end());
    TestParameters::input(param).basis_type = "pw";
    TestParameters::input(param).ks_solver = "lobpcg";
    ModuleIO::InputValues values;
    values["basis_type"] = ModuleIO::InputValue(param.inp.basis_type);
    values["ks_solver"] = ModuleIO::InputValue(param.inp.ks_solver);
    for (const std::string device : {"cpu", "gpu"})
    {
        TestParameters::input(param).device = device;
        values["device"] = ModuleIO::InputValue(device);
        EXPECT_TRUE(item->second.domain.validate_requirements("ks_solver", values).empty());
        item->second.check_value(item->second, param);
    }
}

TEST_F(LobpcgInputTest, LobpcgSubspaceBackends)
{
    Items::iterator item = find_label("diag_subspace", items);
    ASSERT_NE(item, items.end());
    ModuleIO::InputValues values;
    const std::string lobpcg = "lobpcg";
    values["ks_solver"] = ModuleIO::InputValue(lobpcg);
    for (const int backend : {0, 2})
    {
        values["diag_subspace"] = ModuleIO::InputValue(backend);
        EXPECT_TRUE(item->second.domain.validate_requirements("diag_subspace", values).empty());
    }
    values["diag_subspace"] = ModuleIO::InputValue(1);
    EXPECT_FALSE(item->second.domain.validate_requirements("diag_subspace", values).empty());
    const std::string davidson = "dav_subspace";
    values["ks_solver"] = ModuleIO::InputValue(davidson);
    EXPECT_TRUE(item->second.domain.validate_requirements("diag_subspace", values).empty());
}
