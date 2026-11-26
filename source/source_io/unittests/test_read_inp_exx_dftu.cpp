#include "source_io/module_parameter/parameter.h"
#include "source_io/module_parameter/read_input.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

// Match the existing serial INPUT target's friend helper definition.
class TestParameters
{
  public:
    static Input_para &input(Parameter &param) { return param.input; }
    static System_para &sys(Parameter &param) { return param.sys; }
};

static void expect_rejected(const ModuleIO::Input_Item &item, const Parameter &param, const std::string &message)
{
    testing::internal::CaptureStdout();
    EXPECT_EXIT(item.check_value(item, param), ::testing::ExitedWithCode(1), "");
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr(message));
}

TEST(DeviceExxInput, SelectionAndValidation)
{
    ModuleIO::ReadInput readinput(0);
    Parameter param;
    using InputList = std::vector<std::pair<std::string, ModuleIO::Input_Item>>;
    InputList input_lists = readinput.get_input_lists();
    InputList::iterator it = std::find_if(
        input_lists.begin(), input_lists.end(),
        [](const std::pair<std::string, ModuleIO::Input_Item> &item) { return item.first == "device_exx"; });
    ASSERT_NE(it, input_lists.end());
    EXPECT_EQ(param.inp.device_exx, "cpu");
    EXPECT_EQ(it->second.default_value, "cpu");

    Input_para &input = TestParameters::input(param);
    input.basis_type = "lcao";
    input.device = "gpu";
    input.symmetry = "1";
    input.exx_symmetry_realspace = true;
    it->second.check_value(it->second, param);
    it->second.get_final_value(it->second, param);
    EXPECT_EQ(it->second.final_value.str(), "cpu");

    it->second.str_values = {"gpu"};
    it->second.read_value(it->second, param);
    EXPECT_EQ(param.inp.device_exx, "gpu");
    input.basis_type = "pw";
    expect_rejected(it->second, param, "device_exx=gpu requires basis_type=lcao and device=gpu");
    input.basis_type = "lcao";
    input.device = "cpu";
    expect_rejected(it->second, param, "device_exx=gpu requires basis_type=lcao and device=gpu");
    input.device = "gpu";

#ifdef __GPU_RI
    for (const std::string symmetry : {"1", "0", "-1"})
    {
        for (const bool realspace : {false, true})
        {
            input.symmetry = symmetry;
            input.exx_symmetry_realspace = realspace;
            EXPECT_NO_THROW(it->second.check_value(it->second, param));
            EXPECT_EQ(input.device_exx, "gpu");
            EXPECT_EQ(input.symmetry, symmetry);
            EXPECT_EQ(input.exx_symmetry_realspace, realspace);
        }
    }
#else
    expect_rejected(it->second, param, "device_exx=gpu requires a build with ENABLE_LIBRI_GPU=ON");
#endif

    it->second.str_values = {"cpu"};
    it->second.read_value(it->second, param);
    EXPECT_EQ(param.inp.device_exx, "cpu");
    input.symmetry = "1";
    input.exx_symmetry_realspace = true;
    it->second.check_value(it->second, param);
}
