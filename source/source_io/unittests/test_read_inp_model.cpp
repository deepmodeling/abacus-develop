// read_inp_model.cpp: the rVV10 INPUT items. vdw_method accepts rvv10 and declares the unsupported DFPT and
// LR-TDDFT force calculations; vdw_rvv10_b follows dft_functional unless it is set; vdw_rvv10_b and
// vdw_rvv10_c have domains. Complete INPUT files are read for valid input; invalid values are checked on the
// declared domains, because a rejected INPUT ends the process in an MPI build.
#include "source_io/module_parameter/parameter.h"
#include "source_io/module_parameter/read_input.h"

#include "gtest/gtest.h"
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <vector>

/// Parameter grants this class write access to its INPUT values (see parameter.h).
class TestParameters
{
  public:
    static Input_para& input(Parameter& param)
    {
        return param.input;
    }
};

namespace
{

/// The INPUT items by label.
std::map<std::string, ModuleIO::Input_Item> input_items()
{
    ModuleIO::ReadInput reader(0);
    std::map<std::string, ModuleIO::Input_Item> items;
    for (const std::pair<std::string, ModuleIO::Input_Item>& entry: reader.get_input_lists())
    {
        items.emplace(entry.first, entry.second);
    }
    return items;
}

/// Domain and requirement errors of item for the values in param, as read_parameters reports them.
std::string domain_errors(const ModuleIO::Input_Item& item, const Parameter& param)
{
    ModuleIO::ReadInput schema(0);
    ModuleIO::InputValues values;
    for (const std::pair<std::string, ModuleIO::Input_Item>& entry: schema.get_input_lists())
    {
        if (entry.second.domain.bound())
        {
            values[entry.first] = entry.second.domain.value(param);
        }
    }
    std::string errors = item.domain.validate_resolved(item.label, values);
    errors += item.domain.validate_requirements(item.label, values);
    return errors;
}

/// Writes parameters to an INPUT file and reads it.
void read_input(const std::string& parameters, Parameter& param)
{
    const std::string filename = testing::TempDir() + "rvv10_INPUT";
    {
        std::ofstream input(filename.c_str());
        input << "INPUT_PARAMETERS\n" << parameters;
    }
    // the test runs on one process
    ModuleIO::ReadInput reader(0);
    reader.check_ntype_flag = false;
    reader.read_parameters(param, filename);
    std::remove(filename.c_str());
}

} // namespace

TEST(ReadInpModelRvv10, ReadsTheParametersAndPairsBWithTheFunctional)
{
    Parameter param;
    read_input("vdw_method rvv10\ndft_functional scan\nvdw_rvv10_b 15.7\ncalculation cell-relax\n", param);
    EXPECT_EQ(param.inp.vdw_method, "rvv10");
    EXPECT_DOUBLE_EQ(param.inp.vdw_rvv10_b, 15.7);
    EXPECT_DOUBLE_EQ(param.inp.vdw_rvv10_c, 0.0093);
    // b follows dft_functional unless it is set: SCAN 15.7, r2SCAN 11.95, otherwise the rPW86-PBE value 6.3
    const std::vector<std::pair<std::string, double>> b_cases
        = {{"vdw_method rvv10\ndft_functional scan\n", 15.7},
           {"vdw_method rvv10\ndft_functional MGGA_X_SCAN+MGGA_C_SCAN\n", 15.7},
           {"vdw_method rvv10\ndft_functional mgga_x_r2scan+mgga_c_r2scan\n", 11.95},
           {"vdw_method rvv10\ndft_functional scan\nvdw_rvv10_b 9.0\n", 9.0},
           {"vdw_method rvv10\ndft_functional pbe\n", 6.3},
           {"vdw_method rvv10\n", 6.3}};
    for (const std::pair<std::string, double>& b_case: b_cases)
    {
        Parameter param_b;
        read_input(b_case.first, param_b);
        EXPECT_DOUBLE_EQ(param_b.inp.vdw_rvv10_b, b_case.second) << b_case.first;
    }
}

TEST(ReadInpModelRvv10, DeclaresTheUnsupportedCalculationsAndTheDomains)
{
    const std::map<std::string, ModuleIO::Input_Item> items = input_items();
    const ModuleIO::Input_Item& method = items.at("vdw_method");
    EXPECT_TRUE(method.domain.validate_explicit("vdw_method", {"rvv10"}).empty());
    EXPECT_FALSE(method.domain.validate_explicit("vdw_method", {"RVV10"}).empty());

    Parameter param;
    Input_para& input = TestParameters::input(param);
    input.vdw_method = "rvv10";
    EXPECT_TRUE(domain_errors(method, param).empty());
    input.esolver_type = "dfpt";
    EXPECT_NE(domain_errors(method, param).find("Requires esolver_type!=dfpt"), std::string::npos);
    // LR-TDDFT forces; relax and cell-relax set cal_force, which the requirement sees
    for (const std::string esolver: {"lr", "ks-lr"})
    {
        input.esolver_type = esolver;
        input.cal_force = true;
        EXPECT_FALSE(domain_errors(method, param).empty()) << esolver;
        input.cal_force = false;
        EXPECT_TRUE(domain_errors(method, param).empty()) << esolver;
    }

    const ModuleIO::Input_Item& b = items.at("vdw_rvv10_b");
    const ModuleIO::Input_Item& c = items.at("vdw_rvv10_c");
    Parameter rvv10;
    Input_para& rvv10_input = TestParameters::input(rvv10);
    rvv10_input.vdw_method = "rvv10";
    EXPECT_TRUE(domain_errors(b, rvv10).empty());
    EXPECT_TRUE(domain_errors(c, rvv10).empty());
    rvv10_input.vdw_rvv10_b = 0.0;
    rvv10_input.vdw_rvv10_c = -0.1;
    EXPECT_FALSE(domain_errors(b, rvv10).empty());
    EXPECT_FALSE(domain_errors(c, rvv10).empty());
}
