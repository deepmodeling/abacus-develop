#include "source_io/module_parameter/read_input.h"

#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
ModuleIO::Input_Item nac_input_item(const ModuleIO::ReadInput &reader, const std::string &name)
{
    for (const auto &item : reader.get_input_lists())
    {
        if (item.first == name)
        {
            return ModuleIO::Input_Item(item.second);
        }
    }
    throw std::runtime_error("Missing INPUT item: " + name);
}

void read_nac_item(
    const ModuleIO::ReadInput &reader,
    Parameter &parameter,
    const std::string &name,
    const std::vector<std::string> &values)
{
    ModuleIO::Input_Item item = nac_input_item(reader, name);
    item.str_values = values;
    item.read_value(item, parameter);
}

void prepare_nac_input(const ModuleIO::ReadInput &reader, Parameter &parameter)
{
    read_nac_item(
        reader,
        parameter,
        "esolver_type", {"ks-lr"});
    read_nac_item(
        reader,
        parameter,
        "basis_type", {"lcao"});
    read_nac_item(
        reader,
        parameter,
        "gamma_only", {"true"});
    read_nac_item(
        reader,
        parameter,
        "cal_force", {"true"});
    read_nac_item(
        reader,
        parameter,
        "lr_nstates", {"3"});
    read_nac_item(
        reader,
        parameter,
        "lr_nac_pairs",
        {"0",
            "1",
            "1",
            "3"});
}
} // namespace

TEST(LrNacInput, DisabledDefaultAndMetadata)
{
    const ModuleIO::ReadInput reader(0);
    Parameter parameter;
    EXPECT_TRUE(parameter.inp.lr_nac_pairs.empty());
    EXPECT_EQ(parameter.inp.lr_nac_spin, "singlet");
    const ModuleIO::Input_Item pairs = nac_input_item(reader, "lr_nac_pairs");
    EXPECT_EQ(pairs.type, "Vector of Integer");
    EXPECT_TRUE(pairs.default_value.empty());
    EXPECT_TRUE(pairs.unit.empty());
    EXPECT_NE(pairs.description.find("bohr^-1"), std::string::npos);
    EXPECT_NE(pairs.description.find("not full quadratic-response"), std::string::npos);
    EXPECT_NO_THROW(pairs.check_value(pairs, parameter));
    ModuleIO::Input_Item spin = nac_input_item(reader, "lr_nac_spin");
    EXPECT_EQ(spin.default_value, "singlet");
    EXPECT_FALSE(spin.domain.validate_explicit("lr_nac_spin", {"mixed"}).empty());
}

TEST(LrResponseInput, DefaultsAndExplicitPrecision)
{
    const ModuleIO::ReadInput reader(0);
    Parameter parameter;
    EXPECT_DOUBLE_EQ(parameter.inp.lr_grad_thr, 1e-6);
    EXPECT_EQ(parameter.inp.lr_grad_maxiter, 100);
    const ModuleIO::Input_Item tolerance = nac_input_item(reader, "lr_grad_thr");
    const ModuleIO::Input_Item iterations = nac_input_item(reader, "lr_grad_maxiter");
    EXPECT_EQ(tolerance.type, "Real");
    EXPECT_EQ(tolerance.unit, "Ry");
    EXPECT_EQ(tolerance.default_value, "1e-6");
    EXPECT_EQ(iterations.type, "Integer");
    EXPECT_EQ(iterations.default_value, "100");
    read_nac_item(
        reader,
        parameter,
        "lr_grad_thr", {"1e-10"});
    read_nac_item(
        reader,
        parameter,
        "lr_grad_maxiter", {"800"});
    EXPECT_DOUBLE_EQ(parameter.inp.lr_grad_thr, 1e-10);
    EXPECT_EQ(parameter.inp.lr_grad_maxiter, 800);
    EXPECT_NO_THROW(tolerance.check_value(tolerance, parameter));
    EXPECT_NO_THROW(iterations.check_value(iterations, parameter));
}

TEST(LrResponseInput, RejectsNonpositivePrecisionAndIterations)
{
    const ModuleIO::ReadInput reader(0);
    const ModuleIO::Input_Item tolerance = nac_input_item(reader, "lr_grad_thr");
    const ModuleIO::Input_Item iterations = nac_input_item(reader, "lr_grad_maxiter");
    for (const std::string value : {"0",
             "-1e-6",
             "nan",
             "inf"})
    {
        Parameter parameter;
        EXPECT_EXIT(
            {
                read_nac_item(
                    reader,
                    parameter,
                    "lr_grad_thr", {value});
                tolerance.check_value(tolerance, parameter);
            },
            ::testing::ExitedWithCode(1), "");
    }
    for (const std::string value : {"0", "-1"})
    {
        Parameter parameter;
        read_nac_item(
            reader,
            parameter,
            "lr_grad_maxiter", {value});
        EXPECT_EXIT(iterations.check_value(iterations, parameter), ::testing::ExitedWithCode(1), "");
    }
}

TEST(LrNacInput, ReadsSelectedPairsAndClosedShellSpin)
{
    const ModuleIO::ReadInput reader(0);
    Parameter parameter;
    prepare_nac_input(reader, parameter);
    const ModuleIO::Input_Item pairs = nac_input_item(reader, "lr_nac_pairs");
    EXPECT_EQ(
        parameter.inp.lr_nac_pairs, (std::vector<int>{0,
                                        1,
                                        1,
                                        3}));
    EXPECT_NO_THROW(pairs.check_value(pairs, parameter));
    read_nac_item(
        reader,
        parameter,
        "lr_nac_spin", {"triplet"});
    EXPECT_EXIT(pairs.check_value(pairs, parameter), ::testing::ExitedWithCode(1), "");
    read_nac_item(
        reader,
        parameter,
        "nspin", {"2"});
    EXPECT_NO_THROW(pairs.check_value(pairs, parameter));
}

TEST(LrNacInput, RejectsHubbardReferenceWhenNacRequested)
{
    const ModuleIO::ReadInput reader(0);
    const ModuleIO::Input_Item pairs = nac_input_item(reader, "lr_nac_pairs");
    EXPECT_NE(pairs.get_availability().find("dft_plus_u==0"), std::string::npos);
    EXPECT_NE(pairs.description.find("Hubbard response"), std::string::npos);
    for (const std::string mode : {"1", "2"})
    {
        Parameter parameter;
        prepare_nac_input(reader, parameter);
        EXPECT_NO_THROW(pairs.check_value(pairs, parameter));
        read_nac_item(
            reader,
            parameter,
            "dft_plus_u", {mode});
        EXPECT_EXIT(pairs.check_value(pairs, parameter), ::testing::ExitedWithCode(1), "");
    }
}

TEST(LrNacInput, DisabledNacLeavesHubbardInputUntouched)
{
    const ModuleIO::ReadInput reader(0);
    Parameter parameter;
    read_nac_item(
        reader,
        parameter,
        "dft_plus_u", {"2"});
    const ModuleIO::Input_Item pairs = nac_input_item(reader, "lr_nac_pairs");
    ASSERT_TRUE(parameter.inp.lr_nac_pairs.empty());
    EXPECT_NO_THROW(pairs.check_value(pairs, parameter));
    EXPECT_EQ(parameter.inp.dft_plus_u, 2);
}

TEST(LrNacInput, RejectsMalformedIntegerTokens)
{
    const ModuleIO::ReadInput reader(0);
    Parameter parameter;
    for (const std::string token : {"1.5",
             "1junk",
             "1e0",
             "999999999999999999999999999"})
    {
        EXPECT_EXIT(
            read_nac_item(
                reader,
                parameter,
                "lr_nac_pairs", {"0", token}),
            ::testing::ExitedWithCode(1), "");
    }
}

TEST(LrNacInput, RejectsInvalidPairLists)
{
    const ModuleIO::ReadInput reader(0);
    Parameter parameter;
    prepare_nac_input(reader, parameter);
    const ModuleIO::Input_Item pairs = nac_input_item(reader, "lr_nac_pairs");
    const std::vector<std::vector<std::string>> invalid = {{"0"}, {"-1", "1"}, {"1", "1"}, {"2", "1"}, {"0", "4"},
        {"0",
            "1",
            "0",
            "1"}};
    for (const auto &values : invalid)
    {
        read_nac_item(
            reader,
            parameter,
            "lr_nac_pairs",
            values);
        EXPECT_EXIT(pairs.check_value(pairs, parameter), ::testing::ExitedWithCode(1), "");
    }
}

TEST(LrNacInput, RejectsUnsupportedElectronicModels)
{
    const ModuleIO::ReadInput reader(0);
    const ModuleIO::Input_Item pairs = nac_input_item(reader, "lr_nac_pairs");
    const std::vector<std::pair<std::string, std::string>> invalid = {{"esolver_type", "ksdft"}, {"basis_type", "pw"},
        {"gamma_only", "false"}, {"cal_force", "false"}, {"kpar", "2"}, {"lr_unrestricted", "true"}, {"nspin", "4"},
        {"lr_solver", "spectrum"}, {"calculation", "cell-relax"}, {"lr_nac_spin", "mixed"}};
    for (const auto &setting : invalid)
    {
        Parameter parameter;
        prepare_nac_input(reader, parameter);
        read_nac_item(
            reader,
            parameter,
            setting.first, {setting.second});
        EXPECT_EXIT(pairs.check_value(pairs, parameter), ::testing::ExitedWithCode(1), "");
    }
}

TEST(LrNacInput, RejectsRequestedLcaoPoolsAfterInputReset)
{
    const ModuleIO::ReadInput reader(0);
    Parameter parameter;
    prepare_nac_input(reader, parameter);
    read_nac_item(
        reader,
        parameter,
        "kpar", {"2"});
    const ModuleIO::Input_Item pools = nac_input_item(reader, "kpar");
    pools.reset_value(pools, parameter);
#ifdef __LCAO
    ASSERT_EQ(parameter.inp.kpar, 1);
    ASSERT_EQ(parameter.globalv.kpar_lcao, 2);
#else
    ASSERT_EQ(parameter.inp.kpar, 2);
    ASSERT_EQ(parameter.globalv.kpar_lcao, 1);
#endif
    const ModuleIO::Input_Item pairs = nac_input_item(reader, "lr_nac_pairs");
    EXPECT_EXIT(pairs.check_value(pairs, parameter), ::testing::ExitedWithCode(1), "");
}

TEST(LrNacInput, ExplicitUnrestrictedChannelRequiresTwoSpins)
{
    const ModuleIO::ReadInput reader(0);
    Parameter parameter;
    prepare_nac_input(reader, parameter);
    const ModuleIO::Input_Item pairs = nac_input_item(reader, "lr_nac_pairs");
    const ModuleIO::Input_Item spin = nac_input_item(reader, "lr_nac_spin");
    EXPECT_TRUE(spin.domain.validate_explicit("lr_nac_spin", {"updown"}).empty());
    read_nac_item(
        reader,
        parameter,
        "lr_nac_spin", {"updown"});
    EXPECT_EXIT(pairs.check_value(pairs, parameter), ::testing::ExitedWithCode(1), "");
    read_nac_item(
        reader,
        parameter,
        "nspin", {"2"});
    read_nac_item(
        reader,
        parameter,
        "lr_unrestricted", {"true"});
    EXPECT_NO_THROW(pairs.check_value(pairs, parameter));
}

TEST(LrNacInput, RejectsHartreeDisabledOnlyWhenNacIsRequested)
{
    const ModuleIO::ReadInput reader(0);
    Parameter parameter;
    prepare_nac_input(reader, parameter);
    read_nac_item(reader, parameter, "vh_in_h", {"0"});
    const ModuleIO::Input_Item pairs = nac_input_item(reader, "lr_nac_pairs");
    EXPECT_NE(pairs.get_availability().find("vh_in_h==1"), std::string::npos);
    EXPECT_EXIT(pairs.check_value(pairs, parameter), ::testing::ExitedWithCode(1), "");
    read_nac_item(reader, parameter, "lr_nac_pairs", {});
    EXPECT_NO_THROW(pairs.check_value(pairs, parameter));
    EXPECT_EQ(parameter.inp.vh_in_h, 0);
}
