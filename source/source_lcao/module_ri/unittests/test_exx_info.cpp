#include "source_hamilt/module_xc/exx_info.h"
#include "source_io/module_parameter/input_parameter.h"

#include <gtest/gtest.h>
#include <stdexcept>

namespace
{
Input_para hse_input()
{
    Input_para input;
    input.basis_type = "lcao";
    input.esolver_type = "tddft";
    input.dft_functional = "hse";
    input.exx_fock_alpha.clear();
    input.exx_fock_lambda.clear();
    input.exx_erfc_alpha = {"0.25"};
    input.exx_erfc_omega = {"0.11"};
    input.exx_ccp_rmesh_times = "1.5";
    input.exx_real_number = "1";
    input.td_stype = 2;
    return input;
}
}

TEST(ExxInfo, RejectsRealHybridGaugeBeforeRiSetup)
{
    const auto input = hse_input();
    Exx_Info info;
    try
    {
        init_exx_info(info, input);
        FAIL() << "Real RI must be rejected for hybrid-gauge EXX";
    }
    catch (const std::invalid_argument& error)
    {
        EXPECT_STREQ(error.what(), "td_stype=2 EXX requires exx_real_number=0");
    }
}

TEST(ExxInfo, AcceptsComplexHybridGauge)
{
    auto input = hse_input();
    input.exx_real_number = "0";
    Exx_Info info;
    EXPECT_NO_THROW(init_exx_info(info, input));
    EXPECT_FALSE(info.info_ri.real_number);
}

TEST(ExxInfo, PreservesRealRiForOtherGaugesAndGroundState)
{
    auto input = hse_input();
    input.td_stype = 1;
    Exx_Info info;
    EXPECT_NO_THROW(init_exx_info(info, input));
    EXPECT_TRUE(info.info_ri.real_number);
    input.td_stype = 2;
    input.esolver_type = "ksdft";
    EXPECT_NO_THROW(init_exx_info(info, input));
    EXPECT_TRUE(info.info_ri.real_number);
}
