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
    EXPECT_TRUE(info.info_ri.unfold_hybrid_images);
}

TEST(ExxInfo, PreservesRealRiForOtherGaugesAndGroundState)
{
    auto input = hse_input();
    input.td_stype = 1;
    Exx_Info info;
    EXPECT_NO_THROW(init_exx_info(info, input));
    EXPECT_TRUE(info.info_ri.real_number);
    EXPECT_FALSE(info.info_ri.unfold_hybrid_images);
    input.td_stype = 2;
    input.esolver_type = "ksdft";
    EXPECT_NO_THROW(init_exx_info(info, input));
    EXPECT_TRUE(info.info_ri.real_number);
    EXPECT_FALSE(info.info_ri.unfold_hybrid_images);
}

TEST(ExxInfo, RejectsUnsupportedHybridImageKPointPools)
{
    auto input = hse_input();
    input.exx_real_number = "0";
    input.kpar = 2;
    Exx_Info info;
    EXPECT_THROW(init_exx_info(info, input), std::invalid_argument);
}

TEST(ExxInfo, LeavesPlaneWaveExchangeImageDomainUnchanged)
{
    auto input = hse_input();
    input.basis_type = "pw";
    input.exx_real_number = "0";
    Exx_Info info;
    EXPECT_NO_THROW(init_exx_info(info, input));
    EXPECT_FALSE(info.info_ri.unfold_hybrid_images);
}

TEST(ExxInfo, UnfoldsDirectFullRangeKernels)
{
    auto input = hse_input();
    input.exx_real_number = "0";
    input.exx_erfc_alpha.clear();
    input.exx_erfc_omega.clear();
    input.exx_fock_alpha = {"0.25"};
    input.exx_fock_lambda = {"0.3"};
    for (const auto& functional : {"pbe0", "hf"})
    {
        input.dft_functional = functional;
        for (const auto& method : {"spencer", "revised_spencer", "limits"})
        {
            input.exx_singularity_correction = method;
            Exx_Info info;
            EXPECT_NO_THROW(init_exx_info(info, input));
            EXPECT_TRUE(info.info_ri.unfold_hybrid_images);
            input.kpar = 2;
            EXPECT_THROW(init_exx_info(info, input), std::invalid_argument);
            input.kpar = 1;
        }
    }
}

TEST(ExxInfo, PreservesEwaldImagePathIncludingMixedKernels)
{
    auto input = hse_input();
    input.dft_functional = "pbe0";
    input.exx_real_number = "0";
    input.exx_fock_alpha = {"0.25"};
    input.exx_fock_lambda = {"0.3"};
    input.kpar = 2;
    for (const auto& method : {"massidda", "carrier"})
    {
        input.exx_singularity_correction = method;
        for (const bool mixed : {false, true})
        {
            input.exx_erfc_alpha.clear();
            input.exx_erfc_omega.clear();
            if (mixed)
            {
                input.exx_erfc_alpha = {"0.25"};
                input.exx_erfc_omega = {"0.11"};
            }
            Exx_Info info;
            EXPECT_NO_THROW(init_exx_info(info, input));
            EXPECT_FALSE(info.info_ri.unfold_hybrid_images);
        }
    }
}

TEST(ExxInfo, ErfcRemainsDirectRegardlessOfFockSingularityMethod)
{
    auto input = hse_input();
    input.exx_real_number = "0";
    input.exx_singularity_correction = "massidda";
    Exx_Info info;
    EXPECT_NO_THROW(init_exx_info(info, input));
    EXPECT_TRUE(info.info_ri.unfold_hybrid_images);
}
