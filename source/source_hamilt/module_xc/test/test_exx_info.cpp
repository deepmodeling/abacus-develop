#include "../exx_info.h"
#include "source_io/module_parameter/input_parameter.h"
#include <gtest/gtest.h>
#include <stdexcept>

TEST(RpaCutoffConfig, HseDefaultsRemainSupportedButNondefaultCutoffsAreRejected)
{
    Input_para cut_input;
    cut_input.basis_type = "lcao";
    cut_input.rpa = true;
    cut_input.dft_functional = "hse";
    // Prepared HSE defaults, after INPUT reset hooks have resolved "default".
    cut_input.exx_fock_alpha.clear();
    cut_input.exx_fock_lambda.clear();
    cut_input.exx_erfc_alpha = {"0.25"};
    cut_input.exx_erfc_omega = {"0.11"};
    cut_input.exx_singularity_correction = "limits";
    cut_input.exx_real_number = "1";
    cut_input.exx_ccp_rmesh_times = "2";
    Exx_Info cut_default_info;
    EXPECT_NO_THROW(init_exx_info(cut_default_info, cut_input));
    EXPECT_EQ(cut_default_info.info_ri.coulomb_param.count(Conv_Coulomb_Pot_K::Coulomb_Type::Fock), 0);
    EXPECT_EQ(cut_default_info.info_ri.coulomb_param.count(Conv_Coulomb_Pot_K::Coulomb_Type::Erfc), 1);

    cut_input.rpa_cut_coulomb = "log_erfc";
    Exx_Info cut_smooth_info;
    EXPECT_THROW(init_exx_info(cut_smooth_info, cut_input), std::invalid_argument);
    cut_input.rpa_cut_coulomb = "hard";
    cut_input.rpa_cut_radius = "inscribed";
    Exx_Info cut_inscribed_info;
    EXPECT_THROW(init_exx_info(cut_inscribed_info, cut_input), std::invalid_argument);
}

TEST(RpaCutoffConfig, PbeRiInitializationMatchesLibriAvailability)
{
    Input_para cut_input;
    cut_input.basis_type = "lcao";
    cut_input.rpa = true;
    cut_input.dft_functional = "pbe";
    cut_input.exx_fock_lambda.clear();
    cut_input.exx_erfc_omega.clear();
    cut_input.exx_singularity_correction = "limits";
    cut_input.exx_real_number = "1";
    cut_input.exx_ccp_rmesh_times = "2";
    cut_input.rpa_cut_coulomb = "log_erfc";
    cut_input.rpa_cut_radius = "inscribed";
    Exx_Info cut_pbe_info;
    ASSERT_NO_THROW(init_exx_info(cut_pbe_info, cut_input));
#ifdef __EXX
    // ENABLE_LIBRI defines __EXX and enables PBE RPA initialization of info_ri.
    ASSERT_EQ(cut_pbe_info.info_ri.coulomb_param.count(Conv_Coulomb_Pot_K::Coulomb_Type::Fock), 1);
    const auto& cut_fock = cut_pbe_info.info_ri.coulomb_param.at(Conv_Coulomb_Pot_K::Coulomb_Type::Fock);
    ASSERT_EQ(cut_fock.size(), 1);
    EXPECT_EQ(cut_fock.front().at("rpa_cut_coulomb"), "log_erfc");
    EXPECT_EQ(cut_fock.front().at("rpa_cut_radius"), "inscribed");
    EXPECT_DOUBLE_EQ(std::stod(cut_fock.front().at("rpa_cut_width")), cut_input.rpa_cut_width);
#else
    // Keep coverage when LibRI is disabled: PBE must skip RI initialization.
    EXPECT_TRUE(cut_pbe_info.info_ri.coulomb_param.empty());
#endif
}
