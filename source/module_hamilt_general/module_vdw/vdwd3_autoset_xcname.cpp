#ifndef DFTD3_XC_NAME_H
#define DFTD3_XC_NAME_H
/**
 * Intro
 * -----
 * This file stores the mapping from LibXC xcname to the "conventional"
 * 
 * XCNotSupportedError
 * -------------------
 * GGA_X_REVSSB_D
 * GGA_X_SSB_D
 * 
 * in J. Chem. Phys. 131, 094103 2009, a simplified version of PBC (the
 * correlation part of PBE XC) is used as the correlation part, but libXC
 * does not directly support one named as
 * GGA_C_SPBEC.
 * 
 * Certainly, those XC with dispersion correction in form of non-local
 * correlation are not supported. Such as:
 * 
 * vdw-DF family nonlocal dispersion correction included are not supported:
 * GGA_X_OPTB86B_VDW
 * GGA_X_OPTB88_VDW 
 * GGA_X_OPTPBE_VDW 
 * GGA_X_PBEK1_VDW
 * 
 * VV09, VV10 and rVV10 nonlocal correlation included are not supported:
 * GGA_XC_VV10 
 * HYB_GGA_XC_LC_VV10 
 * HYB_MGGA_XC_WB97M_V 
 * HYB_GGA_XC_WB97X_V 
 * MGGA_X_VCML 
 * MGGA_C_REVSCAN_VV10 
 * MGGA_C_SCAN_VV10 
 * MGGA_C_SCANL_VV10 
 * MGGA_XC_B97M_V 
 * MGGA_XC_VCML_RVV10
 * 
 * There is also one quite special, the wB97X-D3BJ functional uses the
 * wB97X-V functionals excluding the VV10 part, then use its own DFT-D3(BJ)
 * parameters. This seems not recorded in simple-dftd3, so it is not supported
 * temporarily:
 * HYB_GGA_XC_WB97X_D3BJ
 * HYB_GGA_XC_WB97X_V
 */
#include <map>
#include <string>
#include <cassert>
#include "module_base/formatter.h"
#include <iostream>
#include <regex>
#include "module_base/tool_quit.h"
#include "module_hamilt_general/module_vdw/vdwd3_parameters.h"

const std::map<std::string, std::string> xcname_libxc_xc_ = {
    {"XC_LDA_XC_TETER93", "teter93"},
    {"XC_LDA_XC_ZLP", "zlp"},
    {"XC_MGGA_XC_OTPSS_D", "otpss_d"}, // DFT-D2
    {"XC_GGA_XC_OPBE_D", "opbe_d"}, // DFT-D2
    {"XC_GGA_XC_OPWLYP_D", "opwlyp_d"}, // DFT-D2
    {"XC_GGA_XC_OBLYP_D", "oblyp_d"}, // DFT-D2
    {"XC_GGA_XC_HCTH_407P", "hcth_407p"},
    {"XC_GGA_XC_HCTH_P76", "hcth_p76"},
    {"XC_GGA_XC_HCTH_P14", "hcth_p14"},
    {"XC_GGA_XC_B97_GGA1", "b97_gga1"},
    {"XC_GGA_XC_KT2", "kt2"},
    {"XC_GGA_XC_TH1", "th1"},
    {"XC_GGA_XC_TH2", "th2"},
    {"XC_GGA_XC_TH3", "th3"},
    {"XC_GGA_XC_TH4", "th4"},
    {"XC_GGA_XC_HCTH_93", "hcth_93"},
    {"XC_GGA_XC_HCTH_120", "hcth_120"},
    {"XC_GGA_XC_HCTH_147", "hcth_147"},
    {"XC_GGA_XC_HCTH_407", "hcth_407"},
    {"XC_GGA_XC_EDF1", "edf1"},
    {"XC_GGA_XC_XLYP", "xlyp"},
    {"XC_GGA_XC_KT1", "kt1"},
    {"XC_GGA_XC_B97_D", "b97_d"}, // DFT-D2?
    {"XC_GGA_XC_PBE1W", "pbe1w"},
    {"XC_GGA_XC_MPWLYP1W", "mpwlyp1w"},
    {"XC_GGA_XC_PBELYP1W", "pbelyp1w"},
    {"XC_HYB_LDA_XC_LDA0", "lda0"},
    {"XC_HYB_LDA_XC_CAM_LDA0", "cam_lda0"},
    {"XC_GGA_XC_NCAP", "ncap"},
    {"XC_GGA_XC_MOHLYP", "mohlyp"},
    {"XC_GGA_XC_MOHLYP2", "mohlyp2"},
    {"XC_GGA_XC_TH_FL", "th_fl"},
    {"XC_GGA_XC_TH_FC", "th_fc"},
    {"XC_GGA_XC_TH_FCFO", "th_fcfo"},
    {"XC_GGA_XC_TH_FCO", "th_fco"},
    {"XC_MGGA_XC_CC06", "cc06"},
    {"XC_MGGA_XC_TPSSLYP1W", "tpsslyp1w"},
    {"XC_MGGA_XC_B97M_V", "b97m_v"},
    {"XC_GGA_XC_VV10", "vv10"},
    {"XC_LDA_XC_KSDT", "ksdt"},
    {"XC_HYB_GGA_XC_B97_1P", "b97_1p"},
    {"XC_HYB_GGA_XC_PBE_MOL0", "pbe_mol0"},
    {"XC_HYB_GGA_XC_PBE_SOL0", "pbe_sol0"},
    {"XC_HYB_GGA_XC_PBEB0", "pbeb0"},
    {"XC_HYB_GGA_XC_PBE_MOLB0", "pbe_molb0"},
    {"XC_GGA_XC_BEEFVDW", "beefvdw"},
    {"XC_MGGA_XC_HLE17", "hle17"},
    {"XC_HYB_GGA_XC_PBE50", "pbe50"},
    {"XC_HYB_GGA_XC_HFLYP", "hflyp"},
    {"XC_HYB_GGA_XC_B3P86_NWCHEM", "b3p86_nwchem"},
    {"XC_LDA_XC_CORRKSDT", "corrksdt"},
    {"XC_HYB_GGA_XC_RELPBE0", "relpbe0"},
    {"XC_GGA_XC_B97_3C", "b97_3c"},
    {"XC_HYB_MGGA_XC_BR3P86", "br3p86"},
    {"XC_HYB_GGA_XC_CASE21", "case21"},
    {"XC_HYB_GGA_XC_PBE_2X", "pbe_2x"},
    {"XC_HYB_GGA_XC_PBE38", "pbe38"},
    {"XC_HYB_GGA_XC_B3LYP3", "b3lyp3"},
    {"XC_HYB_GGA_XC_CAM_O3LYP", "cam_o3lyp"},
    {"XC_HYB_MGGA_XC_TPSS0", "tpss0"},
    {"XC_HYB_MGGA_XC_B94_HYB", "b94_hyb"},
    {"XC_HYB_GGA_XC_WB97X_D3", "wb97x_d3"}, // DFT-D3(0)
    {"XC_HYB_GGA_XC_LC_BLYP", "lc_blyp"},
    {"XC_HYB_GGA_XC_B3PW91", "b3pw91"},
    {"XC_HYB_GGA_XC_B3LYP", "b3lyp"},
    {"XC_HYB_GGA_XC_B3P86", "b3p86"},
    {"XC_HYB_GGA_XC_O3LYP", "o3lyp"},
    {"XC_HYB_GGA_XC_MPW1K", "mpw1k"},
    {"XC_HYB_GGA_XC_PBEH", "pbeh"},
    {"XC_HYB_GGA_XC_B97", "b97"},
    {"XC_HYB_GGA_XC_B97_1", "b97_1"},
    {"XC_HYB_GGA_XC_APF", "apf"},
    {"XC_HYB_GGA_XC_B97_2", "b97_2"},
    {"XC_HYB_GGA_XC_X3LYP", "x3lyp"},
    {"XC_HYB_GGA_XC_B1WC", "b1wc"},
    {"XC_HYB_GGA_XC_B97_K", "b97_k"},
    {"XC_HYB_GGA_XC_B97_3", "b97_3"},
    {"XC_HYB_GGA_XC_MPW3PW", "mpw3pw"},
    {"XC_HYB_GGA_XC_B1LYP", "b1lyp"},
    {"XC_HYB_GGA_XC_B1PW91", "b1pw91"},
    {"XC_HYB_GGA_XC_MPW1PW", "mpw1pw"},
    {"XC_HYB_GGA_XC_MPW3LYP", "mpw3lyp"},
    {"XC_HYB_GGA_XC_SB98_1A", "sb98_1a"},
    {"XC_HYB_GGA_XC_SB98_1B", "sb98_1b"},
    {"XC_HYB_GGA_XC_SB98_1C", "sb98_1c"},
    {"XC_HYB_GGA_XC_SB98_2A", "sb98_2a"},
    {"XC_HYB_GGA_XC_SB98_2B", "sb98_2b"},
    {"XC_HYB_GGA_XC_SB98_2C", "sb98_2c"},
    {"XC_HYB_GGA_XC_HSE03", "hse03"},
    {"XC_HYB_GGA_XC_HSE06", "hse06"},
    {"XC_HYB_GGA_XC_HJS_PBE", "hjs_pbe"},
    {"XC_HYB_GGA_XC_HJS_PBE_SOL", "hjs_pbe_sol"},
    {"XC_HYB_GGA_XC_HJS_B88", "hjs_b88"},
    {"XC_HYB_GGA_XC_HJS_B97X", "hjs_b97x"},
    {"XC_HYB_GGA_XC_CAM_B3LYP", "cam_b3lyp"},
    {"XC_HYB_GGA_XC_TUNED_CAM_B3LYP", "tuned_cam_b3lyp"},
    {"XC_HYB_GGA_XC_BHANDH", "bhandh"},
    {"XC_HYB_GGA_XC_BHANDHLYP", "bhandhlyp"},
    {"XC_HYB_GGA_XC_MB3LYP_RC04", "mb3lyp_rc04"},
    {"XC_HYB_MGGA_XC_B88B95", "b88b95"},
    {"XC_HYB_MGGA_XC_B86B95", "b86b95"},
    {"XC_HYB_MGGA_XC_PW86B95", "pw86b95"},
    {"XC_HYB_MGGA_XC_BB1K", "bb1k"},
    {"XC_HYB_MGGA_XC_MPW1B95", "mpw1b95"},
    {"XC_HYB_MGGA_XC_MPWB1K", "mpwb1k"},
    {"XC_HYB_MGGA_XC_X1B95", "x1b95"},
    {"XC_HYB_MGGA_XC_XB1K", "xb1k"},
    {"XC_HYB_MGGA_XC_PW6B95", "pw6b95"},
    {"XC_HYB_MGGA_XC_PWB6K", "pwb6k"},
    {"XC_HYB_GGA_XC_MPWLYP1M", "mpwlyp1m"},
    {"XC_HYB_GGA_XC_REVB3LYP", "revb3lyp"},
    {"XC_HYB_GGA_XC_CAMY_BLYP", "camy_blyp"},
    {"XC_HYB_GGA_XC_PBE0_13", "pbe0_13"},
    {"XC_HYB_MGGA_XC_TPSSH", "tpssh"},
    {"XC_HYB_MGGA_XC_REVTPSSH", "revtpssh"},
    {"XC_HYB_GGA_XC_B3LYPS", "b3lyps"},
    {"XC_HYB_GGA_XC_QTP17", "qtp17"},
    {"XC_HYB_GGA_XC_B3LYP_MCM1", "b3lyp_mcm1"},
    {"XC_HYB_GGA_XC_B3LYP_MCM2", "b3lyp_mcm2"},
    {"XC_HYB_GGA_XC_WB97", "wb97"},
    {"XC_HYB_GGA_XC_WB97X", "wb97x"},
    {"XC_HYB_GGA_XC_LRC_WPBEH", "lrc_wpbeh"},
    {"XC_HYB_GGA_XC_WB97X_V", "wb97x_v"},
    {"XC_HYB_GGA_XC_LCY_PBE", "lcy_pbe"},
    {"XC_HYB_GGA_XC_LCY_BLYP", "lcy_blyp"},
    {"XC_HYB_GGA_XC_LC_VV10", "lc_vv10"},
    {"XC_HYB_GGA_XC_CAMY_B3LYP", "camy_b3lyp"},
    {"XC_HYB_GGA_XC_WB97X_D", "wb97x_d"}, // DFT-D2
    {"XC_HYB_GGA_XC_HPBEINT", "hpbeint"},
    {"XC_HYB_GGA_XC_LRC_WPBE", "lrc_wpbe"},
    {"XC_HYB_GGA_XC_B3LYP5", "b3lyp5"},
    {"XC_HYB_GGA_XC_EDF2", "edf2"},
    {"XC_HYB_GGA_XC_CAP0", "cap0"},
    {"XC_HYB_GGA_XC_LC_WPBE", "lc_wpbe"},
    {"XC_HYB_GGA_XC_HSE12", "hse12"},
    {"XC_HYB_GGA_XC_HSE12S", "hse12s"},
    {"XC_HYB_GGA_XC_HSE_SOL", "hse_sol"},
    {"XC_HYB_GGA_XC_CAM_QTP_01", "cam_qtp_01"},
    {"XC_HYB_GGA_XC_MPW1LYP", "mpw1lyp"},
    {"XC_HYB_GGA_XC_MPW1PBE", "mpw1pbe"},
    {"XC_HYB_GGA_XC_KMLYP", "kmlyp"},
    {"XC_HYB_GGA_XC_LC_WPBE_WHS", "lc_wpbe_whs"},
    {"XC_HYB_GGA_XC_LC_WPBEH_WHS", "lc_wpbeh_whs"},
    {"XC_HYB_GGA_XC_LC_WPBE08_WHS", "lc_wpbe08_whs"},
    {"XC_HYB_GGA_XC_LC_WPBESOL_WHS", "lc_wpbesol_whs"},
    {"XC_HYB_GGA_XC_CAM_QTP_00", "cam_qtp_00"},
    {"XC_HYB_GGA_XC_CAM_QTP_02", "cam_qtp_02"},
    {"XC_HYB_GGA_XC_LC_QTP", "lc_qtp"},
    {"XC_HYB_GGA_XC_BLYP35", "blyp35"},
    {"XC_HYB_MGGA_XC_WB97M_V", "wb97m_v"},
    {"XC_LDA_XC_1D_EHWLRG_1", "1d_ehwlrg_1"},
    {"XC_LDA_XC_1D_EHWLRG_2", "1d_ehwlrg_2"},
    {"XC_LDA_XC_1D_EHWLRG_3", "1d_ehwlrg_3"},
    {"XC_GGA_XC_HLE16", "hle16"},
    {"XC_LDA_XC_LP_A", "lp_a"},
    {"XC_LDA_XC_LP_B", "lp_b"},
    {"XC_HYB_MGGA_XC_B0KCIS", "b0kcis"},
    {"XC_MGGA_XC_LP90", "lp90"},
    {"XC_HYB_MGGA_XC_MPW1KCIS", "mpw1kcis"},
    {"XC_HYB_MGGA_XC_MPWKCIS1K", "mpwkcis1k"},
    {"XC_HYB_MGGA_XC_PBE1KCIS", "pbe1kcis"},
    {"XC_HYB_MGGA_XC_TPSS1KCIS", "tpss1kcis"},
    {"XC_HYB_GGA_XC_B5050LYP", "b5050lyp"},
    {"XC_LDA_XC_GDSMFB", "gdsmfb"},
    {"XC_GGA_XC_KT3", "kt3"},
    {"XC_HYB_LDA_XC_BN05", "bn05"},
    {"XC_HYB_GGA_XC_LB07", "lb07"},
    {"XC_HYB_MGGA_XC_B98", "b98"},
    {"XC_LDA_XC_TIH", "tih"},
    {"XC_HYB_GGA_XC_APBE0", "apbe0"},
    {"XC_HYB_GGA_XC_HAPBE", "hapbe"},
    {"XC_HYB_GGA_XC_RCAM_B3LYP", "rcam_b3lyp"},
    {"XC_HYB_GGA_XC_WC04", "wc04"},
    {"XC_HYB_GGA_XC_WP04", "wp04"},
    {"XC_HYB_GGA_XC_CAMH_B3LYP", "camh_b3lyp"},
    {"XC_HYB_GGA_XC_WHPBE0", "whpbe0"},
    {"XC_HYB_GGA_XC_LC_BLYP_EA", "lc_blyp_ea"},
    {"XC_HYB_GGA_XC_LC_BOP", "lc_bop"},
    {"XC_HYB_GGA_XC_LC_PBEOP", "lc_pbeop"},
    {"XC_HYB_GGA_XC_LC_BLYPR", "lc_blypr"},
    {"XC_HYB_GGA_XC_MCAM_B3LYP", "mcam_b3lyp"},
    {"XC_MGGA_XC_VCML_RVV10", "vcml_rvv10"},
    {"XC_HYB_MGGA_XC_GAS22", "gas22"},
    {"XC_HYB_MGGA_XC_R2SCANH", "r2scanh"},
    {"XC_HYB_MGGA_XC_R2SCAN0", "r2scan0"},
    {"XC_HYB_MGGA_XC_R2SCAN50", "r2scan50"},
    {"XC_HYB_GGA_XC_CAM_PBEH", "cam_pbeh"},
    {"XC_HYB_GGA_XC_CAMY_PBEH", "camy_pbeh"},
    {"XC_HYB_MGGA_XC_EDMGGAH", "edmggah"},
    {"XC_HYB_MGGA_XC_LC_TMLYP", "lc_tmlyp"},
};
const std::map<std::string, std::string> xcname_libxc_xplusc_ = {
    {"XC_GGA_X_GAM+XC_GGA_C_GAM", "gam"},
    {"XC_GGA_X_HCTH_A+XC_GGA_C_HCTH_A", "hcth_a"},
    {"XC_HYB_MGGA_X_DLDF+XC_MGGA_C_DLDF", "dldf"},
    {"XC_GGA_X_Q2D+XC_GGA_C_Q2D", "q2d"},
    {"XC_GGA_X_PBE_MOL+XC_GGA_C_PBE_MOL", "pbe_mol"},
    {"XC_GGA_X_PBEINT+XC_GGA_C_PBEINT", "pbeint"},
    {"XC_HYB_GGA_X_N12_SX+XC_GGA_C_N12_SX", "n12_sx"},
    {"XC_GGA_X_N12+XC_GGA_C_N12", "n12"},
    {"XC_GGA_X_PBE+XC_GGA_C_PBE", "pbe"},
    {"XC_GGA_X_B88+XC_MGGA_C_B88", "b88"},
    {"XC_GGA_X_PW91+XC_GGA_C_PW91", "pw91"},
    {"XC_GGA_X_PBE_SOL+XC_GGA_C_PBE_SOL", "pbe_sol"},
    {"XC_GGA_X_AM05+XC_GGA_C_AM05", "am05"},
    {"XC_GGA_X_XPBE+XC_GGA_C_XPBE", "xpbe"},
    {"XC_GGA_X_RGE2+XC_GGA_C_RGE2", "rge2"},
    {"XC_GGA_X_SOGGA11+XC_GGA_C_SOGGA11", "sogga11"},
    {"XC_GGA_X_APBE+XC_GGA_C_APBE", "apbe"},
    {"XC_MGGA_X_TPSS+XC_MGGA_C_TPSS", "tpss"},
    {"XC_MGGA_X_M06_L+XC_MGGA_C_M06_L", "m06_l"},
    {"XC_HYB_MGGA_X_TAU_HCTH+XC_GGA_C_TAU_HCTH", "tau_hcth"},
    {"XC_MGGA_X_REVTPSS+XC_MGGA_C_REVTPSS", "revtpss"},
    {"XC_MGGA_X_PKZB+XC_MGGA_C_PKZB", "pkzb"},
    {"XC_MGGA_X_M11_L+XC_MGGA_C_M11_L", "m11_l"},
    {"XC_MGGA_X_MN12_L+XC_MGGA_C_MN12_L", "mn12_l"},
    {"XC_HYB_MGGA_X_MN12_SX+XC_MGGA_C_MN12_SX", "mn12_sx"},
    {"XC_MGGA_X_MN15_L+XC_MGGA_C_MN15_L", "mn15_l"},
    {"XC_MGGA_X_SCAN+XC_MGGA_C_SCAN", "scan"},
    {"XC_GGA_X_PBEFE+XC_GGA_C_PBEFE", "pbefe"},
    {"XC_HYB_MGGA_X_MN15+XC_MGGA_C_MN15", "mn15"},
    {"XC_HYB_MGGA_X_BMK+XC_GGA_C_BMK", "bmk"},
    {"XC_MGGA_X_REVM06_L+XC_MGGA_C_REVM06_L", "revm06_l"},
    {"XC_HYB_MGGA_X_M08_HX+XC_MGGA_C_M08_HX", "m08_hx"},
    {"XC_HYB_MGGA_X_M08_SO+XC_MGGA_C_M08_SO", "m08_so"},
    {"XC_HYB_MGGA_X_M11+XC_MGGA_C_M11", "m11"},
    {"XC_GGA_X_CHACHIYO+XC_GGA_C_CHACHIYO", "chachiyo"},
    {"XC_HYB_MGGA_X_REVM11+XC_MGGA_C_REVM11", "revm11"},
    {"XC_HYB_MGGA_X_REVM06+XC_MGGA_C_REVM06", "revm06"},
    {"XC_HYB_MGGA_X_M06_SX+XC_MGGA_C_M06_SX", "m06_sx"},
    {"XC_GGA_X_PBE_GAUSSIAN+XC_GGA_C_PBE_GAUSSIAN", "pbe_gaussian"},
    {"XC_HYB_GGA_X_SOGGA11_X+XC_GGA_C_SOGGA11_X", "sogga11_x"},
    {"XC_HYB_MGGA_X_M05+XC_MGGA_C_M05", "m05"},
    {"XC_HYB_MGGA_X_M05_2X+XC_MGGA_C_M05_2X", "m05_2x"},
    {"XC_HYB_MGGA_X_M06_HF+XC_MGGA_C_M06_HF", "m06_hf"},
    {"XC_HYB_MGGA_X_M06+XC_MGGA_C_M06", "m06"},
    {"XC_HYB_MGGA_X_M06_2X+XC_MGGA_C_M06_2X", "m06_2x"},
    {"XC_MGGA_X_RSCAN+XC_MGGA_C_RSCAN", "rscan"},
    {"XC_MGGA_X_R2SCAN+XC_MGGA_C_R2SCAN", "r2scan"},
    {"XC_GGA_X_SG4+XC_GGA_C_SG4", "sg4"},
    {"XC_MGGA_X_TM+XC_MGGA_C_TM", "tm"},
    {"XC_MGGA_X_REVSCAN+XC_MGGA_C_REVSCAN", "revscan"},
    {"XC_MGGA_X_REGTPSS+XC_GGA_C_REGTPSS", "regtpss"},
    {"XC_MGGA_X_R2SCAN01+XC_MGGA_C_R2SCAN01", "r2scan01"},
    {"XC_MGGA_X_RPPSCAN+XC_MGGA_C_RPPSCAN", "rppscan"},
    {"XC_MGGA_X_REVTM+XC_MGGA_C_REVTM", "revtm"},
    {"XC_MGGA_X_SCANL+XC_MGGA_C_SCANL", "scanl"},
    {"XC_MGGA_X_MGGAC+XC_GGA_C_MGGAC", "mggac"},
    {"XC_MGGA_X_R2SCANL+XC_MGGA_C_R2SCANL", "r2scanl"},
    {"XC_GGA_X_B88+XC_GGA_C_LYP", "blyp"},
    {"XC_GGA_X_B88+XC_GGA_C_P86", "bp86"},
    {"XC_GGA_X_PW91+XC_GGA_C_PW91", "pw91"},
    {"XC_GGA_X_PBE+XC_GGA_C_PBE", "pbe"},
    {"XC_GGA_X_PBE_SOL+XC_GGA_C_PBE_SOL", "pbesol"},
    {"XC_MGGA_X_PKZB+XC_MGGA_C_PKZB", "pkzb"},
    {"XC_MGGA_X_TPSS+XC_MGGA_C_TPSS", "tpss"},
    {"XC_MGGA_X_REVTPSS+XC_MGGA_C_REVTPSS", "revtpss"},
    {"XC_MGGA_X_SCAN+XC_MGGA_C_SCAN", "scan"},
    {"XC_GGA_X_SOGGA+XC_GGA_C_PBE", "sogga"},
    {"XC_MGGA_X_BLOC+XC_MGGA_C_TPSSLOC", "bloc"},
    {"XC_GGA_X_OPTX+XC_GGA_C_LYP", "olyp"},
    {"XC_GGA_X_RPBE+XC_GGA_C_PBE", "rpbe"},
    {"XC_GGA_X_B88+XC_GGA_C_PBE", "bpbe"},
    {"XC_GGA_X_MPW91+XC_GGA_C_PW91", "mpw91"},
    {"XC_MGGA_X_MS0+XC_GGA_C_REGTPSS", "ms0"},
    {"XC_MGGA_X_MS1+XC_GGA_C_REGTPSS", "ms1"},
    {"XC_MGGA_X_MS2+XC_GGA_C_REGTPSS", "ms2"},
    {"XC_HYB_MGGA_X_MS2H+XC_GGA_C_REGTPSS", "ms2h"},
    {"XC_MGGA_X_MVS+XC_GGA_C_REGTPSS", "mvs"},
    {"XC_HYB_MGGA_X_MVSH+XC_GGA_C_REGTPSS", "mvsh"},
    {"XC_GGA_X_SOGGA11+XC_GGA_C_SOGGA11", "sogga11"},
    {"XC_HYB_GGA_X_SOGGA11_X+XC_GGA_C_SOGGA11_X", "sogga11-x"},
    {"XC_HYB_MGGA_X_DLDF+XC_MGGA_C_DLDF", "dldf"},
    {"XC_GGA_X_GAM+XC_GGA_C_GAM", "gam"},
    {"XC_MGGA_X_M06_L+XC_MGGA_C_M06_L", "m06-l"},
    {"XC_MGGA_X_M11_L+XC_MGGA_C_M11_L", "m11-l"},
    {"XC_MGGA_X_MN12_L+XC_MGGA_C_MN12_L", "mn12-l"},
    {"XC_MGGA_X_MN15_L+XC_MGGA_C_MN15_L", "mn15-l"},
    {"XC_GGA_X_N12+XC_GGA_C_N12", "n12"},
    {"XC_HYB_GGA_X_N12_SX+XC_GGA_C_N12_SX", "n12-sx"},
    {"XC_HYB_MGGA_X_MN12_SX+XC_MGGA_C_MN12_SX", "mn12-sx"},
    {"XC_HYB_MGGA_X_MN15+XC_MGGA_C_MN15", "mn15"},
    {"XC_MGGA_X_MBEEF+XC_GGA_C_PBE_SOL", "mbeef"},
    {"XC_HYB_MGGA_X_SCAN0+XC_MGGA_C_SCAN", "scan0"},
    {"XC_GGA_X_PBE+XC_GGA_C_OP_PBE", "pbeop"},
    {"XC_GGA_X_B88+XC_GGA_C_OP_B88", "bop"}
};

void _xcname_libxc_xplusc(const std::string& xcpattern, std::string& xname)
{
    std::vector<std::string> xc_words = FmtCore::split(xcpattern, "+");
    std::for_each(xc_words.begin(), xc_words.end(), [](std::string& s) {
            s = (FmtCore::startswith(s, "XC_")? s: "XC_" + s); }); // add XC_ if not present
    assert(xc_words.size() == 2);

    std::vector<std::string> words = FmtCore::split(xc_words[0], "_");
    const std::string key = (words[2] == "X")? 
        xc_words[0] + "+" + xc_words[1]: xc_words[1] + "+" + xc_words[0];

    if (xcname_libxc_xplusc_.find(key) != xcname_libxc_xplusc_.end()) {
        xname = xcname_libxc_xplusc_.at(key);
    } else {
        ModuleBase::WARNING_QUIT("ModuleHamiltGeneral::ModuleVDW::DFTD3::xcname_libxc_xplusc",
                                    "XC's LibXC-notation on `" + xcpattern + "` not recognized");
    }
}

void _xcname_libxc_xc(const std::string& xcpattern, std::string& xname)
{
    // add XC_ if not present
    const std::string key = FmtCore::startswith(xcpattern, "XC_")? xcpattern: "XC_" + xcpattern;

    if (xcname_libxc_xc_.find(key) != xcname_libxc_xc_.end()) {
        xname = xcname_libxc_xc_.at(key);
    } else {
        ModuleBase::WARNING_QUIT("ModuleHamiltGeneral::ModuleVDW::DFTD3::xcname_libxc_xc",
                                    "XC's LibXC-notation on `" + xcpattern + "` not recognized");
    }
}

void _xcname_libxc(const std::string& xcpattern, std::string& xname)
{
    if (xcpattern.find("+") != std::string::npos) {
        _xcname_libxc_xplusc(xcpattern, xname);
    } else {
        _xcname_libxc_xc(xcpattern, xname);
    }
}

std::string vdw::Vdwd3Parameters::_vdwd3_xcname(const std::string& xcpattern)
{
    std::string xcname = xcpattern;
    const std::regex pattern("(LDA|GGA|MGGA|HYB|HYB_LDA|HYB_GGA|HYB_MGGA)_(X|C|XC|K)_(.*)");
    // as long as there is piece in xcpattern that can match, we can search for the corresponding name
    if (std::regex_search(xcpattern, pattern)) {
        _xcname_libxc(xcpattern, xcname);
    }
    return xcname;
}

/**
import os
import re
def read_xc_func_h(fn):
    with open(fn) as f:
        lines = f.readlines()
    out = {}
    for line in lines:
        words = line.strip().split()
        xc, xcid = words[1], int(words[2])
        xc_annos = ' '.join(words[4:-1])
        out[xc] = {'id': xcid, 'annos': xc_annos}
    return out

def sort_xc(xc_data):
    '''Sort the xc functionals into x, c, xc, k functionals.
    
    Parameters
    ----------
    xc_data : dict
        from function read_xc_func_h
    
    Returns
    -------
    dict, dict, dict, dict
        The dictionaries of x, c, xc, k functionals, whose keys are the
        like LDA, GGA, MGGA, HYB, HYB_LDA, HYB_GGA, HYB_MGGA, values are
        the dictionaries of the functionals, whose keys are the conventional
        xc name, values include approx, annos, id, full.
    '''
    x, c, xc, k = {}, {}, {}, {}
    dictmap = {'X': x, 'C': c, 'XC': xc, 'K': k}
    xcpat = r'XC_(LDA|GGA|MGGA|HYB|HYB_LDA|HYB_GGA|HYB_MGGA)_(X|C|XC|K)_(.*)'
    for xc_name, data in xc_data.items():
        m = re.match(xcpat, xc_name)
        if m is None:
            print('Warning: cannot match', xc_name)
            continue
        approx, type_, name = m.groups()
        dictmap[type_][name] = {'approx': approx, 'annos': data['annos'],
                                'id': data['id'], 'full': xc_name}
    return x, c, xc, k

def pair_xc(x, c):
    '''
    Pair the x and c functionals.
    
    Parameters
    ----------
    x : dict
        The dictionary of x functionals, whose keys are the conventional
        xc name, values include approx, annos, id, full.
    
    c : dict
        the same as x
    
    Returns
    -------
    dict, dict
        The dictionary of paired and unpaired x and c functionals, whose keys are the
        conventional xc name, values are the dictionary of x and c functionals.
    '''
    paired, unpaired = {}, {}
    for xc_name, data in x.items():
        if xc_name in c:
            paired[xc_name] = {'x': data, 'c': c[xc_name]}
        else:
            unpaired[xc_name] = data
    return paired, unpaired

def xc_to_stdmap(xc, conventional_lower=True):
    '''print the xc in the way of c++ std::map<std::string, std::string>.
    
    Parameters
    ----------
    xc : dict
        The dictionary of xc functionals, whose keys are the conventional
        xc name, values include approx, annos, id, full.
    conventional_lower : bool
        Whether to convert the conventional name to lower case.
    
    Returns
    -------
    str
        The string of c++ code, std::map<std::string, std::string> mapping
        the full name of xc to its conventional name.
    '''
    out = 'const std::map<std::string, std::string> xcname_libxc_xc_ = {\n'
    for name, data in xc.items():
        name = name.lower() if conventional_lower else name
        out += '    {"%s", "%s"},\n' % (data['full'], name)
    out += '};\n'
    return out
    
def paired_xc_to_stdmap(pairs, conventional_lower=True):
    '''print the xc in the way of c++ std::map<std::string, std::string>.
    
    Parameters
    ----------
    pairs : dict
        The dictionary of xc functionals, whose keys are the conventional
        xc name, values include approx, annos, id, full.
    conventional_lower : bool
        Whether to convert the conventional name to lower case.
    
    Returns
    -------
    str
        The string of c++ code, std::map<std::string, std::string> mapping
        the full name of xc to its conventional name.
    '''
    out = 'const std::map<std::string, std::string> xcname_libxc_xplusc_ = {\n'
    for name, data in pairs.items():
        name = name.lower() if conventional_lower else name
        plus = f'{data["x"]["full"]}+{data["c"]["full"]}'
        out += '    {"%s", "%s"},\n' % (plus, name)
        # sulp = f'{data["c"]["full"]}+{data["x"]["full"]}'
        # out += '    {"%s", "%s"},\n' % (sulp, name)
    out += '};\n'
    return out

def special_x_and_c(x, c):
    '''Special pairings of x and c functionals. The following data sheet is 
    from Pyscf: 
    https://github.com/pyscf/pyscf/blob/master/pyscf/dft/xcfun.py
    Thanks for pointing out the bug by @QuantumMiska and the help from wsr (@hebrewsnabla)
    
    
    Parameters
    ----------
    x : dict
        The dictionary of x functionals, whose keys are the conventional
        xc name, values include approx, annos, id, full.
    
    c : dict
        the same as x
    
    Returns
    -------
    dict
        The dictionary of special pairings of x and c functionals.
    '''
    DATA = {
        'BLYP'     : 'B88,LYP',
        'BP86'     : 'B88,P86',
        'PW91'     : 'PW91,PW91',
        'PBE'      : 'PBE,PBE',
        'REVPBE'   : 'REVPBE,PBE',
        'PBESOL'   : 'PBE_SOL,PBE_SOL',
        'PKZB'     : 'PKZB,PKZB',
        'TPSS'     : 'TPSS,TPSS',
        'REVTPSS'  : 'REVTPSS,REVTPSS',
        'SCAN'     : 'SCAN,SCAN',
        'SOGGA'    : 'SOGGA,PBE',
        'BLOC'     : 'BLOC,TPSSLOC',
        'OLYP'     : 'OPTX,LYP',
        'RPBE'     : 'RPBE,PBE',
        'BPBE'     : 'B88,PBE',
        'MPW91'    : 'MPW91,PW91',
        'HFLYP'    : 'HF,LYP',
        'HFPW92'   : 'HF,PWMOD',
        'SPW92'    : 'SLATER,PWMOD',
        'SVWN'     : 'SLATER,VWN',
        'MS0'      : 'MS0,REGTPSS',
        'MS1'      : 'MS1,REGTPSS',
        'MS2'      : 'MS2,REGTPSS',
        'MS2H'     : 'MS2H,REGTPSS',
        'MVS'      : 'MVS,REGTPSS',
        'MVSH'     : 'MVSH,REGTPSS',
        'SOGGA11'  : 'SOGGA11,SOGGA11',
        'SOGGA11-X': 'SOGGA11_X,SOGGA11_X',
        'KT1'      : 'KT1X,VWN',
        'DLDF'     : 'DLDF,DLDF',
        'GAM'      : 'GAM,GAM',
        'M06-L'    : 'M06_L,M06_L',
        'M11-L'    : 'M11_L,M11_L',
        'MN12-L'   : 'MN12_L,MN12_L',
        'MN15-L'   : 'MN15_L,MN15_L',
        'N12'      : 'N12,N12',
        'N12-SX'   : 'N12_SX,N12_SX',
        'MN12-SX'  : 'MN12_SX,MN12_SX',
        'MN15'     : 'MN15,MN15',
        'MBEEF'    : 'MBEEF,PBE_SOL',
        'SCAN0'    : 'SCAN0,SCAN',
        'PBEOP'    : 'PBE,OP_PBE',
        'BOP'      : 'B88,OP_B88',
    }
    paired = {}
    for name, data in DATA.items():
        xname, cname = data.split(',')
        if xname in x and cname in c:
            paired[name] = {'x': x[xname], 'c': c[cname]}
        else:
            print(f'Warning: {name} not found in x or c: {xname}, {cname}')
    return paired

def print_xc(xc):
    print(f'{"Name":20s} {"Full":30s} {"Appr":10s} {"Annos"}')
    for name, data in xc.items():
        print(f'{name:20s} {data["full"]:30s} {data["approx"]:10s} {data["annos"]}')

if __name__ == '__main__':
    libxc = '/root/soft/libxc/libxc-6.2.2'
    f = 'src/xc_funcs.h'
    xc_data = read_xc_func_h(os.path.join(libxc, f))
    x, c, xc, k = sort_xc(xc_data)
    pairs, others = pair_xc(x, c)
    special = special_x_and_c(x, c)
    # print(xc_to_stdmap(xc))
    # print(paired_xc_to_stdmap(pairs))
    # print_xc(others)
    print(paired_xc_to_stdmap(special))
 */
#endif // DFTD3_XCNAME_H_
