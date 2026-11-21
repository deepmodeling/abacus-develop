#include "source_io/module_parameter/dftu_dc_check.h"
#include "source_io/module_parameter/input_parameter.h"
#include <gtest/gtest.h>

#include <functional>
#include <limits>
#include <stdexcept>

namespace
{
Input_para supported_input()
{
    Input_para input;
    input.dft_plus_u = 1;
    input.dftu_dc = "mixed";
    input.dftu_alpha = 0.5;
    input.basis_type = "lcao";
    input.ks_solver = "genelpa";
    input.dft_functional = "PBE";
    input.exx_fock_alpha.clear();
    input.exx_erfc_alpha.clear();
    input.symmetry = "0";
    return input;
}
} // namespace

TEST(DftuDcCheck, LegacyDefaultsHaveNoNewScopeRestrictions)
{
    Input_para input;
    EXPECT_EQ(input.dftu_dc, "fll");
    EXPECT_TRUE(ModuleIO::dftu_dc_error(input, false, 4).empty());
    input.dft_plus_u = 2;
    input.nspin = 4;
    input.lspinorb = true;
    input.cal_force = true;
    input.cal_stress = true;
    input.device = "gpu";
    input.restart_load = true;
    input.mixing_dftu = true;
    EXPECT_TRUE(ModuleIO::dftu_dc_error(input, false, 4).empty());
}

TEST(DftuDcCheck, ExplicitAlphaAndModeValidation)
{
    auto input = supported_input();
    EXPECT_FALSE(ModuleIO::dftu_dc_error(input, false, 1).empty());
    EXPECT_TRUE(ModuleIO::dftu_dc_error(input, true, 1).empty());
    for (double alpha : {0.0, 0.5, 1.0})
    {
        input.dftu_alpha = alpha;
        EXPECT_TRUE(ModuleIO::dftu_dc_error(input, true, 1).empty());
    }
    for (double alpha : {-0.1, 1.1, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()})
    {
        input.dftu_alpha = alpha;
        EXPECT_FALSE(ModuleIO::dftu_dc_error(input, true, 1).empty());
    }
    input.dftu_dc = "amf";
    EXPECT_FALSE(ModuleIO::dftu_dc_error(input, true, 1).empty());
    EXPECT_TRUE(ModuleIO::dftu_dc_error(input, false, 1).empty());
    input.dftu_dc = "fll";
    EXPECT_FALSE(ModuleIO::dftu_dc_error(input, true, 1).empty());
    input.dftu_dc = "AMF";
    EXPECT_FALSE(ModuleIO::dftu_dc_error(input, false, 1).empty());
}

TEST(DftuDcCheck, ParserRejectsTrailingGarbageAndNonfiniteValues)
{
    EXPECT_DOUBLE_EQ(ModuleIO::parse_dftu_alpha("0"), 0.0);
    EXPECT_DOUBLE_EQ(ModuleIO::parse_dftu_alpha("1"), 1.0);
    EXPECT_DOUBLE_EQ(ModuleIO::parse_dftu_alpha("2.5e-1"), 0.25);
    for (const std::string token : {"", "nan", "inf", "-inf", "0.5junk", "0.5 0.2", "-0.01", "1.01", "1e9999"})
    {
        EXPECT_THROW(ModuleIO::parse_dftu_alpha(token), std::exception) << token;
    }
}

TEST(DftuDcCheck, ScopeGuardsApplyToMixedEndpointsButNeverLegacyFll)
{
    const std::vector<std::function<void(Input_para&)>> forbidden = {
        [](Input_para& p) { p.dft_plus_u = 0; },
        [](Input_para& p) { p.dft_plus_u = 2; },
        [](Input_para& p) { p.dft_plus_u = 3; },
        [](Input_para& p) { p.basis_type = "pw"; },
        [](Input_para& p) { p.device = "gpu"; },
        [](Input_para& p) { p.precision = "single"; },
        [](Input_para& p) { p.calculation = "relax"; },
        [](Input_para& p) { p.calculation = "md"; },
        [](Input_para& p) { p.calculation = "nscf"; },
        [](Input_para& p) { p.cal_force = true; },
        [](Input_para& p) { p.cal_stress = true; },
        [](Input_para& p) { p.esolver_type = "tddft"; },
        [](Input_para& p) { p.socket_driver = true; },
        [](Input_para& p) { p.nspin = 4; },
        [](Input_para& p) { p.lspinorb = true; },
        [](Input_para& p) { p.noncolin = true; },
        [](Input_para& p) { p.sc_mag_switch = true; },
        [](Input_para& p) { p.nupdown = 12.0; },
        [](Input_para& p) { p.ks_solver = "pexsi"; },
        [](Input_para& p) { p.kpar = 2; },
        [](Input_para& p) { p.bndpar = 2; },
        [](Input_para& p) { p.symmetry = "1"; },
        [](Input_para& p) { p.init_chg = "file"; },
        [](Input_para& p) { p.init_wfc = "file"; },
        [](Input_para& p) { p.init_occ_mat = 1; },
        [](Input_para& p) { p.omc = 2; },
        [](Input_para& p) { p.restart_load = true; },
        [](Input_para& p) { p.restart_save = true; },
        [](Input_para& p) { p.dm_to_rho = true; },
        [](Input_para& p) { p.uramping_eV = 0.1; },
        [](Input_para& p) { p.yukawa_potential = true; },
        [](Input_para& p) { p.mixing_dftu = true; },
        [](Input_para& p) { p.mixing_dmr = true; },
        [](Input_para& p) { p.dft_functional = "hse"; },
        [](Input_para& p) { p.dft_functional = "default"; },
        [](Input_para& p) { p.exx_fock_alpha = {"0.25"}; },
        [](Input_para& p) { p.exx_erfc_alpha = {"0.25"}; },
        [](Input_para& p) { p.deepks_scf = true; },
        [](Input_para& p) { p.dft_plus_dmft = true; },
        [](Input_para& p) { p.rpa = true; },
        [](Input_para& p) { p.out_hsk[0] = 1; },
        [](Input_para& p) { p.out_hsr[0] = 1; },
        [](Input_para& p) { p.out_mat_hs2[0] = 1; },
        [](Input_para& p) { p.out_dmr[0] = 1; },
        [](Input_para& p) { p.out_eband = true; },
        [](Input_para& p) { p.out_eband_terms_compat = true; },
        [](Input_para& p) { p.out_mat_tk[0] = 1; },
        [](Input_para& p) { p.out_mat_vnlk[0] = 1; },
        [](Input_para& p) { p.out_mat_vlk[0] = 1; },
        [](Input_para& p) { p.out_mat_vhk[0] = 1; },
        [](Input_para& p) { p.out_mat_vxck[0] = 1; },
        [](Input_para& p) { p.out_mat_exxk[0] = 1; },
        [](Input_para& p) { p.out_mat_dh_vxc[0] = 1; },
        [](Input_para& p) { p.towannier90 = true; },
        [](Input_para& p) { p.qo_switch = true; }};
    for (std::size_t index = 0; index < forbidden.size(); ++index)
    {
        SCOPED_TRACE(index);
        auto input = supported_input();
        forbidden[index](input);
        for (double alpha : {0.0, 0.5, 1.0})
        {
            input.dftu_alpha = alpha;
            EXPECT_FALSE(ModuleIO::dftu_dc_error(input, true, 1).empty());
        }
        input.dftu_dc = "amf";
        EXPECT_FALSE(ModuleIO::dftu_dc_error(input, false, 1).empty());
        input.dftu_dc = "fll";
        EXPECT_TRUE(ModuleIO::dftu_dc_error(input, false, 1).empty());
    }
    auto input = supported_input();
    EXPECT_FALSE(ModuleIO::dftu_dc_error(input, true, 2).empty());
}

TEST(DftuDcCheck, CollinearSpinAndTimeReversalOptions)
{
    auto input = supported_input();
    for (int spin : {1, 2})
    {
        input.nspin = spin;
        for (const std::string symmetry : {"0", "-1"})
        {
            input.symmetry = symmetry;
            EXPECT_TRUE(ModuleIO::dftu_dc_error(input, true, 1).empty());
        }
    }
}
