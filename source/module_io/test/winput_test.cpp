#include "gtest/gtest.h"
#include "gmock/gmock.h"
/************************************************
 *  unit test of winput.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Read()
 *     - Read in parameters about Wannier functions
 *   - Print()
 *     - Print out parameters about Wannier functions
 */

#define private public
#include "../winput.h"

class WInputTest : public testing::Test
{
protected:
	std::ifstream ifs;
	std::string output;
};

TEST_F(WInputTest, Read)
{
	winput::Default();
	std::stringstream ss1;
	ss1 << "./support/WINPUT";
	winput::Read(ss1.str());
	EXPECT_EQ(winput::target,"test");
	EXPECT_EQ(winput::wlmr_dir        ,"./");
	EXPECT_EQ(winput::rcut            ,10);
	EXPECT_EQ(winput::before_iter     ,0);
	EXPECT_EQ(winput::after_iter      ,0);
	EXPECT_EQ(winput::begin_stop_flag ,0);
	EXPECT_EQ(winput::end_flag        ,0);
	EXPECT_EQ(winput::wf_type         ,"V");
	EXPECT_EQ(winput::build_wf        ,0);
	EXPECT_EQ(winput::imp_pao         ,0);
	EXPECT_EQ(winput::b_out_wf        ,0);
	EXPECT_EQ(winput::b_fftwan        ,0);
	EXPECT_EQ(winput::b_plot_build    ,0);
	EXPECT_EQ(winput::b_plot_atomic   ,0);
	EXPECT_EQ(winput::trial           ,"atomic");
	EXPECT_DOUBLE_EQ(winput::bs       ,2.5);
	EXPECT_EQ(winput::bp              ,2);
	EXPECT_EQ(winput::px              ,2);
	EXPECT_EQ(winput::g1              ,3);
	EXPECT_EQ(winput::g2              ,3);
	EXPECT_EQ(winput::bloch_begin     ,0);
	EXPECT_EQ(winput::bloch_end       ,0);
	EXPECT_EQ(winput::no_center       ,0);
	EXPECT_EQ(winput::sph_proj        ,0);
	EXPECT_EQ(winput::sph_type        ,0);
	EXPECT_EQ(winput::b_recon         ,0);
	EXPECT_EQ(winput::speed_mode      ,1);
	EXPECT_EQ(winput::recon_wanq      ,0);
	EXPECT_EQ(winput::b_mix_wf        ,0);
	EXPECT_EQ(winput::mix_wf          ,0);
	EXPECT_EQ(winput::b_near_atom     ,0);
	EXPECT_EQ(winput::range0          ,0);
	EXPECT_EQ(winput::range1          ,0);
	EXPECT_EQ(winput::L_start         ,0);
	EXPECT_EQ(winput::L_end           ,2);
	EXPECT_EQ(winput::atom_start      ,0);
	EXPECT_EQ(winput::atom_end        ,1);
	EXPECT_EQ(winput::trunc_ao        ,6);
	EXPECT_EQ(winput::trunc_wlmr      ,14);
	EXPECT_EQ(winput::trunc_wan       ,6);
	EXPECT_EQ(winput::fermi_t         ,1);
	EXPECT_DOUBLE_EQ(winput::clm2_lowest,1e-07);
	EXPECT_EQ(winput::plot_wanq       ,0);
	EXPECT_EQ(winput::plot_option     ,"(110)");
	EXPECT_EQ(winput::n_unitcell      ,2);
	EXPECT_EQ(winput::out_all         ,0);
	EXPECT_EQ(winput::out_chg         ,0);
	EXPECT_EQ(winput::compare_atomic  ,0);
	EXPECT_EQ(winput::cal_bands       ,0);
	EXPECT_EQ(winput::cal_bands2      ,0);
	EXPECT_EQ(winput::charge_type     ,"planewave");
	EXPECT_EQ(winput::cal_dos         ,0);
	EXPECT_EQ(winput::out_spillage    ,0);
	EXPECT_EQ(winput::spillage_outdir ,"./");
	EXPECT_EQ(winput::mesh            ,999);
	EXPECT_DOUBLE_EQ(winput::dr       ,0.01);
	EXPECT_EQ(winput::sum_lm          ,0);
}

TEST_F(WInputTest, Print)
{
	winput::Default();
	std::stringstream ss1;
	ss1 << "WINPUT_out";
	winput::Print(ss1.str());
	ifs.open("WINPUT_out");
	getline(ifs,output);
	// test output in warning.log file
	EXPECT_THAT(output,testing::HasSubstr("WANNIER_PARAMETERS"));
	ifs.close();
	remove("WINPUT_out");
}
