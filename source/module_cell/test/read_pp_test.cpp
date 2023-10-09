#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include<memory>
/************************************************
 *  unit test of read_pp
 ***********************************************/

/**
 * - Tested Functions:
 *   - ReadUPF100
 *     - read_pseudo_upf
 *       - read 1.0.0 version of upf pseudopotential file
 *     - read_pseudo_header
 *     - read_pseudo_mesh
 *     - read_pseudo_nl
 *     - read_pseudo_local
 *     - read_pseudo_pswfc
 *     - read_pseudo_rhoatom
 *     - read_pseudo_so
 *     - print_pseudo_upf
 *   - ReadUPF201
 *   - ReadUSPPUPF201
 *     - read_pseudo_upf201
 *       - read 2.0.1 version of upf pseudopotential file
 *     - getnameval
 *     - read_pseudo_upf201_header
 *     - read_pseudo_upf201_mesh
 *     - read_pseudo_upf201_nonlocal
 *     - read_pseudo_upf201_pswfc
 *     - void read_pseudo_upf201_so
 *   - ReadUSppErr100
 *     - read_pseudo_nl
 *     - read_pseudo_nlcc
 *     - ultrasoft is not supported
 *   - HeaderErr201
 *     - coulomb and paw pp are not supported
 *     - no pp_header
 *     - arbitrary is not read in
 *     - semi-local pp is not supported
 *   - ReadUPF201MESH2
 *     - a different "<PP_MESH" header in UPF file
 *   - ReadUPF201FR
 *     - read a full-relativistic pp
 *   - VWR
 *     - read vwr type of pp
 *   - VWR
 *     - read vwr type of pp
 *   - BLPS
 *     - read blps type of pp
 *   - SetPseudoType
 *     - find the type of pp, upf201 or upf
 *   - Trim
 *     - trim: an iterative function to delete all tab and space in a string
 *     - trimend: trim tab and space and space at two ends of a string
 *   - InitReader
 *     - init_pseudo_reader: actual pseudo reader
 *   - SetEmptyElement
 *     - set_empty_element: used in BSSE calculation
 *   - SetPseudoType
 *     - set_pseudo_type: tell it is upf or upf201 from the pp file's header line
 *   - Average
 *     - average_p: modulate the soc effect in pseudopotential
 */

#define private public
#include "module_cell/read_pp.h"

class ReadPPTest : public testing::Test
{
protected:
	std::string output;
	std::unique_ptr<Pseudopot_upf> upf{new Pseudopot_upf};
};

TEST_F(ReadPPTest, ReadUPF100)
{
	std::ifstream ifs;
	ifs.open("./support/Te.pbe-rrkj.UPF");
	upf->read_pseudo_upf(ifs);
	EXPECT_TRUE(upf->has_so); // has soc info
	EXPECT_EQ(upf->nv,0); // number of version
	EXPECT_EQ(upf->psd,"Te"); // element label
	EXPECT_EQ(upf->pp_type,"NC"); // pp_type
	EXPECT_FALSE(upf->tvanp); // not ultrasoft
	EXPECT_FALSE(upf->nlcc); // no Nonlinear core correction
	EXPECT_EQ(upf->xc_func,"PBE"); // Exchange-Correlation functional
	EXPECT_EQ(upf->zp,6); // Z valence
	EXPECT_DOUBLE_EQ(upf->etotps,-15.54533017755); // total energy
	EXPECT_DOUBLE_EQ(upf->ecutwfc,0.0); // suggested cutoff for wfc
	EXPECT_DOUBLE_EQ(upf->ecutrho,0.0); // suggested cutoff for rho
	EXPECT_EQ(upf->lmax,2); // max angular momentum component
	EXPECT_EQ(upf->mesh,1191); // Number of points in mesh
	EXPECT_EQ(upf->nwfc,3); // Number of wavefunctions
	EXPECT_EQ(upf->nbeta,3); // Number of projectors
	EXPECT_EQ(upf->els[0],"5S"); // label for i-th atomic orbital
	EXPECT_EQ(upf->els[1],"5P"); // label for i-th atomic orbital
	EXPECT_EQ(upf->els[2],"5D"); // label for i-th atomic orbital
	EXPECT_EQ(upf->lchi[0],0); // angluar momentum of each atomic orbital
	EXPECT_EQ(upf->lchi[1],1); // angluar momentum of each atomic orbital
	EXPECT_EQ(upf->lchi[2],2); // angluar momentum of each atomic orbital
	EXPECT_DOUBLE_EQ(upf->oc[0],2.0); // occupation of each atomic orbital
	EXPECT_DOUBLE_EQ(upf->oc[1],3.0); // occupation of each atomic orbital
	EXPECT_DOUBLE_EQ(upf->oc[2],1.0); // occupation of each atomic orbital
	EXPECT_DOUBLE_EQ(upf->r[0],1.75361916453E-05); // r
	EXPECT_DOUBLE_EQ(upf->r[upf->mesh-1],5.05901190442E+01); // r
	EXPECT_DOUBLE_EQ(upf->rab[0],2.19202395566E-07); // rab
	EXPECT_DOUBLE_EQ(upf->rab[upf->mesh-1],6.32376488053E-01); // rab
	EXPECT_DOUBLE_EQ(upf->vloc[0],-5.00890143222E+00); // vloc
	EXPECT_DOUBLE_EQ(upf->vloc[upf->mesh-1],-2.37200471955E-01); // vloc
	EXPECT_EQ(upf->lll[0],0); // BETA
    EXPECT_EQ(upf->kbeta[0], 957);
    EXPECT_DOUBLE_EQ(upf->beta(0, 0), -1.82560984478E-03);
    EXPECT_DOUBLE_EQ(upf->beta(0, upf->kbeta[0] - 1), -1.61398366674E-03);
    EXPECT_EQ(upf->lll[1], 0); // BETA
    EXPECT_EQ(upf->kbeta[1], 957);
    EXPECT_DOUBLE_EQ(upf->beta(1, 0), 1.51692136792E-03);
    EXPECT_DOUBLE_EQ(upf->beta(1, upf->kbeta[1] - 1), 1.51458412218E-03);
    EXPECT_EQ(upf->lll[2], 2); // BETA
    EXPECT_EQ(upf->kbeta[2], 957);
    EXPECT_DOUBLE_EQ(upf->beta(2, 0), -3.10582746893E-13);
    EXPECT_DOUBLE_EQ(upf->beta(2, upf->kbeta[2] - 1), -4.17131335030E-04);
    EXPECT_EQ(upf->nd,4); // DIJ
	EXPECT_DOUBLE_EQ(upf->dion(0,0),-1.70394647943E-01);
	EXPECT_DOUBLE_EQ(upf->dion(0,1),-1.76521654672E-01);
	EXPECT_DOUBLE_EQ(upf->dion(1,1),-1.80323263809E-01);
	EXPECT_DOUBLE_EQ(upf->dion(2,2),1.16612440320E-02);
	EXPECT_DOUBLE_EQ(upf->chi(0,0),1.40252610787E-06);  // PSWFC
	EXPECT_DOUBLE_EQ(upf->chi(0,upf->mesh-1),6.15962544650E-25);
	EXPECT_DOUBLE_EQ(upf->chi(1,0),7.65306256201E-11);
	EXPECT_DOUBLE_EQ(upf->chi(1,upf->mesh-1),1.44320361049E-17);
	EXPECT_DOUBLE_EQ(upf->chi(2,0),4.37015997370E-16);
	EXPECT_DOUBLE_EQ(upf->chi(2,upf->mesh-1),6.20093850585E-05);
	EXPECT_DOUBLE_EQ(upf->rho_at[0],3.93415898407E-12); // RhoAtom
	EXPECT_DOUBLE_EQ(upf->rho_at[upf->mesh-1],3.84516383534E-09);
	EXPECT_EQ(upf->nn[0],1); // nn
	EXPECT_EQ(upf->nn[1],2);
	EXPECT_EQ(upf->nn[2],3);
	EXPECT_DOUBLE_EQ(upf->jchi[0],0.0); // jchi
	EXPECT_DOUBLE_EQ(upf->jchi[1],0.0);
	EXPECT_NEAR(upf->jchi[3],0.0,1e-20);
	EXPECT_DOUBLE_EQ(upf->jjj[0],0.0); // jjj
	EXPECT_DOUBLE_EQ(upf->jjj[1],0.0);
	EXPECT_DOUBLE_EQ(upf->jjj[2],0.0);
	//EXPECT_EQ
	ifs.close();
	std::ofstream ofs;
	ofs.open("tmp");
	upf->print_pseudo_upf(ofs);
	ofs.close();
	ifs.open("tmp");
	getline(ifs, output);
	EXPECT_THAT(output, testing::HasSubstr("==== read_pseudo_upf ==="));
	ifs.close();
}

TEST_F(ReadPPTest, ReadUSppErr100)
{
	std::ifstream ifs;
	ifs.open("./support/Zn.pw91-n-van.UPF");
	//upf->read_pseudo_upf(ifs);
	testing::internal::CaptureStdout();
	EXPECT_EXIT(upf->read_pseudo_upf(ifs),
			::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	// test output on screening
	// EXPECT_THAT(output,testing::HasSubstr("this function is called")); // read_pseudo_nlcc
	EXPECT_THAT(output,testing::HasSubstr("Ultra Soft Pseudopotential not available yet."));
	ifs.close();
}

TEST_F(ReadPPTest, ReadUPF201)
{
	std::ifstream ifs;
	ifs.open("./support/Cu_ONCV_PBE-1.0.upf");
	upf->read_pseudo_upf201(ifs);
	EXPECT_EQ(upf->psd,"Cu");
	EXPECT_EQ(upf->pp_type,"NC");
	EXPECT_FALSE(upf->has_so);
	EXPECT_FALSE(upf->nlcc);
	EXPECT_EQ(upf->xc_func,"PBE");
	EXPECT_EQ(upf->zp,19);
	EXPECT_DOUBLE_EQ(upf->etotps,-1.82394100797E+02);
	EXPECT_EQ(upf->lmax,2);
	EXPECT_EQ(upf->mesh,601); // mesh -= 1 at line 388 (why? Let's see)
	EXPECT_EQ(upf->nwfc,0);
	EXPECT_EQ(upf->nbeta,6);
	EXPECT_DOUBLE_EQ(upf->r[0],0.0);
	EXPECT_DOUBLE_EQ(upf->r[601],6.01);
	EXPECT_DOUBLE_EQ(upf->rab[0],0.01);
	EXPECT_DOUBLE_EQ(upf->rab[601],0.01);
	EXPECT_DOUBLE_EQ(upf->vloc[0],-5.3426582174E+01);
	EXPECT_DOUBLE_EQ(upf->vloc[601],-6.3227960060E+00);
	EXPECT_EQ(upf->lll[0],0);
    EXPECT_EQ(upf->kbeta[0], 196);
    EXPECT_DOUBLE_EQ(upf->beta(0,0),0.0);
	EXPECT_DOUBLE_EQ(upf->beta(0,601),0.0);
	EXPECT_EQ(upf->lll[1],0);
    EXPECT_EQ(upf->kbeta[1], 196);
    EXPECT_DOUBLE_EQ(upf->beta(1,0),0.0);
	EXPECT_DOUBLE_EQ(upf->beta(1,601),0.0);
	EXPECT_EQ(upf->lll[2],1);
    EXPECT_EQ(upf->kbeta[2], 196);
    EXPECT_DOUBLE_EQ(upf->beta(2,0),0.0);
	EXPECT_DOUBLE_EQ(upf->beta(2,601),0.0);
	EXPECT_EQ(upf->lll[3],1);
    EXPECT_EQ(upf->kbeta[3], 196);
    EXPECT_DOUBLE_EQ(upf->beta(3,0),0.0);
	EXPECT_DOUBLE_EQ(upf->beta(3,601),0.0);
	EXPECT_EQ(upf->lll[4],2);
    EXPECT_EQ(upf->kbeta[4], 196);
    EXPECT_DOUBLE_EQ(upf->beta(4,0),0.0);
	EXPECT_DOUBLE_EQ(upf->beta(4,601),0.0);
	EXPECT_EQ(upf->lll[5],2);
    EXPECT_EQ(upf->kbeta[5], 196);
    EXPECT_DOUBLE_EQ(upf->beta(5,0),0.0);
	EXPECT_DOUBLE_EQ(upf->beta(5,601),0.0);
	EXPECT_DOUBLE_EQ(upf->dion(0,0),-6.6178420255E+00);
	EXPECT_DOUBLE_EQ(upf->dion(5,5),-7.0938557228E+00);
	EXPECT_DOUBLE_EQ(upf->rho_at[0],0.0);
	EXPECT_DOUBLE_EQ(upf->rho_at[601],3.1742307110E-02);
	ifs.close();
}

TEST_F(ReadPPTest, ReadUSPPUPF201)
{
    std::ifstream ifs;
    ifs.open("./support/Al.pbe-sp-van.UPF");
    upf->read_pseudo_upf201(ifs);
    EXPECT_EQ(upf->psd, "Al");
    EXPECT_EQ(upf->pp_type, "US");
    EXPECT_EQ(upf->relativistic, "no");
    EXPECT_TRUE(upf->tvanp);
    EXPECT_FALSE(upf->has_so);
    EXPECT_FALSE(upf->nlcc);
    EXPECT_EQ(upf->xc_func, "PBE");
    EXPECT_EQ(upf->zp, 11);
    EXPECT_EQ(upf->nv, 0);
    EXPECT_DOUBLE_EQ(upf->etotps, -1.596432307730e2);
    EXPECT_DOUBLE_EQ(upf->ecutwfc, 0.0);
    EXPECT_DOUBLE_EQ(upf->ecutrho, 0.0);
    EXPECT_EQ(upf->lmax, 2);
    EXPECT_EQ(upf->lmax_rho, 0);
    EXPECT_EQ(upf->lloc, 0);
    EXPECT_EQ(upf->mesh, 893);
    EXPECT_EQ(upf->nwfc, 4);
    EXPECT_EQ(upf->nbeta, 5);
    EXPECT_EQ(upf->kkbeta, 617);
    EXPECT_DOUBLE_EQ(upf->rmax, 2.006810756590e2);
    EXPECT_DOUBLE_EQ(upf->zmesh, 13.0);
    EXPECT_FALSE(upf->q_with_l);
    EXPECT_EQ(upf->nqf, 8);
    EXPECT_EQ(upf->nqlc, 5);
    EXPECT_DOUBLE_EQ(upf->r[0], 0.0);
    EXPECT_DOUBLE_EQ(upf->r[892], 2.006810756590000e2);
    EXPECT_DOUBLE_EQ(upf->rab[0], 1.169079443020000e-6);
    EXPECT_DOUBLE_EQ(upf->rab[892], 3.344685763390000e0);
    EXPECT_DOUBLE_EQ(upf->vloc[0], 3.456089057550000e0);
    EXPECT_DOUBLE_EQ(upf->vloc[892], -1.096266796840000e-1);
    EXPECT_EQ(upf->rho_atc, nullptr);
    EXPECT_EQ(upf->lll[0], 0);
    EXPECT_EQ(upf->kbeta[0], 617);
    EXPECT_EQ(upf->els_beta[0], "2S");
    EXPECT_DOUBLE_EQ(upf->rcut[0], 0.0);
    EXPECT_DOUBLE_EQ(upf->rcutus[0], 1.4);
    EXPECT_DOUBLE_EQ(upf->beta(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(upf->beta(0, 892), 0.0);
    EXPECT_EQ(upf->lll[1], 0);
    EXPECT_EQ(upf->kbeta[1], 617);
    EXPECT_EQ(upf->els_beta[1], "2P");
    EXPECT_DOUBLE_EQ(upf->rcut[1], 0.0);
    EXPECT_DOUBLE_EQ(upf->rcutus[1], 1.42);
    EXPECT_DOUBLE_EQ(upf->beta(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(upf->beta(1, 892), 0.0);
    EXPECT_DOUBLE_EQ(upf->dion(0, 0), -2.182408428460000e2);
    EXPECT_DOUBLE_EQ(upf->dion(4, 4), -3.087562171130000e1);
    EXPECT_DOUBLE_EQ(upf->qqq(0, 0), 3.896280866700000e0);
    EXPECT_DOUBLE_EQ(upf->qqq(4, 4), 7.218068659650000e-1);
    EXPECT_DOUBLE_EQ(upf->qfcoef(0, 0, 0, 0), 8.705252055130002e1);
    EXPECT_DOUBLE_EQ(upf->qfcoef(4, 4, 4, 7), 9.910935792140002e1);
    EXPECT_DOUBLE_EQ(upf->rinner[0], 1.1);
    EXPECT_DOUBLE_EQ(upf->rinner[4], 1.1);
    EXPECT_DOUBLE_EQ(upf->qfunc(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(upf->qfunc(0, 892), 0.0);
    EXPECT_EQ(upf->els[0], "2S");
    EXPECT_EQ(upf->lchi[0], 0);
    EXPECT_EQ(upf->nchi[0], 0);
    EXPECT_DOUBLE_EQ(upf->oc[0], 2.0);
    EXPECT_DOUBLE_EQ(upf->epseu[0], 0.0);
    EXPECT_DOUBLE_EQ(upf->rcut_chi[0], 0.0);
    EXPECT_DOUBLE_EQ(upf->rcutus_chi[0], 1.4);
    EXPECT_DOUBLE_EQ(upf->chi(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(upf->chi(0, 892), 0.0);
    EXPECT_DOUBLE_EQ(upf->rho_at[0], 0.0);
    EXPECT_DOUBLE_EQ(upf->rho_at[892], 0.0);
    EXPECT_EQ(upf->jchi, nullptr);
    EXPECT_EQ(upf->jjj, nullptr);
    EXPECT_EQ(upf->nn, nullptr);
    ifs.close();
}

TEST_F(ReadPPTest, HeaderErr2011)
{
	std::ifstream ifs;
	// 1st
	ifs.open("./support/HeaderError1");
	//upf->read_pseudo_upf201(ifs);
	testing::internal::CaptureStdout();
	EXPECT_EXIT(upf->read_pseudo_upf201(ifs),
			::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Found no PP_HEADER"));
    ifs.close();
}

TEST_F(ReadPPTest, HeaderErr2012)
{
	std::ifstream ifs;
	// 2nd
	ifs.open("./support/HeaderError2");
	//upf->read_pseudo_upf201(ifs);
	testing::internal::CaptureStdout();
	EXPECT_EXIT(upf->read_pseudo_upf201(ifs),
			::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("SEMI-LOCAL PSEUDOPOTENTIAL IS NOT SUPPORTED"));
    ifs.close();
}

TEST_F(ReadPPTest, HeaderErr2013)
{
	std::ifstream ifs;
	// 3rd
	ifs.open("./support/HeaderError3");
	//upf->read_pseudo_upf201(ifs);
	testing::internal::CaptureStdout();
	EXPECT_EXIT(upf->read_pseudo_upf201(ifs),
			::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("PAW POTENTIAL IS NOT SUPPORTED"));
    ifs.close();
}

TEST_F(ReadPPTest, HeaderErr2014)
{
	std::ifstream ifs;
    // 3rd
    ifs.open("./support/HeaderError4");
    // upf->read_pseudo_upf201(ifs);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(upf->read_pseudo_upf201(ifs), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("COULOMB POTENTIAL IS NOT SUPPORTED"));
    ifs.close();
}

TEST_F(ReadPPTest, HeaderErr2015)
{
    std::ifstream ifs;
    // 4th
    GlobalV::ofs_warning.open("warning.log");
    ifs.open("./support/HeaderError5");
    upf->read_pseudo_upf201(ifs);
    GlobalV::ofs_warning.close();
    ifs.close();
    ifs.open("warning.log");
    getline(ifs, output);
    EXPECT_THAT(
        output,
        testing::HasSubstr("arbitrary is not read in. Please add this parameter in read_pp_upf201.cpp if needed."));
    ifs.close();
    remove("warning.log");
}

TEST_F(ReadPPTest, ReadUPF201FR)
{
	std::ifstream ifs;
	// this is a dojo full-relativisitic pp
	ifs.open("./support/C.upf");
	upf->read_pseudo_upf201(ifs);
	EXPECT_EQ(upf->psd,"C");
	EXPECT_TRUE(upf->has_so);
	EXPECT_TRUE(upf->nlcc);
	EXPECT_EQ(upf->mesh,1247);
	//RELBETA
	EXPECT_EQ(upf->nbeta,6);
	EXPECT_EQ(upf->lll[0],0);
	EXPECT_EQ(upf->lll[1],0);
	EXPECT_EQ(upf->lll[2],1);
	EXPECT_EQ(upf->lll[3],1);
	EXPECT_EQ(upf->lll[4],1);
	EXPECT_EQ(upf->lll[5],1);
	EXPECT_DOUBLE_EQ(upf->jjj[0],0.5);
	EXPECT_DOUBLE_EQ(upf->jjj[1],0.5);
	EXPECT_DOUBLE_EQ(upf->jjj[2],0.5);
	EXPECT_DOUBLE_EQ(upf->jjj[3],1.5);
	EXPECT_DOUBLE_EQ(upf->jjj[4],0.5);
	EXPECT_DOUBLE_EQ(upf->jjj[5],1.5);
	//RELWFC
	EXPECT_EQ(upf->nwfc,3);
	EXPECT_EQ(upf->nn[0],1);
	EXPECT_EQ(upf->nn[1],2);
	EXPECT_EQ(upf->nn[2],2);
	EXPECT_EQ(upf->lchi[0],0);
	EXPECT_EQ(upf->lchi[1],1);
	EXPECT_EQ(upf->lchi[2],1);
	EXPECT_DOUBLE_EQ(upf->jchi[0],0.5);
	EXPECT_DOUBLE_EQ(upf->jchi[1],1.5);
	EXPECT_DOUBLE_EQ(upf->jchi[2],0.5);
	//PSWFC
	EXPECT_EQ(upf->els[0],"2S");
	EXPECT_EQ(upf->lchi[0],0);
	EXPECT_DOUBLE_EQ(upf->oc[0],2.0);
	EXPECT_EQ(upf->els[1],"2P");
	EXPECT_EQ(upf->lchi[1],1);
	EXPECT_DOUBLE_EQ(upf->oc[1],1.333);
	EXPECT_EQ(upf->els[2],"2P");
	EXPECT_EQ(upf->lchi[2],1);
	EXPECT_DOUBLE_EQ(upf->oc[2],0.667);
	EXPECT_DOUBLE_EQ(upf->chi(0,0),2.0715339166E-12);
	EXPECT_DOUBLE_EQ(upf->chi(2,upf->mesh-1),1.1201306967E-03);
	//NLCC
	EXPECT_DOUBLE_EQ(upf->rho_atc[0],8.7234550809E-01);
	EXPECT_DOUBLE_EQ(upf->rho_atc[upf->mesh-1],0.0);
	ifs.close();
}

TEST_F(ReadPPTest, ReadUPF201MESH2)
{
	std::ifstream ifs;
	// this pp file has gipaw, thus a different header
	ifs.open("./support/Fe.pbe-sp-mt_gipaw.UPF");
	upf->read_pseudo_upf201(ifs);
	EXPECT_EQ(upf->psd,"Fe");
	ifs.close();
}

TEST_F(ReadPPTest, VWR)
{
	std::ifstream ifs;
	// this pp file is a vwr type of pp
	ifs.open("./support/vwr.Si");
	upf->read_pseudo_vwr(ifs);
	EXPECT_EQ(upf->xc_func,"PZ");
	EXPECT_EQ(upf->pp_type,"NC");
	EXPECT_FALSE(upf->tvanp);
	EXPECT_EQ(upf->mesh,1073);
	EXPECT_FALSE(upf->nlcc);
	EXPECT_EQ(upf->psd,"14");
	EXPECT_EQ(upf->zp,4);
	EXPECT_EQ(upf->spd_loc,2);
	EXPECT_FALSE(upf->has_so);
	EXPECT_EQ(upf->iTB_s,1);
	EXPECT_EQ(upf->iTB_p,1);
	EXPECT_EQ(upf->iTB_d,0);
	EXPECT_EQ(upf->nwfc,2);
	EXPECT_DOUBLE_EQ(upf->oc[0],2);
	EXPECT_DOUBLE_EQ(upf->oc[1],2);
	EXPECT_EQ(upf->lchi[0],0);
	EXPECT_EQ(upf->lchi[1],1);
	EXPECT_EQ(upf->els[0],"S");
	EXPECT_EQ(upf->els[1],"P");
	EXPECT_DOUBLE_EQ(upf->r[0],.22270617E-05);
	EXPECT_DOUBLE_EQ(upf->r[upf->mesh-1],.11832572E+03);
	EXPECT_NEAR(upf->rho_at[0],6.18479e-13,1.0e-17);
	EXPECT_NEAR(upf->rho_at[upf->mesh-1],3.46232e-56,1.0e-60);
	EXPECT_EQ(upf->nbeta,1);
	EXPECT_NEAR(upf->beta(0,2),2.67501e-05,1.0e-9);
	EXPECT_EQ(upf->lll[0],0);
	ifs.close();
}

TEST_F(ReadPPTest, BLPS)
{
	std::ifstream ifs;
	// this pp file is a vwr type of pp
	ifs.open("./support/si.lda.lps");
	GlobalV::DFT_FUNCTIONAL="default";
	upf->read_pseudo_blps(ifs);
	EXPECT_FALSE(upf->nlcc);
	EXPECT_FALSE(upf->tvanp);
	EXPECT_FALSE(upf->has_so);
	EXPECT_EQ(upf->nbeta,0);
	EXPECT_EQ(upf->psd,"Silicon");
	EXPECT_EQ(upf->zp,4);
	EXPECT_EQ(upf->lmax,0);
	EXPECT_EQ(upf->mesh,1601);
	EXPECT_EQ(upf->xc_func,"PZ");
	EXPECT_DOUBLE_EQ(upf->r[0],0.0);
	EXPECT_DOUBLE_EQ(upf->r[upf->mesh-1],16.0);
	EXPECT_DOUBLE_EQ(upf->vloc[0],2.4189229665506291*2.0);
	EXPECT_DOUBLE_EQ(upf->vloc[upf->mesh-1],-0.25*2.0);
	EXPECT_DOUBLE_EQ(upf->rho_at[0],0.25);
	EXPECT_DOUBLE_EQ(upf->rho_at[upf->mesh-1],0.25);
	ifs.close();
}

TEST_F(ReadPPTest, SetPseudoType)
{
	std::string pp_address = "./support/Cu_ONCV_PBE-1.0.upf";
	std::string type = "auto";
	upf->set_pseudo_type(pp_address,type);
	EXPECT_EQ(type,"upf201");
	pp_address = "./support/Te.pbe-rrkj.UPF";
	upf->set_pseudo_type(pp_address,type);
	EXPECT_EQ(type,"upf");
}

TEST_F(ReadPPTest, Trim)
{
	std::string tmp_string = "   aaa   \t  bbb\t  ";
	output = upf->trim(tmp_string);
	EXPECT_EQ(output,"aaabbb");
	tmp_string = "   \taaa\tbbb\t   ";
	output = upf->trimend(tmp_string);
	EXPECT_EQ(output,"aaa\tbbb");
}

TEST_F(ReadPPTest, SetEmptyElement)
{
	upf->mesh = 10;
	upf->nbeta = 10;
	upf->vloc = new double[upf->mesh];
	upf->rho_at = new double[upf->mesh];
	upf->dion.create(upf->nbeta,upf->nbeta);
	upf->set_empty_element();
	for(int ir=0;ir<upf->mesh;++ir)
	{
		EXPECT_DOUBLE_EQ(upf->vloc[ir],0.0);
		EXPECT_DOUBLE_EQ(upf->rho_at[ir],0.0);
	}
	for(int i=0;i<upf->nbeta;++i)
	{
		for(int j=0;j<upf->nbeta;++j)
		{
			EXPECT_DOUBLE_EQ(upf->dion(i,j),0.0);
		}
	}
}

TEST_F(ReadPPTest, SetUpfQ)
{
    std::ifstream ifs;
    ifs.open("./support/Al.pbe-sp-van.UPF");
    upf->read_pseudo_upf201(ifs);
    upf->set_upf_q();
    EXPECT_DOUBLE_EQ(upf->qfuncl(0, 0, 0), 0.0);
    EXPECT_DOUBLE_EQ(upf->qfuncl(0, 0, 100), 7.8994151918886213e-06);
    EXPECT_DOUBLE_EQ(upf->qfuncl(0, 1, 0), 0.0);
    EXPECT_DOUBLE_EQ(upf->qfuncl(0, 1, 100), -2.1915710970869145e-05);
    EXPECT_DOUBLE_EQ(upf->qfuncl(0, 2, 0), 0.0);
    EXPECT_DOUBLE_EQ(upf->qfuncl(0, 2, 100), 5.9614487166963409e-05);
    EXPECT_DOUBLE_EQ(upf->qfuncl(1, 0, 0), 0.0);
    EXPECT_DOUBLE_EQ(upf->qfuncl(1, 0, 100), 0.0);
    EXPECT_DOUBLE_EQ(upf->qfuncl(1, 1, 0), 0.0);
    EXPECT_DOUBLE_EQ(upf->qfuncl(1, 1, 100), 0.0);
    EXPECT_DOUBLE_EQ(upf->qfuncl(1, 2, 0), 0.0);
    EXPECT_DOUBLE_EQ(upf->qfuncl(1, 2, 100), 0.0);

    ifs.close();
}

TEST_F(ReadPPTest, SetQfNew)
{
    // Set up input data
    int nqf = 3;
    int mesh = 10;
    int l = 2;
    int n = 1;

    double qfcoef[] = {1.0, 2.0, 3.0};
    double r[mesh];
    double rho[mesh];

    for (int i = 0; i < mesh; ++i)
    {
        r[i] = i + 1; // Assuming some values for r
    }

    // Call the function under test
    upf->setqfnew(nqf, mesh, l, n, qfcoef, r, rho);

    // Validate the output
    for (int ir = 0; ir < mesh; ++ir)
    {
        double rr = r[ir] * r[ir];
        double expectedValue = qfcoef[0];
        for (int iq = 1; iq < nqf; ++iq)
        {
            expectedValue += qfcoef[iq] * pow(rr, iq);
        }
        expectedValue *= pow(r[ir], l + n);
        EXPECT_DOUBLE_EQ(expectedValue, rho[ir]);
    }
}

TEST_F(ReadPPTest, InitReader)
{
	std::string pp_file = "arbitrary";
	std::string type = "auto";
	int info = upf->init_pseudo_reader(pp_file,type);
	EXPECT_EQ(info,1);
	pp_file = "./support/Te.pbe-rrkj.UPF";
	info = upf->init_pseudo_reader(pp_file,type);
	EXPECT_EQ(type,"upf");
	EXPECT_EQ(info,0);
	pp_file = "./support/Cu_ONCV_PBE-1.0.upf";
	info = upf->init_pseudo_reader(pp_file,type);
	EXPECT_EQ(info,2);
	pp_file = "./support/Cu_ONCV_PBE-1.0.upf";
	type = "auto";
	info = upf->init_pseudo_reader(pp_file,type);
	EXPECT_EQ(type,"upf201");
	EXPECT_EQ(info,0);
	pp_file = "./support/vwr.Si";
	type = "vwr";
	info = upf->init_pseudo_reader(pp_file,type);
	EXPECT_EQ(info,0);
	pp_file = "./support/si.lda.lps";
	type = "blps";
	info = upf->init_pseudo_reader(pp_file,type);
	EXPECT_EQ(info,0);
}

TEST_F(ReadPPTest, AverageSimpleReturns)
{
	int ierr;
	double lambda = 1.0;
	// first return
	GlobalV::LSPINORB = 1;
	upf->has_so = 0;
	ierr = upf->average_p(lambda);
	EXPECT_EQ(ierr,1);
	// second return
	upf->has_so = 1;
	ierr = upf->average_p(lambda);
	EXPECT_EQ(ierr,0);
}

TEST_F(ReadPPTest, AverageErrReturns)
{
	int ierr;
	double lambda = 1.0;
	// LSPINORB = 0
	std::ifstream ifs;
	ifs.open("./support/Te.pbe-rrkj.UPF");
	upf->read_pseudo_upf(ifs);
	EXPECT_TRUE(upf->has_so); // has soc info
	GlobalV::LSPINORB = 0;
	ierr = upf->average_p(lambda);
	EXPECT_EQ(upf->nbeta,3);
	EXPECT_EQ(ierr,1);
	// LSPINORB = 1
	ierr = upf->average_p(lambda);
	EXPECT_EQ(ierr,1);
	ifs.close();
}

TEST_F(ReadPPTest, AverageLSPINORB0)
{
	std::ifstream ifs;
	// this is a dojo full-relativisitic pp
	ifs.open("./support/C.upf");
	upf->read_pseudo_upf201(ifs);
	EXPECT_TRUE(upf->has_so); // has soc info
	int ierr;
	double lambda = 1.0;
	// LSPINORB = 0
	GlobalV::LSPINORB = 0;
	ierr = upf->average_p(lambda);
	EXPECT_EQ(ierr,0);
	EXPECT_EQ(upf->nbeta,4);
	EXPECT_FALSE(upf->has_so); // has not soc info,why?
}

TEST_F(ReadPPTest, AverageLSPINORB1)
{
	std::ifstream ifs;
	// this is a dojo full-relativisitic pp
	ifs.open("./support/C.upf");
	upf->read_pseudo_upf201(ifs);
	EXPECT_TRUE(upf->has_so); // has soc info
	int ierr;
	double lambda = 1.1;
	// LSPINORB = 0
	GlobalV::LSPINORB = 1;
	ierr = upf->average_p(lambda);
	EXPECT_EQ(ierr,0);
	EXPECT_EQ(upf->nbeta,6);
	EXPECT_TRUE(upf->has_so); // has soc info
}
#undef private
