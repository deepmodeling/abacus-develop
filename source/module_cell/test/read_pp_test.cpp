#include "gtest/gtest.h"
#include "gmock/gmock.h"
#define private public
#include "module_parameter/parameter.h"
#undef private
#include<memory>
/************************************************
 *  unit test of read_pp
 ***********************************************/

/**
 * - Tested Functions:
 *   - ReadUPF100
 *   - ReadUPF100USPP
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
#include "module_cell/atom_pseudo.h"
#undef private
class ReadPPTest : public testing::Test
{
protected:
	std::string output;
	std::unique_ptr<Pseudopot_upf> read_pp{new Pseudopot_upf};
	std::unique_ptr<Atom_pseudo> upf{new Atom_pseudo};
};

TEST_F(ReadPPTest, ReadUPF100_Coulomb)
{
	std::ifstream ifs;
	ifs.open("./support/Te.pbe-coulomb.UPF");
	read_pp->read_pseudo_upf(ifs, *upf);
	EXPECT_TRUE(upf->vloc_at.empty());
	EXPECT_EQ(read_pp->coulomb_potential, true);
	EXPECT_EQ(upf->tvanp, false);
	EXPECT_EQ(upf->nbeta, 0);
	EXPECT_EQ(upf->lmax, 0);
	EXPECT_EQ(read_pp->lloc, 0);
	ifs.close();
}

TEST_F(ReadPPTest, ReadUPF100)
{
	std::ifstream ifs;
	ifs.open("./support/Te.pbe-rrkj.UPF");
	read_pp->read_pseudo_upf(ifs, *upf);
	EXPECT_TRUE(upf->has_so); // has soc info
	EXPECT_EQ(upf->nv,0); // number of version
	EXPECT_EQ(upf->psd,"Te"); // element label
	EXPECT_EQ(upf->pp_type,"NC"); // pp_type
	EXPECT_FALSE(upf->tvanp); // not ultrasoft
	EXPECT_FALSE(upf->nlcc); // no Nonlinear core correction
	EXPECT_EQ(upf->xc_func,"PBE"); // Exchange-Correlation functional
	EXPECT_EQ(upf->zv,6); // Z valence
	EXPECT_DOUBLE_EQ(upf->etotps,-15.54533017755); // total energy
	EXPECT_DOUBLE_EQ(upf->ecutwfc,0.0); // suggested cutoff for wfc
	EXPECT_DOUBLE_EQ(upf->ecutrho,0.0); // suggested cutoff for rho
	EXPECT_EQ(upf->lmax,2); // max angular momentum component
	EXPECT_EQ(upf->mesh,1191); // Number of points in mesh
	EXPECT_EQ(upf->nchi,3); // Number of wavefunctions
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
	EXPECT_DOUBLE_EQ(upf->vloc_at[0],-5.00890143222E+00); // vloc
	EXPECT_DOUBLE_EQ(upf->vloc_at[upf->mesh-1],-2.37200471955E-01); // vloc
	EXPECT_EQ(upf->lll[0],0); // BETA
    EXPECT_EQ(read_pp->kbeta[0], 957);
    EXPECT_DOUBLE_EQ(upf->betar(0, 0), -1.82560984478E-03);
    EXPECT_DOUBLE_EQ(upf->betar(0, read_pp->kbeta[0] - 1), -1.61398366674E-03);
    EXPECT_EQ(upf->lll[1], 0); // BETA
    EXPECT_EQ(read_pp->kbeta[1], 957);
    EXPECT_DOUBLE_EQ(upf->betar(1, 0), 1.51692136792E-03);
    EXPECT_DOUBLE_EQ(upf->betar(1, read_pp->kbeta[1] - 1), 1.51458412218E-03);
    EXPECT_EQ(upf->lll[2], 2); // BETA
    EXPECT_EQ(read_pp->kbeta[2], 957);
    EXPECT_DOUBLE_EQ(upf->betar(2, 0), -3.10582746893E-13);
    EXPECT_DOUBLE_EQ(upf->betar(2, read_pp->kbeta[2] - 1), -4.17131335030E-04);
    EXPECT_EQ(read_pp->nd,4); // DIJ
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
	EXPECT_DOUBLE_EQ(upf->jchi[2],0.0);
	EXPECT_DOUBLE_EQ(upf->jjj[0],0.0); // jjj
	EXPECT_DOUBLE_EQ(upf->jjj[1],0.0);
	EXPECT_DOUBLE_EQ(upf->jjj[2],0.0);
	//EXPECT_EQ
	ifs.close();
	std::ofstream ofs;
	ofs.open("tmp");
	read_pp->print_pseudo_upf(ofs, *upf);
	ofs.close();
	ifs.open("tmp");
	getline(ifs, output);
	EXPECT_THAT(output, testing::HasSubstr("==== read_pseudo_upf ==="));
	ifs.close();
}

TEST_F(ReadPPTest, ReadUPF100USPP)
{
    std::ifstream ifs;
    ifs.open("./support/fe_pbe_v1.5.uspp.F.UPF");
	read_pp->read_pseudo_upf(ifs, *upf);
    EXPECT_FALSE(upf->has_so);                                        // has soc info
    EXPECT_FALSE(read_pp->q_with_l);                                      // q_with_l
    EXPECT_EQ(upf->nv, 0);                                            // number of version
    EXPECT_EQ(upf->psd, "Fe");                                        // element label
    EXPECT_EQ(upf->pp_type, "US");                                    // pp_type
    EXPECT_TRUE(upf->tvanp);                                          // not ultrasoft
    EXPECT_TRUE(upf->nlcc);                                           // no Nonlinear core correction
    EXPECT_EQ(upf->xc_func, "PBE");                                   // Exchange-Correlation functional
    EXPECT_EQ(upf->zv, 16);                                           // Z valence
    EXPECT_DOUBLE_EQ(upf->etotps, -248.63387366200);                  // total energy
    EXPECT_DOUBLE_EQ(upf->ecutwfc, 0.0);                              // suggested cutoff for wfc
    EXPECT_DOUBLE_EQ(upf->ecutrho, 0.0);                              // suggested cutoff for rho
    EXPECT_EQ(upf->lmax, 2);                                          // max angular momentum component
    EXPECT_EQ(upf->mesh, 861);                                        // Number of points in mesh
    EXPECT_EQ(upf->nchi, 5);                                          // Number of wavefunctions
    EXPECT_EQ(upf->nbeta, 6);                                         // Number of projectors
    EXPECT_EQ(upf->els[0], "3S");                                     // label for i-th atomic orbital
    EXPECT_EQ(upf->els[1], "3P");                                     // label for i-th atomic orbital
    EXPECT_EQ(upf->els[2], "3D");                                     // label for i-th atomic orbital
    EXPECT_EQ(upf->els[3], "4S");                                     // label for i-th atomic orbital
    EXPECT_EQ(upf->els[4], "4P");                                     // label for i-th atomic orbital
    EXPECT_EQ(upf->lchi[0], 0);                                       // angluar momentum of each atomic orbital
    EXPECT_EQ(upf->lchi[1], 1);                                       // angluar momentum of each atomic orbital
    EXPECT_EQ(upf->lchi[2], 2);                                       // angluar momentum of each atomic orbital
    EXPECT_EQ(upf->lchi[3], 0);                                       // angluar momentum of each atomic orbital
    EXPECT_EQ(upf->lchi[4], 1);                                       // angluar momentum of each atomic orbital
    EXPECT_DOUBLE_EQ(upf->oc[0], 2.0);                                // occupation of each atomic orbital
    EXPECT_DOUBLE_EQ(upf->oc[1], 6.0);                                // occupation of each atomic orbital
    EXPECT_DOUBLE_EQ(upf->oc[2], 5.0);                                // occupation of each atomic orbital
    EXPECT_DOUBLE_EQ(upf->oc[3], 2.0);                                // occupation of each atomic orbital
    EXPECT_DOUBLE_EQ(upf->oc[4], 0.0);                                // occupation of each atomic orbital
    EXPECT_DOUBLE_EQ(upf->r[0], 0.00000000000E+00);                   // r
    EXPECT_DOUBLE_EQ(upf->r[upf->mesh - 1], 2.04011054501E+02);       // r
    EXPECT_DOUBLE_EQ(upf->rab[0], 1.61587495219E-06);                 // rab
    EXPECT_DOUBLE_EQ(upf->rab[upf->mesh - 1], 3.45781609895E+00);     // rab
    EXPECT_DOUBLE_EQ(upf->rho_atc[1], 3.48284864470E+00);             // nlcc
    EXPECT_DOUBLE_EQ(upf->rho_atc[upf->mesh - 1], 0.00000000000E+00); // nlcc
    EXPECT_DOUBLE_EQ(upf->vloc_at[0], -5.13189996435E+01);               // vloc
    EXPECT_DOUBLE_EQ(upf->vloc_at[upf->mesh - 1], -1.56854245366E-01);   // vloc
    EXPECT_EQ(upf->lll[0], 0);                                        // BETA
    EXPECT_EQ(read_pp->kbeta[0], 607);
    EXPECT_DOUBLE_EQ(upf->betar(0, 0), 0.00000000000E+00);
    EXPECT_DOUBLE_EQ(upf->betar(0, read_pp->kbeta[0] - 1), 0.00000000000E+00);
    EXPECT_EQ(upf->lll[1], 0); // BETA
    EXPECT_EQ(read_pp->kbeta[1], 607);
    EXPECT_DOUBLE_EQ(upf->betar(1, 0), 0.00000000000E+00);
    EXPECT_DOUBLE_EQ(upf->betar(1, read_pp->kbeta[0] - 1), 0.00000000000E+00);
    EXPECT_EQ(upf->lll[2], 1); // BETA
    EXPECT_EQ(read_pp->kbeta[2], 607);
    EXPECT_DOUBLE_EQ(upf->betar(2, 1), 5.76970159520E-12);
    EXPECT_DOUBLE_EQ(upf->betar(2, read_pp->kbeta[2] - 1), 0.00000000000E+00);
    EXPECT_EQ(upf->lll[5], 2); // BETA
    EXPECT_EQ(read_pp->kbeta[5], 607);
    EXPECT_DOUBLE_EQ(upf->betar(5, 2), -3.30692076673E-11);
    EXPECT_DOUBLE_EQ(upf->betar(5, read_pp->kbeta[2] - 1), 0.00000000000E+00);
    EXPECT_EQ(read_pp->nd, 9); // DIJ
    EXPECT_DOUBLE_EQ(upf->dion(0, 0), -2.44502386412E-02);
    EXPECT_DOUBLE_EQ(upf->dion(0, 1), 1.88646719481E+00);
    EXPECT_DOUBLE_EQ(upf->dion(1, 1), 2.82162386984E+00);
    EXPECT_DOUBLE_EQ(upf->dion(2, 2), 9.10650114165E+00);
    EXPECT_DOUBLE_EQ(upf->dion(2, 3), -1.66542402638E+01);
    EXPECT_DOUBLE_EQ(upf->dion(3, 3), 2.82741263018E+01);
    EXPECT_DOUBLE_EQ(upf->dion(4, 4), 5.48893904730E+01);
    EXPECT_DOUBLE_EQ(upf->dion(4, 5), 6.28094728901E+01);
    EXPECT_DOUBLE_EQ(upf->dion(5, 5), 7.17648258086E+01);
    EXPECT_EQ(read_pp->nqf, 8);                // QIJ
    EXPECT_EQ(upf->nqlc, 5);               // nqlc
    EXPECT_DOUBLE_EQ(read_pp->rinner[0], 1.0); // rinner
    EXPECT_DOUBLE_EQ(read_pp->rinner[1], 1.0);
    EXPECT_DOUBLE_EQ(read_pp->rinner[2], 1.0);
    EXPECT_DOUBLE_EQ(read_pp->rinner[3], 1.0);
    EXPECT_DOUBLE_EQ(read_pp->rinner[4], 1.0);
    EXPECT_DOUBLE_EQ(upf->qqq(0, 0), 9.61156723771E-02);
    EXPECT_DOUBLE_EQ(read_pp->qfunc(0, 1), -2.94880294140E-11);
    EXPECT_DOUBLE_EQ(read_pp->qfunc(0, upf->mesh - 1), 0.00000000000E+00);
    EXPECT_DOUBLE_EQ(read_pp->qfcoef(0, 0, 0, 0), -1.11034753554E+01);
    EXPECT_DOUBLE_EQ(read_pp->qfcoef(0, 0, 4, 7), 0.00000000000E+00);
    EXPECT_DOUBLE_EQ(upf->qqq(0, 1), 6.30706989525E-02);
    EXPECT_DOUBLE_EQ(read_pp->qfunc(1, 1), -9.73254487126E-12);
    EXPECT_DOUBLE_EQ(read_pp->qfunc(0, upf->mesh - 1), 0.00000000000E+00);
    EXPECT_DOUBLE_EQ(read_pp->qfcoef(0, 1, 0, 0), -3.66470985929E+00);
    EXPECT_DOUBLE_EQ(read_pp->qfcoef(0, 1, 4, 7), 0.00000000000E+00);
    EXPECT_DOUBLE_EQ(upf->qqq(5, 5), 4.88172232559E+00);
    EXPECT_DOUBLE_EQ(read_pp->qfunc(20, 1), 1.69461625558E-10);
    EXPECT_DOUBLE_EQ(read_pp->qfunc(0, upf->mesh - 1), 0.00000000000E+00);
    EXPECT_DOUBLE_EQ(read_pp->qfcoef(5, 5, 0, 0), 6.38093836870E+01);
    EXPECT_DOUBLE_EQ(read_pp->qfcoef(5, 5, 4, 7), 8.40128670914E+03);
    EXPECT_DOUBLE_EQ(upf->chi(0, 1), 4.36429934825E-06); // PSWFC
    EXPECT_DOUBLE_EQ(upf->chi(0, upf->mesh - 1), 0.00000000000E+00);
    EXPECT_DOUBLE_EQ(upf->chi(1, 1), 6.46349114028E-12);
    EXPECT_DOUBLE_EQ(upf->chi(1, upf->mesh - 1), 0.00000000000E+00);
    EXPECT_DOUBLE_EQ(upf->chi(2, 1), 7.06492930658E-18);
    EXPECT_DOUBLE_EQ(upf->chi(2, upf->mesh - 1), 0.00000000000E+00);
    EXPECT_DOUBLE_EQ(upf->chi(3, 1), 1.45872860322E-06);
    EXPECT_DOUBLE_EQ(upf->chi(3, upf->mesh - 1), 0.00000000000E+00);
    EXPECT_DOUBLE_EQ(upf->chi(4, 1), 1.27728623256E-12);
    EXPECT_DOUBLE_EQ(upf->chi(4, upf->mesh - 1), 0.00000000000E+00);
    EXPECT_DOUBLE_EQ(upf->rho_at[1], 7.94937989763E-11); // RhoAtom
    EXPECT_DOUBLE_EQ(upf->rho_at[upf->mesh - 1], 0.00000000000E+00);
    // EXPECT_EQ
    ifs.close();
    std::ofstream ofs;
    ofs.open("tmp");
    read_pp->print_pseudo_upf(ofs, *upf);
    ofs.close();
    ifs.open("tmp");
    getline(ifs, output);
    EXPECT_THAT(output, testing::HasSubstr("==== read_pseudo_upf ==="));
    ifs.close();
}

TEST_F(ReadPPTest, ReadUPF201_Coulomb)
{
	std::ifstream ifs;
	ifs.open("./support/Al.pbe-coulomb.UPF");
	read_pp->read_pseudo_upf201(ifs, *upf);
	EXPECT_TRUE(upf->vloc_at.empty());
	EXPECT_EQ(read_pp->coulomb_potential, true);
	EXPECT_EQ(upf->nbeta, 0);
	EXPECT_EQ(upf->lmax, 0);
	EXPECT_EQ(read_pp->lloc, 0);
	ifs.close();
}

TEST_F(ReadPPTest, ReadUPF201)
{
	std::ifstream ifs;
	ifs.open("./support/Cu_ONCV_PBE-1.0.upf");
	read_pp->read_pseudo_upf201(ifs, *upf);
	EXPECT_EQ(upf->psd,"Cu");
	EXPECT_EQ(upf->pp_type,"NC");
	EXPECT_FALSE(upf->has_so);
	EXPECT_FALSE(upf->nlcc);
	EXPECT_EQ(upf->xc_func,"PBE");
	EXPECT_EQ(upf->zv,19);
	EXPECT_DOUBLE_EQ(upf->etotps,-1.82394100797E+02);
	EXPECT_EQ(upf->lmax,2);
	EXPECT_EQ(upf->mesh,601); // mesh -= 1 at line 388 (why? Let's see)
	EXPECT_EQ(upf->nchi,0);
	EXPECT_EQ(upf->nbeta,6);
	EXPECT_DOUBLE_EQ(upf->r[0],0.0);
	EXPECT_DOUBLE_EQ(upf->r[600],6.00);
	EXPECT_DOUBLE_EQ(upf->rab[0],0.01);
	EXPECT_DOUBLE_EQ(upf->rab[600],0.01);
	EXPECT_DOUBLE_EQ(upf->vloc_at[0],-5.3426582174E+01);
	EXPECT_DOUBLE_EQ(upf->vloc_at[600],-6.3333339776E+00);
	EXPECT_EQ(upf->lll[0],0);
    EXPECT_EQ(read_pp->kbeta[0], 196);
    EXPECT_DOUBLE_EQ(upf->betar(0,0),0.0);
	EXPECT_DOUBLE_EQ(upf->betar(0,600),0.0);
	EXPECT_EQ(upf->lll[1],0);
    EXPECT_EQ(read_pp->kbeta[1], 196);
    EXPECT_DOUBLE_EQ(upf->betar(1,0),0.0);
	EXPECT_DOUBLE_EQ(upf->betar(1,600),0.0);
	EXPECT_EQ(upf->lll[2],1);
    EXPECT_EQ(read_pp->kbeta[2], 196);
    EXPECT_DOUBLE_EQ(upf->betar(2,0),0.0);
	EXPECT_DOUBLE_EQ(upf->betar(2,600),0.0);
	EXPECT_EQ(upf->lll[3],1);
    EXPECT_EQ(read_pp->kbeta[3], 196);
    EXPECT_DOUBLE_EQ(upf->betar(3,0),0.0);
	EXPECT_DOUBLE_EQ(upf->betar(3,600),0.0);
	EXPECT_EQ(upf->lll[4],2);
    EXPECT_EQ(read_pp->kbeta[4], 196);
    EXPECT_DOUBLE_EQ(upf->betar(4,0),0.0);
	EXPECT_DOUBLE_EQ(upf->betar(4,600),0.0);
	EXPECT_EQ(upf->lll[5],2);
    EXPECT_EQ(read_pp->kbeta[5], 196);
    EXPECT_DOUBLE_EQ(upf->betar(5,0),0.0);
	EXPECT_DOUBLE_EQ(upf->betar(5,600),0.0);
	EXPECT_DOUBLE_EQ(upf->dion(0,0),-6.6178420255E+00);
	EXPECT_DOUBLE_EQ(upf->dion(5,5),-7.0938557228E+00);
	EXPECT_DOUBLE_EQ(upf->rho_at[0],0.0);
	EXPECT_DOUBLE_EQ(upf->rho_at[600],3.2115793029E-02);
	ifs.close();
}

TEST_F(ReadPPTest, ReadUSPPUPF201)
{
    std::ifstream ifs;
    ifs.open("./support/Al.pbe-sp-van.UPF");
    read_pp->read_pseudo_upf201(ifs, *upf);
    EXPECT_EQ(upf->psd, "Al");
    EXPECT_EQ(upf->pp_type, "US");
    EXPECT_EQ(read_pp->relativistic, "no");
    EXPECT_TRUE(upf->tvanp);
    EXPECT_FALSE(upf->has_so);
    EXPECT_FALSE(upf->nlcc);
    EXPECT_EQ(upf->xc_func, "PBE");
    EXPECT_EQ(upf->zv, 11);
    EXPECT_EQ(upf->nv, 0);
    EXPECT_DOUBLE_EQ(upf->etotps, -1.596432307730e2);
    EXPECT_DOUBLE_EQ(upf->ecutwfc, 0.0);
    EXPECT_DOUBLE_EQ(upf->ecutrho, 0.0);
    EXPECT_EQ(upf->lmax, 2);
    EXPECT_EQ(read_pp->lmax_rho, 0);
    EXPECT_EQ(read_pp->lloc, 0);
    EXPECT_EQ(upf->mesh, 893);
    EXPECT_EQ(upf->nchi, 4);
    EXPECT_EQ(upf->nbeta, 5);
    EXPECT_EQ(upf->kkbeta, 617);
    EXPECT_DOUBLE_EQ(read_pp->rmax, 2.006810756590e2);
    EXPECT_DOUBLE_EQ(read_pp->zmesh, 13.0);
    EXPECT_FALSE(read_pp->q_with_l);
    EXPECT_EQ(read_pp->nqf, 8);
    EXPECT_EQ(upf->nqlc, 5);
    EXPECT_DOUBLE_EQ(upf->r[0], 0.0);
    EXPECT_DOUBLE_EQ(upf->r[892], 2.006810756590000e2);
    EXPECT_DOUBLE_EQ(upf->rab[0], 1.169079443020000e-6);
    EXPECT_DOUBLE_EQ(upf->rab[892], 3.344685763390000e0);
    EXPECT_DOUBLE_EQ(upf->vloc_at[0], 3.456089057550000e0);
    EXPECT_DOUBLE_EQ(upf->vloc_at[892], -1.096266796840000e-1);
	EXPECT_TRUE(upf->rho_atc.empty());
    EXPECT_EQ(upf->lll[0], 0);
    EXPECT_EQ(read_pp->kbeta[0], 617);
    EXPECT_EQ(read_pp->els_beta[0], "2S");
    EXPECT_DOUBLE_EQ(read_pp->rcut[0], 0.0);
    EXPECT_DOUBLE_EQ(read_pp->rcutus[0], 1.4);
    EXPECT_DOUBLE_EQ(upf->betar(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(upf->betar(0, 892), 0.0);
    EXPECT_EQ(upf->lll[1], 0);
    EXPECT_EQ(read_pp->kbeta[1], 617);
    EXPECT_EQ(read_pp->els_beta[1], "2P");
    EXPECT_DOUBLE_EQ(read_pp->rcut[1], 0.0);
    EXPECT_DOUBLE_EQ(read_pp->rcutus[1], 1.42);
    EXPECT_DOUBLE_EQ(upf->betar(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(upf->betar(1, 892), 0.0);
    EXPECT_DOUBLE_EQ(upf->dion(0, 0), -2.182408428460000e2);
    EXPECT_DOUBLE_EQ(upf->dion(4, 4), -3.087562171130000e1);
    EXPECT_DOUBLE_EQ(upf->qqq(0, 0), 3.896280866700000e0);
    EXPECT_DOUBLE_EQ(upf->qqq(4, 4), 7.218068659650000e-1);
    EXPECT_DOUBLE_EQ(read_pp->qfcoef(0, 0, 0, 0), 8.705252055130002e1);
    EXPECT_DOUBLE_EQ(read_pp->qfcoef(4, 4, 4, 7), 9.910935792140002e1);
    EXPECT_DOUBLE_EQ(read_pp->rinner[0], 1.1);
    EXPECT_DOUBLE_EQ(read_pp->rinner[4], 1.1);
    EXPECT_DOUBLE_EQ(read_pp->qfunc(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(read_pp->qfunc(0, 892), 0.0);
    EXPECT_EQ(upf->els[0], "2S");
    EXPECT_EQ(upf->lchi[0], 0);
    EXPECT_EQ(read_pp->nchi[0], 0);
    EXPECT_DOUBLE_EQ(upf->oc[0], 2.0);
    EXPECT_DOUBLE_EQ(read_pp->epseu[0], 0.0);
    EXPECT_DOUBLE_EQ(read_pp->rcut_chi[0], 0.0);
    EXPECT_DOUBLE_EQ(read_pp->rcutus_chi[0], 1.4);
    EXPECT_DOUBLE_EQ(upf->chi(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(upf->chi(0, 892), 0.0);
    EXPECT_DOUBLE_EQ(upf->rho_at[0], 0.0);
    EXPECT_DOUBLE_EQ(upf->rho_at[892], 0.0);
    EXPECT_TRUE(upf->jchi.empty());
	EXPECT_TRUE(upf->jjj.empty());
	EXPECT_TRUE(upf->nn.empty());
    ifs.close();
}

TEST_F(ReadPPTest, HeaderErr2011)
{
	std::ifstream ifs;
	// 1st
	ifs.open("./support/HeaderError1");
	//read_pp->read_pseudo_upf201(ifs, *upf);
	testing::internal::CaptureStdout();
	EXPECT_EXIT(read_pp->read_pseudo_upf201(ifs, *upf),
			::testing::ExitedWithCode(1),"");
	output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Found no PP_HEADER"));
    ifs.close();
}

TEST_F(ReadPPTest, HeaderErr2012)
{
	std::ifstream ifs;
	// 2nd
	ifs.open("./support/HeaderError2");
	//read_pp->read_pseudo_upf201(ifs, *upf);
	testing::internal::CaptureStdout();
	EXPECT_EXIT(read_pp->read_pseudo_upf201(ifs, *upf),
			::testing::ExitedWithCode(1),"");
	output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("SEMI-LOCAL PSEUDOPOTENTIAL IS NOT SUPPORTED"));
    ifs.close();
}

TEST_F(ReadPPTest, HeaderErr2013)
{
	std::ifstream ifs;
	// 3rd
	ifs.open("./support/HeaderError3");
	//read_pp->read_pseudo_upf201(ifs, *upf);
	testing::internal::CaptureStdout();
	EXPECT_EXIT(read_pp->read_pseudo_upf201(ifs, *upf),
			::testing::ExitedWithCode(1),"");
	output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("PAW POTENTIAL IS NOT SUPPORTED"));
    ifs.close();
}

TEST_F(ReadPPTest, HeaderErr2015)
{
    std::ifstream ifs;
    // 4th
    GlobalV::ofs_warning.open("warning.log");
    ifs.open("./support/HeaderError5");
	upf->mesh = 1; // avoid assert(pp.mesh > 0) in line 406 of read_pp_upf201.cpp
    read_pp->read_pseudo_upf201(ifs, *upf);
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
	read_pp->read_pseudo_upf201(ifs, *upf);
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
	EXPECT_EQ(upf->nchi,3);
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
	read_pp->read_pseudo_upf201(ifs, *upf);
	EXPECT_EQ(upf->psd,"Fe");
	ifs.close();
}

TEST_F(ReadPPTest, VWR)
{
	std::ifstream ifs;
	// this pp file is a vwr type of pp
	ifs.open("./support/vwr.Si");
	read_pp->read_pseudo_vwr(ifs, *upf);
	EXPECT_EQ(upf->xc_func,"PZ");
	EXPECT_EQ(upf->pp_type,"NC");
	EXPECT_FALSE(upf->tvanp);
	EXPECT_EQ(upf->mesh,1073);
	EXPECT_FALSE(upf->nlcc);
	EXPECT_EQ(upf->psd,"14");
	EXPECT_EQ(upf->zv,4);
	EXPECT_EQ(read_pp->spd_loc,2);
	EXPECT_FALSE(upf->has_so);
	EXPECT_EQ(read_pp->iTB_s,1);
	EXPECT_EQ(read_pp->iTB_p,1);
	EXPECT_EQ(read_pp->iTB_d,0);
	EXPECT_EQ(upf->nchi,2);
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
	EXPECT_NEAR(upf->betar(0,2),2.67501e-05,1.0e-9);
	EXPECT_EQ(upf->lll[0],0);
	ifs.close();
}

TEST_F(ReadPPTest, BLPS)
{
	std::ifstream ifs;
	// this pp file is a vwr type of pp
	ifs.open("./support/si.lda.lps");
	PARAM.input.dft_functional="default";
	read_pp->read_pseudo_blps(ifs, *upf);
	EXPECT_FALSE(upf->nlcc);
	EXPECT_FALSE(upf->tvanp);
	EXPECT_FALSE(upf->has_so);
	EXPECT_EQ(upf->nbeta,0);
	EXPECT_EQ(upf->psd,"Si");
	EXPECT_EQ(upf->zv,4);
	EXPECT_EQ(upf->lmax,0);
	EXPECT_EQ(upf->mesh,1601);
	EXPECT_EQ(upf->xc_func,"PZ");
	EXPECT_DOUBLE_EQ(upf->r[0],0.0);
	EXPECT_DOUBLE_EQ(upf->r[upf->mesh-1],16.0);
	EXPECT_DOUBLE_EQ(upf->vloc_at[0],2.4189229665506291*2.0);
	EXPECT_DOUBLE_EQ(upf->vloc_at[upf->mesh-1],-0.25*2.0);
	EXPECT_DOUBLE_EQ(upf->rho_at[0],0.25);
	EXPECT_DOUBLE_EQ(upf->rho_at[upf->mesh-1],0.25);
	ifs.close();
}

TEST_F(ReadPPTest, SetPseudoType)
{
	std::string pp_address = "./support/Cu_ONCV_PBE-1.0.upf";
	std::string type = "auto";
	read_pp->set_pseudo_type(pp_address,type);
	EXPECT_EQ(type,"upf201");
	pp_address = "./support/Te.pbe-rrkj.UPF";
	read_pp->set_pseudo_type(pp_address,type);
	EXPECT_EQ(type,"upf");
}

TEST_F(ReadPPTest, Trim)
{
	std::string tmp_string = "   aaa   \t  bbb\t  ";
	output = read_pp->trim(tmp_string);
	EXPECT_EQ(output,"aaabbb");
	tmp_string = "   \taaa\tbbb\t   ";
	output = read_pp->trimend(tmp_string);
	EXPECT_EQ(output,"aaa\tbbb");
}

TEST_F(ReadPPTest, SetEmptyElement)
{
	upf->mesh = 10;
	upf->nbeta = 10;
	upf->vloc_at = std::vector<double>(upf->mesh, 0.0);
	upf->rho_at = std::vector<double>(upf->mesh, 0.0);
	upf->dion.create(upf->nbeta,upf->nbeta);
	read_pp->set_empty_element(*upf);
	for(int ir=0;ir<upf->mesh;++ir)
	{
		EXPECT_DOUBLE_EQ(upf->vloc_at[ir],0.0);
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
    read_pp->read_pseudo_upf201(ifs, *upf);
    read_pp->set_upf_q(*upf);
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
    read_pp->setqfnew(nqf, mesh, l, n, qfcoef, r, rho);

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
	int info = read_pp->init_pseudo_reader(pp_file,type,*upf);
	EXPECT_EQ(info,1);
	pp_file = "./support/Te.pbe-rrkj.UPF";
	info = read_pp->init_pseudo_reader(pp_file,type,*upf);
	EXPECT_EQ(type,"upf");
	EXPECT_EQ(info,0);
	pp_file = "./support/Cu_ONCV_PBE-1.0.upf";
	info = read_pp->init_pseudo_reader(pp_file,type,*upf);
	EXPECT_EQ(info,2);
	pp_file = "./support/Cu_ONCV_PBE-1.0.upf";
	type = "auto";
	info = read_pp->init_pseudo_reader(pp_file,type,*upf);
	EXPECT_EQ(type,"upf201");
	EXPECT_EQ(info,0);
	pp_file = "./support/vwr.Si";
	type = "vwr";
	info = read_pp->init_pseudo_reader(pp_file,type,*upf);
	EXPECT_EQ(info,0);
	pp_file = "./support/si.lda.lps";
	type = "blps";
	info = read_pp->init_pseudo_reader(pp_file,type,*upf);
	EXPECT_EQ(info,0);
}

TEST_F(ReadPPTest, AverageSimpleReturns)
{
	int ierr;
	double lambda = 1.0;
	// first return
	PARAM.input.lspinorb = 1;
	upf->has_so = 0;
	ierr = read_pp->average_p(lambda, *upf);
	EXPECT_EQ(ierr,1);
	// second return
	upf->has_so = 1;
	ierr = read_pp->average_p(lambda, *upf);
	EXPECT_EQ(ierr,0);
    upf->has_so = 1;
    upf->tvanp = 1;
    ierr = read_pp->average_p(lambda, *upf);
    EXPECT_EQ(ierr, 1);
}

TEST_F(ReadPPTest, AverageErrReturns)
{
	int ierr;
	double lambda = 1.0;
	// LSPINORB = 0
	std::ifstream ifs;
	ifs.open("./support/Te.pbe-rrkj.UPF");
	read_pp->read_pseudo_upf(ifs, *upf);
	EXPECT_TRUE(upf->has_so); // has soc info
	PARAM.input.lspinorb = 0;
	ierr = read_pp->average_p(lambda, *upf);
	EXPECT_EQ(upf->nbeta,3);
	EXPECT_EQ(ierr,1);
	// LSPINORB = 1
	ierr = read_pp->average_p(lambda, *upf);
	EXPECT_EQ(ierr,1);
	ifs.close();
}

TEST_F(ReadPPTest, AverageLSPINORB0)
{
	std::ifstream ifs;
	// this is a dojo full-relativisitic pp
	ifs.open("./support/C.upf");
	read_pp->read_pseudo_upf201(ifs, *upf);
	EXPECT_TRUE(upf->has_so); // has soc info
	int ierr;
	double lambda = 1.0;
	// LSPINORB = 0
	PARAM.input.lspinorb = 0;
	ierr = read_pp->average_p(lambda, *upf);
	EXPECT_EQ(ierr,0);
	EXPECT_EQ(upf->nbeta,4);
	EXPECT_FALSE(upf->has_so); // has not soc info,why?
}

TEST_F(ReadPPTest, AverageLSPINORB1)
{
	std::ifstream ifs;
	// this is a dojo full-relativisitic pp
	ifs.open("./support/C.upf");
	read_pp->read_pseudo_upf201(ifs, *upf);
	EXPECT_TRUE(upf->has_so); // has soc info
	int ierr;
	double lambda = 1.1;
	// LSPINORB = 0
	PARAM.input.lspinorb = 1;
	ierr = read_pp->average_p(lambda, *upf);
	EXPECT_EQ(ierr,0);
	EXPECT_EQ(upf->nbeta,6);
	EXPECT_TRUE(upf->has_so); // has soc info
}
