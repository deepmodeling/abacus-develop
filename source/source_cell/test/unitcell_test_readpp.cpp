#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "memory"
#include "source_base/global_variable.h"
#include "source_base/mathzone.h"
#include "source_cell/check_atomic_stru.h"
#include "source_cell/unitcell.h"
#include "source_cell/cal_nelec_nband.h"
#include "source_estate/read_pseudo.h"
#include <valarray>
#include <vector>
#include "string.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include "prepare_unitcell.h"


Magnetism::Magnetism() {
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}
Magnetism::~Magnetism() { }
/************************************************
 *  unit test of class UnitCell
 ***********************************************/

/**
 * - Tested Functions:
 *   - ReadCellPPWarning1
 *     - read_cell_pseudopots(): error when average the pseudopotential:
 * error_ap
 *   - ReadCellPPWarning2
 *     - read_cell_pseudopots(): Couldn't find pseudopotential file:: error == 1
 *   - ReadCellPPWarning3
 *     - read_cell_pseudopots(): Pseudopotential data do not match: error ==2
 *     - error==3 is currently difficult to reach in read_pseudo_vwr
 *   - ReadCellPPWarning4
 *     - read_cell_pseudopots(): dft_functional from INPUT does not match that
 * in pseudopot file
 *   - ReadCellPPWarning5
 *     - read_cell_pseudopots(): Unknown pseudopotential type
 *   - ReadCellPP
 *     - read_cell_pseudopots(): read pp files with flag_empty_element set
 *   - CalMeshx
 *     - cal_meshx(): calculate max mesh info from atomic pseudo potential file
 *   - CalNatomwfc1
 *     - cal_natomwfc(): calculate total number of atomic orbitals in pseudo
 * potential file
 *     - NSPIN != 4
 *     - this corresponds to number_of_wfc, PP_CHI in pp file, and
 * atoms[it].ncpp.lchi[ncpp.nchi]
 *     - setup the total number of PAOs: pseudopotential atomic orbitals
 *   - CalNatomwfc2
 *     - cal_natomwfc(): calculate total number of atomic orbitals in pseudo
 * potential file
 *     - NSPIN ==4, has_so = false
 *   - CalNatomwfc3
 *     - cal_natomwfc(): calculate total number of atomic orbitals in pseudo
 * potential file
 *     - NSPIN ==4, has_so = true
 *   - CalNwfc1
 *     - cal_nwfc(): calcuate the total number of local basis: NSPIN != 4
 *     - this corresponds to number_of_proj, PP_BETA in pp file, and
 * atoms[it].l_nchi[nw], nw from orb file
 *     - setup nlocal parameter
 *     - interfaces initialed in this function:
 *       - itia2iat
 *       - iat2iwt
 *       - itiaiw2iwt
 *       - iwt2iat
 *       - iwt2iw
 *   - CalNwfc2
 *     - cal_nwfc(): calcuate the total number of local basis: NSPIN == 4
 *   - CheckStructure
 *     - check_atomic_stru(): check if too atoms are two close
 *   - ReadPseudoWarning1
 *     - read_pseudo(): All DFT functional must consistent.
 *   - ReadPseudoWarning2
 *     - read_pseudo(): number valence electrons > corresponding minimum
 * possible of an element
 *   - CalNelec: UnitCell::cal_nelec
 *     - calculate the total number of valence electrons from psp files
 *   - CalNbands: unitcell::cal_nbands()
 *     - calculate the number of bands
 */

class UcellTest : public ::testing::Test {
  protected:
    UcellTestPrepare utp = UcellTestLib["C1H2-Read"];
    std::unique_ptr<UnitCell> ucell;
    std::ofstream ofs;
    std::string pp_dir;
    std::string output;
    void SetUp() {
        ofs.open("running.log");
        ucell = utp.SetUcellInfo();
        pp_dir = "./support/";
    }
    void TearDown() { ofs.close(); }
};

using UcellDeathTest = UcellTest;

TEST_F(UcellDeathTest, ReadCellPPWarning1) {
    const bool lspinorb = true;
    const double pseudo_rcut = 15.0;
    const double soc_lambda = 0.0;
    ucell->pseudo_fn[1] = "H_sr.upf";
    testing::internal::CaptureStdout();
    const std::string global_out_dir = "./";
    const std::string dft_functional = "default";
    EXPECT_EXIT(elecstate::read_cell_pseudopots(pp_dir, ofs, *ucell, global_out_dir, dft_functional, lspinorb, pseudo_rcut, soc_lambda),
                ::testing::ExitedWithCode(1),
                "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,
                testing::HasSubstr("error when average the pseudopotential."));
}

TEST_F(UcellDeathTest, ReadCellPPWarning2) {
    const bool lspinorb = false;
    const double pseudo_rcut = 15.0;
    const double soc_lambda = 0.0;
    pp_dir = "./arbitrary/";
    testing::internal::CaptureStdout();
    const std::string global_out_dir = "./";
    const std::string dft_functional = "default";
    EXPECT_EXIT(elecstate::read_cell_pseudopots(pp_dir, ofs, *ucell, global_out_dir, dft_functional, lspinorb, pseudo_rcut, soc_lambda),
                ::testing::ExitedWithCode(1),
                "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,
                testing::HasSubstr("Couldn't find pseudopotential file"));
}

TEST_F(UcellDeathTest, ReadCellPPWarning3) {
    const bool lspinorb = false;
    const double pseudo_rcut = 15.0;
    const double soc_lambda = 0.0;
    ucell->pseudo_type[0] = "upf";
    testing::internal::CaptureStdout();
    const std::string global_out_dir = "./";
    const std::string dft_functional = "default";
    EXPECT_EXIT(elecstate::read_cell_pseudopots(pp_dir, ofs, *ucell, global_out_dir, dft_functional, lspinorb, pseudo_rcut, soc_lambda),
                ::testing::ExitedWithCode(1),
                "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,
                testing::HasSubstr("Pseudopotential data do not match."));
}

TEST_F(UcellDeathTest, ReadCellPPWarning4) {
    const bool lspinorb = false;
    const double pseudo_rcut = 15.0;
    const double soc_lambda = 0.0;
    const std::string dft_functional = "LDA";
    testing::internal::CaptureStdout();
    const std::string global_out_dir = "./";
    EXPECT_NO_THROW(elecstate::read_cell_pseudopots(pp_dir, ofs, *ucell, global_out_dir, dft_functional, lspinorb, pseudo_rcut, soc_lambda));
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("DFT FUNC. (PSEUDO)   : PBE"));
    EXPECT_THAT(output, testing::HasSubstr("DFT FUNC. (SET TO)   : LDA"));
}

TEST_F(UcellDeathTest, ReadCellPPWarning5) {
    const bool lspinorb = false;
    const double pseudo_rcut = 15.0;
    const double soc_lambda = 0.0;
    ucell->pseudo_type[0] = "upf0000";
    testing::internal::CaptureStdout();
    const std::string global_out_dir = "./";
    const std::string dft_functional = "default";
    EXPECT_EXIT(elecstate::read_cell_pseudopots(pp_dir, ofs, *ucell, global_out_dir, dft_functional, lspinorb, pseudo_rcut, soc_lambda),
                ::testing::ExitedWithCode(1),
                "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Unknown pseudopotential type."));
}

TEST_F(UcellTest, ReadCellPP) {
    const bool lspinorb = false;
    const double pseudo_rcut = 15.0;
    const double soc_lambda = 0.0;
    ucell->atoms[1].flag_empty_element = true;
    const std::string global_out_dir = "./";
    const std::string dft_functional = "default";
    elecstate::read_cell_pseudopots(pp_dir, ofs, *ucell, global_out_dir, dft_functional, lspinorb, pseudo_rcut, soc_lambda);
    EXPECT_EQ(ucell->atoms[0].ncpp.pp_type, "NC");
    EXPECT_FALSE(ucell->atoms[0].ncpp.has_so); // becomes false in average_p
    EXPECT_FALSE(ucell->atoms[1].ncpp.has_so);
    EXPECT_EQ(ucell->atoms[0].ncpp.nchi, 2); // 3=>2 in average_p
    EXPECT_EQ(ucell->atoms[1].ncpp.nchi, 1);
    ofs.close();
    std::ifstream ifs;
    ifs.open("running.log");
    std::string str((std::istreambuf_iterator<char>(ifs)),
                    std::istreambuf_iterator<char>());
    EXPECT_THAT(str,
                testing::HasSubstr("Pseudopotential file = C.upf"));
    EXPECT_THAT(str, testing::HasSubstr("Pseudopotential type = NC"));
    EXPECT_THAT(str,
                testing::HasSubstr("Exchange-correlation functional = PBE"));
    EXPECT_THAT(str, testing::HasSubstr("Valence electrons = 4"));
    EXPECT_THAT(str,
                testing::HasSubstr("Pseudopotential file = H.upf"));
    EXPECT_THAT(str, testing::HasSubstr("Valence electrons = 0"));
}

TEST_F(UcellTest, CalMeshx) {
    const bool lspinorb = false;
    const double pseudo_rcut = 15.0;
    const double soc_lambda = 0.0;
    const std::string global_out_dir = "./";
    const std::string dft_functional = "default";
    elecstate::read_cell_pseudopots(pp_dir, ofs, *ucell, global_out_dir, dft_functional, lspinorb, pseudo_rcut, soc_lambda);
    elecstate::cal_meshx(ucell->meshx,ucell->atoms,ucell->ntype);
    EXPECT_EQ(ucell->atoms[0].ncpp.msh, 1247);
    EXPECT_EQ(ucell->atoms[1].ncpp.msh, 1165);
    EXPECT_EQ(ucell->meshx, 1247);
}

TEST_F(UcellTest, CalNatomwfc1) {
    const bool lspinorb = false;
    const double pseudo_rcut = 15.0;
    const double soc_lambda = 0.0;
    const std::string global_out_dir = "./";
    const std::string dft_functional = "default";
    const int nspin = 1;
    elecstate::read_cell_pseudopots(pp_dir, ofs, *ucell, global_out_dir, dft_functional, lspinorb, pseudo_rcut, soc_lambda);
    EXPECT_FALSE(ucell->atoms[0].ncpp.has_so);
    EXPECT_FALSE(ucell->atoms[1].ncpp.has_so);
    elecstate::cal_natomwfc(ofs,ucell->natomwfc,ucell->ntype,ucell->atoms,nspin);
    EXPECT_EQ(ucell->atoms[0].ncpp.nchi, 2);
    EXPECT_EQ(ucell->atoms[1].ncpp.nchi, 1);
    EXPECT_EQ(ucell->atoms[0].na, 1);
    EXPECT_EQ(ucell->atoms[1].na, 2);
    EXPECT_EQ(ucell->natomwfc, (1 + 3) * 1 + 1 * 2);
}

TEST_F(UcellTest, CalNatomwfc2) {
    const bool lspinorb = false;
    const double pseudo_rcut = 15.0;
    const double soc_lambda = 0.0;
    const int nspin = 4;
    const std::string global_out_dir = "./";
    const std::string dft_functional = "default";
    elecstate::read_cell_pseudopots(pp_dir, ofs, *ucell, global_out_dir, dft_functional, lspinorb, pseudo_rcut, soc_lambda);
    EXPECT_FALSE(ucell->atoms[0].ncpp.has_so);
    EXPECT_FALSE(ucell->atoms[1].ncpp.has_so);
    elecstate::cal_natomwfc(ofs,ucell->natomwfc,ucell->ntype,ucell->atoms,nspin);
    EXPECT_EQ(ucell->atoms[0].ncpp.nchi, 2);
    EXPECT_EQ(ucell->atoms[1].ncpp.nchi, 1);
    EXPECT_EQ(ucell->atoms[0].na, 1);
    EXPECT_EQ(ucell->atoms[1].na, 2);
    EXPECT_EQ(ucell->natomwfc, ((1 + 3) * 1 + 1 * 2) * 2);
}

TEST_F(UcellTest, CalNatomwfc3) {
    const bool lspinorb = true;
    const double pseudo_rcut = 15.0;
    const double soc_lambda = 0.0;
    const int nspin = 4;
    const std::string global_out_dir = "./";
    const std::string dft_functional = "default";
    elecstate::read_cell_pseudopots(pp_dir, ofs, *ucell, global_out_dir, dft_functional, lspinorb, pseudo_rcut, soc_lambda);
    EXPECT_TRUE(ucell->atoms[0].ncpp.has_so);
    EXPECT_TRUE(ucell->atoms[1].ncpp.has_so);
    elecstate::cal_natomwfc(ofs,ucell->natomwfc,ucell->ntype,ucell->atoms,nspin);
    EXPECT_EQ(ucell->atoms[0].ncpp.nchi, 3);
    EXPECT_EQ(ucell->atoms[1].ncpp.nchi, 1);
    EXPECT_EQ(ucell->atoms[0].na, 1);
    EXPECT_EQ(ucell->atoms[1].na, 2);
    EXPECT_EQ(ucell->natomwfc,
              ((2 * 0 + 2) + (2 * 1 + 2) + (2 * 1)) * 1 + (2 * 0 + 2) * 2);
}

TEST_F(UcellTest, CalNwfc1) {
    const bool lspinorb = false;
    const double pseudo_rcut = 15.0;
    const double soc_lambda = 0.0;
    const std::string global_out_dir = "./";
    const std::string dft_functional = "default";
    const int nspin = 1;
    const int nlocal = 27;
    const int npol = 1;
    const std::string basis_type = "pw";
    const std::string esolver_type = "ksdft";
    const std::string init_wfc = "";
    const int nbands = 6;
    elecstate::read_cell_pseudopots(pp_dir, ofs, *ucell, global_out_dir, dft_functional, lspinorb, pseudo_rcut, soc_lambda);
    EXPECT_FALSE(ucell->atoms[0].ncpp.has_so);
    EXPECT_FALSE(ucell->atoms[1].ncpp.has_so);
    elecstate::cal_nwfc(ofs,*ucell,ucell->atoms, nspin, nlocal, npol, basis_type, esolver_type, init_wfc, nbands);
    EXPECT_EQ(ucell->atoms[0].iw2l[8], 2);
    EXPECT_EQ(ucell->atoms[0].iw2n[8], 0);
    EXPECT_EQ(ucell->atoms[0].iw2m[8], 4);
    EXPECT_EQ(ucell->atoms[1].iw2l[8], 2);
    EXPECT_EQ(ucell->atoms[1].iw2n[8], 0);
    EXPECT_EQ(ucell->atoms[1].iw2m[8], 4);
    EXPECT_EQ(ucell->atoms[1].iw2_ylm[8], 8);
    // here is the default table for pw basis calculation
    //  nw = 1*1 + 3*1 + 5*1 = 9
    //     L N m  L*L+m
    //  0  0 0 0    0
    //  1  1 0 0    0
    //  2  1 0 1    2
    //  3  1 0 2    2
    //  4  2 0 0    4
    //  5  2 0 1    5
    //  6  2 0 2    6
    //  7  2 0 3    7
    //  8  2 0 4    8
    EXPECT_EQ(ucell->atoms[0].na, 1);
    EXPECT_EQ(ucell->atoms[1].na, 2);
    EXPECT_EQ(ucell->namax, 2);
    EXPECT_EQ(ucell->atoms[0].nw, 9);
    EXPECT_EQ(ucell->atoms[1].nw, 9);
    EXPECT_EQ(ucell->nwmax, 9);
    // check itia2iat
    EXPECT_EQ(ucell->itia2iat.getSize(), 4);
    EXPECT_EQ(ucell->itia2iat(0, 0), 0);
    EXPECT_EQ(ucell->itia2iat(0, 1), 0);
    EXPECT_EQ(ucell->itia2iat(1, 0), 1);
    EXPECT_EQ(ucell->itia2iat(1, 1), 2);
    // check iat2iwt
    EXPECT_EQ(ucell->get_npol(), 1);
    EXPECT_EQ(ucell->get_iat2iwt()[0], 0);
    EXPECT_EQ(ucell->get_iat2iwt()[1], 9);
    EXPECT_EQ(ucell->get_iat2iwt()[2], 18);
    // check itiaiw2iwt
    EXPECT_EQ(ucell->itiaiw2iwt(0, 0, 0), 0);
    EXPECT_EQ(ucell->itiaiw2iwt(0, 0, 1), 1);
    EXPECT_EQ(ucell->itiaiw2iwt(0, 0, 8), 8);
    EXPECT_EQ(ucell->itiaiw2iwt(1, 0, 0), 9);
    EXPECT_EQ(ucell->itiaiw2iwt(1, 1, 0), 18);
    // check itia2iat
    EXPECT_EQ(ucell->itia2iat.getSize(), 4);
    EXPECT_EQ(ucell->itia2iat(0, 0), 0);
    EXPECT_EQ(ucell->itia2iat(0, 1), 0);
    EXPECT_EQ(ucell->itia2iat(1, 0), 1);
    EXPECT_EQ(ucell->itia2iat(1, 1), 2);
    // check iwt2iat
    EXPECT_EQ(ucell->iwt2iat[0], 0);
    EXPECT_EQ(ucell->iwt2iat[10], 1);
    EXPECT_EQ(ucell->iwt2iat[20], 2);
    // check iwt2iw
    EXPECT_EQ(ucell->iwt2iw[0], 0);
    EXPECT_EQ(ucell->iwt2iw[10], 1);
    EXPECT_EQ(ucell->iwt2iw[20], 2);
}

TEST_F(UcellTest, CalNwfc2) {
    const bool lspinorb = false;
    const double pseudo_rcut = 15.0;
    const double soc_lambda = 0.0;
    const int nspin = 4;
    const int nlocal = 54;
    const int npol = 2;
    const std::string basis_type = "lcao";
    const std::string esolver_type = "ksdft";
    const std::string init_wfc = "";
    const int nbands = 6;
    const std::string global_out_dir = "./";
    const std::string dft_functional = "default";
    elecstate::read_cell_pseudopots(pp_dir, ofs, *ucell, global_out_dir, dft_functional, lspinorb, pseudo_rcut, soc_lambda);
    EXPECT_FALSE(ucell->atoms[0].ncpp.has_so);
    EXPECT_FALSE(ucell->atoms[1].ncpp.has_so);
    EXPECT_NO_THROW(elecstate::cal_nwfc(ofs,*ucell,ucell->atoms, nspin, nlocal, npol, basis_type, esolver_type, init_wfc, nbands));
}

TEST_F(UcellDeathTest, CheckStructure) {
    const bool lspinorb = false;
    const double pseudo_rcut = 15.0;
    const double soc_lambda = 0.0;
    const std::string global_out_dir = "./";
    const std::string dft_functional = "default";
    elecstate::read_cell_pseudopots(pp_dir, ofs, *ucell, global_out_dir, dft_functional, lspinorb, pseudo_rcut, soc_lambda);
    EXPECT_FALSE(ucell->atoms[0].ncpp.has_so);
    EXPECT_FALSE(ucell->atoms[1].ncpp.has_so);
    // trial 1
    
    testing::internal::CaptureStdout();
    double factor = 0.2;
    ucell->set_iat2itia();
    EXPECT_NO_THROW(unitcell::check_atomic_stru(*ucell, factor));
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,testing::HasSubstr("WARNING: Some atoms are too close!!!"));
    // trial 2
    GlobalV::ofs_running.open("CheckStructure2.txt");
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    factor = 0.4;
    EXPECT_EXIT(unitcell::check_atomic_stru(*ucell, factor),
                ::testing::ExitedWithCode(1),
                "");
    std::ifstream ifs("CheckStructure2.txt");
    if (ifs.is_open()) 
    {
        std::string line;
        while (std::getline(ifs, line)) {
            output+=line;
        }
    }
    EXPECT_THAT(output, testing::HasSubstr("The structure is unreasonable!"));
    GlobalV::ofs_running.open("running.log");
    // trial 3
    ucell->atoms[0].label = "arbitrary";
    testing::internal::CaptureStdout();
    factor = 0.2;
    EXPECT_NO_THROW(unitcell::check_atomic_stru(*ucell, factor));
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,testing::HasSubstr("Notice: symbol 'arbitrary' is not an element "
                           "symbol!!!! set the covalent radius to be 0."));
    // trial 4
    ucell->atoms[0].label = "Fe1";
    testing::internal::CaptureStdout();
    factor = 0.2;
    EXPECT_NO_THROW(unitcell::check_atomic_stru(*ucell, factor));
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,testing::HasSubstr("WARNING: Some atoms are too close!!!"));
}

TEST_F(UcellDeathTest, ReadPseudoWarning1) {
    const std::string pseudo_dir = pp_dir;
    const std::string global_out_dir = "./";
    const bool out_element_info = true;
    const std::string dft_functional = "default";
    const bool lspinorb = false;
    const double pseudo_rcut = 15.0;
    const double soc_lambda = 0.0;
    const int nspin = 1;
    const int npol = 1;
    const std::string basis_type = "pw";
    const std::string esolver_type = "ksdft";
    const std::string init_wfc = "";
    const int nbands = 6;
    const bool two_fermi = false;
    const double nelec_delta = 0.0;
    const std::string smearing_method = "none";
    const std::string ks_solver = "genelpa";
    const int bndpar = 1;
    ucell->pseudo_fn[1] = "H_sr_lda.upf";
    testing::internal::CaptureStdout();
    const double nelec = 0.0;
    EXPECT_EXIT(elecstate::read_pseudo(ofs, *ucell, pseudo_dir, global_out_dir, out_element_info, dft_functional, lspinorb, pseudo_rcut, soc_lambda, nspin, npol, basis_type, esolver_type, init_wfc, nbands, two_fermi, nelec_delta, smearing_method, ks_solver, bndpar, nelec), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,
                testing::HasSubstr("All DFT functional must consistent."));
}

TEST_F(UcellDeathTest, ReadPseudoWarning2) {
    const std::string pseudo_dir = pp_dir;
    const std::string global_out_dir = "./";
    const bool out_element_info = true;
    const std::string dft_functional = "default";
    const bool lspinorb = false;
    const double pseudo_rcut = 15.0;
    const double soc_lambda = 0.0;
    const int nspin = 1;
    const int npol = 1;
    const std::string basis_type = "pw";
    const std::string esolver_type = "ksdft";
    const std::string init_wfc = "";
    const int nbands = 6;
    const bool two_fermi = false;
    const double nelec_delta = 0.0;
    const std::string smearing_method = "none";
    const std::string ks_solver = "genelpa";
    const int bndpar = 1;
    ucell->pseudo_fn[0] = "Al_ONCV_PBE-1.0.upf";
    testing::internal::CaptureStdout();
    const double nelec = 0.0;
    EXPECT_NO_THROW(elecstate::read_pseudo(ofs, *ucell, pseudo_dir, global_out_dir, out_element_info, dft_functional, lspinorb, pseudo_rcut, soc_lambda, nspin, npol, basis_type, esolver_type, init_wfc, nbands, two_fermi, nelec_delta, smearing_method, ks_solver, bndpar, nelec));
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(
        output,
        testing::HasSubstr("Warning: the number of valence electrons in "
                           "pseudopotential > 3 for Al: [Ne] 3s2 3p1"));
}

TEST_F(UcellTest, CalNelec) {
    const bool lspinorb = false;
    const double pseudo_rcut = 15.0;
    const double soc_lambda = 0.0;
    const std::string global_out_dir = "./";
    const std::string dft_functional = "default";
    elecstate::read_cell_pseudopots(pp_dir, ofs, *ucell, global_out_dir, dft_functional, lspinorb, pseudo_rcut, soc_lambda);
    EXPECT_EQ(4, ucell->atoms[0].ncpp.zv);
    EXPECT_EQ(1, ucell->atoms[1].ncpp.zv);
    EXPECT_EQ(1, ucell->atoms[0].na);
    EXPECT_EQ(2, ucell->atoms[1].na);
    double nelec = 0;
    const double nelec_delta = 0.0;
    unitcell::cal_nelec(ucell->atoms, ucell->ntype, nelec, nelec_delta);
    EXPECT_DOUBLE_EQ(6, nelec);
}

TEST_F(UcellTest, CalNbands)
{
    const int nelec = 10;
    const int nlocal = 6;
    const int nbands_in = 6;
    int nbands = nbands_in;
    const std::string esolver_type = "ksdft";
    const bool lspinorb = false;
    const int nspin = 1;
    const std::string basis_type = "pw";
    const std::string smearing_method = "fixed";
    std::vector<double> nelec_spin(2, 5.0);
    unitcell::cal_nbands(nelec, nlocal, nelec_spin, nbands, esolver_type, lspinorb, nspin, basis_type, smearing_method);
    EXPECT_EQ(nbands, 6);
}

TEST_F(UcellTest, CalNbandsFractionElec)
{
    const int nelec = 9;
    const int nlocal = 6;
    const int nbands_in = 6;
    int nbands = nbands_in;
    const std::string esolver_type = "ksdft";
    const bool lspinorb = false;
    const int nspin = 1;
    const std::string basis_type = "pw";
    const std::string smearing_method = "fixed";
    std::vector<double> nelec_spin(2, 5.0);
    unitcell::cal_nbands(nelec, nlocal, nelec_spin, nbands, esolver_type, lspinorb, nspin, basis_type, smearing_method);
    EXPECT_EQ(nbands, 6);
}

TEST_F(UcellTest, CalNbandsSOC)
{
    const int nelec = 10;
    const int nlocal = 6;
    int nbands = 0;
    const std::string esolver_type = "ksdft";
    const bool lspinorb = true;
    const int nspin = 1;
    const std::string basis_type = "pw";
    const std::string smearing_method = "fixed";
    std::vector<double> nelec_spin(2, 5.0);
    unitcell::cal_nbands(nelec, nlocal, nelec_spin, nbands, esolver_type, lspinorb, nspin, basis_type, smearing_method);
    EXPECT_EQ(nbands, 20);
}

TEST_F(UcellTest, CalNbandsSDFT)
{
    const int nelec = 10;
    const int nlocal = 6;
    const int nbands_in = 6;
    int nbands = nbands_in;
    const std::string esolver_type = "sdft";
    const bool lspinorb = false;
    const int nspin = 1;
    const std::string basis_type = "pw";
    const std::string smearing_method = "fixed";
    std::vector<double> nelec_spin(2, 5.0);
    EXPECT_NO_THROW(unitcell::cal_nbands(nelec, nlocal, nelec_spin, nbands, esolver_type, lspinorb, nspin, basis_type, smearing_method));
}

TEST_F(UcellTest, CalNbandsLCAO)
{
    const int nelec = 10;
    const int nlocal = 6;
    const int nbands_in = 6;
    int nbands = nbands_in;
    const std::string esolver_type = "ksdft";
    const bool lspinorb = false;
    const int nspin = 1;
    const std::string basis_type = "lcao";
    const std::string smearing_method = "fixed";
    std::vector<double> nelec_spin(2, 5.0);
    EXPECT_NO_THROW(unitcell::cal_nbands(nelec, nlocal, nelec_spin, nbands, esolver_type, lspinorb, nspin, basis_type, smearing_method));
}

TEST_F(UcellTest, CalNbandsLCAOINPW)
{
    const int nelec = 10;
    const int nlocal = 5;
    const int nbands_in = 6;
    int nbands = nbands_in;
    const std::string esolver_type = "ksdft";
    const bool lspinorb = false;
    const int nspin = 1;
    const std::string basis_type = "lcao_in_pw";
    const std::string smearing_method = "fixed";
    std::vector<double> nelec_spin(2, 5.0);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::cal_nbands(nelec, nlocal, nelec_spin, nbands, esolver_type, lspinorb, nspin, basis_type, smearing_method), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Number of basis (NLOCAL) < Number of electronic states (NBANDS)"));
}

TEST_F(UcellTest, CalNbandsWarning1)
{
    const int nelec = 10;
    const int nlocal = 6;
    int nbands = 4;
    const std::string esolver_type = "ksdft";
    const bool lspinorb = false;
    const int nspin = 1;
    const std::string basis_type = "pw";
    const std::string smearing_method = "fixed";
    std::vector<double> nelec_spin(2, 5.0);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::cal_nbands(nelec, nlocal, nelec_spin, nbands, esolver_type, lspinorb, nspin, basis_type, smearing_method), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Too few bands!"));
}

TEST_F(UcellTest, CalNbandsWarning2)
{
    const int nelec = 10;
    const int nlocal = 6;
    const int nbands_in = 6;
    int nbands = nbands_in;
    const std::string esolver_type = "ksdft";
    const bool lspinorb = false;
    const int nspin = 2;
    const std::string basis_type = "pw";
    const std::string smearing_method = "fixed";
    std::vector<double> nelec_spin(2);
    nelec_spin[0] = 7.0;
    nelec_spin[1] = 3.0;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::cal_nbands(nelec, nlocal, nelec_spin, nbands, esolver_type, lspinorb, nspin, basis_type, smearing_method), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Too few spin up bands!"));
}

TEST_F(UcellTest, CalNbandsWarning3)
{
    const int nelec = 10;
    const int nlocal = 6;
    const int nbands_in = 6;
    int nbands = nbands_in;
    const std::string esolver_type = "ksdft";
    const bool lspinorb = false;
    const int nspin = 2;
    const std::string basis_type = "pw";
    const std::string smearing_method = "fixed";
    std::vector<double> nelec_spin(2);
    nelec_spin[0] = 3.0;
    nelec_spin[1] = 7.0;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::cal_nbands(nelec, nlocal, nelec_spin, nbands, esolver_type, lspinorb, nspin, basis_type, smearing_method), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Too few spin down bands!"));
}

TEST_F(UcellTest, CalNbandsSpin1)
{
    const int nelec = 10;
    const int nlocal = 6;
    int nbands = 0;
    const std::string esolver_type = "ksdft";
    const bool lspinorb = false;
    const int nspin = 1;
    const std::string basis_type = "pw";
    const std::string smearing_method = "fixed";
    std::vector<double> nelec_spin(2, 5.0);
    unitcell::cal_nbands(nelec, nlocal, nelec_spin, nbands, esolver_type, lspinorb, nspin, basis_type, smearing_method);
    EXPECT_EQ(nbands, 15);
}

TEST_F(UcellTest, CalNbandsSpin1LCAO)
{
    const int nelec = 10;
    const int nlocal = 6;
    int nbands = 0;
    const std::string esolver_type = "ksdft";
    const bool lspinorb = false;
    const int nspin = 1;
    const std::string basis_type = "lcao";
    const std::string smearing_method = "fixed";
    std::vector<double> nelec_spin(2, 5.0);
    unitcell::cal_nbands(nelec, nlocal, nelec_spin, nbands, esolver_type, lspinorb, nspin, basis_type, smearing_method);
    EXPECT_EQ(nbands, 6);
}

TEST_F(UcellTest, CalNbandsSpin4)
{
    const int nelec = 10;
    const int nlocal = 6;
    int nbands = 0;
    const std::string esolver_type = "ksdft";
    const bool lspinorb = false;
    const int nspin = 4;
    const std::string basis_type = "pw";
    const std::string smearing_method = "fixed";
    std::vector<double> nelec_spin(2, 5.0);
    unitcell::cal_nbands(nelec, nlocal, nelec_spin, nbands, esolver_type, lspinorb, nspin, basis_type, smearing_method);
    EXPECT_EQ(nbands, 30);
}

TEST_F(UcellTest, CalNbandsSpin4LCAO)
{
    const int nelec = 10;
    const int nlocal = 6;
    int nbands = 0;
    const std::string esolver_type = "ksdft";
    const bool lspinorb = false;
    const int nspin = 4;
    const std::string basis_type = "lcao";
    const std::string smearing_method = "fixed";
    std::vector<double> nelec_spin(2, 5.0);
    unitcell::cal_nbands(nelec, nlocal, nelec_spin, nbands, esolver_type, lspinorb, nspin, basis_type, smearing_method);
    EXPECT_EQ(nbands, 6);
}

TEST_F(UcellTest, CalNbandsSpin2)
{
    const int nelec = 10;
    const int nlocal = 6;
    int nbands = 0;
    const std::string esolver_type = "ksdft";
    const bool lspinorb = false;
    const int nspin = 2;
    const std::string basis_type = "pw";
    const std::string smearing_method = "fixed";
    std::vector<double> nelec_spin(2, 5.0);
    unitcell::cal_nbands(nelec, nlocal, nelec_spin, nbands, esolver_type, lspinorb, nspin, basis_type, smearing_method);
    EXPECT_EQ(nbands, 16);
}

TEST_F(UcellTest, CalNbandsSpin2LCAO)
{
    const int nelec = 10;
    const int nlocal = 6;
    int nbands = 0;
    const std::string esolver_type = "ksdft";
    const bool lspinorb = false;
    const int nspin = 2;
    const std::string basis_type = "lcao";
    const std::string smearing_method = "fixed";
    std::vector<double> nelec_spin(2, 5.0);
    unitcell::cal_nbands(nelec, nlocal, nelec_spin, nbands, esolver_type, lspinorb, nspin, basis_type, smearing_method);
    EXPECT_EQ(nbands, 6);
}

TEST_F(UcellTest, CalNbandsGaussWarning)
{
    const int nelec = 10;
    const int nlocal = 6;
    int nbands = 5;
    const std::string esolver_type = "ksdft";
    const bool lspinorb = false;
    const int nspin = 1;
    const std::string basis_type = "pw";
    const std::string smearing_method = "gaussian";
    std::vector<double> nelec_spin(2, 5.0);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::cal_nbands(nelec, nlocal, nelec_spin, nbands, esolver_type, lspinorb, nspin, basis_type, smearing_method), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("for smearing, num. of bands > num. of occupied bands"));
}

#ifdef __MPI
#include "mpi.h"
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);

    MPI_Comm_size(MPI_COMM_WORLD, &GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);

    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
#endif
