#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#define private public
#include "source_cell/klist.h"
#include "source_estate/elecstate.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_pot/efield.h"
#include "source_estate/module_pot/gatefield.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/parameter.h"
#include "source_estate/elecstate_print.h"
#undef private 
/***************************************************************
 *  mock functions
 ****************************************************************/
namespace elecstate
{
const double* ElecState::getRho(int spin) const
{
    return &(this->eferm.ef);
} // just for mock
double Efield::etotefield = 1.1;
double elecstate::Gatefield::etotgatefield = 2.2;

} // namespace elecstate
UnitCell::UnitCell(){}
UnitCell::~UnitCell(){}
Magnetism::Magnetism(){}
Magnetism::~Magnetism(){}
InfoNonlocal::InfoNonlocal(){}
InfoNonlocal::~InfoNonlocal(){}
Charge::Charge()
{
}
Charge::~Charge()
{
}

int XC_Functional::func_type = 0;
bool XC_Functional::ked_flag = false;

/***************************************************************
 *  unit test of functions in elecstate_print.cpp
 ****************************************************************/

/**
 * - Tested functions:
 *   - ElecState::print_format()
 */

class ElecStatePrintTest : public ::testing::Test
{
  protected:
    elecstate::ElecState elecstate;
    UnitCell ucell;
    std::string output;
    std::ifstream ifs;
    std::ofstream ofs;
    K_Vectors* p_klist;
    void SetUp()
    {
        p_klist = new K_Vectors;
        p_klist->set_nks(2);
        p_klist->set_nkstot(2);
        p_klist->isk = {0, 1};
        p_klist->ngk = {100, 101};
        p_klist->kvec_c.resize(2);
        p_klist->kvec_c[0].set(0.1, 0.11, 0.111);
        p_klist->kvec_c[1].set(0.2, 0.22, 0.222);
        p_klist->ik2iktot.resize(2);
        p_klist->ik2iktot[0] = 0;
        p_klist->ik2iktot[1] = 1;
        // initialize klist of elecstate
        elecstate.klist = p_klist;
        // initialize ekb of elecstate
        elecstate.ekb.create(2, 2);
        elecstate.ekb(0, 0) = 1.0;
        elecstate.ekb(0, 1) = 2.0;
        elecstate.ekb(1, 0) = 3.0;
        elecstate.ekb(1, 1) = 4.0;
        // initialize wg of elecstate
        elecstate.wg.create(2, 2);
        elecstate.wg(0, 0) = 0.1;
        elecstate.wg(0, 1) = 0.2;
        elecstate.wg(1, 0) = 0.3;
        elecstate.wg(1, 1) = 0.4;
        ucell.magnet.tot_mag = 1.1;
        ucell.magnet.abs_mag = 2.2;
        ucell.magnet.tot_mag_nc[0] = 3.3;
        ucell.magnet.tot_mag_nc[1] = 4.4;
        ucell.magnet.tot_mag_nc[2] = 5.5;
        PARAM.input.ks_solver = "dav";
        PARAM.sys.log_file = "test.dat";
    }
    void TearDown()
    {
        delete p_klist;
    }
};

TEST_F(ElecStatePrintTest, PrintFormat)
{
    GlobalV::ofs_running.open("test.dat", std::ios::out);
    elecstate::print_format("test", 0.1);
    GlobalV::ofs_running.close();
    ifs.open("test.dat", std::ios::in);
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("test                          +0.1                      +1.36057"));
    ifs.close();
    std::remove("test.dat");
}

TEST_F(ElecStatePrintTest, PrintEtot)
{
    GlobalV::ofs_running.open("test.dat", std::ios::out);
    bool converged = false;
    int iter = 1;
    double scf_thr = 0.1;
    double scf_thr_kin = 0.0;
    double duration = 2.0;
    double pw_diag_thr = 0.1;
    int avg_iter = 2;
    bool print = true;
    elecstate.charge = new Charge;
    elecstate.charge->nrxx = 100;
    elecstate.charge->nxyz = 1000;
    PARAM.input.out_freq_elec = 1;
    PARAM.input.imp_sol = true;
    PARAM.input.efield_flag = true;
    PARAM.input.gate_flag = true;
    PARAM.sys.two_fermi = true;
    PARAM.input.out_bandgap = true;
    GlobalV::MY_RANK = 0;
    PARAM.input.basis_type = "pw";
    PARAM.input.nspin = 2;

    // iteration of different vdw_method
    std::vector<std::string> vdw_methods = {"d2", "d3_0", "d3_bj"};
    for (int i = 0; i < vdw_methods.size(); i++)
    {
        PARAM.input.vdw_method = vdw_methods[i];
        elecstate::print_etot(ucell.magnet,elecstate, converged, iter, scf_thr, 
        scf_thr_kin, duration, pw_diag_thr, avg_iter, false);
    }

    // iteration of different ks_solver
    std::vector<std::string> ks_solvers = {"cg", "lapack", "genelpa", "dav", "scalapack_gvx", "cusolver"};
    for (int i = 0; i < ks_solvers.size(); i++)
    {
        PARAM.input.ks_solver = ks_solvers[i];
        testing::internal::CaptureStdout();

        elecstate::print_etot(ucell.magnet,elecstate,converged, iter, scf_thr, 
        scf_thr_kin, duration, pw_diag_thr, avg_iter, print);

        output = testing::internal::GetCapturedStdout();
        if (PARAM.input.ks_solver == "cg")
        {
            EXPECT_THAT(output, testing::HasSubstr("CG"));
        }
        else if (PARAM.input.ks_solver == "lapack")
        {
            EXPECT_THAT(output, testing::HasSubstr("LA"));
        }
        else if (PARAM.input.ks_solver == "genelpa")
        {
            EXPECT_THAT(output, testing::HasSubstr("GE"));
        }
        else if (PARAM.input.ks_solver == "dav")
        {
            EXPECT_THAT(output, testing::HasSubstr("DA"));
        }
        else if (PARAM.input.ks_solver == "scalapack_gvx")
        {
            EXPECT_THAT(output, testing::HasSubstr("GV"));
        }
        else if (PARAM.input.ks_solver == "cusolver")
        {
            EXPECT_THAT(output, testing::HasSubstr("CU"));
        }
    }
    GlobalV::ofs_running.close();
    ifs.open("test.dat", std::ios::in);
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("Electron density deviation 0.1"));
    EXPECT_THAT(str, testing::HasSubstr("Diago Threshold = 0.1"));
    EXPECT_THAT(str, testing::HasSubstr("E_KohnSham"));
    EXPECT_THAT(str, testing::HasSubstr("E_vdwD2"));
    EXPECT_THAT(str, testing::HasSubstr("E_vdwD3"));
    EXPECT_THAT(str, testing::HasSubstr("E_sol_el"));
    EXPECT_THAT(str, testing::HasSubstr("E_sol_cav"));
    EXPECT_THAT(str, testing::HasSubstr("E_efield"));
    EXPECT_THAT(str, testing::HasSubstr("E_gatefield"));
    ifs.close();
    delete elecstate.charge;
    std::remove("test.dat");
}

TEST_F(ElecStatePrintTest, PrintEtotColorS2)
{
    bool converged = false;
    int iter = 1;
    double scf_thr = 2.0;
    double scf_thr_kin = 0.0;
    double duration = 2.0;
    double pw_diag_thr = 0.1;
    int avg_iter = 2;
    bool print = true;
    elecstate.charge = new Charge;
    elecstate.charge->nrxx = 100;
    elecstate.charge->nxyz = 1000;

    PARAM.input.out_freq_elec = 1;
    PARAM.input.imp_sol = true;
    PARAM.input.efield_flag = true;
    PARAM.input.gate_flag = true;
    PARAM.sys.two_fermi = true;
    PARAM.input.out_bandgap = true;
    PARAM.input.nspin = 2;
    GlobalV::MY_RANK = 0;

    elecstate::print_etot(ucell.magnet,elecstate,converged, iter, scf_thr, 
    scf_thr_kin, duration, pw_diag_thr, avg_iter, print);

    delete elecstate.charge;
}


TEST_F(ElecStatePrintTest, PrintEtotColorS4)
{
    bool converged = false;
    int iter = 1;
    double scf_thr = 0.1;
    double scf_thr_kin = 0.0;
    double duration = 2.0;
    double pw_diag_thr = 0.1;
    int avg_iter = 2;
    bool print = true;
    elecstate.charge = new Charge;
    elecstate.charge->nrxx = 100;
    elecstate.charge->nxyz = 1000;

    PARAM.input.out_freq_elec = 1;
    PARAM.input.imp_sol = true;
    PARAM.input.efield_flag = true;
    PARAM.input.gate_flag = true;
    PARAM.sys.two_fermi = true;
    PARAM.input.out_bandgap = true;
    PARAM.input.nspin = 4;
    PARAM.input.noncolin = true;
    GlobalV::MY_RANK = 0;

    elecstate::print_etot(ucell.magnet,elecstate, converged, iter, scf_thr, scf_thr_kin, 
    duration, pw_diag_thr, avg_iter, print);

    delete elecstate.charge;
}
