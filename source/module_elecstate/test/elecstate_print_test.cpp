#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_elecstate/elecstate.h"
#include "module_elecstate/elecstate_getters.h"
#include "module_elecstate/potentials/efield.h"
#include "module_elecstate/potentials/gatefield.h"
#include "module_elecstate/module_charge/charge.h"

/***************************************************************
 *  mock functions
 ****************************************************************/
namespace elecstate
{
const double* ElecState::getRho(int spin) const
{
    return &(this->eferm.ef);
} // just for mock
void ElecState::calculate_weights()
{
    return;
} // just for mock
double Efield::etotefield = 1.1;
double elecstate::Gatefield::etotgatefield = 2.2;
std::string tmp_vdw_method = "d2";
std::string get_input_vdw_method()
{
    return tmp_vdw_method;
}
double get_ucell_tot_magnetization()
{
    return 1.1;
}
double get_ucell_abs_magnetization()
{
    return 2.2;
}
double get_ucell_tot_magnetization_nc_x()
{
    return 3.3;
}
double get_ucell_tot_magnetization_nc_y()
{
    return 4.4;
}
double get_ucell_tot_magnetization_nc_z()
{
    return 5.5;
}
std::string tmp_ks_solver = "dav";
std::string get_ks_solver_type()
{
    return tmp_ks_solver;
}
} // namespace elecstate

K_Vectors::K_Vectors()
{
}
K_Vectors::~K_Vectors()
{
}
Charge::Charge()
{
}
Charge::~Charge()
{
}

/***************************************************************
 *  unit test of functions in elecstate_print.cpp
 ****************************************************************/

/**
 * - Tested functions:
 *   - ElecState::print_format()
 *   - ElecState::print_eigenvalue()
 */

class ElecStatePrintTest : public ::testing::Test
{
  protected:
    elecstate::ElecState elecstate;
    std::string output;
    std::ifstream ifs;
    std::ofstream ofs;
    K_Vectors* p_klist;
    void SetUp()
    {
        p_klist = new K_Vectors;
        p_klist->nks = 2;
        p_klist->isk = {0, 1};
        p_klist->ngk = {100, 101};
        p_klist->kvec_c.resize(2);
        p_klist->kvec_c[0].set(0.1, 0.11, 0.111);
        p_klist->kvec_c[1].set(0.2, 0.22, 0.222);
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
    }
    void TearDown()
    {
        delete p_klist;
    }
};

TEST_F(ElecStatePrintTest, PrintFormat)
{
    GlobalV::ofs_running.open("test.dat", std::ios::out);
    elecstate.print_format("test", 0.1);
    GlobalV::ofs_running.close();
    ifs.open("test.dat", std::ios::in);
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("test                          +0.1                      +1.36057"));
    ifs.close();
    std::remove("test.dat");
}

TEST_F(ElecStatePrintTest, PrintEigenvalueS2)
{
    GlobalV::NSPIN = 2;
    GlobalV::ofs_running.open("test.dat", std::ios::out);
    // print eigenvalue
    elecstate.print_eigenvalue(GlobalV::ofs_running);
    GlobalV::ofs_running.close();
    ifs.open("test.dat", std::ios::in);
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("STATE ENERGY(eV) AND OCCUPATIONS"));
    EXPECT_THAT(str, testing::HasSubstr("NSPIN == 2"));
    EXPECT_THAT(str, testing::HasSubstr("SPIN UP :"));
    EXPECT_THAT(str, testing::HasSubstr("1/1 kpoint (Cartesian) = 0.10000 0.11000 0.11100 (100 pws)"));
    EXPECT_THAT(str, testing::HasSubstr("1        13.6057       0.100000"));
    EXPECT_THAT(str, testing::HasSubstr("2        27.2114       0.200000"));
    EXPECT_THAT(str, testing::HasSubstr("SPIN DOWN :"));
    EXPECT_THAT(str, testing::HasSubstr("1/1 kpoint (Cartesian) = 0.20000 0.22000 0.22200 (101 pws)"));
    EXPECT_THAT(str, testing::HasSubstr("1        40.8171       0.300000"));
    EXPECT_THAT(str, testing::HasSubstr("2        54.4228       0.400000"));
    ifs.close();
    std::remove("test.dat");
}

TEST_F(ElecStatePrintTest, PrintEigenvalueS4)
{
    GlobalV::NSPIN = 4;
    GlobalV::ofs_running.open("test.dat", std::ios::out);
    // print eigenvalue
    elecstate.print_eigenvalue(GlobalV::ofs_running);
    GlobalV::ofs_running.close();
    ifs.open("test.dat", std::ios::in);
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("STATE ENERGY(eV) AND OCCUPATIONS"));
    EXPECT_THAT(str, testing::HasSubstr("NSPIN == 4"));
    EXPECT_THAT(str, testing::HasSubstr("1/2 kpoint (Cartesian) = 0.10000 0.11000 0.11100 (100 pws)"));
    EXPECT_THAT(str, testing::HasSubstr("1        13.6057       0.100000"));
    EXPECT_THAT(str, testing::HasSubstr("2        27.2114       0.200000"));
    EXPECT_THAT(str, testing::HasSubstr("2/2 kpoint (Cartesian) = 0.20000 0.22000 0.22200 (101 pws)"));
    EXPECT_THAT(str, testing::HasSubstr("1        40.8171       0.300000"));
    EXPECT_THAT(str, testing::HasSubstr("2        54.4228       0.400000"));
    ifs.close();
    std::remove("test.dat");
}

TEST_F(ElecStatePrintTest, PrintBand)
{
    GlobalV::NSPIN = 1;
    GlobalV::CURRENT_SPIN = 0;
    GlobalV::NBANDS = 2;
    GlobalV::MY_RANK = 0;
    GlobalV::ofs_running.open("test.dat", std::ios::out);
    // print eigenvalue
    elecstate.print_band(0, 1, 0);
    GlobalV::ofs_running.close();
    ifs.open("test.dat", std::ios::in);
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("Energy (eV) & Occupations  for spin=1 K-point=1"));
    EXPECT_THAT(str, testing::HasSubstr("1        13.6057       0.100000"));
    EXPECT_THAT(str, testing::HasSubstr("2        27.2114       0.200000"));
    ifs.close();
    std::remove("test.dat");
}

TEST_F(ElecStatePrintTest, PrintEigenvalueWarning)
{
    elecstate.ekb(0, 0) = 1.0e11;
    GlobalV::NSPIN = 4;
    GlobalV::ofs_running.open("test.dat", std::ios::out);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(elecstate.print_eigenvalue(GlobalV::ofs_running), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Eigenvalues are too large!"));
    GlobalV::ofs_running.close();
    std::remove("test.dat");
}

TEST_F(ElecStatePrintTest, PrintBandWarning)
{
    elecstate.ekb(0, 0) = 1.0e11;
    GlobalV::NSPIN = 4;
    GlobalV::ofs_running.open("test.dat", std::ios::out);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(elecstate.print_band(0, 1, 0), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Eigenvalues are too large!"));
    GlobalV::ofs_running.close();
    std::remove("test.dat");
}

TEST_F(ElecStatePrintTest, PrintEtot)
{
    GlobalV::ofs_running.open("test.dat", std::ios::out);
    bool converged = false;
    int iter = 1;
    double scf_thr = 0.1;
    double duration = 2.0;
    int printe = 1;
    double pw_diag_thr = 0.1;
    int avg_iter = 2;
    bool print = true;
    elecstate.charge = new Charge;
    elecstate.charge->nrxx = 100;
    elecstate.charge->nxyz = 1000;
    GlobalV::imp_sol = true;
    GlobalV::EFIELD_FLAG = true;
    GlobalV::GATE_FLAG = true;
    GlobalV::TWO_EFERMI = true;
    GlobalV::out_bandgap = true;
    GlobalV::COLOUR = false;
    GlobalV::MY_RANK = 0;
    GlobalV::BASIS_TYPE = "pw";
    GlobalV::NSPIN = 2;
    // iteration of different vdw_method
    std::vector<std::string> vdw_methods = {"d2", "d3_0", "d3_bj"};
    for (int i = 0; i < vdw_methods.size(); i++)
    {
        elecstate::tmp_vdw_method = vdw_methods[i];
        elecstate.print_etot(converged, iter, scf_thr, duration, printe, pw_diag_thr, avg_iter, false);
    }
    // iteration of different ks_solver
    std::vector<std::string> ks_solvers = {"cg", "lapack", "genelpa", "dav", "scalapack_gvx", "cusolver"};
    for (int i = 0; i < ks_solvers.size(); i++)
    {
        elecstate::tmp_ks_solver = ks_solvers[i];
        testing::internal::CaptureStdout();
        elecstate.print_etot(converged, iter, scf_thr, duration, printe, pw_diag_thr, avg_iter, print);
        output = testing::internal::GetCapturedStdout();
        if (elecstate::tmp_ks_solver == "cg")
        {
            EXPECT_THAT(output, testing::HasSubstr("CG"));
        }
        else if (elecstate::tmp_ks_solver == "lapack")
        {
            EXPECT_THAT(output, testing::HasSubstr("LA"));
        }
        else if (elecstate::tmp_ks_solver == "genelpa")
        {
            EXPECT_THAT(output, testing::HasSubstr("GE"));
        }
        else if (elecstate::tmp_ks_solver == "dav")
        {
            EXPECT_THAT(output, testing::HasSubstr("DA"));
        }
        else if (elecstate::tmp_ks_solver == "scalapack_gvx")
        {
            EXPECT_THAT(output, testing::HasSubstr("GV"));
        }
        else if (elecstate::tmp_ks_solver == "cusolver")
        {
            EXPECT_THAT(output, testing::HasSubstr("CU"));
        }
    }
    GlobalV::ofs_running.close();
    ifs.open("test.dat", std::ios::in);
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("Density error is 0.1"));
    EXPECT_THAT(str, testing::HasSubstr("Error Threshold = 0.1"));
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

TEST_F(ElecStatePrintTest, PrintEtot2)
{
    GlobalV::ofs_running.open("test.dat", std::ios::out);
    bool converged = false;
    int iter = 1;
    double scf_thr = 0.1;
    double duration = 2.0;
    int printe = 0;
    double pw_diag_thr = 0.1;
    int avg_iter = 2;
    bool print = true;
    elecstate.charge = new Charge;
    elecstate.charge->nrxx = 100;
    elecstate.charge->nxyz = 1000;
    GlobalV::imp_sol = true;
    GlobalV::EFIELD_FLAG = true;
    GlobalV::GATE_FLAG = true;
    GlobalV::TWO_EFERMI = false;
    GlobalV::out_bandgap = true;
    GlobalV::COLOUR = false;
    GlobalV::MY_RANK = 0;
    GlobalV::BASIS_TYPE = "pw";
    GlobalV::SCF_NMAX = 100;
    elecstate.print_etot(converged, iter, scf_thr, duration, printe, pw_diag_thr, avg_iter, print);
    GlobalV::ofs_running.close();
    ifs.open("test.dat", std::ios::in);
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("Density error is 0.1"));
    EXPECT_THAT(str, testing::HasSubstr("Error Threshold = 0.1"));
    EXPECT_THAT(str, testing::HasSubstr("E_KohnSham"));
    EXPECT_THAT(str, testing::HasSubstr("E_Harris"));
    EXPECT_THAT(str, testing::HasSubstr("E_Fermi"));
    EXPECT_THAT(str, testing::HasSubstr("E_bandgap"));
    ifs.close();
    delete elecstate.charge;
    std::remove("test.dat");
}

TEST_F(ElecStatePrintTest, PrintEtotColorS2)
{
    bool converged = false;
    int iter = 1;
    double scf_thr = 2.0;
    double duration = 2.0;
    int printe = 1;
    double pw_diag_thr = 0.1;
    int avg_iter = 2;
    bool print = true;
    elecstate.charge = new Charge;
    elecstate.charge->nrxx = 100;
    elecstate.charge->nxyz = 1000;
    GlobalV::imp_sol = true;
    GlobalV::EFIELD_FLAG = true;
    GlobalV::GATE_FLAG = true;
    GlobalV::TWO_EFERMI = true;
    GlobalV::out_bandgap = true;
    GlobalV::COLOUR = true;
    GlobalV::NSPIN = 2;
    GlobalV::MY_RANK = 0;
    elecstate.print_etot(converged, iter, scf_thr, duration, printe, pw_diag_thr, avg_iter, print);
}

TEST_F(ElecStatePrintTest, PrintEtotColorS4)
{
    bool converged = false;
    int iter = 1;
    double scf_thr = 0.1;
    double duration = 2.0;
    int printe = 1;
    double pw_diag_thr = 0.1;
    int avg_iter = 2;
    bool print = true;
    elecstate.charge = new Charge;
    elecstate.charge->nrxx = 100;
    elecstate.charge->nxyz = 1000;
    GlobalV::imp_sol = true;
    GlobalV::EFIELD_FLAG = true;
    GlobalV::GATE_FLAG = true;
    GlobalV::TWO_EFERMI = true;
    GlobalV::out_bandgap = true;
    GlobalV::COLOUR = true;
    GlobalV::NSPIN = 4;
    GlobalV::NONCOLIN = true;
    GlobalV::MY_RANK = 0;
    elecstate.print_etot(converged, iter, scf_thr, duration, printe, pw_diag_thr, avg_iter, print);
}

TEST_F(ElecStatePrintTest, PrintEtotWarning)
{
    GlobalV::ofs_running.open("test.dat", std::ios::out);
    bool converged = false;
    int iter = 1;
    double scf_thr = 0.1;
    double duration = 2.0;
    int printe = 0;
    double pw_diag_thr = 0.1;
    int avg_iter = 2;
    bool print = true;
    elecstate.charge = new Charge;
    elecstate.charge->nrxx = 100;
    elecstate.charge->nxyz = 1000;
    GlobalV::imp_sol = true;
    GlobalV::EFIELD_FLAG = true;
    GlobalV::GATE_FLAG = true;
    GlobalV::TWO_EFERMI = false;
    GlobalV::out_bandgap = true;
    GlobalV::COLOUR = false;
    GlobalV::MY_RANK = 0;
    GlobalV::BASIS_TYPE = "pw";
    GlobalV::SCF_NMAX = 100;
    elecstate::tmp_ks_solver = "unknown";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(elecstate.print_etot(converged, iter, scf_thr, duration, printe, pw_diag_thr, avg_iter, print), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("print_etot found unknown ks_solver_type"));
    GlobalV::ofs_running.close();
    delete elecstate.charge;
    std::remove("test.dat");
}

#include "mpi.h"
int main(int argc, char **argv)
{
        MPI_Init(&argc, &argv);
        testing::InitGoogleTest(&argc, argv);

        int nproc;
        int my_rank;

        MPI_Comm_size(MPI_COMM_WORLD,&nproc);
        MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

        int result = RUN_ALL_TESTS();

        MPI_Finalize();

        return result;
}