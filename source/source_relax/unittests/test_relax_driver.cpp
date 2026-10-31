#include "source_cell/cell_tools.h"
#include "source_cell/print_cell.h"
#include "source_cell/update_cell.h"
#include "source_esolver/esolver.h"
#include "for_test.h"
#include "source_relax/relax_driver.h"
#include "source_relax/relax_stru_io.h"
#include "source_relax/socket_driver.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <fstream>
#include <iterator>
#include <string>
#include <unistd.h>

namespace unitcell
{
bool allow_cell_change = false;

bool if_atoms_can_move(const Atom*, const int)
{
    return true;
}

bool if_cell_can_change(const std::vector<int>&)
{
    return allow_cell_change;
}

void update_pos_tau(const Lattice& lat, const double* pos, const int ntype, const int nat, Atom* atoms)
{
    int iat = 0;
    for (int it = 0; it < ntype; ++it)
    {
        for (int ia = 0; ia < atoms[it].na; ++ia)
        {
            for (int ik = 0; ik < 3; ++ik)
            {
                atoms[it].tau[ia][ik] = pos[3 * iat + ik] / lat.lat0;
            }
            ++iat;
        }
    }
    EXPECT_EQ(iat, nat);
}

void update_pos_taud(const Lattice&, const double*, const int, const int, Atom*)
{
}

void setup_cell_after_vc(UnitCell& ucell, std::ofstream&, const int)
{
    ucell.omega = std::abs(ucell.latvec.Det()) * ucell.lat0 * ucell.lat0 * ucell.lat0;
}

void remake_cell(Lattice&)
{
}

void print_tau(Atom*, const std::string&, const int, const double, std::ofstream&)
{
}
} // namespace unitcell

namespace ModuleIO
{
void print_screen(const int&, const int&, const int&)
{
}

void output_after_relax(const bool, const bool, std::ofstream&)
{
}
} // namespace ModuleIO

namespace relax_stru_io
{
std::string build_stru_header(const int,
                              const double,
                              const ModuleBase::matrix&,
                              const Input_para&,
                              const bool,
                              const bool)
{
    return std::string();
}

bool need_orbital(const Input_para&)
{
    return false;
}

void write_stru(UnitCell&,
                const Input_para&,
                const std::string&,
                const std::string&,
                const ModuleBase::matrix&,
                const bool,
                const bool,
                const int,
                const bool)
{
}
} // namespace relax_stru_io

void Socket_Driver::socket_driver(ModuleESolver::ESolver*, UnitCell&, const Input_para&, std::ofstream&)
{
}

namespace
{

class DriverTestESolver : public ModuleESolver::ESolver
{
  public:
    explicit DriverTestESolver(const double force_x) : force_x_(force_x)
    {
    }

    void before_all_runners(BaseCell&, const Input_para&) override
    {
    }

    void runner(BaseCell& cell, const int) override
    {
        ++runner_count;
        runner_omega = cell.omega();
        runner_latvec_e11 = cell.latvec().e11;
    }

    void after_all_runners(BaseCell& cell) override
    {
        ++after_all_count;
        after_all_omega = cell.omega();
        after_all_latvec_e11 = cell.latvec().e11;
    }

    double cal_energy() override
    {
        return 0.0;
    }

    void cal_force(BaseCell&, ModuleBase::matrix& force) override
    {
        force.zero_out();
        force(0, 0) = force_x_;
    }

    void cal_stress(BaseCell&, ModuleBase::matrix& stress) override
    {
        stress.zero_out();
        stress(0, 0) = stress_x_;
        stress(1, 1) = stress_x_;
        stress(2, 2) = stress_x_;
    }

    void set_stress(const double stress_x)
    {
        stress_x_ = stress_x;
    }

    int runner_count = 0;
    int after_all_count = 0;
    double runner_omega = 0.0;
    double after_all_omega = 0.0;
    double runner_latvec_e11 = 0.0;
    double after_all_latvec_e11 = 0.0;

  private:
    double force_x_;
    double stress_x_ = 0.0;
};

class RelaxDriverTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        log_file_ = testing::TempDir() + "running.log";
    }

    static void make_cell(UnitCell& ucell)
    {
        ucell.ntype = 1;
        ucell.nat = 1;
        ucell.lat0 = 1.0;
        ucell.omega = 1.0;
        ucell.lat.lat0 = 1.0;
        ucell.lat.latvec.Identity();
        ucell.latvec.Identity();
        ucell.atoms = new Atom[1];
        ucell.atoms[0].na = 1;
        ucell.atoms[0].tau.resize(1);
        ucell.atoms[0].taud.resize(1);
        ucell.atoms[0].mbl.resize(1, {1, 1, 1});
        ucell.iat2it.resize(1);
        ucell.iat2ia.resize(1);
        ucell.iat2it[0] = 0;
        ucell.iat2ia[0] = 0;
    }

    std::string read_running_log() const
    {
        std::ifstream ifs(log_file_);
        return std::string((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    }

    std::string log_file_;
};

TEST_F(RelaxDriverTest, RelaxNmaxZeroDoesNotRunESolver)
{
    Input_para inp;
    inp.calculation = "relax";
    inp.relax_method = "lbfgs";
    inp.relax_nmax = 0;
    inp.out_level = "m";

    UnitCell ucell;
    make_cell(ucell);
    DriverTestESolver esolver(1.0);
    std::ofstream ofs(log_file_);
    Relax_Driver driver;
    testing::internal::CaptureStdout();
    driver.relax_driver(&esolver, ucell, inp, ofs);
    const std::string output = testing::internal::GetCapturedStdout();
    ofs.close();

    EXPECT_EQ(esolver.runner_count, 0);
    EXPECT_EQ(esolver.after_all_count, 1);
    EXPECT_THAT(output, testing::HasSubstr("relax_nmax = 0, DRY RUN TEST SUCCEEDS"));
}

TEST_F(RelaxDriverTest, FirstStepConvergesOnNonSimultaneousPath)
{
    Input_para inp;
    inp.calculation = "relax";
    inp.relax_method = "lbfgs";
    inp.relax_nmax = 3;
    inp.force_thr = 0.1;
    inp.force_thr_ev = 0.1 * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A;
    inp.cal_force = 1;
    inp.out_level = "m";

    UnitCell ucell;
    make_cell(ucell);
    DriverTestESolver esolver(1.0e-4);
    std::ofstream ofs(log_file_);
    Relax_Driver driver;
    testing::internal::CaptureStdout();
    driver.relax_driver(&esolver, ucell, inp, ofs);
    const std::string output = testing::internal::GetCapturedStdout();
    ofs.close();

    EXPECT_EQ(esolver.runner_count, 1);
    EXPECT_EQ(esolver.after_all_count, 1);
    EXPECT_THAT(read_running_log(), testing::HasSubstr("Relaxation is converged!"));
    EXPECT_THAT(output, testing::HasSubstr("Geometry relaxation thresholds are reached within 1 steps."));
}

TEST_F(RelaxDriverTest, FirstStepConvergesOnSimultaneousPath)
{
    Input_para inp;
    inp.calculation = "relax";
    inp.relax_method = "cg";
    inp.relax_nmax = 3;
    inp.force_thr = 0.1;
    inp.force_thr_ev = 0.1 * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A;
    inp.cal_force = 1;
    inp.out_level = "m";

    UnitCell ucell;
    make_cell(ucell);
    DriverTestESolver esolver(1.0e-4);
    std::ofstream ofs(log_file_);
    Relax_Driver driver;
    testing::internal::CaptureStdout();
    driver.relax_driver(&esolver, ucell, inp, ofs);
    const std::string output = testing::internal::GetCapturedStdout();
    ofs.close();

    EXPECT_EQ(esolver.runner_count, 1);
    EXPECT_EQ(esolver.after_all_count, 1);
    EXPECT_THAT(read_running_log(), testing::HasSubstr("Relaxation is converged!"));
    EXPECT_THAT(output, testing::HasSubstr("Geometry relaxation thresholds are reached within 1 steps."));
}

TEST_F(RelaxDriverTest, MaxStepsMovesFinalNonSimultaneousGeometry)
{
    Input_para inp;
    inp.calculation = "relax";
    inp.relax_method = "lbfgs";
    inp.relax_nmax = 1;
    inp.force_thr = 0.0;
    inp.force_thr_ev = 0.0;
    inp.cal_force = 1;
    inp.relax_bfgs_rmax = 0.2;
    inp.out_level = "m";

    UnitCell ucell;
    make_cell(ucell);
    DriverTestESolver esolver(1.0);
    std::ofstream ofs(log_file_);
    Relax_Driver driver;
    driver.relax_driver(&esolver, ucell, inp, ofs);
    ofs.close();

    EXPECT_EQ(esolver.runner_count, 1);
    EXPECT_EQ(esolver.after_all_count, 1);
    EXPECT_THAT(read_running_log(), testing::HasSubstr("reaching relax_nmax"));
    EXPECT_THAT(read_running_log(), testing::HasSubstr("stopped after 1 ionic step(s) (relax_nmax = 1 reached)"));
}

TEST_F(RelaxDriverTest, MaxStepsStopsOnSimultaneousPath)
{
    Input_para inp;
    inp.calculation = "relax";
    inp.relax_method = "cg";
    inp.relax_nmax = 1;
    inp.force_thr = -1.0;
    inp.force_thr_ev = -1.0;
    inp.cal_force = 1;
    inp.out_level = "m";

    UnitCell ucell;
    make_cell(ucell);
    DriverTestESolver esolver(1.0);
    std::ofstream ofs(log_file_);
    Relax_Driver driver;
    driver.relax_driver(&esolver, ucell, inp, ofs);
    ofs.close();

    EXPECT_EQ(esolver.runner_count, 1);
    EXPECT_EQ(esolver.after_all_count, 1);
    EXPECT_THAT(read_running_log(), testing::HasSubstr("reaching relax_nmax"));
    EXPECT_THAT(read_running_log(), testing::HasSubstr("stopped after 1 ionic step(s) (relax_nmax = 1 reached)"));
}

TEST_F(RelaxDriverTest, MaxStepsAfterAllUsesEvaluatedCellRelaxGeometry)
{
    unitcell::allow_cell_change = true;

    Input_para inp;
    inp.calculation = "cell-relax";
    inp.relax_method = "cg";
    inp.relax_nmax = 1;
    inp.force_thr = -1.0;
    inp.force_thr_ev = -1.0;
    inp.stress_thr = -1.0;
    inp.cal_force = 1;
    inp.cal_stress = 1;
    inp.fixed_axes = "None";
    inp.out_level = "m";

    UnitCell ucell;
    make_cell(ucell);
    ucell.lat_axis_free.assign(3, 1);
    DriverTestESolver esolver(0.0);
    esolver.set_stress(1.0);
    std::ofstream ofs(log_file_);
    Relax_Driver driver;
    driver.relax_driver(&esolver, ucell, inp, ofs);
    ofs.close();

    unitcell::allow_cell_change = false;

    EXPECT_EQ(esolver.runner_count, 1);
    EXPECT_EQ(esolver.after_all_count, 1);
    EXPECT_DOUBLE_EQ(esolver.after_all_omega, esolver.runner_omega);
    EXPECT_DOUBLE_EQ(esolver.after_all_latvec_e11, esolver.runner_latvec_e11);
    EXPECT_NE(ucell.omega, esolver.runner_omega);
    EXPECT_THAT(read_running_log(), testing::HasSubstr("reaching relax_nmax"));
}

TEST_F(RelaxDriverTest, FinalAllowedStepConvergesAndCallsAfterAllOnce)
{
    Input_para inp;
    inp.calculation = "relax";
    inp.relax_method = "lbfgs";
    inp.relax_nmax = 1;
    inp.force_thr = 0.1;
    inp.force_thr_ev = 0.1 * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A;
    inp.cal_force = 1;
    inp.out_level = "m";

    UnitCell ucell;
    make_cell(ucell);
    DriverTestESolver esolver(1.0e-4);
    std::ofstream ofs(log_file_);
    Relax_Driver driver;
    driver.relax_driver(&esolver, ucell, inp, ofs);
    ofs.close();

    EXPECT_EQ(esolver.runner_count, 1);
    EXPECT_EQ(esolver.after_all_count, 1);
    EXPECT_THAT(read_running_log(), testing::HasSubstr("Relaxation is converged!"));
}

TEST_F(RelaxDriverTest, ExitStopsBeforeNonSimultaneousRelaxStep)
{
    char old_cwd[4096];
    ASSERT_NE(::getcwd(old_cwd, sizeof(old_cwd)), nullptr);
    ASSERT_EQ(::chdir(testing::TempDir().c_str()), 0);
    {
        std::ofstream exit_file("EXIT");
        exit_file << "stop_ion true\n";
    }

    Input_para inp;
    inp.calculation = "relax";
    inp.relax_method = "lbfgs";
    inp.relax_nmax = 3;
    inp.force_thr = 0.0;
    inp.force_thr_ev = 0.0;
    inp.cal_force = 1;
    inp.out_level = "m";

    UnitCell ucell;
    make_cell(ucell);
    DriverTestESolver esolver(1.0);
    std::ofstream ofs(log_file_);
    Relax_Driver driver;
    driver.relax_driver(&esolver, ucell, inp, ofs);
    ofs.close();

    ASSERT_EQ(::chdir(old_cwd), 0);
    EXPECT_EQ(esolver.runner_count, 1);
    EXPECT_EQ(esolver.after_all_count, 1);
    EXPECT_DOUBLE_EQ(ucell.atoms[0].tau[0][0], 0.0);
    EXPECT_THAT(read_running_log(), testing::HasSubstr("stopped by EXIT"));
    EXPECT_THAT(read_running_log(), testing::HasSubstr("stopped by EXIT after 1 ionic step(s)"));
}

TEST_F(RelaxDriverTest, ExitStopsBeforeSimultaneousRelaxStep)
{
    char old_cwd[4096];
    ASSERT_NE(::getcwd(old_cwd, sizeof(old_cwd)), nullptr);
    ASSERT_EQ(::chdir(testing::TempDir().c_str()), 0);
    {
        std::ofstream exit_file("EXIT");
        exit_file << "stop_ion true\n";
    }

    Input_para inp;
    inp.calculation = "relax";
    inp.relax_method = "cg";
    inp.relax_nmax = 3;
    inp.force_thr = -1.0;
    inp.force_thr_ev = -1.0;
    inp.cal_force = 1;
    inp.out_level = "m";

    UnitCell ucell;
    make_cell(ucell);
    DriverTestESolver esolver(1.0);
    std::ofstream ofs(log_file_);
    Relax_Driver driver;
    driver.relax_driver(&esolver, ucell, inp, ofs);
    ofs.close();

    ASSERT_EQ(::chdir(old_cwd), 0);
    EXPECT_EQ(esolver.runner_count, 1);
    EXPECT_EQ(esolver.after_all_count, 1);
    EXPECT_THAT(read_running_log(), testing::HasSubstr("stopped by EXIT"));
    EXPECT_THAT(read_running_log(), testing::HasSubstr("stopped by EXIT after 1 ionic step(s)"));
}

} // namespace
