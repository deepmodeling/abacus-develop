#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_base/global_variable.h"
#include "source_cell/unitcell.h"
#include "source_io/module_parameter/parameter.h"
#include "source_relax/lattice_change_basic.h"
#include "source_relax/relax_nsync.h"

/************************************************
 *  unit tests of class IonCellOptimizer (relax_nsync.cpp)
 ***********************************************/

/**
 * The running-log convergence reporting is centralized in
 * IonCellOptimizer::relax_step. These tests lock the unified format shared by
 * every relax method (sd/cg/bfgs/bfgs_trad/lbfgs) and by cell-relax:
 *   Largest force is ... eV/Angstrom while threshold is ... eV/Angstrom
 *   Largest stress is ... kbar while threshold is ... kbar
 *   Relaxation is converged!
 *   Relaxation is not converged yet!
 *   Relaxation is not converged after reaching relax_nmax!
 * The markers are what the ASE interface greps for (issue #6546).
 */

namespace unitcell
{
// Stub: the real update_pos_tau performs lattice-vector algebra that needs a
// fully built UnitCell; these tests only care about the log output.
void update_pos_tau(const Lattice&, const double*, const int, const int, Atom*)
{
}
// Stub: the real setup_cell_after_vc recomputes reciprocal lattice etc.
void setup_cell_after_vc(UnitCell&, std::ofstream&, const int)
{
}
} // namespace unitcell

UnitCell::UnitCell()
{
    Coordinate = "Direct";
    latName = "none";
    lat0 = 1.0;
    latvec.Identity();

    ntype = 1;
    nat = 1;
    itia2iat.create(1, 1);

    atoms = new Atom[ntype];
    set_atom_flag = true;
    atoms[0].label = "Si";
}
UnitCell::~UnitCell()
{
}
Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}
Atom::Atom()
{
    na = 1;
    tau.resize(na);
    dis.resize(na);
    mag.resize(na);
    mbl.resize(na);
    vel.resize(na);
    taud.resize(na);
    mbl[0] = {1, 1, 1};
}
Atom::~Atom()
{
}
#ifdef __MPI
void Atom::bcast_atom()
{
}
void Atom::bcast_atom2()
{
}
#endif
Atom_pseudo::Atom_pseudo()
{
}
Atom_pseudo::~Atom_pseudo()
{
}
pseudo::pseudo()
{
}
pseudo::~pseudo()
{
}
SepPot::SepPot()
{
}
SepPot::~SepPot()
{
}
Sep_Cell::Sep_Cell() noexcept
{
}
Sep_Cell::~Sep_Cell() noexcept
{
}
int ModuleSymmetry::Symmetry::symm_flag = 0;
void ModuleSymmetry::Symmetry::symmetrize_mat3(ModuleBase::matrix& sigma, const Lattice& lat) const
{
}
void ModuleSymmetry::Symmetry::symmetrize_vec3_nat(double* v) const
{
}

// Friend of Parameter; the only sanctioned way to write PARAM.input in tests.
class TestParameters
{
  public:
    static void set_relax(const std::string& out_level,
                          const double relax_bfgs_rmax,
                          const double relax_bfgs_rmin,
                          const double relax_bfgs_init,
                          const std::string& fixed_axes)
    {
        PARAM.input.out_level = out_level;
        PARAM.input.relax_bfgs_rmax = relax_bfgs_rmax;
        PARAM.input.relax_bfgs_rmin = relax_bfgs_rmin;
        PARAM.input.relax_bfgs_init = relax_bfgs_init;
        PARAM.input.fixed_axes = fixed_axes;
    }
};

class IonCellOptimizerTest : public ::testing::Test
{
  protected:
    IonCellOptimizer optimizer;
    Input_para inp;
    const int natom = 1;
    const double etot = 0.0;
    const std::string log_file = "relax_nsync_test.log";
    const std::string warning_file = "relax_nsync_test_warning.log";

    void SetUp() override
    {
        // keep stdout quiet (out_level "m") and give bfgs sane trust radii
        TestParameters::set_relax("m", 0.2, 1.0e-5, 0.5, "None");

        inp.calculation = "relax";
        inp.relax_method = {"lbfgs", "1"};
        inp.relax_nmax = 50;
        inp.force_thr = 0.1;                          // Ry/Bohr
        inp.force_thr_ev = 0.1 * 13.6058 / 0.529177;  // eV/Angstrom
        inp.stress_thr = 0.5;                         // kbar
        inp.cal_force = 1;
        inp.cal_stress = 1;
        inp.out_level = "m";
        inp.test_relax_method = 0;
        inp.fixed_ibrav = false;
        inp.nspin = 1;
    }

    void TearDown() override
    {
        std::remove(log_file.c_str());
        std::remove(warning_file.c_str());
    }

    // UnitCell is non-copyable (unique_ptr member), so fill it in place.
    void setup_ucell(UnitCell& ucell)
    {
        ucell.lat_axis_free[0] = 1;
        ucell.lat_axis_free[1] = 1;
        ucell.lat_axis_free[2] = 1;
        // The L-BFGS update loop dereferences iat2it/iat2ia from the second
        // iteration on, so give the single atom a valid type/index.
        ucell.iat2it.resize(natom);
        ucell.iat2ia.resize(natom);
        ucell.iat2it[0] = 0;
        ucell.iat2ia[0] = 0;
    }

    std::string read_log()
    {
        std::ifstream ifs(log_file);
        std::string content((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
        ifs.close();
        return content;
    }

    std::string read_warning_log()
    {
        std::ifstream ifs(warning_file);
        std::string content((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
        ifs.close();
        return content;
    }
};

// relax + small force: the log must carry the force line and the ASE marker.
TEST_F(IonCellOptimizerTest, RelaxConverged)
{
    UnitCell ucell;
    setup_ucell(ucell);
    optimizer.init_relax(natom, inp);
    ModuleBase::matrix force(natom, 3);
    force(0, 0) = 1.0e-4; // well below force_thr
    ModuleBase::matrix stress(3, 3);
    int force_step = 1;
    int stress_step = 1;
    std::ofstream ofs(log_file);

    const bool done = optimizer.relax_step(1, etot, ucell, force, stress, force_step, stress_step, ofs);
    ofs.close();

    EXPECT_TRUE(done);
    const std::string log = read_log();
    EXPECT_THAT(log, testing::HasSubstr(" Largest force is "));
    EXPECT_THAT(log, testing::HasSubstr(" eV/Angstrom while threshold is "));
    EXPECT_THAT(log, testing::HasSubstr(" Relaxation method: lbfgs"));
    EXPECT_THAT(log, testing::HasSubstr(" Relaxation converged in 1 step(s)."));
    EXPECT_THAT(log, testing::HasSubstr(" Largest force per step (eV/Angstrom):"));
    EXPECT_THAT(log, testing::HasSubstr("\n Relaxation is converged!"));
    EXPECT_THAT(log, testing::Not(testing::HasSubstr("not converged")));
}

// relax + large force: not converged, and no ASE converged marker may appear.
TEST_F(IonCellOptimizerTest, RelaxNotConverged)
{
    UnitCell ucell;
    setup_ucell(ucell);
    optimizer.init_relax(natom, inp);
    ModuleBase::matrix force(natom, 3);
    force(0, 0) = 1.0; // above force_thr
    ModuleBase::matrix stress(3, 3);
    int force_step = 1;
    int stress_step = 1;
    std::ofstream ofs(log_file);

    const bool done = optimizer.relax_step(1, etot, ucell, force, stress, force_step, stress_step, ofs);
    ofs.close();

    EXPECT_FALSE(done);
    const std::string log = read_log();
    EXPECT_THAT(log, testing::HasSubstr(" Largest force is "));
    EXPECT_THAT(log, testing::HasSubstr("\n Relaxation is not converged yet!"));
    EXPECT_THAT(log, testing::Not(testing::HasSubstr("Relaxation is converged!")));
}

// Reaching relax_nmax must terminate with an explicit not-converged marker.
TEST_F(IonCellOptimizerTest, RelaxNmaxReached)
{
    UnitCell ucell;
    setup_ucell(ucell);
    optimizer.init_relax(natom, inp);
    ModuleBase::matrix force(natom, 3);
    ModuleBase::matrix stress(3, 3);
    int force_step = 1;
    int stress_step = 1;
    std::ofstream ofs(log_file);

    const bool done = optimizer.relax_step(inp.relax_nmax, etot, ucell, force, stress, force_step, stress_step, ofs);
    ofs.close();

    EXPECT_TRUE(done);
    // Key invariant: when relax_nmax terminates the loop, cal_movement is
    // skipped and the geometry is NOT touched. The force/stress captured by
    // the driver before relax_step are still consistent with ucell, so the
    // driver's geometry_evaluated=true is honest.
    EXPECT_FALSE(ucell.ionic_position_updated);
    const std::string log = read_log();
    EXPECT_THAT(log, testing::HasSubstr(" Relaxation method: lbfgs"));
    EXPECT_THAT(log, testing::HasSubstr(" ionic step(s) (relax_nmax = 50 reached)."));
    EXPECT_THAT(log, testing::HasSubstr(" Relaxation is not converged after reaching relax_nmax!"));
    EXPECT_THAT(log, testing::Not(testing::HasSubstr("Relaxation is converged!")));
}

// cell-relax + small force + small stress: both force and stress lines are
// reported and the converged marker is printed once the cell step converges.
TEST_F(IonCellOptimizerTest, CellRelaxConverged)
{
    inp.calculation = "cell-relax";
    UnitCell ucell;
    setup_ucell(ucell);
    optimizer.init_relax(natom, inp);
    ModuleBase::matrix force(natom, 3);
    force(0, 0) = 1.0e-4;
    ModuleBase::matrix stress(3, 3); // zero stress, below stress_thr
    int force_step = 1;
    int stress_step = 1;
    std::ofstream ofs(log_file);

    const bool done = optimizer.relax_step(1, etot, ucell, force, stress, force_step, stress_step, ofs);
    ofs.close();

    EXPECT_TRUE(done);
    const std::string log = read_log();
    EXPECT_THAT(log, testing::HasSubstr(" Largest force is "));
    EXPECT_THAT(log, testing::HasSubstr(" Largest stress is "));
    EXPECT_THAT(log, testing::HasSubstr(" kbar while threshold is "));
    EXPECT_THAT(log, testing::HasSubstr(" Relaxation method: lbfgs"));
    EXPECT_THAT(log, testing::HasSubstr(" Relaxation converged in 1 step(s)."));
    EXPECT_THAT(log, testing::HasSubstr(" Largest force per step (eV/Angstrom):"));
    EXPECT_THAT(log, testing::HasSubstr(" Largest stress per step (kbar):"));
    EXPECT_THAT(log, testing::HasSubstr("\n Relaxation is converged!"));
}

// cell-relax + small force + large stress: the ionic part converged but the
// cell part did not, so the log shows not-converged and no ASE marker.
TEST_F(IonCellOptimizerTest, CellRelaxNotConverged)
{
    inp.calculation = "cell-relax";
    UnitCell ucell;
    setup_ucell(ucell);
    optimizer.init_relax(natom, inp);
    ModuleBase::matrix force(natom, 3);
    force(0, 0) = 1.0e-4;
    ModuleBase::matrix stress(3, 3);
    stress(0, 0) = 1.0; // ~2.9e5 kbar, above stress_thr
    int force_step = 1;
    int stress_step = 1;
    std::ofstream ofs(log_file);

    const bool done = optimizer.relax_step(1, etot, ucell, force, stress, force_step, stress_step, ofs);
    ofs.close();

    EXPECT_FALSE(done);
    const std::string log = read_log();
    EXPECT_THAT(log, testing::HasSubstr(" Largest stress is "));
    EXPECT_THAT(log, testing::HasSubstr("\n Relaxation is not converged yet!"));
    EXPECT_THAT(log, testing::Not(testing::HasSubstr("Relaxation is converged!")));
}

// Decision: the force history survives a cell change (reset_after_cell_change
// resets force_step and the ionic optimizer, but max_force_history_ must keep
// every ionic step from the whole run). After one unconverged cell step, a
// second ionic step must append, not restart, the history.
TEST_F(IonCellOptimizerTest, CellRelaxForceHistoryNotClearedAcrossCellChange)
{
    inp.calculation = "cell-relax";
    UnitCell ucell;
    setup_ucell(ucell);
    optimizer.init_relax(natom, inp);

    std::ofstream ofs(log_file);
    int force_step = 1;
    int stress_step = 1;

    // Step 1: force converged but stress not -> triggers cell change, which
    // resets force_step to 1 and calls IMM.reset_after_cell_change.
    ModuleBase::matrix force(natom, 3);
    force(0, 0) = 1.0e-4;
    ModuleBase::matrix stress(3, 3);
    stress(0, 0) = 1.0; // ~2.9e5 kbar, above stress_thr -> cell change
    const bool done1 = optimizer.relax_step(1, etot, ucell, force, stress, force_step, stress_step, ofs);
    EXPECT_FALSE(done1);
    EXPECT_EQ(force_step, 1); // reset by the cell change
    ASSERT_EQ(optimizer.get_max_force_history().size(), 1u);

    // Step 2: another ionic step after the reset. The history must accumulate.
    const bool done2 = optimizer.relax_step(2, etot, ucell, force, stress, force_step, stress_step, ofs);
    ofs.close();
    EXPECT_FALSE(done2);
    EXPECT_EQ(optimizer.get_max_force_history().size(), 2u);
}

// Regression test for the review on PR #8067: relax_nmax = 0 is a valid
// dry-run mode. The not-converged branch must stay silent so the running log
// never contains "Relaxation is not converged after reaching relax_nmax!".
TEST_F(IonCellOptimizerTest, DryRunRelaxNmaxZero)
{
    inp.calculation = "relax";
    inp.relax_nmax = 0;
    UnitCell ucell;
    setup_ucell(ucell);
    optimizer.init_relax(natom, inp);
    ModuleBase::matrix force(natom, 3);
    force(0, 0) = 1.0; // would be above threshold, but no step should run
    ModuleBase::matrix stress(3, 3);
    int force_step = 1;
    int stress_step = 1;
    std::ofstream ofs(log_file);

    // istep == relax_nmax == 0 hits the former not-converged branch.
    const bool done = optimizer.relax_step(0, etot, ucell, force, stress, force_step, stress_step, ofs);
    ofs.close();

    EXPECT_TRUE(done);
    const std::string log = read_log();
    EXPECT_THAT(log, testing::Not(testing::HasSubstr("not converged")));
    EXPECT_THAT(log, testing::Not(testing::HasSubstr("relax_nmax = 0 reached")));
}

// cell-relax + fixed_axes = abc: the lattice is fully fixed but the atoms can
// still move. Once the forces converge there is no cell step to run, and the
// log must still carry the unified converged summary (PR #8067 review).
TEST_F(IonCellOptimizerTest, CellRelaxFixedLatticeConverged)
{
    inp.calculation = "cell-relax";
    UnitCell ucell;
    setup_ucell(ucell);
    // fixed_axes = abc: no lattice vector may change.
    ucell.lat_axis_free[0] = 0;
    ucell.lat_axis_free[1] = 0;
    ucell.lat_axis_free[2] = 0;
    optimizer.init_relax(natom, inp);
    ModuleBase::matrix force(natom, 3);
    force(0, 0) = 1.0e-4; // below force_thr
    ModuleBase::matrix stress(3, 3);
    int force_step = 1;
    int stress_step = 1;
    std::ofstream ofs(log_file);

    const bool done = optimizer.relax_step(1, etot, ucell, force, stress, force_step, stress_step, ofs);
    ofs.close();

    EXPECT_TRUE(done);
    const std::string log = read_log();
    EXPECT_THAT(log, testing::HasSubstr(" Largest force is "));
    EXPECT_THAT(log, testing::HasSubstr(" Relaxation method: lbfgs"));
    EXPECT_THAT(log, testing::HasSubstr(" Relaxation converged in 1 step(s)."));
    EXPECT_THAT(log, testing::HasSubstr(" Largest force per step (eV/Angstrom):"));
    EXPECT_THAT(log, testing::HasSubstr("\n Relaxation is converged!"));
    // No cell step ran, so no stress history line and no not-converged marker.
    EXPECT_THAT(log, testing::Not(testing::HasSubstr(" Largest stress per step (kbar):")));
    EXPECT_THAT(log, testing::Not(testing::HasSubstr("not converged")));
}

// Regression test for the no-atom-can-move branch of relax_step: with every
// atom fixed (mbl = 0), relax mode is a valid no-op relaxation. The running
// log must still carry the unified converged summary so that downstream tools
// (e.g. the ASE interface) find an explicit success marker.
TEST_F(IonCellOptimizerTest, RelaxAllAtomsFixedConverged)
{
    inp.calculation = "relax";
    UnitCell ucell;
    setup_ucell(ucell);
    // Fix the only atom in all directions: no atoms are allowed to move.
    ucell.atoms[0].mbl[0] = {0, 0, 0};
    optimizer.init_relax(natom, inp);
    ModuleBase::matrix force(natom, 3);
    ModuleBase::matrix stress(3, 3);
    int force_step = 1;
    int stress_step = 1;
    std::ofstream ofs(log_file);

    // Redirect warnings to a file so the no-op warning can be inspected.
    GlobalV::ofs_warning.open(warning_file);
    const bool done = optimizer.relax_step(1, etot, ucell, force, stress, force_step, stress_step, ofs);
    GlobalV::ofs_warning.close();
    ofs.close();

    EXPECT_TRUE(done);
    const std::string log = read_log();
    EXPECT_THAT(log, testing::HasSubstr(" Relaxation method: lbfgs"));
    EXPECT_THAT(log, testing::HasSubstr(" Relaxation converged in 1 step(s)."));
    EXPECT_THAT(log, testing::HasSubstr("\n Relaxation is converged!"));
    // No ionic step ran, so no force history and no not-converged marker.
    EXPECT_THAT(log, testing::Not(testing::HasSubstr(" Largest force per step (eV/Angstrom):")));
    EXPECT_THAT(log, testing::Not(testing::HasSubstr("not converged")));
    // The no-op branch must warn that no atoms are allowed to move.
    const std::string warning_log = read_warning_log();
    EXPECT_THAT(warning_log, testing::HasSubstr("No atoms are allowed to move!"));
}

// Same fixed-lattice setup, but the forces are above threshold: relaxation is
// genuinely in progress, so the log must not print the converged marker.
TEST_F(IonCellOptimizerTest, CellRelaxFixedLatticeNotConverged)
{
    inp.calculation = "cell-relax";
    UnitCell ucell;
    setup_ucell(ucell);
    ucell.lat_axis_free[0] = 0;
    ucell.lat_axis_free[1] = 0;
    ucell.lat_axis_free[2] = 0;
    optimizer.init_relax(natom, inp);
    ModuleBase::matrix force(natom, 3);
    force(0, 0) = 1.0; // above force_thr
    ModuleBase::matrix stress(3, 3);
    int force_step = 1;
    int stress_step = 1;
    std::ofstream ofs(log_file);

    const bool done = optimizer.relax_step(1, etot, ucell, force, stress, force_step, stress_step, ofs);
    ofs.close();

    EXPECT_FALSE(done);
    const std::string log = read_log();
    EXPECT_THAT(log, testing::HasSubstr("\n Relaxation is not converged yet!"));
    EXPECT_THAT(log, testing::Not(testing::HasSubstr("Relaxation is converged!")));
}
