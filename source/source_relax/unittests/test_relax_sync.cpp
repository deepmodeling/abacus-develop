#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <iomanip>
#include "../relax_sync.h"
#include "source_cell/unitcell.h"
#include "source_io/module_parameter/input_parameter.h"
#include "relax_test.h"
#include <fstream>

// Fixed mock SCF outputs (forces, stresses, energies) feeding relax_step,
// previously read from support/*.txt. Inlining them keeps the test
// independent of the working directory. Values are the first 3 steps of the
// original files, so result_ref is unchanged.
namespace {
constexpr int kRelaxSteps = 3;
constexpr int kRelaxNat = 5;
constexpr int kRelaxNforce = kRelaxSteps * kRelaxNat * 3;  // 45
constexpr int kRelaxNstress = kRelaxSteps * 3 * 3;           // 27

const double kRelaxForces[kRelaxNforce] = {
    4.56e-08, -8.74e-08, 4.0671950e-03,
    -5.59e-08, 1.479e-07, 1.81171376e-02,
    2.942e-07, -5.56e-08, -1.7979905e-03,
    -1.691e-07, 1.612e-07, -1.7977699e-03,
    -1.149e-07, -1.660e-07, -1.85885722e-02,
    6.73e-08, 5.90e-08, 4.8608620e-03,
    -1.632e-07, -1.833e-07, -2.24081615e-02,
    4.93e-08, 1.246e-07, -1.8938742e-03,
    9.84e-08, 4.34e-08, -1.8939930e-03,
    -5.19e-08, -4.38e-08, 2.13351667e-02,
    8.25e-08, 2.095e-07, 4.5297162e-03,
    -5.19e-08, 3.57e-08, -3.6659292e-03,
    4.95e-08, -1.654e-07, -1.8863808e-03,
    4.14e-08, 2.440e-07, -1.8870376e-03,
    -1.215e-07, -3.239e-07, 2.9096315e-03,
};

const double kRelaxStresses[kRelaxNstress] = {
    4.79987e-05, -2.0e-10, 6.0e-10,
    -2.0e-10, 4.80038e-05, -4.0e-10,
    6.0e-10, -4.0e-10, -2.751687e-04,
    1.040725e-04, 2.0e-10, -1.6e-09,
    2.0e-10, 1.040682e-04, -7.0e-10,
    -1.6e-09, -7.0e-10, 1.982021e-04,
    8.00586e-05, 1.0e-10, -1.0e-09,
    1.0e-10, 8.00526e-05, 2.8e-09,
    -1.0e-09, 2.8e-09, -6.4010e-06,
};

const double kRelaxEnergies[kRelaxSteps] = {
    4.79987e-05,
    -2.0e-10,
    6.0e-10,
};
} // namespace

class Test_SETGRAD : public testing::Test
{
    protected:
        Relax rl;
        std::vector<double> result;
        Input_para inp;
        UnitCell ucell;
        std::ofstream ofs;

        void SetUp()
        {
            inp.force_thr  = 0.001;
            inp.calculation = "cell-relax";

            ModuleBase::matrix force_in, stress_in;
            int nat = 3;

            force_in.create(nat,3);
            stress_in.create(3,3);

            force_in(0,0) = 0.1; force_in(0,1) = 0.1; force_in(0,2)= 0.1;
            force_in(1,0) = 0; force_in(1,1) = 0.1; force_in(1,2)= 0.1;
            force_in(2,0) = 0; force_in(2,1) = 0; force_in(2,2)= 0.1;

            stress_in(0,0) = 1; stress_in(0,1) = 1; stress_in(0,2)= 1;
            stress_in(1,0) = 0; stress_in(1,1) = 1; stress_in(1,2)= 1;
            stress_in(2,0) = 0; stress_in(2,1) = 0; stress_in(2,2)= 1;

            ucell.ntype = 1;
            ucell.nat = nat;
            ucell.atoms = new Atom[1];
            ucell.atoms[0].na = nat;
            ucell.omega = 1.0;
            ucell.lat0 = 1.0;
            
            ucell.iat2it = new int[nat];
            ucell.iat2ia = new int[nat];
            ucell.atoms[0].mbl.resize(nat);
            ucell.atoms[0].taud.resize(nat);
            ucell.atoms[0].tau.resize(nat);
            ucell.atoms[0].dis.resize(nat);
            ucell.atoms[0].mag.resize(nat);
            ucell.atoms[0].vel.resize(nat);

            ucell.iat2it[0] = 0;
            ucell.iat2it[1] = 0;
            ucell.iat2it[2] = 0;

            ucell.iat2ia[0] = 0;
            ucell.iat2ia[1] = 1;
            ucell.iat2ia[2] = 2;

            ucell.atoms[0].mbl[0].x = 0;
            ucell.atoms[0].mbl[0].y = 0;
            ucell.atoms[0].mbl[0].z = 1;

            ucell.atoms[0].mbl[1].x = 0;
            ucell.atoms[0].mbl[1].y = 1;
            ucell.atoms[0].mbl[1].z = 0;

            ucell.atoms[0].mbl[2].x = 1;
            ucell.atoms[0].mbl[2].y = 0;
            ucell.atoms[0].mbl[2].z = 0;

            ucell.atoms[0].taud[0] = 0.0;
            ucell.atoms[0].taud[1] = 0.0;
            ucell.atoms[0].taud[2] = 0.0;

            ucell.atoms[0].tau.resize(nat);

            ucell.lat_axis_free[0] = 1;
            ucell.lat_axis_free[1] = 1;
            ucell.lat_axis_free[2] = 1;

            rl.init_relax(nat, inp);
            rl.relax_step(ucell,force_in,stress_in,0.0, ofs);

            for(int i=0;i<3;i++)
            {
                result.push_back(ucell.atoms[0].taud[i].x);
                result.push_back(ucell.atoms[0].taud[i].y);
                result.push_back(ucell.atoms[0].taud[i].z);
            }
            push_result();

            //reset lattice vector
            ucell.latvec.Identity();
            inp.fixed_axes = "shape";
            rl.init_relax(nat, inp);
            rl.relax_step(ucell,force_in,stress_in,0.0, ofs);
            push_result();

            //reset lattice vector
            ucell.latvec.Identity();
            inp.fixed_axes = "volume";
            rl.init_relax(nat, inp);
            rl.relax_step(ucell,force_in,stress_in,0.0, ofs);
            push_result();

            //reset lattice vector
            ucell.latvec.Identity();
            inp.fixed_axes = "a"; //anything other than "None"
            inp.fixed_ibrav = true;
            ucell.latName = "sc";
            ucell.lat_axis_free[0] = 0;
            ucell.lat_axis_free[1] = 0;
            ucell.lat_axis_free[2] = 0;
            rl.init_relax(nat, inp);
            rl.relax_step(ucell,force_in,stress_in,0.0, ofs);
            
            push_result();
        }

        void push_result()
        {
            result.push_back(ucell.latvec.e11);
            result.push_back(ucell.latvec.e12);
            result.push_back(ucell.latvec.e13);
            result.push_back(ucell.latvec.e21);
            result.push_back(ucell.latvec.e22);
            result.push_back(ucell.latvec.e23);
            result.push_back(ucell.latvec.e31);
            result.push_back(ucell.latvec.e32);
            result.push_back(ucell.latvec.e33);             
        }

};

TEST_F(Test_SETGRAD, relax_new)
{
    std::vector<double> result_ref = 
    {
        0,           0,            0.24293434145, 
        0,           0.242934341453,           0,
        0,           0,           0,
        //paramter for taud
        1.2267616333,0.2267616333,0.22676163333, 
        0,            1.2267616333 ,0.2267616333,
        0,            0,           1.22676163333,
        // paramter for fisrt time after relaxation
        1.3677603495, 0,            0,
        0,            1.36776034956, 0,
        0,            0,            1.36776034956,
        // paramter for second time after relaxation
        1.3677603495  ,0.3633367476,0.36333674766,
        0,            1.3677603495 ,0.36333674766,
        0,            0,            1.3677603495 ,
        // paramter for third time after relaxation
        1,0,0,0,1,0,0,0,1
        // paramter for fourth time after relaxation
    };
    for(int i=0;i<result.size();i++)
    {
        EXPECT_NEAR(result[i],result_ref[i],1e-8);
    }
}

class Test_RELAX : public testing::Test
{
    protected:
        Relax rl;
        std::vector<double> result;
        Input_para inp;
        UnitCell ucell;
        std::ofstream ofs;

        void SetUp()
        {
            int nstep = 3;
            int nat = 5;
            inp.calculation = "cell-relax";
            inp.force_thr = 0.001;
            inp.fixed_axes = "a";
            inp.stress_thr = 0.01;
            inp.fixed_ibrav = false;

            this->setup_cell();

            ModuleBase::matrix force_in, stress_in;
            force_in.create(nat,3);
            stress_in.create(3,3);

            rl.init_relax(nat, inp);

            for(int istep=0;istep<nstep;istep++)
            {
                const int fbase = istep * nat * 3;
                for(int i=0;i<nat;i++)
                {
                    for(int j=0;j<3;j++)
                    {
                        force_in(i,j) = kRelaxForces[fbase + i*3 + j];
                    }
                }
                const int sbase = istep * 9;
                for(int i=0;i<3;i++)
                {
                    for(int j=0;j<3;j++)
                    {
                        stress_in(i,j) = kRelaxStresses[sbase + i*3 + j];
                    }
                }

                const double energy = kRelaxEnergies[istep];

                inp.fixed_ibrav = false;
                rl.relax_step(ucell,force_in,stress_in,energy, ofs);
                result.push_back(ucell.atoms[0].taud[0].x);
                result.push_back(ucell.atoms[0].taud[0].y);
                result.push_back(ucell.atoms[0].taud[0].z);
                result.push_back(ucell.atoms[1].taud[0].x);
                result.push_back(ucell.atoms[1].taud[0].y);
                result.push_back(ucell.atoms[1].taud[0].z);
                result.push_back(ucell.atoms[2].taud[0].x);
                result.push_back(ucell.atoms[2].taud[0].y);
                result.push_back(ucell.atoms[2].taud[0].z);
                result.push_back(ucell.atoms[2].taud[1].x);
                result.push_back(ucell.atoms[2].taud[1].y);
                result.push_back(ucell.atoms[2].taud[1].z);
                result.push_back(ucell.atoms[2].taud[2].x);
                result.push_back(ucell.atoms[2].taud[2].y);
                result.push_back(ucell.atoms[2].taud[2].z);
                result.push_back(ucell.latvec.e11);
                result.push_back(ucell.latvec.e12);
                result.push_back(ucell.latvec.e13);
                result.push_back(ucell.latvec.e21);
                result.push_back(ucell.latvec.e22);
                result.push_back(ucell.latvec.e23);
                result.push_back(ucell.latvec.e31);
                result.push_back(ucell.latvec.e32);
                result.push_back(ucell.latvec.e33);
            }
        }

        void setup_cell()
        {
            int ntype = 3, nat = 5;
            ucell.ntype = ntype;
            ucell.nat = nat;

            ucell.omega = 452.590903143121;
            ucell.lat0 = 1.8897259886;
            ucell.iat2it = new int[nat];
            ucell.iat2ia = new int[nat];
            ucell.iat2it[0] = 0;
            ucell.iat2it[1] = 1;
            ucell.iat2it[2] = 2;
            ucell.iat2it[3] = 2;
            ucell.iat2it[4] = 2;

            ucell.iat2ia[0] = 0;
            ucell.iat2ia[1] = 0;
            ucell.iat2ia[2] = 0;
            ucell.iat2ia[3] = 1;
            ucell.iat2ia[4] = 2;

            ucell.atoms = new Atom[ntype];
            ucell.atoms[0].na = 1;
            ucell.atoms[1].na = 1;
            ucell.atoms[2].na = 3;
            
            for(int i=0;i<ntype;i++)
            {
                int na = ucell.atoms[i].na;
                ucell.atoms[i].label="test";
                ucell.atoms[i].mbl.resize(na);
                ucell.atoms[i].taud.resize(na);
                ucell.atoms[i].tau.resize(na);
                ucell.atoms[i].dis.resize(na);
                ucell.atoms[i].mag.resize(na);
                ucell.atoms[i].vel.resize(na);
                for (int j=0;j<na;j++)
                {
                    ucell.atoms[i].mbl[j] = {1,1,1};
                }
            }
            ucell.atoms[0].taud[0] = {0.5,0.5,0.00413599999956205};
            ucell.atoms[1].taud[0] = {0  ,0  ,0.524312999999893  };
            ucell.atoms[2].taud[0] = {0  ,0.5,0.479348999999274  };
            ucell.atoms[2].taud[1] = {0.5,0  ,0.479348999999274  };
            ucell.atoms[2].taud[2] = {0  ,0  ,0.958854000000429  };
            
            ucell.lat_axis_free[0] = 1;
            ucell.lat_axis_free[1] = 1;
            ucell.lat_axis_free[2] = 1;

            ucell.latvec.e11 = 3.96;
            ucell.latvec.e12 = 0;
            ucell.latvec.e13 = 0;
            ucell.latvec.e21 = 0;
            ucell.latvec.e22 = 3.96;
            ucell.latvec.e23 = 0;
            ucell.latvec.e31 = 0;
            ucell.latvec.e32 = 0;
            ucell.latvec.e33 = 4.2768;
        }
};

TEST_F(Test_RELAX, relax_new)
{
    int size = 72;
    double tmp;
    std::vector<double> result_ref=
    {
        0.5000000586,0.4999998876,0.009364595811,
        0.9999999281,1.901333279e-07,0.5476035454,
        3.782097706e-07,0.4999999285,0.4770375874,
        0.4999997826,2.072311863e-07,0.477037871,
        0.9999998523,0.9999997866,0.9349574003,
        // paramter for taud after first relaxation
        4.006349654,-1.93128788e-07,5.793863639e-07,
        -1.93128788e-07,4.006354579,-3.86257576e-07,
        6.757962549e-07,-4.505308366e-07,3.966870038,
        // paramter for latvec after first relaxation
        0.5000000566,0.4999998916,0.009177239183,
        0.9999999308,1.832935626e-07,0.5467689737,
        3.647769323e-07,0.4999999311,0.4771204124,
        0.4999997903,1.998879545e-07,0.4771206859,
        0.9999998574,0.9999997943,0.9358136888,
        // paramter for taud after second relaxation
        3.999761277,-1.656764727e-07,4.97029418e-07,
        -1.656764727e-07,3.999765501,-3.313529453e-07,
        5.797351131e-07,-3.864900754e-07,4.010925071,
        // paramter for latvec after second relaxation
        0.500000082,0.4999999574,0.01057784352,
        0.9999999149,1.939640249e-07,0.5455830599,
        3.795967781e-07,0.4999998795,0.4765373919,
        0.4999998037,2.756298268e-07,0.4765374602,
        0.9999998196,0.9999996936,0.9367652445,
        // paramter for taud after third relaxation
        4.017733155,-1.420363309e-07,2.637046077e-07,
        -1.420364243e-07,4.017735987,3.126225134e-07,
        3.479123171e-07,2.578467568e-07,4.011674933
        // paramter for latvec after third relaxation
    };
    for(int i=0;i<size;i++)
    {
        EXPECT_NEAR(result_ref[i],result[i],1e-8);
    }
}

// Drive the simultaneous CG path (relax_sync) to convergence in one step and
// check the unified summary printed before "Relaxation is converged!".
TEST(RelaxSyncSummary, ConvergedPrintsSummary)
{
    const int nat = 1;
    Input_para inp;
    inp.calculation = "relax";
    inp.relax_method = {"cg", "2"};
    inp.force_thr = 0.001;   // force_thr_eva ~ 0.0257 eV/Angstrom
    inp.force_thr_ev = inp.force_thr * 13.6058 / 0.529177;

    // iat2it/iat2ia are owned by UnitCell's internal Statistics member, whose
    // destructor releases them; do not delete them again (mirror the other
    // tests in this file).
    UnitCell ucell;
    ucell.ntype = 1;
    ucell.nat = nat;
    ucell.atoms = new Atom[1];
    ucell.atoms[0].na = nat;
    ucell.atoms[0].label = "Si";
    ucell.omega = 1.0;
    ucell.lat0 = 1.0;
    ucell.iat2it = new int[nat];
    ucell.iat2ia = new int[nat];
    ucell.iat2it[0] = 0;
    ucell.iat2ia[0] = 0;
    ucell.atoms[0].mbl.resize(nat);
    ucell.atoms[0].taud.resize(nat);
    ucell.atoms[0].tau.resize(nat);
    ucell.atoms[0].dis.resize(nat);
    ucell.atoms[0].mag.resize(nat);
    ucell.atoms[0].vel.resize(nat);
    ucell.atoms[0].mbl[0] = {1, 1, 1};
    ucell.atoms[0].taud[0] = {0.0, 0.0, 0.0};
    ucell.latvec.Identity();

    // Well below the threshold -> converged immediately.
    ModuleBase::matrix force_in(nat, 3);
    ModuleBase::matrix stress_in(3, 3);
    force_in(0, 0) = 1.0e-4;

    Relax rl;
    rl.init_relax(nat, inp);

    std::ofstream ofs("./running_relax_sync_test.log");
    const bool done = rl.relax_step(ucell, force_in, stress_in, 0.0, ofs);
    ofs.close();

    EXPECT_TRUE(done);

    std::ifstream ifs("./running_relax_sync_test.log");
    const std::string log((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("./running_relax_sync_test.log");

    EXPECT_THAT(log, testing::HasSubstr(" Relaxation method: cg"));
    EXPECT_THAT(log, testing::HasSubstr(" Relaxation converged in 1 step(s)."));
    EXPECT_THAT(log, testing::HasSubstr(" Largest force per step (eV/Angstrom):"));
    EXPECT_THAT(log, testing::HasSubstr("\n Relaxation is converged!"));

    delete[] ucell.atoms;
}

// First step above threshold (not converged), second step below (converged).
// Locks the converged_step = istep+1 accounting and the per-step history.
TEST(RelaxSyncSummary, TwoStepConvergence)
{
    const int nat = 1;
    Input_para inp;
    inp.calculation = "relax";
    inp.relax_method = {"cg", "2"};
    inp.force_thr = 0.001;   // force_thr_eva ~ 0.0257 eV/Angstrom
    inp.force_thr_ev = inp.force_thr * 13.6058 / 0.529177;

    UnitCell ucell;
    ucell.ntype = 1;
    ucell.nat = nat;
    ucell.atoms = new Atom[1];
    ucell.atoms[0].na = nat;
    ucell.atoms[0].label = "Si";
    ucell.omega = 1.0;
    ucell.lat0 = 1.0;
    ucell.iat2it = new int[nat];
    ucell.iat2ia = new int[nat];
    ucell.iat2it[0] = 0;
    ucell.iat2ia[0] = 0;
    ucell.atoms[0].mbl.resize(nat);
    ucell.atoms[0].taud.resize(nat);
    ucell.atoms[0].tau.resize(nat);
    ucell.atoms[0].dis.resize(nat);
    ucell.atoms[0].mag.resize(nat);
    ucell.atoms[0].vel.resize(nat);
    ucell.atoms[0].mbl[0] = {1, 1, 1};
    ucell.atoms[0].taud[0] = {0.0, 0.0, 0.0};
    ucell.latvec.Identity();

    ModuleBase::matrix stress_in(3, 3);

    Relax rl;
    rl.init_relax(nat, inp);

    // Step 1: force above threshold -> not converged. The magnitude must stay
    // small enough that the CG step keeps the atom inside the unit cell
    // (latvec is the identity), otherwise ABACUS aborts with "Movement of
    // atom is larger than the cell length".
    ModuleBase::matrix force_big(nat, 3);
    force_big(0, 0) = 0.01; // ~0.257 eV/Angstrom, above threshold
    std::ofstream ofs("./running_relax_sync_test.log");
    const bool done1 = rl.relax_step(ucell, force_big, stress_in, 0.0, ofs);
    EXPECT_FALSE(done1);

    // Step 2: small force -> converged.
    ModuleBase::matrix force_small(nat, 3);
    force_small(0, 0) = 1.0e-4;
    const bool done2 = rl.relax_step(ucell, force_small, stress_in, 0.0, ofs);
    ofs.close();
    EXPECT_TRUE(done2);

    std::ifstream ifs("./running_relax_sync_test.log");
    const std::string log((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("./running_relax_sync_test.log");

    EXPECT_THAT(log, testing::HasSubstr("\n Relaxation is not converged yet!"));
    EXPECT_THAT(log, testing::HasSubstr(" Relaxation converged in 2 step(s)."));
    EXPECT_THAT(log, testing::HasSubstr("\n Relaxation is converged!"));

    // Each relax_step pushes exactly one force entry.
    ASSERT_EQ(rl.get_max_force_history().size(), 2u);

    delete[] ucell.atoms;
}

// cell-relax converged: the summary must also carry the stress history line.
TEST(RelaxSyncSummary, CellRelaxConvergedPrintsStressHistory)
{
    const int nat = 1;
    Input_para inp;
    inp.calculation = "cell-relax";
    inp.relax_method = {"cg", "2"};
    inp.force_thr = 0.001;
    inp.force_thr_ev = inp.force_thr * 13.6058 / 0.529177;
    inp.stress_thr = 0.5; // kbar

    UnitCell ucell;
    ucell.ntype = 1;
    ucell.nat = nat;
    ucell.atoms = new Atom[1];
    ucell.atoms[0].na = nat;
    ucell.atoms[0].label = "Si";
    ucell.omega = 1.0;
    ucell.lat0 = 1.0;
    ucell.iat2it = new int[nat];
    ucell.iat2ia = new int[nat];
    ucell.iat2it[0] = 0;
    ucell.iat2ia[0] = 0;
    ucell.atoms[0].mbl.resize(nat);
    ucell.atoms[0].taud.resize(nat);
    ucell.atoms[0].tau.resize(nat);
    ucell.atoms[0].dis.resize(nat);
    ucell.atoms[0].mag.resize(nat);
    ucell.atoms[0].vel.resize(nat);
    ucell.atoms[0].mbl[0] = {1, 1, 1};
    ucell.atoms[0].taud[0] = {0.0, 0.0, 0.0};
    ucell.latvec.Identity();
    ucell.lat_axis_free[0] = 1;
    ucell.lat_axis_free[1] = 1;
    ucell.lat_axis_free[2] = 1;

    ModuleBase::matrix force_in(nat, 3);
    ModuleBase::matrix stress_in(3, 3);
    force_in(0, 0) = 1.0e-4; // force converged
    // zero stress -> below stress_thr -> cell converged too

    Relax rl;
    rl.init_relax(nat, inp);

    std::ofstream ofs("./running_relax_sync_test.log");
    const bool done = rl.relax_step(ucell, force_in, stress_in, 0.0, ofs);
    ofs.close();
    EXPECT_TRUE(done);

    std::ifstream ifs("./running_relax_sync_test.log");
    const std::string log((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("./running_relax_sync_test.log");

    EXPECT_THAT(log, testing::HasSubstr(" Largest stress is "));
    EXPECT_THAT(log, testing::HasSubstr(" Largest stress per step (kbar):"));
    EXPECT_THAT(log, testing::HasSubstr("\n Relaxation is converged!"));

    delete[] ucell.atoms;
}

// ---------------------------------------------------------------------------
// Behavior tests locking move_cell_ions branches before its refactor.
// All drive the public relax_step and assert only externally observable
// UnitCell state (latvec / taud / omega), so the private-method split must
// keep these numbers bit-stable.
// ---------------------------------------------------------------------------

namespace
{
// Build a minimal 2-type, 3-atom cell with the given per-atom move flags.
// Atom 0 is type 0 (1 atom), atoms 1-2 are type 1 (2 atoms). taud start at 0.
// UnitCell is non-copyable, so the caller supplies the object to fill.
void make_two_type_cell(UnitCell& ucell, const int mbl_flat[9])
{
    ucell.ntype = 2;
    ucell.nat = 3;
    ucell.omega = 8.0;   // 2x2x2 cell
    ucell.lat0 = 1.0;
    ucell.atoms = new Atom[2];
    ucell.atoms[0].na = 1;
    ucell.atoms[1].na = 2;
    for (int t = 0; t < 2; t++)
    {
        const int na = ucell.atoms[t].na;
        ucell.atoms[t].label = "X";
        ucell.atoms[t].mbl.resize(na);
        ucell.atoms[t].taud.resize(na);
        ucell.atoms[t].tau.resize(na);
        ucell.atoms[t].dis.resize(na);
        ucell.atoms[t].mag.resize(na);
        ucell.atoms[t].vel.resize(na);
    }
    // iat: 0 -> (type0, ia0); 1 -> (type1, ia0); 2 -> (type1, ia1)
    ucell.iat2it = new int[3];
    ucell.iat2ia = new int[3];
    ucell.iat2it[0] = 0; ucell.iat2it[1] = 1; ucell.iat2it[2] = 1;
    ucell.iat2ia[0] = 0; ucell.iat2ia[1] = 0; ucell.iat2ia[2] = 1;

    // mbl flags, laid out per atom
    ucell.atoms[0].mbl[0] = {mbl_flat[0], mbl_flat[1], mbl_flat[2]};
    ucell.atoms[1].mbl[0] = {mbl_flat[3], mbl_flat[4], mbl_flat[5]};
    ucell.atoms[1].mbl[1] = {mbl_flat[6], mbl_flat[7], mbl_flat[8]};

    ucell.atoms[0].taud[0] = {0.0, 0.0, 0.0};
    ucell.atoms[1].taud[0] = {0.0, 0.0, 0.0};
    ucell.atoms[1].taud[1] = {0.0, 0.0, 0.0};

    ucell.latvec.Identity();
    ucell.latvec.e11 = 2.0;
    ucell.latvec.e22 = 2.0;
    ucell.latvec.e33 = 2.0;
    ucell.lat_axis_free[0] = 1;
    ucell.lat_axis_free[1] = 1;
    ucell.lat_axis_free[2] = 1;
}

void free_two_type_cell(UnitCell& ucell)
{
    delete[] ucell.atoms;
    // iat2it / iat2ia are owned by the mock; release explicitly.
    delete[] ucell.iat2it;
    delete[] ucell.iat2ia;
    ucell.iat2it = nullptr;
    ucell.iat2ia = nullptr;
}
} // namespace

// fixed_axes="shape": only the hydrostatic part of the stress survives, so the
// lattice scales isotropically (off-diagonals stay zero, diagonal scales equal).
TEST(RelaxSyncCellMove, FixedShapeIsotropicScale)
{
    const int mbl[9] = {1,1,1, 1,1,1, 1,1,1};
    UnitCell ucell;
    make_two_type_cell(ucell, mbl);

    Input_para inp;
    inp.calculation = "cell-relax";
    inp.force_thr = 0.0;      // force never converged -> take the CG move
    inp.force_thr_ev = 0.0;
    inp.stress_thr = 0.0;     // stress never converged either
    inp.fixed_axes = "shape";
    inp.fixed_ibrav = false;

    ModuleBase::matrix force_in(3, 3);
    force_in(0, 0) = 0.01;    // nonzero -> not converged -> performs the move
    ModuleBase::matrix stress_in(3, 3);
    stress_in(0, 0) = 1.0; stress_in(1, 1) = 2.0; stress_in(2, 2) = 3.0;
    stress_in(0, 1) = 0.5; // anisotropic + off-diagonal parts must be dropped

    Relax rl;
    rl.init_relax(3, inp);
    std::ofstream ofs("./running_relax_sync_shape.log");
    rl.relax_step(ucell, force_in, stress_in, 0.0, ofs);
    ofs.close();
    std::remove("./running_relax_sync_shape.log");

    // Off-diagonal lattice components unchanged (remain zero).
    EXPECT_NEAR(ucell.latvec.e12, 0.0, 1e-12);
    EXPECT_NEAR(ucell.latvec.e13, 0.0, 1e-12);
    EXPECT_NEAR(ucell.latvec.e21, 0.0, 1e-12);
    EXPECT_NEAR(ucell.latvec.e23, 0.0, 1e-12);
    EXPECT_NEAR(ucell.latvec.e31, 0.0, 1e-12);
    EXPECT_NEAR(ucell.latvec.e32, 0.0, 1e-12);
    // Isotropic scaling: all three diagonals shifted by the same amount.
    EXPECT_NEAR(ucell.latvec.e11 - 2.0, ucell.latvec.e22 - 2.0, 1e-12);
    EXPECT_NEAR(ucell.latvec.e22 - 2.0, ucell.latvec.e33 - 2.0, 1e-12);

    free_two_type_cell(ucell);
}

// ---------------------------------------------------------------------------
// Max-step termination on the Relax (simultaneous) path.
//
// Unlike IonCellOptimizer::relax_step (which checks istep == relax_nmax and
// returns true without moving), Relax::relax_step has no max-step guard.
// When the driver loop hits relax_nmax, the last relax_step still calls
// move_cell_ions, so the geometry IS updated and the force/stress in the
// driver become stale. The driver then sets geometry_evaluated=false.
//
// This test locks that observable behavior: relax_step moves taud and returns
// false (not converged).
// ---------------------------------------------------------------------------
TEST(RelaxSyncMaxStep, GeometryMovedAndNotConverged)
{
    const int nat = 1;
    Input_para inp;
    inp.calculation = "relax";
    inp.relax_method = {"cg", "2"};
    inp.force_thr = 0.0;      // never converged -> always take the move
    inp.force_thr_ev = 0.0;
    inp.fixed_axes = "None";
    inp.fixed_ibrav = false;

    UnitCell ucell;
    ucell.ntype = 1;
    ucell.nat = nat;
    ucell.atoms = new Atom[1];
    ucell.atoms[0].na = nat;
    ucell.atoms[0].label = "Si";
    ucell.omega = 1.0;
    ucell.lat0 = 1.0;
    ucell.iat2it = new int[nat];
    ucell.iat2ia = new int[nat];
    ucell.iat2it[0] = 0;
    ucell.iat2ia[0] = 0;
    ucell.atoms[0].mbl.resize(nat);
    ucell.atoms[0].taud.resize(nat);
    ucell.atoms[0].tau.resize(nat);
    ucell.atoms[0].dis.resize(nat);
    ucell.atoms[0].mag.resize(nat);
    ucell.atoms[0].vel.resize(nat);
    ucell.atoms[0].mbl[0] = {1, 1, 1};
    ucell.atoms[0].taud[0] = {0.0, 0.0, 0.0};
    ucell.latvec.Identity();

    ModuleBase::matrix force_in(nat, 3);
    ModuleBase::matrix stress_in(3, 3);
    force_in(0, 0) = 0.01;   // nonzero -> not converged -> move

    Relax rl;
    rl.init_relax(nat, inp);
    std::ofstream ofs("./running_relax_sync_maxstep.log");
    const bool done = rl.relax_step(ucell, force_in, stress_in, 0.0, ofs);
    ofs.close();
    std::remove("./running_relax_sync_maxstep.log");

    EXPECT_FALSE(done);
    // The geometry was moved: taud is no longer the initial value.
    EXPECT_NE(ucell.atoms[0].taud[0].x, 0.0);
    EXPECT_TRUE(ucell.ionic_position_updated);

    delete[] ucell.atoms;
}

// fixed_axes="volume": the cell may change shape but its volume is re-scaled
// back to the original omega after each move (relax_sync.cpp Step 1).
TEST(RelaxSyncCellMove, FixedVolumePreservesOmega)
{
    const int mbl[9] = {1,1,1, 1,1,1, 1,1,1};
    UnitCell ucell;
    make_two_type_cell(ucell, mbl);
    const double omega0 = ucell.omega;

    Input_para inp;
    inp.calculation = "cell-relax";
    inp.force_thr = 0.0;
    inp.force_thr_ev = 0.0;
    inp.stress_thr = 0.0;
    inp.fixed_axes = "volume";
    inp.fixed_ibrav = false;

    ModuleBase::matrix force_in(3, 3);
    force_in(0, 0) = 0.01;    // nonzero -> not converged
    ModuleBase::matrix stress_in(3, 3);
    stress_in(0, 0) = 1.0; stress_in(1, 1) = -0.5; stress_in(0, 1) = 0.3;

    Relax rl;
    rl.init_relax(3, inp);
    std::ofstream ofs("./running_relax_sync_vol.log");
    rl.relax_step(ucell, force_in, stress_in, 0.0, ofs);
    ofs.close();
    std::remove("./running_relax_sync_vol.log");

    const double omega_new = std::abs(ucell.latvec.Det()) * ucell.lat0 * ucell.lat0 * ucell.lat0;
    EXPECT_NEAR(omega_new, omega0, 1e-8);

    free_two_type_cell(ucell);
}

// lat_axis_free = {1,0,0}: only the first lattice row may move; rows 2 and 3
// must remain exactly at their initial values (relax_sync.cpp 572-589).
TEST(RelaxSyncCellMove, PartialAxisFreeKeepsLockedRows)
{
    const int mbl[9] = {1,1,1, 1,1,1, 1,1,1};
    UnitCell ucell;
    make_two_type_cell(ucell, mbl);
    ucell.lat_axis_free[0] = 1;
    ucell.lat_axis_free[1] = 0; // lock row 2
    ucell.lat_axis_free[2] = 0; // lock row 3

    Input_para inp;
    inp.calculation = "cell-relax";
    inp.force_thr = 0.0;
    inp.force_thr_ev = 0.0;
    inp.stress_thr = 0.0;
    inp.fixed_axes = "None";
    inp.fixed_ibrav = false;

    ModuleBase::matrix force_in(3, 3);
    force_in(0, 0) = 0.01;    // nonzero -> not converged
    ModuleBase::matrix stress_in(3, 3);
    stress_in(0, 0) = 1.0; stress_in(1, 1) = 1.0; stress_in(2, 2) = 1.0;
    stress_in(0, 1) = 0.2; stress_in(1, 0) = 0.2;

    Relax rl;
    rl.init_relax(3, inp);
    std::ofstream ofs("./running_relax_sync_axis.log");
    rl.relax_step(ucell, force_in, stress_in, 0.0, ofs);
    ofs.close();
    std::remove("./running_relax_sync_axis.log");

    // Locked rows unchanged.
    EXPECT_NEAR(ucell.latvec.e21, 0.0, 1e-12);
    EXPECT_NEAR(ucell.latvec.e22, 2.0, 1e-12);
    EXPECT_NEAR(ucell.latvec.e23, 0.0, 1e-12);
    EXPECT_NEAR(ucell.latvec.e31, 0.0, 1e-12);
    EXPECT_NEAR(ucell.latvec.e32, 0.0, 1e-12);
    EXPECT_NEAR(ucell.latvec.e33, 2.0, 1e-12);

    free_two_type_cell(ucell);
}

// Mixed mbl flags across two atom types: each atom's direct-coordinate
// displacement must be zero along the constrained (mbl==0) directions and
// nonzero only along free ones (relax_sync.cpp Step 2&3, iat2it/iat2ia + mbl).
TEST(RelaxSyncCellMove, MixedMblMovesOnlyFreeComponents)
{
    // atom0: free x only; atom1(type1,ia0): free y only; atom2(type1,ia1): free z only
    const int mbl[9] = {1,0,0, 0,1,0, 0,0,1};
    UnitCell ucell;
    make_two_type_cell(ucell, mbl);

    Input_para inp;
    inp.calculation = "relax";  // ions only, keep lattice fixed for clarity
    inp.force_thr = 0.0;
    inp.force_thr_ev = 0.0;
    inp.fixed_axes = "None";
    inp.fixed_ibrav = false;

    // Uniform force on every atom/component; mbl must mask the movement.
    ModuleBase::matrix force_in(3, 3);
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            force_in(i, j) = 0.05;
        }
    }
    ModuleBase::matrix stress_in(3, 3);

    Relax rl;
    rl.init_relax(3, inp);
    std::ofstream ofs("./running_relax_sync_mbl.log");
    rl.relax_step(ucell, force_in, stress_in, 0.0, ofs);
    ofs.close();
    std::remove("./running_relax_sync_mbl.log");

    // atom0: only x moved
    EXPECT_NE(ucell.atoms[0].taud[0].x, 0.0);
    EXPECT_NEAR(ucell.atoms[0].taud[0].y, 0.0, 1e-12);
    EXPECT_NEAR(ucell.atoms[0].taud[0].z, 0.0, 1e-12);
    // atom1: only y moved
    EXPECT_NEAR(ucell.atoms[1].taud[0].x, 0.0, 1e-12);
    EXPECT_NE(ucell.atoms[1].taud[0].y, 0.0);
    EXPECT_NEAR(ucell.atoms[1].taud[0].z, 0.0, 1e-12);
    // atom2: only z moved
    EXPECT_NEAR(ucell.atoms[1].taud[1].x, 0.0, 1e-12);
    EXPECT_NEAR(ucell.atoms[1].taud[1].y, 0.0, 1e-12);
    EXPECT_NE(ucell.atoms[1].taud[1].z, 0.0);

    free_two_type_cell(ucell);
}