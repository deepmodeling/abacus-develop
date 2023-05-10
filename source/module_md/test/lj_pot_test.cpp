#include "gtest/gtest.h"
#include "module_esolver/esolver_lj.h"
#include "module_md/md_func.h"
#include "setcell.h"

#define doublethreshold 1e-12

/************************************************
 *  unit test of functions in esolver_lj.h
 ***********************************************/

/**
 * - Tested Function
 *   - ESolver_LJ::Run
 *     - calculate energy, force, virial for lj pot
 */

class LJ_pot_test : public testing::Test
{
  protected:
    ModuleBase::Vector3<double>* force;
    ModuleBase::matrix stress;
    double potential;
    int natom;

    void SetUp()
    {
        UnitCell ucell;
        Setcell::setupcell(ucell);

        natom = ucell.nat;
        force = new ModuleBase::Vector3<double>[natom];
        stress.create(3, 3);

        Setcell::parameters();

        ModuleESolver::ESolver* p_esolver = new ModuleESolver::ESolver_LJ();
        p_esolver->Init(INPUT, ucell);
        MD_func::force_virial(p_esolver, 0, ucell, potential, force, true, stress);
    }

    void TearDown()
    {
        delete[] force;
    }
};

TEST_F(LJ_pot_test, potential)
{
    EXPECT_NEAR(potential, -0.011957818623534381, doublethreshold);
}

TEST_F(LJ_pot_test, force)
{
    EXPECT_NEAR(force[0].x, 0.00049817733089377704, doublethreshold);
    EXPECT_NEAR(force[0].y, 0.00082237246837022328, doublethreshold);
    EXPECT_NEAR(force[0].z, -3.0493186101154812e-20, doublethreshold);
    EXPECT_NEAR(force[1].x, -0.00064758615201580339, doublethreshold);
    EXPECT_NEAR(force[1].y, -0.00066924999462089304, doublethreshold);
    EXPECT_NEAR(force[1].z, -1.8634724839594607e-20, doublethreshold);
    EXPECT_NEAR(force[2].x, -0.00035411224839165616, doublethreshold);
    EXPECT_NEAR(force[2].y, 0.0008091080910885112, doublethreshold);
    EXPECT_NEAR(force[2].z, -1.1858461261560205e-20, doublethreshold);
    EXPECT_NEAR(force[3].x, 0.00050352106951368229, doublethreshold);
    EXPECT_NEAR(force[3].y, -0.00096223056483784122, doublethreshold);
    EXPECT_NEAR(force[3].z, 2.0328790734103208e-20, doublethreshold);
}

TEST_F(LJ_pot_test, stress)
{
    EXPECT_NEAR(stress(0, 0), 8.0360222227631859e-07, doublethreshold);
    EXPECT_NEAR(stress(0, 1), 1.7207745586539077e-07, doublethreshold);
    EXPECT_NEAR(stress(0, 2), 0, doublethreshold);
    EXPECT_NEAR(stress(1, 0), 1.7207745586539077e-07, doublethreshold);
    EXPECT_NEAR(stress(1, 1), 1.0630708613186662e-06, doublethreshold);
    EXPECT_NEAR(stress(1, 2), -1.1858461261560206e-22, doublethreshold);
    EXPECT_NEAR(stress(2, 0), 0, doublethreshold);
    EXPECT_NEAR(stress(2, 1), -1.1858461261560206e-22, doublethreshold);
    EXPECT_NEAR(stress(2, 2), 6.4275429572682057e-07, doublethreshold);
}