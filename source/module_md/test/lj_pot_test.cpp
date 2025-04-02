#include "gtest/gtest.h"
#define private public
#include "module_parameter/parameter.h"
#include "module_esolver/esolver_lj.h"
#include "module_md/md_func.h"
#include "setcell.h"
#undef private
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
    UnitCell ucell;
    Input_para input;

    void SetUp()
    {
        Setcell::setupcell(ucell);

        natom = ucell.nat;
        force = new ModuleBase::Vector3<double>[natom];
        stress.create(3, 3);

        Setcell::parameters(input);
    }

    void TearDown()
    {
        delete[] force;
    }
};

TEST_F(LJ_pot_test, potential)
{
    ModuleESolver::ESolver* p_esolver = new ModuleESolver::ESolver_LJ();
    p_esolver->before_all_runners(ucell, input);
    MD_func::force_virial(p_esolver, 0, ucell, potential, force, true, stress);
    EXPECT_NEAR(potential, -0.011957818623534381, doublethreshold);
}

TEST_F(LJ_pot_test, force)
{
    ModuleESolver::ESolver* p_esolver = new ModuleESolver::ESolver_LJ();
    p_esolver->before_all_runners(ucell, input);
    MD_func::force_virial(p_esolver, 0, ucell, potential, force, true, stress);
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
    ModuleESolver::ESolver* p_esolver = new ModuleESolver::ESolver_LJ();
    p_esolver->before_all_runners(ucell, input);
    MD_func::force_virial(p_esolver, 0, ucell, potential, force, true, stress);
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

TEST_F(LJ_pot_test, RcutSearchRadius)
{
    ModuleESolver::ESolver_LJ* p_esolver = new ModuleESolver::ESolver_LJ();
    ucell.ntype = 2;
    std::vector<double> rcut = {3.0};
    p_esolver->rcut_search_radius(ucell.ntype, rcut);

    for (int i = 0; i < ucell.ntype; i++)
    {
        for (int j = 0; j < ucell.ntype; j++)
        {
            EXPECT_NEAR(p_esolver->lj_rcut(i, j), 3.0 * ModuleBase::ANGSTROM_AU, doublethreshold);
        }
    }
    EXPECT_NEAR(p_esolver->search_radius, 3.0 * ModuleBase::ANGSTROM_AU + 0.01, doublethreshold);

    rcut = {3.0, 4.0, 5.0};
    p_esolver->rcut_search_radius(ucell.ntype, rcut);
    EXPECT_NEAR(p_esolver->lj_rcut(0, 0), 3.0 * ModuleBase::ANGSTROM_AU, doublethreshold);
    EXPECT_NEAR(p_esolver->lj_rcut(0, 1), 4.0 * ModuleBase::ANGSTROM_AU, doublethreshold);
    EXPECT_NEAR(p_esolver->lj_rcut(1, 0), 4.0 * ModuleBase::ANGSTROM_AU, doublethreshold);
    EXPECT_NEAR(p_esolver->lj_rcut(1, 1), 5.0 * ModuleBase::ANGSTROM_AU, doublethreshold);
    EXPECT_NEAR(p_esolver->search_radius, 5.0 * ModuleBase::ANGSTROM_AU + 0.01, doublethreshold);
}

TEST_F(LJ_pot_test, SetC6C12)
{
    ModuleESolver::ESolver_LJ* p_esolver = new ModuleESolver::ESolver_LJ();
    ucell.ntype = 2;

    // no rule
    int rule = 1;
    std::vector<double> lj_epsilon = {0.1, 0.2, 0.3};
    std::vector<double> lj_sigma = {0.2, 0.4, 0.6};

    p_esolver->set_c6_c12(ucell.ntype, rule, lj_epsilon, lj_sigma);

    for (int i = 0; i < ucell.ntype; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            int k = i * (i + 1) / 2 + j;
            double temp = pow(lj_sigma[k] * ModuleBase::ANGSTROM_AU, 6);
            EXPECT_NEAR(p_esolver->lj_c6(i, j), 4.0 * lj_epsilon[k] * temp / ModuleBase::Ry_to_eV, doublethreshold);
            EXPECT_NEAR(p_esolver->lj_c12(i, j), p_esolver->lj_c6(i, j) * temp, doublethreshold);
            EXPECT_DOUBLE_EQ(p_esolver->lj_c6(i, j), p_esolver->lj_c6(j, i));
            EXPECT_DOUBLE_EQ(p_esolver->lj_c12(i, j), p_esolver->lj_c12(j, i));
        }
    }

    // rule 1
    rule = 1;
    lj_epsilon = {0.1, 0.2};
    lj_sigma = {0.2, 0.4};

    p_esolver->set_c6_c12(ucell.ntype, rule, lj_epsilon, lj_sigma);

    for (int i = 0; i < ucell.ntype; i++)
    {
        double temp = pow(lj_sigma[i] * ModuleBase::ANGSTROM_AU, 6);
        EXPECT_NEAR(p_esolver->lj_c6(i, i), 4.0 * lj_epsilon[i] * temp / ModuleBase::Ry_to_eV, doublethreshold);
        EXPECT_NEAR(p_esolver->lj_c12(i, i), p_esolver->lj_c6(i, i) * temp, doublethreshold);

        for (int j = 0; j < i; j++)
        {
            EXPECT_NEAR(p_esolver->lj_c6(i, j),
                        std::sqrt(p_esolver->lj_c6(i, i) * p_esolver->lj_c6(j, j)),
                        doublethreshold);
            EXPECT_NEAR(p_esolver->lj_c12(i, j),
                        std::sqrt(p_esolver->lj_c12(i, i) * p_esolver->lj_c12(j, j)),
                        doublethreshold);
            EXPECT_DOUBLE_EQ(p_esolver->lj_c6(i, j), p_esolver->lj_c6(j, i));
            EXPECT_DOUBLE_EQ(p_esolver->lj_c12(i, j), p_esolver->lj_c12(j, i));
        }
    }

    // rule 2
    rule = 2;
    lj_epsilon = {0.1, 0.2};
    lj_sigma = {0.2, 0.4};

    p_esolver->set_c6_c12(ucell.ntype, rule, lj_epsilon, lj_sigma);

    for (int i = 0; i < ucell.ntype; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            double temp = pow((lj_sigma[i] + lj_sigma[j]) / 2 * ModuleBase::ANGSTROM_AU, 6);
            EXPECT_NEAR(p_esolver->lj_c6(i, j),
                        4.0 * std::sqrt(lj_epsilon[i] * lj_epsilon[j]) * temp / ModuleBase::Ry_to_eV,
                        doublethreshold);
            EXPECT_NEAR(p_esolver->lj_c12(i, j), p_esolver->lj_c6(i, j) * temp, doublethreshold);
            EXPECT_DOUBLE_EQ(p_esolver->lj_c6(i, j), p_esolver->lj_c6(j, i));
            EXPECT_DOUBLE_EQ(p_esolver->lj_c12(i, j), p_esolver->lj_c12(j, i));
        }
    }
}

TEST_F(LJ_pot_test, CalEnShift)
{
    ModuleESolver::ESolver_LJ* p_esolver = new ModuleESolver::ESolver_LJ();
    ucell.ntype = 2;

    std::vector<double> rcut = {3.0};
    p_esolver->rcut_search_radius(ucell.ntype, rcut);

    int rule = 1;
    std::vector<double> lj_epsilon = {0.1, 0.2, 0.3};
    std::vector<double> lj_sigma = {0.2, 0.4, 0.6};
    p_esolver->set_c6_c12(ucell.ntype, rule, lj_epsilon, lj_sigma);

    // false
    p_esolver->cal_en_shift(ucell.ntype, false);
    for (int i = 0; i < ucell.ntype; i++)
    {
        for (int j = 0; j < ucell.ntype; j++)
        {
            EXPECT_DOUBLE_EQ(p_esolver->en_shift(i, j), 0.0);
        }
    }

    // true
    p_esolver->cal_en_shift(ucell.ntype, true);
    EXPECT_NEAR(p_esolver->en_shift(0, 0), -2.5810212013100967e-09, doublethreshold);
    EXPECT_NEAR(p_esolver->en_shift(0, 1), -3.303688865319793e-07, doublethreshold);
    EXPECT_NEAR(p_esolver->en_shift(1, 0), -3.303688865319793e-07, doublethreshold);
    EXPECT_NEAR(p_esolver->en_shift(1, 1), -5.6443326024140752e-06, doublethreshold);
}