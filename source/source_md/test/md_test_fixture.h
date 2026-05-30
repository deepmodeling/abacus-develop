#ifndef MD_TEST_FIXTURE_H
#define MD_TEST_FIXTURE_H

#include "gtest/gtest.h"
#include "source_esolver/esolver_lj.h"
#include "source_md/md_base.h"
#include "setcell.h"

#include <memory>
#include <vector>

class MdTestBase : public testing::Test
{
  protected:
    UnitCell ucell;
    Parameter param_in;
    std::unique_ptr<ModuleESolver::ESolver> p_esolver;

    void SetUp() override
    {
        Setcell::setupcell(ucell);
        Setcell::parameters(param_in.input);

        p_esolver.reset(new ModuleESolver::ESolver_LJ());
        p_esolver->before_all_runners(ucell, param_in.inp);
    }
};

template <class Integrator>
class MdIntegratorFixture : public MdTestBase
{
  protected:
    std::unique_ptr<MD_base> mdrun;

    void SetUp() override
    {
        MdTestBase::SetUp();
        mdrun.reset(new Integrator(param_in, ucell));
        mdrun->setup(p_esolver.get(), PARAM.sys.global_readin_dir);
    }
};

class MdFuncTestFixture : public testing::Test
{
  protected:
    UnitCell ucell;
    std::vector<double> allmass_store;
    std::vector<ModuleBase::Vector3<double>> pos_store;
    std::vector<ModuleBase::Vector3<double>> vel_store;
    std::vector<ModuleBase::Vector3<int>> ionmbl_store;
    std::vector<ModuleBase::Vector3<double>> force_store;
    double* allmass = nullptr;
    ModuleBase::Vector3<double>* pos = nullptr;
    ModuleBase::Vector3<double>* vel = nullptr;
    ModuleBase::Vector3<int>* ionmbl = nullptr;
    ModuleBase::Vector3<double>* force = nullptr;
    ModuleBase::matrix virial;
    ModuleBase::matrix stress;
    double potential = 0.0;
    int natom = 0;
    double temperature = 0.0;
    int frozen_freedom = 0;
    Parameter param_in;

    void SetUp() override
    {
        Setcell::setupcell(ucell);
        Setcell::parameters(param_in.input);
        natom = ucell.nat;

        allmass_store.resize(natom);
        pos_store.resize(natom);
        vel_store.resize(natom);
        ionmbl_store.resize(natom);
        force_store.resize(natom);
        allmass = allmass_store.data();
        pos = pos_store.data();
        vel = vel_store.data();
        ionmbl = ionmbl_store.data();
        force = force_store.data();
        stress.create(3, 3);
        virial.create(3, 3);
    }
};

class LjPotTestFixture : public testing::Test
{
  protected:
    std::vector<ModuleBase::Vector3<double>> force_store;
    ModuleBase::Vector3<double>* force = nullptr;
    ModuleBase::matrix stress;
    double potential = 0.0;
    int natom = 0;
    UnitCell ucell;
    Input_para input;

    void SetUp() override
    {
        Setcell::setupcell(ucell);

        natom = ucell.nat;
        force_store.resize(natom);
        force = force_store.data();
        stress.create(3, 3);

        Setcell::parameters(input);
    }
};

#endif // MD_TEST_FIXTURE_H
