#include "gtest/gtest.h"

#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "source_estate/elecstate.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_charge/chg_mix.h"
#include "source_estate/module_charge/chg_routine.h"
#include "source_pw/module_pwdft/dftu_base.h"

#include <complex>
#include <vector>

// charge.cpp references Magnetism; provide a lightweight stub.
Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}
Magnetism::~Magnetism()
{
}

/************************************************
 *  unit test of module_charge/chg_routine.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - chgmixing_ks_pw: iter==1 init path sets mixing_restart_step.
 *   - chgmixing_ks_lcao: iter==1 mix_reset path sets mixing_restart_step.
 *   - chgmixing_ks: convergence branches (conv_esolver true skips mixing;
 *     drho < hsolver_error skips mixing).
 */

namespace
{

MixingConfig make_plain_cfg(int nspin)
{
    return {"plain", 0.5, 8, 0.0, false, 0.5, 0.0, 0.0, -1.0,
            false, nspin, 2, false, false, false, false, 100};
}

} // namespace

class ChgRoutineTest : public ::testing::Test
{
  protected:
    ModulePW::PW_Basis pw_basis;
    Charge charge;
    UnitCell ucell;

    void SetUp() override
    {
        pw_basis.initgrids(4, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 20);
        pw_basis.initparameters(false, 20);
        pw_basis.setuptransform();
        pw_basis.collect_local_pw();
        charge.set_rhopw(&pw_basis);
        charge.allocate(1, false, false, 0);
        ucell.omega = pw_basis.omega;
        ucell.tpiba = 1.0;
    }
};

TEST_F(ChgRoutineTest, ChgmixingKsPwIter1SetsRestartStep)
{
    Charge_Mixing cm;
    cm.set_mixing(make_plain_cfg(1), &pw_basis, &pw_basis, ucell.omega, ucell.tpiba);
    Plus_U_Base dftu;
    Input_para inp;
    inp.scf_nmax = 50;
    inp.mixing_restart = 0.0;
    inp.dft_plus_u = false;

    module_charge::chgmixing_ks_pw(1, &cm, dftu, true, inp);

    EXPECT_EQ(cm.mixing_restart_step, inp.scf_nmax + 1);
}

TEST_F(ChgRoutineTest, ChgmixingKsLcaoIter1SetsRestartStep)
{
    Charge_Mixing cm;
    cm.set_mixing(make_plain_cfg(1), &pw_basis, &pw_basis, ucell.omega, ucell.tpiba);
    Plus_U_Base dftu;
    Input_para inp;
    inp.scf_nmax = 50;
    inp.mixing_restart = 0.0;
    inp.dft_plus_u = false;

    module_charge::chgmixing_ks_lcao(1, &cm, dftu, 0, inp);

    EXPECT_EQ(cm.mixing_restart_step, inp.scf_nmax + 1);
}

TEST_F(ChgRoutineTest, ChgmixingKsConvergedSkipsMixing)
{
    Charge_Mixing cm;
    cm.set_mixing(make_plain_cfg(1), &pw_basis, &pw_basis, ucell.omega, ucell.tpiba);
    Input_para inp;
    inp.mixing_restart = 0.0;
    inp.scf_os_stop = false;
    inp.scf_thr_type = 2;
    inp.calculation = "scf";
    inp.nelec = 1.0;

    module_charge::ScfMixingCtx ctx;
    ctx.hsolver_error = 1e-6;
    ctx.scf_thr = 1e-4;
    ctx.scf_ene_thr = 0.0;
    ctx.converged_u = true;
    ctx.ks_run = true;
    ctx.drho = 1e-6; // < scf_thr => converged

    std::vector<double> rho_before(charge.rho[0], charge.rho[0] + pw_basis.nrxx);

    module_charge::chgmixing_ks(2, ucell, nullptr, charge, pw_basis, &cm, ctx, inp);

    EXPECT_TRUE(ctx.conv_esolver);
    // rho must be unchanged because conv_esolver is true
    for (int ir = 0; ir < pw_basis.nrxx; ++ir)
    {
        EXPECT_EQ(charge.rho[0][ir], rho_before[ir]);
    }
}

TEST_F(ChgRoutineTest, ChgmixingKsDrhoBelowHsolverSkipsMixing)
{
    Charge_Mixing cm;
    cm.set_mixing(make_plain_cfg(1), &pw_basis, &pw_basis, ucell.omega, ucell.tpiba);
    Input_para inp;
    inp.mixing_restart = 0.0;
    inp.scf_os_stop = false;
    inp.scf_thr_type = 2;
    inp.calculation = "scf";
    inp.nelec = 1.0;

    module_charge::ScfMixingCtx ctx;
    ctx.hsolver_error = 1e-3;
    ctx.scf_thr = 1e-4;
    ctx.scf_ene_thr = 0.0;
    ctx.converged_u = true;
    ctx.ks_run = true;
    ctx.drho = 1e-5; // < hsolver_error

    std::vector<double> rho_before(charge.rho[0], charge.rho[0] + pw_basis.nrxx);

    module_charge::chgmixing_ks(2, ucell, nullptr, charge, pw_basis, &cm, ctx, inp);

    // rho unchanged because drho < hsolver_error
    for (int ir = 0; ir < pw_basis.nrxx; ++ir)
    {
        EXPECT_EQ(charge.rho[0][ir], rho_before[ir]);
    }
}
