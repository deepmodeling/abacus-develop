// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#include "esolver_dfpt_pw.h"

#include "source_estate/module_charge/charge.h"
#include "source_estate/module_pot/pot_xc_fdm.h"
#include "source_base/macros.h"
#include "source_base/tool_quit.h"
#include "source_io/module_parameter/parameter.h"
#include "source_pw/module_dfpt/dfpt_pw.h"
#include "source_pw/module_dfpt/dfpt_rho.h"

#include <complex>
#include <vector>

namespace {

/**
 * @brief Re/Im-split finite-difference adapter over elecstate::PotXC_FDM.
 *
 * The q-shifted complex density amplitude is applied twice (real and
 * imaginary part around the ground-state density); the two real
 * finite-difference responses delta V_xc recombine linearly into the
 * complex first-order kernel (exact up to O(|drho|^2)).
 */
class XC_First_Order_FDM : public ModuleDFPT::XC_First_Order
{
  public:
    XC_First_Order_FDM(ModulePW::PW_Basis* rho_basis,
                       const Charge* chg0,
                       const UnitCell* ucell)
        : ucell_(ucell)
    {
        fdm_ = new elecstate::PotXC_FDM(rho_basis, chg0, ucell);
        chg1_ = new Charge();
        chg1_->set_rhopw(rho_basis);
        chg1_->allocate(chg0->nspin, false);
        veff_1_.create(chg0->nspin, chg0->nrxx);
    }

    ~XC_First_Order_FDM()
    {
        delete fdm_;
        delete chg1_;
    }

    void apply(const std::vector<std::complex<double>>& drho_r,
               std::vector<std::complex<double>>& dvxc_r) const override
    {
        const int nrxx = veff_1_.nc;
        // real part: V_xc[rho0 + Re drho] - V_xc[rho0]
        for (int ir = 0; ir < nrxx; ++ir)
        {
            chg1_->rho[0][ir] = drho_r[ir].real();
        }
        veff_1_.zero_out();
        fdm_->cal_v_eff(chg1_, ucell_, veff_1_);
        for (int ir = 0; ir < nrxx; ++ir)
        {
            dvxc_r[ir] = veff_1_(0, ir);
        }
        // imaginary part: V_xc[rho0 + Im drho] - V_xc[rho0]
        for (int ir = 0; ir < nrxx; ++ir)
        {
            chg1_->rho[0][ir] = drho_r[ir].imag();
        }
        veff_1_.zero_out();
        fdm_->cal_v_eff(chg1_, ucell_, veff_1_);
        for (int ir = 0; ir < nrxx; ++ir)
        {
            dvxc_r[ir] += std::complex<double>(0.0, 1.0) * veff_1_(0, ir);
        }
    }

  private:
    elecstate::PotXC_FDM* fdm_ = nullptr;
    Charge* chg1_ = nullptr;
    mutable ModuleBase::matrix veff_1_;
    const UnitCell* ucell_ = nullptr;
};

} // namespace

namespace ModuleESolver
{

ESolver_DFPT_PW::ESolver_DFPT_PW()
{
    this->classname = "ESolver_DFPT_PW";
    this->basisname = "PW";
    gs_done_ = false;
    dfpt_wired_ = false;
    dfpt_ = nullptr;
    xc_adapter_ = nullptr;
}

ESolver_DFPT_PW::~ESolver_DFPT_PW()
{
    if (dfpt_ != nullptr)
    {
        delete dfpt_;
        dfpt_ = nullptr;
    }
    if (xc_adapter_ != nullptr)
    {
        delete xc_adapter_;
        xc_adapter_ = nullptr;
    }
}

void ESolver_DFPT_PW::before_all_runners(BaseCell& basecell, const Input_para& inp)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    ModuleBase::TITLE("ESolver_DFPT_PW", "before_all_runners");

    ESolver_KS_PW<std::complex<double>, base_device::DEVICE_CPU>::before_all_runners(ucell, inp);

    // capture the (possibly autoset) ground-state scalars once; inp aliases
    // the global input record and read_pseudo/ParamUpdater have run inside
    // the base call
    nspin_ = inp.nspin;
    nelec_ = inp.nelec;
    ecutwfc_ = inp.ecutwfc;
    dft_plus_u_ = inp.dft_plus_u;

    // static DFPT configuration; the ground-state data wiring happens in
    // init_dfpt after the SCF has converged
    dfpt_ = new ModuleDFPT::DFPT_PW();
    dfpt_->set_parameters("dfpt.in");
    dfpt_->set_qmesh(1, 1, 1);
    dfpt_->set_conv_thr(1e-8);
    dfpt_->set_max_iter(100);
}

void ESolver_DFPT_PW::runner(BaseCell& basecell, const int istep)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    ModuleBase::TITLE("ESolver_DFPT_PW", "runner");

    if (!gs_done_)
    {
        run_gs(ucell);
        gs_done_ = true;
    }

    if (dfpt_ != nullptr)
    {
        if (!dfpt_wired_)
        {
            init_dfpt(ucell);
            dfpt_wired_ = true;
        }
        dfpt_->run();
    }

    run_post_process(ucell);
}

void ESolver_DFPT_PW::after_all_runners(BaseCell& basecell)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    ModuleBase::TITLE("ESolver_DFPT_PW", "after_all_runners");

    ESolver_KS_PW<std::complex<double>, base_device::DEVICE_CPU>::after_all_runners(ucell);
}

void ESolver_DFPT_PW::run_gs(UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_DFPT_PW", "run_gs");

    ESolver_KS_PW<std::complex<double>, base_device::DEVICE_CPU>::runner(ucell, 0);
}

void ESolver_DFPT_PW::init_dfpt(UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_DFPT_PW", "init_dfpt");

    if (nspin_ != 1)
    {
        ModuleBase::WARNING_QUIT("ESolver_DFPT_PW::init_dfpt", "DFPT currently supports nspin = 1 only");
    }
    if (this->pelec == nullptr || this->pelec->charge == nullptr || this->pelec->pot == nullptr
        || this->stp.psi_cpu == nullptr || this->pw_rho == nullptr || this->pw_wfc == nullptr)
    {
        ModuleBase::WARNING_QUIT("ESolver_DFPT_PW::init_dfpt", "ground state is not ready");
    }
    if (this->pelec->charge->nrxx != this->pw_rho->nrxx)
    {
        ModuleBase::WARNING_QUIT("ESolver_DFPT_PW::init_dfpt",
                                 "charge is not on the rho grid (DFPT supports NCPP)");
    }

    // converged effective potential on the shared rho grid (row 0: nspin = 1)
    const ModuleBase::matrix& veff_smooth = this->pelec->pot->get_veff_smooth();
    if (veff_smooth.nc != this->pw_rho->nrxx)
    {
        ModuleBase::WARNING_QUIT("ESolver_DFPT_PW::init_dfpt", "veff_smooth is not on the rho grid");
    }
    std::vector<double> veff_r(veff_smooth.nc, 0.0);
    for (int ir = 0; ir < veff_smooth.nc; ++ir)
    {
        veff_r[ir] = veff_smooth(0, ir);
    }

    // first-order XC kernel adapter around the converged ground-state density
    xc_adapter_ = new XC_First_Order_FDM(this->pw_rho, this->pelec->charge, &ucell);

    dfpt_->init(ucell, *this->stp.psi_cpu, this->pw_rho, this->pw_wfc, &this->sf, veff_r,
                this->pelec->wg, this->pelec->ekb, xc_adapter_, nelec_, ecutwfc_,
                dft_plus_u_ ? &this->dftu : nullptr);
}

void ESolver_DFPT_PW::run_post_process(UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_DFPT_PW", "run_post_process");
}

} // namespace ModuleESolver
