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
#include "source_base/math_polyint.h"
#include "source_base/tool_quit.h"
#include "source_io/module_parameter/parameter.h"
#include "source_pw/module_dfpt/dfpt_pert.h"
#include "source_pw/module_dfpt/dfpt_pw.h"
#include "source_pw/module_dfpt/dfpt_rho.h"

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
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
        if (static_cast<int>(drho_r.size()) != nrxx)
        {
            ModuleBase::WARNING_QUIT("XC_First_Order_FDM", "drho_r is not on the rho grid");
        }
        if (dvxc_r.size() != drho_r.size())
        {
            dvxc_r.resize(drho_r.size());
        }
        // central difference with a small probe amplitude: a forward
        // difference Vxc[rho0 + drho] - Vxc[rho0] carries the curvature
        // term ~ Vxc'' * drho^2 / 2, which is quadratic in the T2 response
        // and leaks a spurious A1 component into the screened potential
        const double eta = 1.0e-6;
        std::vector<double> v_plus(nrxx);
        // real part: (Vxc[rho0 + eta Re drho] - Vxc[rho0 - eta Re drho]) / 2
        for (int ir = 0; ir < nrxx; ++ir)
        {
            chg1_->rho[0][ir] = eta * drho_r[ir].real();
        }
        veff_1_.zero_out();
        fdm_->cal_v_eff(chg1_, ucell_, veff_1_);
        for (int ir = 0; ir < nrxx; ++ir)
        {
            v_plus[ir] = veff_1_(0, ir);
        }
        for (int ir = 0; ir < nrxx; ++ir)
        {
            chg1_->rho[0][ir] = -eta * drho_r[ir].real();
        }
        veff_1_.zero_out();
        fdm_->cal_v_eff(chg1_, ucell_, veff_1_);
        for (int ir = 0; ir < nrxx; ++ir)
        {
            dvxc_r[ir] = (v_plus[ir] - veff_1_(0, ir)) / (2.0 * eta);
        }
        if (getenv("DFPT_XCDBG") != nullptr)
        {
            static int call = 0;
            if (call < 4)
            {
                double dp_dv = 0.0;
                double dp_dp = 0.0;
                double dv_dv = 0.0;
                for (int ir = 0; ir < nrxx; ++ir)
                {
                    dp_dv += drho_r[ir].real() * dvxc_r[ir].real();
                    dp_dp += drho_r[ir].real() * drho_r[ir].real();
                    dv_dv += dvxc_r[ir].real() * dvxc_r[ir].real();
                }
                int imax = 0;
                double dmax = 0.0;
                for (int ir = 0; ir < nrxx; ++ir)
                {
                    if (std::abs(drho_r[ir].real()) > dmax)
                    {
                        dmax = std::abs(drho_r[ir].real());
                        imax = ir;
                    }
                }
                std::cout << "XCDBG call=" << call
                          << " <drho|dvxc>/|drho|^2=" << (dp_dv / dp_dp)
                          << " |dvxc|/|drho|=" << (std::sqrt(dv_dv) / std::sqrt(dp_dp))
                          << " ratio@max=" << (dvxc_r[imax].real() / drho_r[imax].real())
                          << " drho@max=" << drho_r[imax].real()
                          << std::endl;
            }
            ++call;
        }
        // imaginary part: same central difference on Im drho
        for (int ir = 0; ir < nrxx; ++ir)
        {
            chg1_->rho[0][ir] = eta * drho_r[ir].imag();
        }
        veff_1_.zero_out();
        fdm_->cal_v_eff(chg1_, ucell_, veff_1_);
        for (int ir = 0; ir < nrxx; ++ir)
        {
            v_plus[ir] = veff_1_(0, ir);
        }
        for (int ir = 0; ir < nrxx; ++ir)
        {
            chg1_->rho[0][ir] = -eta * drho_r[ir].imag();
        }
        veff_1_.zero_out();
        fdm_->cal_v_eff(chg1_, ucell_, veff_1_);
        for (int ir = 0; ir < nrxx; ++ir)
        {
            dvxc_r[ir] += std::complex<double>(0.0, 1.0)
                          * (v_plus[ir] - veff_1_(0, ir)) / (2.0 * eta);
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

    // static DFPT configuration from INPUT (explicit passing, rule 1); the
    // ground-state data wiring happens in init_dfpt after the SCF converges
    dfpt_ = new ModuleDFPT::DFPT_PW();
    dfpt_->set_qmesh(inp.dfpt_qmesh[0], inp.dfpt_qmesh[1], inp.dfpt_qmesh[2]);
    dfpt_->set_qfile(inp.dfpt_qfile);
    dfpt_->set_conv_thr(inp.dfpt_conv_thr);
    dfpt_->set_max_iter(inp.dfpt_max_iter);
    dfpt_->set_mix_beta(inp.dfpt_mix_beta);
    dfpt_->set_compute_q0(inp.dfpt_compute_q0);
    dfpt_->set_loto(inp.dfpt_loto);
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

    if (getenv("DFPT_VKB") != nullptr)
    {
        // design-phase validation: elementwise comparison of the module's
        // NC projectors against the ground-state ppcell vkb at the last k
        const int ik = 0;
        const int npw = this->pw_wfc->npwk[ik];
        std::vector<ModuleBase::Vector3<double>> gk(npw);
        for (int ig = 0; ig < npw; ++ig)
        {
            gk[ig] = this->pw_wfc->getgpluskcar(ik, ig);
        }
        ModuleDFPT::DFPT_Pert pert;
        pert.init(ucell, this->pw_rho, this->pw_wfc, this->sf);
        std::vector<std::vector<std::complex<double>>> vkb;
        pert.build_vkb(0, 0, gk, vkb);
        std::vector<std::vector<std::complex<double>>> vkb1;
        pert.build_vkb(0, 1, gk, vkb1);
        const int nh = ucell.atoms[0].ncpp.nh;
        const std::complex<double>* gsvkb = this->ppcell.get_vkb_data<double>();
        std::cout << "VKBCHK npw=" << npw << " nh=" << nh
                  << " nkb=" << this->ppcell.nkb << std::endl;
        for (int mu = 0; mu < nh; ++mu)
        {
            std::complex<double> dot(0.0, 0.0);
            double nrm_mine = 0.0;
            double nrm_gs = 0.0;
            for (int ig = 0; ig < npw; ++ig)
            {
                const std::complex<double> gs = gsvkb[mu * this->pw_wfc->npwk_max + ig];
                dot += std::conj(vkb[mu][ig]) * gs;
                nrm_mine += std::norm(vkb[mu][ig]);
                nrm_gs += std::norm(gs);
            }
            std::cout << "VKBCHK mu=" << mu << " <mine|gs>=" << dot
                      << " |mine|^2=" << nrm_mine << " |gs|^2=" << nrm_gs << std::endl;
        }
        for (int mu = 0; mu < nh; ++mu)
        {
            std::complex<double> dot(0.0, 0.0);
            double nrm_mine = 0.0;
            double nrm_gs = 0.0;
            for (int ig = 0; ig < npw; ++ig)
            {
                const std::complex<double> gs = gsvkb[(nh + mu) * this->pw_wfc->npwk_max + ig];
                dot += std::conj(vkb1[mu][ig]) * gs;
                nrm_mine += std::norm(vkb1[mu][ig]);
                nrm_gs += std::norm(gs);
            }
            std::cout << "VKBCHK a1 mu=" << mu << " <mine|gs>=" << dot
                      << " |mine|^2=" << nrm_mine << " |gs|^2=" << nrm_gs << std::endl;
        }
        for (int ig = 0; ig < 6; ++ig)
        {
            std::cout << "VKBEL a1 mu=0 ig=" << ig << " mine=" << vkb1[0][ig]
                      << " gs=" << gsvkb[4 * this->pw_wfc->npwk_max + ig]
                      << " gcar=" << gk[ig].x << "," << gk[ig].y << "," << gk[ig].z << std::endl;
        }
    }

    dfpt_->init(ucell, *this->stp.psi_cpu, this->pw_rho, this->pw_wfc, &this->sf, veff_r,
                this->pelec->wg, this->pelec->ekb, xc_adapter_, nelec_, ecutwfc_,
                dft_plus_u_ ? &this->dftu : nullptr);
}

void ESolver_DFPT_PW::run_post_process(UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_DFPT_PW", "run_post_process");

    if (dfpt_ == nullptr)
    {
        return;
    }
    // multi-q frequency report (one block per q of the list, plus the LO-TO
    // corrected Gamma block along the data-layer direction when enabled);
    // the tensor blocks below stay design-phase std::cout until the io
    // layer integration of the data-layer consolidation stage
    const int nq = dfpt_->get_nq();
    for (int q_idx = 0; q_idx < nq; ++q_idx)
    {
        std::cout << dfpt_->format_q_report(q_idx);
        if (q_idx == 0)
        {
            std::cout << dfpt_->format_loto_report();
        }
    }
    const ModuleBase::matrix& eps = dfpt_->get_dielectric_tensor();
    if (eps.nr == 3 && eps.nc == 3)
    {
        std::cout << " DFPT dielectric tensor (epsilon_inf):" << std::endl;
        for (int a = 0; a < eps.nr; ++a)
        {
            std::cout << "   ";
            for (int b = 0; b < eps.nc; ++b)
            {
                std::cout << eps(a, b) << " ";
            }
            std::cout << std::endl;
        }
    }
    for (int iat = 0; iat < ucell.nat; ++iat)
    {
        const ModuleBase::matrix& zstar = dfpt_->get_born_charges(iat);
        if (zstar.nr != 3 || zstar.nc != 3)
        {
            continue;
        }
        std::cout << " DFPT Born effective charge atom " << iat << ":" << std::endl;
        for (int a = 0; a < zstar.nr; ++a)
        {
            std::cout << "   ";
            for (int b = 0; b < zstar.nc; ++b)
            {
                std::cout << zstar(a, b) << " ";
            }
            std::cout << std::endl;
        }
    }
}

} // namespace ModuleESolver
