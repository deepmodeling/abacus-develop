#include "chg_symm.h"

#include "chg_symm_detail.h"
#include "source_estate/module_charge/charge.h"

namespace module_charge
{

// TODO: callers currently pass different grids for LCAO-based paths:
// esolver_ks_lcao (and TDDFT/SDFT LCAO) pass the smooth pw_rho, while
// esolver_ks_lcaopw (LIP) and get_pchg_pw pass the dense pw_rhod. The two
// coincide only because LCAO rejects USPP, so double_grid is always false there
// (see uspp_support.cpp). If LCAO is extended to USPP, verify which grid the
// charge symmetrization must use before relaxing that restriction.
void symmetrize_rho(const int nspin,
                    const Charge& chr,
                    const ModulePW::PW_Basis* pw,
                    ModuleSymmetry::Symmetry& symm)
{
    if (nspin == 4)
    {
        // nspin=4 (non-collinear/SOC): rho[0] is the charge density rho^0 (scalar, symmetrized
        // spatially like nspin=1); rho[1,2,3] are the spin density (rho^x, rho^y, rho^z) which
        // must be symmetrized TOGETHER with the per-operation spin rotation W(g).
        cal_rhog_symm(0, chr, pw, symm);
        cal_rhog_symm_soc(chr.rho, chr.rhog, pw, symm);
        return;
    }
    for (int is = 0; is < nspin; is++)
    {
        cal_rhog_symm(is, chr, pw, symm);
    }
}

void cal_rhog_symm(const int& spin_now,
                   const Charge& chr,
                   const ModulePW::PW_Basis* rho_basis,
                   ModuleSymmetry::Symmetry& symm)
{
    assert(spin_now < 4); // added by zhengdy-soc

    if (ModuleSymmetry::Symmetry::symm_flag != 1)
    {
        return;
    }

    ModuleBase::TITLE("module_charge", "cal_rhog_symm");
    ModuleBase::timer::start("module_charge", "cal_rhog_symm");

    rho_basis->real2recip(chr.rho[spin_now], chr.rhog[spin_now]);

    detail::psymmg(chr.rhog[spin_now], rho_basis, symm);

    rho_basis->recip2real(chr.rhog[spin_now], chr.rho[spin_now]);

    if (chr.kin_r != nullptr)
    {
        // Use std::vector to manage kin_g instead of raw pointer
        std::vector<std::complex<double>> kin_g(rho_basis->npw);
        rho_basis->real2recip(chr.kin_r[spin_now], kin_g.data());
        detail::psymmg(kin_g.data(), rho_basis, symm);
        rho_basis->recip2real(kin_g.data(), chr.kin_r[spin_now]);
    }

    ModuleBase::timer::end("module_charge", "cal_rhog_symm");
    return;
}

void cal_rhog_symm(const int& spin_now,
                   double** rho,
                   std::complex<double>** rhog,
                   int ngmc,
                   double** kin_r,
                   const ModulePW::PW_Basis* rho_basis,
                   ModuleSymmetry::Symmetry& symm)
{
    assert(spin_now < 4); // added by zhengdy-soc

    if (ModuleSymmetry::Symmetry::symm_flag != 1)
    {
        return;
    }

    ModuleBase::TITLE("module_charge", "cal_rhog_symm");
    ModuleBase::timer::start("module_charge", "cal_rhog_symm");

    {
        rho_basis->real2recip(rho[spin_now], rhog[spin_now]);
        detail::psymmg(rhog[spin_now], rho_basis, symm);
        rho_basis->recip2real(rhog[spin_now], rho[spin_now]);

        if (kin_r != nullptr)
        {
            std::vector<std::complex<double>> kin_g(ngmc);
            rho_basis->real2recip(kin_r[spin_now], kin_g.data());
            detail::psymmg(kin_g.data(), rho_basis, symm);
            rho_basis->recip2real(kin_g.data(), kin_r[spin_now]);
        }
    }

    ModuleBase::timer::end("module_charge", "cal_rhog_symm");
    return;
}

void cal_rhog_symm_soc(double** rho,
                       std::complex<double>** rhog,
                       const ModulePW::PW_Basis* rho_basis,
                       ModuleSymmetry::Symmetry& symm)
{
    if (ModuleSymmetry::Symmetry::symm_flag != 1)
    {
        return;
    }

    ModuleBase::TITLE("module_charge", "cal_rhog_symm_soc");
    ModuleBase::timer::start("module_charge", "cal_rhog_symm_soc");

    // the three spin components are coupled by the spin rotation, so they are transformed to
    // reciprocal space and symmetrized together (rho[1]=rho^x, rho[2]=rho^y, rho[3]=rho^z).
    for (int is = 1; is < 4; ++is)
    {
        rho_basis->real2recip(rho[is], rhog[is]);
    }

    detail::psymmg_soc(rhog[1], rhog[2], rhog[3], rho_basis, symm);

    for (int is = 1; is < 4; ++is)
    {
        rho_basis->recip2real(rhog[is], rho[is]);
    }

    ModuleBase::timer::end("module_charge", "cal_rhog_symm_soc");
    return;
}

} // namespace module_charge
