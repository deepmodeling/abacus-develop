#include "evolve_phi.h"

#include "source_io/module_parameter/parameter.h"
#include <iostream>

#include "source_base/parallel_reduce.h"

void EVOLVE_PHI::get_Hpsi(elecstate::ElecState* pelec, const Charge& chr, UnitCell& ucell, const std::complex<double>* const* psi_, ModulePW::PW_Basis* pw_rho, std::complex<double>** Hpsi)
{
    // update rho
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            chr.rho[is][ir] = abs(psi_[is][ir])*abs(psi_[is][ir]);
        }
    }

    pelec->pot->update_from_charge(&chr, &ucell); // Hartree + XC + external
    this->get_tf_potential(chr.rho,pw_rho ,pelec->pot->get_effective_v()); // TF potential
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        const double* vr_eff = pelec->pot->get_effective_v(is);
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            Hpsi[is][ir] = vr_eff[ir]*psi_[is][ir];
        }
    }
    this->get_vw_potential_phi(psi_, pw_rho, Hpsi);
}

void EVOLVE_PHI::get_tf_potential(const double* const* prho, ModulePW::PW_Basis* pw_rho, ModuleBase::matrix& rpot)
{
    if (PARAM.inp.nspin == 1)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            rpot(0, ir) += 5.0 / 3.0 * this->c_tf_ * std::pow(prho[0][ir], 2. / 3.);
        }
    }
    else if (PARAM.inp.nspin == 2)
    {
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            for (int ir = 0; ir < pw_rho->nrxx; ++ir)
            {
                rpot(is, ir) += 5.0 / 3.0 * this->c_tf_ * std::pow(2. * prho[is][ir], 2. / 3.);
            }
        }
    }
}

void EVOLVE_PHI::get_vw_potential_phi(const std::complex<double>* const* pphi, ModulePW::PW_Basis* pw_rho, std::complex<double>** Hpsi)
{
    std::complex<double>** rLapPhi = new std::complex<double>*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is) {
        rLapPhi[is] = new std::complex<double>[pw_rho->nrxx];
    }
    std::complex<double>** recipPhi = new std::complex<double>*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        recipPhi[is] = new std::complex<double>[pw_rho->npw];

        pw_rho->real2recip(pphi[is], recipPhi[is]);
        for (int ik = 0; ik < pw_rho->npw; ++ik)
        {
            recipPhi[is][ik] *= pw_rho->gg[ik] * pw_rho->tpiba2;
        }
        pw_rho->recip2real(recipPhi[is], rLapPhi[is]);
        for (int ik = 0; ik < pw_rho->npw; ++ik)
        {
            Hpsi[is][ik]+=rLapPhi[is][ik];
        }
    }

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        delete[] recipPhi[is];
        delete[] rLapPhi[is];
    }
    delete[] recipPhi;
    delete[] rLapPhi;
}

void EVOLVE_PHI::get_CD_potential(const std::complex<double>* const* psi_, ModulePW::PW_Basis* pw_rho, ModuleBase::matrix& rpot)
{
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        //recipCurrent = new std::complex<double>[pw_rho->npw];
        //delete[] recipCurrent;
    }
}

void EVOLVE_PHI::propagate_psi(elecstate::ElecState* pelec, const Charge& chr, UnitCell& ucell, std::complex<double>** pphi_, ModulePW::PW_Basis* pw_rho)
{
    ModuleBase::timer::tick("ESolver_OF_TDDFT", "propagte_psi");

    std::complex<double> imag(0.0,1.0);
    double dt=PARAM.inp.mdp.md_dt;
    std::complex<double>** K1 = new std::complex<double>*[PARAM.inp.nspin];
    std::complex<double>** K2 = new std::complex<double>*[PARAM.inp.nspin];
    std::complex<double>** K3 = new std::complex<double>*[PARAM.inp.nspin];
    std::complex<double>** K4 = new std::complex<double>*[PARAM.inp.nspin];
    std::complex<double>** psi1 = new std::complex<double>*[PARAM.inp.nspin];
    std::complex<double>** psi2 = new std::complex<double>*[PARAM.inp.nspin];
    std::complex<double>** psi3 = new std::complex<double>*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is) {
        K1[is] = new std::complex<double>[pw_rho->nrxx];
        K2[is] = new std::complex<double>[pw_rho->nrxx];
        K3[is] = new std::complex<double>[pw_rho->nrxx];
        K4[is] = new std::complex<double>[pw_rho->nrxx];
        psi1[is] = new std::complex<double>[pw_rho->nrxx];
        psi2[is] = new std::complex<double>[pw_rho->nrxx];
        psi3[is] = new std::complex<double>[pw_rho->nrxx];
    }
    get_Hpsi(pelec,chr,ucell,pphi_,pw_rho,K1);
    for (int is = 0; is < PARAM.inp.nspin; ++is){
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            K1[is][ir]=-1.0*K1[is][ir]*dt*imag;
            psi1[is][ir]=pphi_[is][ir]+0.5*K1[is][ir];
        }
    }
    get_Hpsi(pelec,chr,ucell,psi1,pw_rho,K2);
    for (int is = 0; is < PARAM.inp.nspin; ++is){
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            K2[is][ir]=-1.0*K2[is][ir]*dt*imag;
            psi2[is][ir]=pphi_[is][ir]+0.5*K2[is][ir];
        }
    }
    get_Hpsi(pelec,chr,ucell,psi2,pw_rho,K3);
    for (int is = 0; is < PARAM.inp.nspin; ++is){
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            K3[is][ir]=-1.0*K3[is][ir]*dt*imag;
            psi3[is][ir]=pphi_[is][ir]+K3[is][ir];
        }
    }
    get_Hpsi(pelec,chr,ucell,psi3,pw_rho,K4);
    for (int is = 0; is < PARAM.inp.nspin; ++is){
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            K4[is][ir]=-1.0*K4[is][ir]*dt*imag;
            pphi_[is][ir]+=1.0/6.0*(K1[is][ir]+2.0*K2[is][ir]+2.0*K3[is][ir]+K4[is][ir]);
        }
    }

    for (int is = 0; is < PARAM.inp.nspin; ++is) {
        delete[] K1[is];
        delete[] K2[is];
        delete[] K3[is];
        delete[] K4[is];
        delete[] psi1[is];
        delete[] psi2[is];
        delete[] psi3[is];
    }
    delete[] K1;
    delete[] K2;
    delete[] K3;
    delete[] K4;
    delete[] psi1;
    delete[] psi2;
    delete[] psi3;
    ModuleBase::timer::tick("ESolver_OF_TDDFT", "propagte_psi");
}