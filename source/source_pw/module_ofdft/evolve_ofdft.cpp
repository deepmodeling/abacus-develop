#include "evolve_ofdft.h"

#include "source_io/module_parameter/parameter.h"
#include <iostream>

#include "source_base/parallel_reduce.h"

void Evolve_OFDFT::cal_Hpsi(elecstate::ElecState* pelec, 
                            const Charge& chr, 
                            UnitCell& ucell, 
                            std::vector<std::complex<double>> psi_, 
                            ModulePW::PW_Basis* pw_rho, 
                            std::vector<std::complex<double>> Hpsi)
{
    // update rho
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            chr.rho[is][ir] = abs(psi_[is * pw_rho->nrxx + ir])*abs(psi_[is * pw_rho->nrxx + ir]);
        }
    }

    pelec->pot->update_from_charge(&chr, &ucell); // Hartree + XC + external
    this->cal_tf_potential(chr.rho,pw_rho ,pelec->pot->get_effective_v()); // TF potential

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        const double* vr_eff = pelec->pot->get_effective_v(is);
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            Hpsi[is * pw_rho->nrxx + ir] = vr_eff[ir]*psi_[is * pw_rho->nrxx + ir];
        }
    }
    this->cal_vw_potential_phi(psi_, pw_rho, Hpsi);
}

void Evolve_OFDFT::cal_tf_potential(const double* const* prho, ModulePW::PW_Basis* pw_rho, ModuleBase::matrix& rpot)
{
    if (PARAM.inp.nspin == 1)
    {
#ifdef _OPENMP
#pragma omp parallel for 
#endif
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            rpot(0, ir) += 5.0 / 3.0 * this->c_tf_ * std::pow(prho[0][ir], 2. / 3.);
        }
    }
    else if (PARAM.inp.nspin == 2)
    {
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            for (int ir = 0; ir < pw_rho->nrxx; ++ir)
            {
                rpot(is, ir) += 5.0 / 3.0 * this->c_tf_ * std::pow(2. * prho[is][ir], 2. / 3.);
            }
        }
    }
}

void Evolve_OFDFT::cal_vw_potential_phi(std::vector<std::complex<double>> pphi, 
                                        ModulePW::PW_Basis* pw_rho, 
                                        std::vector<std::complex<double>> Hpsi)
{
    std::complex<double>** rLapPhi = new std::complex<double>*[PARAM.inp.nspin];
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int is = 0; is < PARAM.inp.nspin; ++is) {
        rLapPhi[is] = new std::complex<double>[pw_rho->nrxx];
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            rLapPhi[is][ir]=pphi[is * pw_rho->nrxx + ir];
        }
    }
    std::complex<double>** recipPhi = new std::complex<double>*[PARAM.inp.nspin];

#ifdef _OPENMP
#pragma omp parallel for
#endif  
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        recipPhi[is] = new std::complex<double>[pw_rho->npw];

        pw_rho->real2recip(rLapPhi[is], recipPhi[is]);
        for (int ik = 0; ik < pw_rho->npw; ++ik)
        {
            recipPhi[is][ik] *= pw_rho->gg[ik] * pw_rho->tpiba2;
        }
        pw_rho->recip2real(recipPhi[is], rLapPhi[is]);
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            Hpsi[is * pw_rho->nrxx + ir]+=rLapPhi[is][ir];
        }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        delete[] recipPhi[is];
        delete[] rLapPhi[is];
    }
    delete[] recipPhi;
    delete[] rLapPhi;
}

void Evolve_OFDFT::cal_CD_potential(std::vector<std::complex<double>> psi_, 
                                    ModulePW::PW_Basis* pw_rho, 
                                    ModuleBase::matrix& rpot)
{
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        //recipCurrent = new std::complex<double>[pw_rho->npw];
        //delete[] recipCurrent;
    }
}

void Evolve_OFDFT::propagate_psi(elecstate::ElecState* pelec, 
                                 const Charge& chr, UnitCell& ucell, 
                                 std::vector<std::complex<double>> pphi_, 
                                 ModulePW::PW_Basis* pw_rho)
{
    ModuleBase::timer::tick("ESolver_OF_TDDFT", "propagte_psi");

    std::complex<double> imag(0.0,1.0);
    double dt=PARAM.inp.mdp.md_dt;
    const int nspin = PARAM.inp.nspin;
    const int nrxx = pw_rho->nrxx;
    const int total_size = nspin * nrxx;
    std::vector<std::complex<double>> K1(total_size);
    std::vector<std::complex<double>> K2(total_size);
    std::vector<std::complex<double>> K3(total_size);
    std::vector<std::complex<double>> K4(total_size);
    std::vector<std::complex<double>> psi1(total_size);
    std::vector<std::complex<double>> psi2(total_size);
    std::vector<std::complex<double>> psi3(total_size);

    cal_Hpsi(pelec,chr,ucell,pphi_,pw_rho,K1);
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int is = 0; is < PARAM.inp.nspin; ++is){
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            K1[is * nrxx + ir]=-1.0*K1[is * nrxx + ir]*dt*imag;
            psi1[is * nrxx + ir]=pphi_[is * nrxx + ir]+0.5*K1[is * nrxx + ir];
        }
    }
    cal_Hpsi(pelec,chr,ucell,psi1,pw_rho,K2);
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int is = 0; is < PARAM.inp.nspin; ++is){
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            K2[is * nrxx + ir]=-1.0*K2[is * nrxx + ir]*dt*imag;
            psi2[is * nrxx + ir]=pphi_[is * nrxx + ir]+0.5*K2[is * nrxx + ir];
        }
    }
    cal_Hpsi(pelec,chr,ucell,psi2,pw_rho,K3);
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int is = 0; is < PARAM.inp.nspin; ++is){
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            K3[is * nrxx + ir]=-1.0*K3[is * nrxx + ir]*dt*imag;
            psi3[is * nrxx + ir]=pphi_[is * nrxx + ir]+K3[is * nrxx + ir];
        }
    }
    cal_Hpsi(pelec,chr,ucell,psi3,pw_rho,K4);
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int is = 0; is < PARAM.inp.nspin; ++is){
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            K4[is * nrxx + ir]=-1.0*K4[is * nrxx + ir]*dt*imag;
            pphi_[is * nrxx + ir]+=1.0/6.0*(K1[is * nrxx + ir]+2.0*K2[is * nrxx + ir]+2.0*K3[is * nrxx + ir]+K4[is * nrxx + ir]);
        }
    }

    ModuleBase::timer::tick("ESolver_OF_TDDFT", "propagte_psi");
}