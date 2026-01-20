#include "evolve_ofdft.h"

#include "source_io/module_parameter/parameter.h"
#include <iostream>

#include "source_base/parallel_reduce.h"

void Evolve_OFDFT::cal_Hpsi(elecstate::ElecState* pelec, 
                            Charge& chr, 
                            UnitCell& ucell, 
                            std::vector<std::complex<double>>& psi_, 
                            ModulePW::PW_Basis* pw_rho, 
                            std::vector<std::complex<double>>& Hpsi)
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
    this->renormalize_psi(chr, pw_rho, psi_);

    pelec->pot->update_from_charge(&chr, &ucell); // Hartree + XC + external
    this->cal_tf_potential(chr.rho, pw_rho, pelec->pot->get_eff_v()); // TF potential
    if (PARAM.inp.of_cd)
    {
        this->cal_CD_potential(psi_, pw_rho, pelec->pot->get_eff_v(), PARAM.inp.of_mCD_alpha); // CD potential
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        const double* vr_eff = pelec->pot->get_eff_v(is);
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            Hpsi[is * pw_rho->nrxx + ir] = vr_eff[ir]*psi_[is * pw_rho->nrxx + ir];
        }
    }
    this->cal_vw_potential_phi(psi_, pw_rho, Hpsi);
}

void Evolve_OFDFT::renormalize_psi(Charge& chr, ModulePW::PW_Basis* pw_rho, std::vector<std::complex<double>>& pphi_)
{
    const double sr = chr.sum_rho();
    const double normalize_factor = PARAM.inp.nelec / sr;

    std::cout<<"sr="<<sr<<" nelec="<<PARAM.inp.nelec<<" normalize_factor="<<normalize_factor<<std::endl;
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ir++)
        {
            pphi_[is * pw_rho->nrxx + ir] *= sqrt(normalize_factor);
            chr.rho[is][ir] *= normalize_factor;
        }
    }
    return;
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

void Evolve_OFDFT::cal_vw_potential_phi(std::vector<std::complex<double>>& pphi, 
                                        ModulePW::PW_Basis* pw_rho, 
                                        std::vector<std::complex<double>>& Hpsi)
{
    if (PARAM.inp.nspin <= 0) {
        ModuleBase::WARNING_QUIT("Evolve_OFDFT","nspin must be positive");
    }
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
            Hpsi[is * pw_rho->nrxx + ir] += rLapPhi[is][ir];
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

void Evolve_OFDFT::cal_CD_potential(std::vector<std::complex<double>>& psi_, 
                                    ModulePW::PW_Basis* pw_rho, 
                                    ModuleBase::matrix& rpot,
                                    double mCD_para)
{
    std::complex<double> imag(0.0,1.0);

    if (PARAM.inp.nspin <= 0) {
        ModuleBase::WARNING_QUIT("Evolve_OFDFT","nspin must be positive");
    }
    std::complex<double>** recipPhi = new std::complex<double>*[PARAM.inp.nspin];
    std::complex<double>** rPhi = new std::complex<double>*[PARAM.inp.nspin];
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int is = 0; is < PARAM.inp.nspin; ++is) {
        rPhi[is] = new std::complex<double>[pw_rho->nrxx];
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            rPhi[is][ir]=psi_[is * pw_rho->nrxx + ir];
        }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        std::vector<std::complex<double>> recipCurrent_x(pw_rho->npw);
        std::vector<std::complex<double>> recipCurrent_y(pw_rho->npw);
        std::vector<std::complex<double>> recipCurrent_z(pw_rho->npw);
        std::vector<std::complex<double>> recipCDPotential(pw_rho->npw);
        std::vector<std::complex<double>> rCurrent_x(pw_rho->nrxx);
        std::vector<std::complex<double>> rCurrent_y(pw_rho->nrxx);
        std::vector<std::complex<double>> rCurrent_z(pw_rho->nrxx);
        std::vector<std::complex<double>> kF_r(pw_rho->nrxx);
        std::vector<std::complex<double>> rCDPotential(pw_rho->nrxx);
        recipPhi[is] = new std::complex<double>[pw_rho->npw];

        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            kF_r[ir]=std::pow(3*std::pow(ModuleBase::PI*std::abs(rPhi[is][ir]),2),1.0/3.0);
        }

        pw_rho->real2recip(rPhi[is], recipPhi[is]);
        for (int ik = 0; ik < pw_rho->npw; ++ik)
        {
            recipCurrent_x[ik]=imag*pw_rho->gcar[ik].x*recipPhi[is][ik]* pw_rho->tpiba;
            recipCurrent_y[ik]=imag*pw_rho->gcar[ik].y*recipPhi[is][ik]* pw_rho->tpiba;
            recipCurrent_z[ik]=imag*pw_rho->gcar[ik].z*recipPhi[is][ik]* pw_rho->tpiba;
        }
        pw_rho->recip2real(recipCurrent_x.data(),rCurrent_x.data());
        pw_rho->recip2real(recipCurrent_y.data(),rCurrent_y.data());
        pw_rho->recip2real(recipCurrent_z.data(),rCurrent_z.data());
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            rCurrent_x[ir]=std::imag(rCurrent_x[ir]*std::conj(rPhi[is][ir]));
            rCurrent_y[ir]=std::imag(rCurrent_y[ir]*std::conj(rPhi[is][ir]));
            rCurrent_z[ir]=std::imag(rCurrent_z[ir]*std::conj(rPhi[is][ir]));
        }
        pw_rho->real2recip(rCurrent_x.data(),recipCurrent_x.data());
        pw_rho->real2recip(rCurrent_y.data(),recipCurrent_y.data());
        pw_rho->real2recip(rCurrent_z.data(),recipCurrent_z.data());
        for (int ik = 0; ik < pw_rho->npw; ++ik)
        {
            recipCDPotential[ik]=recipCurrent_x[ik]*pw_rho->gcar[ik].x+recipCurrent_y[ik]*pw_rho->gcar[ik].y+recipCurrent_z[ik]*pw_rho->gcar[ik].z;
            if (pw_rho->gg[ik]==0) 
            {
                recipCDPotential[ik]=0.0;
            }
            else
            {
                recipCDPotential[ik]*=imag/sqrt(pw_rho->gg[ik]);
            }
        }
        pw_rho->recip2real(recipCDPotential.data(),rCDPotential.data());

        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            rpot(0, ir) -= mCD_para*2.0*std::real(rCDPotential[ir])*std::pow(ModuleBase::PI,3) 
                        / (2.0*std::pow(std::real(kF_r[ir]),2));
            if (std::isnan(rpot(0, ir))) 
            {
                rpot(0, ir)=0.0;
            }
        }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        delete[] recipPhi[is];
        delete[] rPhi[is];
    }
    delete[] recipPhi;
    delete[] rPhi;
}

void Evolve_OFDFT::propagate_psi_RK4(elecstate::ElecState* pelec, 
                                 Charge& chr, 
                                 UnitCell& ucell, 
                                 std::vector<std::complex<double>>& pphi_, 
                                 ModulePW::PW_Basis* pw_rho)
{
    ModuleBase::timer::tick("ESolver_OF_TDDFT", "propagate_psi_RK4");

    std::complex<double> imag(0.0,1.0);
    double dt=PARAM.inp.mdp.md_dt / ModuleBase::AU_to_FS;
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
            K1[is * nrxx + ir]=-0.5*K1[is * nrxx + ir]*dt*imag;   // 0.5 convert Ry to Hartree
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
            K2[is * nrxx + ir]=-0.5*K2[is * nrxx + ir]*dt*imag;
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
            K3[is * nrxx + ir]=-0.5*K3[is * nrxx + ir]*dt*imag;
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
            K4[is * nrxx + ir]=-0.5*K4[is * nrxx + ir]*dt*imag;
            pphi_[is * nrxx + ir]+=1.0/6.0*(K1[is * nrxx + ir]+2.0*K2[is * nrxx + ir]+2.0*K3[is * nrxx + ir]+K4[is * nrxx + ir]);
        }
    }

#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            chr.rho[is][ir] = abs(pphi_[is * pw_rho->nrxx + ir])*abs(pphi_[is * pw_rho->nrxx + ir]);
        }
    }
    this->renormalize_psi(chr, pw_rho, pphi_);

    ModuleBase::timer::tick("ESolver_OF_TDDFT", "propagate_psi_RK4");
}

void Evolve_OFDFT::propagate_psi_RK2(elecstate::ElecState* pelec, 
                                 Charge& chr, 
                                 UnitCell& ucell, 
                                 std::vector<std::complex<double>>& pphi_, 
                                 ModulePW::PW_Basis* pw_rho)
{
    ModuleBase::timer::tick("ESolver_OF_TDDFT", "propagate_psi_RK2");

    const std::complex<double> imag(0.0, 1.0);
    double dt=PARAM.inp.mdp.md_dt / ModuleBase::AU_to_FS;
    const int nspin = PARAM.inp.nspin;
    const int nrxx = pw_rho->nrxx;
    const int total_size = nspin * nrxx;

    std::vector<std::complex<double>> K1(total_size);
    std::vector<std::complex<double>> K2(total_size);
    std::vector<std::complex<double>> psi_mid(total_size);

    cal_Hpsi(pelec, chr, ucell, pphi_, pw_rho, K1);

#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int is = 0; is < nspin; ++is) {
        for (int ir = 0; ir < nrxx; ++ir) {
            const int idx = is * nrxx + ir;
            K1[idx] = -0.5 * K1[idx] * dt * imag;
            psi_mid[idx] = pphi_[idx] + 0.5 * K1[idx];
        }
    }

    cal_Hpsi(pelec, chr, ucell, psi_mid, pw_rho, K2);

#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int is = 0; is < nspin; ++is) {
        for (int ir = 0; ir < nrxx; ++ir) {
            const int idx = is * nrxx + ir;
            K2[idx] = -0.5 * K2[idx] * dt * imag;
            pphi_[idx] += K2[idx];
        }
    }

#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int is = 0; is < nspin; ++is) {
        for (int ir = 0; ir < nrxx; ++ir) {
            chr.rho[is][ir] = std::norm(pphi_[is * nrxx + ir]); 
        }
    }

    this->renormalize_psi(chr, pw_rho, pphi_);

    ModuleBase::timer::tick("ESolver_OF_TDDFT", "propagate_psi_RK2");
}
