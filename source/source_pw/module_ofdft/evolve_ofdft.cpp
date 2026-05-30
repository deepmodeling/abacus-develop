#include "evolve_ofdft.h"

#include "source_io/module_parameter/parameter.h"
#include <complex>

#include "source_base/parallel_reduce.h"

// Data race fix: Cache shared member variables to local const variables before OpenMP parallel regions.
// ThreadSanitizer detected race conditions when accessing member variables through pointers
// in parallel loops, even though these values are logically read-only. By caching them
// to local const variables, the compiler can guarantee they won't change within the parallel region,
// eliminating false positive data race warnings from ThreadSanitizer.

void Evolve_OFDFT::cal_Hpsi(elecstate::ElecState* pelec, 
                            Charge& chr, 
                            UnitCell& ucell, 
                            std::vector<std::complex<double>>& psi_, 
                            ModulePW::PW_Basis* pw_rho, 
                            std::vector<std::complex<double>>& Hpsi)
{
    const int nspin = PARAM.inp.nspin;
    const int nrxx = pw_rho->nrxx;

    // update rho
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < nrxx; ++ir)
        {
            chr.rho[is][ir] = std::norm(psi_[is * nrxx + ir]);
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
    for (int is = 0; is < nspin; ++is)
    {
        const double* vr_eff = pelec->pot->get_eff_v(is);
        for (int ir = 0; ir < nrxx; ++ir)
        {
            Hpsi[is * nrxx + ir] = vr_eff[ir]*psi_[is * nrxx + ir];
        }
    }
    this->cal_vw_potential_phi(psi_, pw_rho, Hpsi);
}

void Evolve_OFDFT::renormalize_psi(Charge& chr, ModulePW::PW_Basis* pw_rho, std::vector<std::complex<double>>& pphi_)
{
    const double sr = chr.sum_rho();
    const double normalize_factor = PARAM.inp.nelec / sr;
    const int nspin = PARAM.inp.nspin;
    const int nrxx = pw_rho->nrxx;


#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int is = 0; is < nspin; is++)
    {
        for (int ir = 0; ir < nrxx; ir++)
        {
            pphi_[is * nrxx + ir] *= sqrt(normalize_factor);
            chr.rho[is][ir] *= normalize_factor;
        }
    }
    return;
}

void Evolve_OFDFT::cal_tf_potential(const double* const* prho, ModulePW::PW_Basis* pw_rho, ModuleBase::matrix& rpot)
{
    const int nspin = PARAM.inp.nspin;
    const int nrxx = pw_rho->nrxx;
    const double c_tf_local = this->c_tf_;

    if (nspin == 1)
    {
#ifdef _OPENMP
#pragma omp parallel for 
#endif
        for (int ir = 0; ir < nrxx; ++ir)
        {
            rpot(0, ir) += 5.0 / 3.0 * c_tf_local * std::pow(prho[0][ir], 2. / 3.);
        }
    }
    else if (nspin == 2)
    {
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
        for (int is = 0; is < nspin; ++is)
        {
            for (int ir = 0; ir < nrxx; ++ir)
            {
                rpot(is, ir) += 5.0 / 3.0 * c_tf_local * std::pow(2. * prho[is][ir], 2. / 3.);
            }
        }
    }
}

void Evolve_OFDFT::cal_vw_potential_phi(std::vector<std::complex<double>>& pphi, 
                                        ModulePW::PW_Basis* pw_rho, 
                                        std::vector<std::complex<double>>& Hpsi)
{
    const int nspin = PARAM.inp.nspin;
    const int nrxx = pw_rho->nrxx;
    const int npw = pw_rho->npw;
    const double* gg = pw_rho->gg;
    const double tpiba2 = pw_rho->tpiba2;

    std::vector<std::vector<std::complex<double>>> rLapPhi(nspin, std::vector<std::complex<double>>(nrxx));
    std::vector<std::vector<std::complex<double>>> recipPhi(nspin, std::vector<std::complex<double>>(npw));

    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < nrxx; ++ir)
        {
            rLapPhi[is][ir] = pphi[is * nrxx + ir];
        }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int is = 0; is < nspin; ++is)
    {
        pw_rho->real2recip(rLapPhi[is].data(), recipPhi[is].data());
        for (int ik = 0; ik < npw; ++ik)
        {
            recipPhi[is][ik] *= gg[ik] * tpiba2;
        }
        pw_rho->recip2real(recipPhi[is].data(), rLapPhi[is].data());
        for (int ir = 0; ir < nrxx; ++ir)
        {
            Hpsi[is * nrxx + ir] += rLapPhi[is][ir];
        }
    }
}

void Evolve_OFDFT::cal_CD_potential(std::vector<std::complex<double>>& psi_, 
                                    ModulePW::PW_Basis* pw_rho, 
                                    ModuleBase::matrix& rpot,
                                    double mCD_para)
{
    const std::complex<double> imag(0.0,1.0);
    const int nspin = PARAM.inp.nspin;
    const int nrxx = pw_rho->nrxx;
    const int npw = pw_rho->npw;
    const double tpiba = pw_rho->tpiba;
    const double* gg = pw_rho->gg;
    const ModuleBase::Vector3<double>* gcar = pw_rho->gcar;

    std::vector<std::vector<std::complex<double>>> recipPhi(nspin, std::vector<std::complex<double>>(npw));
    std::vector<std::vector<std::complex<double>>> rPhi(nspin, std::vector<std::complex<double>>(nrxx));

    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < nrxx; ++ir)
        {
            rPhi[is][ir] = psi_[is * nrxx + ir];
        }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int is = 0; is < nspin; ++is)
    {
        std::vector<std::complex<double>> recipCurrent_x(npw);
        std::vector<std::complex<double>> recipCurrent_y(npw);
        std::vector<std::complex<double>> recipCurrent_z(npw);
        std::vector<std::complex<double>> recipCDPotential(npw);
        std::vector<std::complex<double>> rCurrent_x(nrxx);
        std::vector<std::complex<double>> rCurrent_y(nrxx);
        std::vector<std::complex<double>> rCurrent_z(nrxx);
        std::vector<std::complex<double>> kF_r(nrxx);
        std::vector<std::complex<double>> rCDPotential(nrxx);

        for (int ir = 0; ir < nrxx; ++ir)
        {
            kF_r[ir]=std::pow(3*std::pow(ModuleBase::PI*std::abs(rPhi[is][ir]),2),1.0/3.0);
        }

        pw_rho->real2recip(rPhi[is].data(), recipPhi[is].data());
        for (int ik = 0; ik < npw; ++ik)
        {
            recipCurrent_x[ik]=imag*gcar[ik].x*recipPhi[is][ik]*tpiba;
            recipCurrent_y[ik]=imag*gcar[ik].y*recipPhi[is][ik]*tpiba;
            recipCurrent_z[ik]=imag*gcar[ik].z*recipPhi[is][ik]*tpiba;
        }
        pw_rho->recip2real(recipCurrent_x.data(),rCurrent_x.data());
        pw_rho->recip2real(recipCurrent_y.data(),rCurrent_y.data());
        pw_rho->recip2real(recipCurrent_z.data(),rCurrent_z.data());
        for (int ir = 0; ir < nrxx; ++ir)
        {
            rCurrent_x[ir]=std::imag(rCurrent_x[ir]*std::conj(rPhi[is][ir]));
            rCurrent_y[ir]=std::imag(rCurrent_y[ir]*std::conj(rPhi[is][ir]));
            rCurrent_z[ir]=std::imag(rCurrent_z[ir]*std::conj(rPhi[is][ir]));
        }
        pw_rho->real2recip(rCurrent_x.data(),recipCurrent_x.data());
        pw_rho->real2recip(rCurrent_y.data(),recipCurrent_y.data());
        pw_rho->real2recip(rCurrent_z.data(),recipCurrent_z.data());
        for (int ik = 0; ik < npw; ++ik)
        {
            recipCDPotential[ik]=recipCurrent_x[ik]*gcar[ik].x+recipCurrent_y[ik]*gcar[ik].y+recipCurrent_z[ik]*gcar[ik].z;
            if (gg[ik]==0) 
            {
                recipCDPotential[ik]=0.0;
            }
            else
            {
                recipCDPotential[ik]*=imag/sqrt(gg[ik]);
            }
        }
        pw_rho->recip2real(recipCDPotential.data(),rCDPotential.data());

        for (int ir = 0; ir < nrxx; ++ir)
        {
            rpot(is, ir) -= mCD_para*2.0*std::real(rCDPotential[ir])*std::pow(ModuleBase::PI,3)
                        / (2.0*std::pow(std::real(kF_r[ir]),2));
            if (std::isnan(rpot(is, ir)))
            {
                rpot(is, ir)=0.0;
            }
        }
    }
}

void Evolve_OFDFT::propagate_psi_RK4(elecstate::ElecState* pelec, 
                                 Charge& chr, 
                                 UnitCell& ucell, 
                                 std::vector<std::complex<double>>& pphi_, 
                                 ModulePW::PW_Basis* pw_rho)
{
    ModuleBase::timer::start("ESolver_OF_TDDFT", "propagate_psi_RK4");

    const std::complex<double> imag(0.0,1.0);
    const double dt=PARAM.inp.mdp.md_dt / ModuleBase::AU_to_FS;
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
    for (int is = 0; is < nspin; ++is){
        for (int ir = 0; ir < nrxx; ++ir)
        {
            K1[is * nrxx + ir]=-0.5*K1[is * nrxx + ir]*dt*imag;   // 0.5 convert Ry to Hartree
            psi1[is * nrxx + ir]=pphi_[is * nrxx + ir]+0.5*K1[is * nrxx + ir];
        }
    }
    cal_Hpsi(pelec,chr,ucell,psi1,pw_rho,K2);
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int is = 0; is < nspin; ++is){
        for (int ir = 0; ir < nrxx; ++ir)
        {
            K2[is * nrxx + ir]=-0.5*K2[is * nrxx + ir]*dt*imag;
            psi2[is * nrxx + ir]=pphi_[is * nrxx + ir]+0.5*K2[is * nrxx + ir];
        }
    }
    cal_Hpsi(pelec,chr,ucell,psi2,pw_rho,K3);
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int is = 0; is < nspin; ++is){
        for (int ir = 0; ir < nrxx; ++ir)
        {
            K3[is * nrxx + ir]=-0.5*K3[is * nrxx + ir]*dt*imag;
            psi3[is * nrxx + ir]=pphi_[is * nrxx + ir]+K3[is * nrxx + ir];
        }
    }
    cal_Hpsi(pelec,chr,ucell,psi3,pw_rho,K4);
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int is = 0; is < nspin; ++is){
        for (int ir = 0; ir < nrxx; ++ir)
        {
            K4[is * nrxx + ir]=-0.5*K4[is * nrxx + ir]*dt*imag;
            pphi_[is * nrxx + ir]+=1.0/6.0*(K1[is * nrxx + ir]+2.0*K2[is * nrxx + ir]+2.0*K3[is * nrxx + ir]+K4[is * nrxx + ir]);
        }
    }

#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < nrxx; ++ir)
        {
            chr.rho[is][ir] = abs(pphi_[is * nrxx + ir])*abs(pphi_[is * nrxx + ir]);
        }
    }
    this->renormalize_psi(chr, pw_rho, pphi_);

    ModuleBase::timer::end("ESolver_OF_TDDFT", "propagate_psi_RK4");
}

void Evolve_OFDFT::propagate_psi_RK2(elecstate::ElecState* pelec, 
                                 Charge& chr, 
                                 UnitCell& ucell, 
                                 std::vector<std::complex<double>>& pphi_, 
                                 ModulePW::PW_Basis* pw_rho)
{
    ModuleBase::timer::start("ESolver_OF_TDDFT", "propagate_psi_RK2");

    const std::complex<double> imag(0.0, 1.0);
    const double dt=PARAM.inp.mdp.md_dt / ModuleBase::AU_to_FS;
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

    ModuleBase::timer::end("ESolver_OF_TDDFT", "propagate_psi_RK2");
}
