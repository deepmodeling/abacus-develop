#include "charge_mixing.h"

#include "module_parameter/parameter.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

void Charge_Mixing::Kerker_screen_recip(std::complex<double>* drhog)
{
    if (this->mixing_gg0 <= 0.0 || this->mixing_beta <= 0.1) {
        return;
}
    double fac = 0.0;
    double gg0 = 0.0;
    double amin = 0.0;

    /// consider a resize for mixing_angle
    int resize_tmp = 1;
    if (PARAM.inp.nspin == 4 && this->mixing_angle > 0) { resize_tmp = 2;
}

    /// implement Kerker for density and magnetization separately
    for (int is = 0; is < PARAM.inp.nspin / resize_tmp; ++is)
    {
        /// new mixing method only support nspin=2 not nspin=4
        if (is >= 1)
        {
            if (this->mixing_gg0_mag <= 0.0001 || this->mixing_beta_mag <= 0.1)
            {
#ifdef __DEBUG
                assert(is == 1); // make sure break works
#endif
                double is_mag = PARAM.inp.nspin - 1;
                //for (int ig = 0; ig < this->rhopw->npw * is_mag; ig++)
                //{
                //    drhog[is * this->rhopw->npw + ig] *= 1;
                //}
                break;
            }
            fac = this->mixing_gg0_mag;
            amin = this->mixing_beta_mag;
        }
        else
        {
            fac = this->mixing_gg0;
            amin = this->mixing_beta;
        }

        gg0 = std::pow(fac * 0.529177 / *this->tpiba, 2);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
        for (int ig = 0; ig < this->rhopw->npw; ++ig)
        {
            double gg = this->rhopw->gg[ig];
            double filter_g = std::max(gg / (gg + gg0), this->mixing_gg0_min / amin);
            drhog[is * this->rhopw->npw + ig] *= filter_g;
        }
    }
    return;
}

void Charge_Mixing::Kerker_screen_real(double* drhor)
{
    if (this->mixing_gg0 <= 0.0001 || this->mixing_beta <= 0.1)
        return;

    const int spin_channels = PARAM.inp.nspin == 4 && this->mixing_angle > 0 ? 2 : PARAM.inp.nspin;
    const bool screen_magnetization = this->mixing_gg0_mag > 0.0001 && this->mixing_beta_mag > 0.1;
    const int screened_channels = screen_magnetization ? spin_channels : 1;
    std::vector<std::complex<double>> drhog(this->rhopw->npw);
    std::vector<double> drhor_filter(this->rhopw->nrxx);
    for (int spin = 0; spin < screened_channels; ++spin)
    {
        double* residual = drhor + spin * this->rhopw->nrxx;
        const double screening = spin == 0 ? this->mixing_gg0 : this->mixing_gg0_mag;
        const double beta = spin == 0 ? this->mixing_beta : this->mixing_beta_mag;
        const double gg0 = std::pow(screening * 0.529177 / *this->tpiba, 2);
        this->rhopw->real2recip(residual, drhog.data());
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
        for (int ig = 0; ig < this->rhopw->npw; ++ig)
        {
            const double gg = this->rhopw->gg[ig];
            const double filter = std::max(gg / (gg + gg0), this->mixing_gg0_min / beta);
            drhog[ig] *= 1 - filter;
        }
        // Subtracting the low-frequency correction preserves residual components above the basis cutoff.
        // A disabled magnetic filter is the identity and requires no transform.
        this->rhopw->recip2real(drhog.data(), drhor_filter.data());
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
        for (int ir = 0; ir < this->rhopw->nrxx; ++ir)
        {
            residual[ir] -= drhor_filter[ir];
        }
    }
}
