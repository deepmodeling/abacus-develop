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

        gg0 = std::pow(fac * 0.529177 / GlobalC::ucell.tpiba, 2);
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
    if (this->mixing_gg0 <= 0.0001 || this->mixing_beta <= 0.1) {
        return;
}
    /// consider a resize for mixing_angle
    int resize_tmp = 1;
    if (PARAM.inp.nspin == 4 && this->mixing_angle > 0) { resize_tmp = 2;
}
    
    std::vector<std::complex<double>> drhog(this->rhopw->npw * PARAM.inp.nspin / resize_tmp);
    std::vector<double> drhor_filter(this->rhopw->nrxx * PARAM.inp.nspin / resize_tmp);
    for (int is = 0; is < PARAM.inp.nspin / resize_tmp; ++is)
    {
        // Note after this process some G which is higher than Gmax will be filtered.
        // Thus we cannot use Kerker_screen_recip(drhog.data()) directly after it.
        this->rhopw->real2recip(drhor + is * this->rhopw->nrxx, drhog.data() + is * this->rhopw->npw);
    }
    /// implement Kerker for density and magnetization separately
    double fac = 0.0;
    double gg0 = 0.0;
    double amin = 0.0;

    for (int is = 0; is < PARAM.inp.nspin / resize_tmp; is++)
    {

        if (is >= 1)
        {
            if (this->mixing_gg0_mag <= 0.0001 || this->mixing_beta_mag <= 0.1)
            {
#ifdef __DEBUG
                assert(is == 1); /// make sure break works
#endif
                double is_mag = PARAM.inp.nspin - 1;
                if (PARAM.inp.nspin == 4 && this->mixing_angle > 0) { is_mag = 1;
}
                for (int ig = 0; ig < this->rhopw->npw * is_mag; ig++)
                {
                    drhog[is * this->rhopw->npw + ig] = 0;
                }
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
        
        gg0 = std::pow(fac * 0.529177 / GlobalC::ucell.tpiba, 2);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
        for (int ig = 0; ig < this->rhopw->npw; ig++)
        {
            double gg = this->rhopw->gg[ig];
            // I have not decided how to handle gg=0 part, will be changed in future
            //if (gg == 0)
            //{
            //    drhog[is * this->rhopw->npw + ig] *= 0;
            //    continue;
            //}
            double filter_g = std::max(gg / (gg + gg0), this->mixing_gg0_min / amin);
            drhog[is * this->rhopw->npw + ig] *= (1 - filter_g);
        }
    }
    /// inverse FT
    for (int is = 0; is < PARAM.inp.nspin / resize_tmp; ++is)
    {
        this->rhopw->recip2real(drhog.data() + is * this->rhopw->npw, drhor_filter.data() + is * this->rhopw->nrxx);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
    for (int ir = 0; ir < this->rhopw->nrxx * PARAM.inp.nspin / resize_tmp; ir++)
    {
        drhor[ir] -= drhor_filter[ir];
    }
}
