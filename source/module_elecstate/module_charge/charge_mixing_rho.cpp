#include "charge_mixing.h"

#include "module_parameter/parameter.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

void Charge_Mixing::mix_rho_recip(Charge* chr)
{
    std::complex<double>* rhog_in = nullptr;
    std::complex<double>* rhog_out = nullptr;
    // for smooth part
    std::complex<double>* rhogs_in = chr->rhog_save[0];
    std::complex<double>* rhogs_out = chr->rhog[0];
    // for high_frequency part
    std::complex<double>* rhoghf_in = nullptr;
    std::complex<double>* rhoghf_out = nullptr;

    if ( PARAM.globalv.double_grid)
    {
        // divide into smooth part and high_frequency part
        divide_data(chr->rhog_save[0], rhogs_in, rhoghf_in);
        divide_data(chr->rhog[0], rhogs_out, rhoghf_out);
    }

    //  inner_product_recip_hartree is a hartree-like sum, unit is Ry
    auto inner_product
        = std::bind(&Charge_Mixing::inner_product_recip_hartree, this, std::placeholders::_1, std::placeholders::_2);

    // DIIS Mixing Only for smooth part, while high_frequency part is mixed by plain mixing method.
    if (PARAM.inp.nspin == 1)
    {
        rhog_in = rhogs_in;
        rhog_out = rhogs_out;
        auto screen = std::bind(&Charge_Mixing::Kerker_screen_recip, this, std::placeholders::_1);
        this->mixing->push_data(this->rho_mdata, rhog_in, rhog_out, screen, true);
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhog_out);
    }
    else if (PARAM.inp.nspin == 2)
    {
        // magnetic density
        std::complex<double> *rhog_mag = nullptr;
        std::complex<double> *rhog_mag_save = nullptr;
        const int npw = this->rhopw->npw;
        // allocate rhog_mag[is*ngmc] and rhog_mag_save[is*ngmc]
        rhog_mag = new std::complex<double>[npw * PARAM.inp.nspin];
        rhog_mag_save = new std::complex<double>[npw * PARAM.inp.nspin];
        ModuleBase::GlobalFunc::ZEROS(rhog_mag, npw * PARAM.inp.nspin);
        ModuleBase::GlobalFunc::ZEROS(rhog_mag_save, npw * PARAM.inp.nspin);
        // get rhog_mag[is*ngmc] and rhog_mag_save[is*ngmc]
        for (int ig = 0; ig < npw; ig++)
        {
            rhog_mag[ig] = chr->rhog[0][ig] + chr->rhog[1][ig];
            rhog_mag_save[ig] = chr->rhog_save[0][ig] + chr->rhog_save[1][ig];
        }
        for (int ig = 0; ig < npw; ig++)
        {
            rhog_mag[ig + npw] = chr->rhog[0][ig] - chr->rhog[1][ig];
            rhog_mag_save[ig + npw] = chr->rhog_save[0][ig] - chr->rhog_save[1][ig];
        }
        //
        rhog_in = rhog_mag_save;
        rhog_out = rhog_mag;
        //
        auto screen = std::bind(&Charge_Mixing::Kerker_screen_recip, this, std::placeholders::_1);
        auto twobeta_mix
            = [this, npw](std::complex<double>* out, const std::complex<double>* in, const std::complex<double>* sres) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
                  for (int i = 0; i < npw; ++i)
                  {
                      out[i] = in[i] + this->mixing_beta * sres[i];
                  }
            // magnetism
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
                  for (int i = npw; i < 2 * npw; ++i)
                  {
                      out[i] = in[i] + this->mixing_beta_mag * sres[i];
                  }
              };
        this->mixing->push_data(this->rho_mdata, rhog_in, rhog_out, screen, twobeta_mix, true);
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhog_out);
        // get rhog[is][ngmc] from rhog_mag[is*ngmc]
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(chr->rhog[is], npw);
        }
        for (int ig = 0; ig < npw; ig++)
        {
            chr->rhog[0][ig] = 0.5 * (rhog_mag[ig] + rhog_mag[ig+npw]);
            chr->rhog[1][ig] = 0.5 * (rhog_mag[ig] - rhog_mag[ig+npw]);
        }
        // delete
        delete[] rhog_mag;
        delete[] rhog_mag_save;
        // get rhogs_out for combine_data()
        if ( PARAM.globalv.double_grid)
        {
            for (int ig = 0; ig < npw; ig++)
            {
                rhogs_out[ig] = chr->rhog[0][ig];
                rhogs_out[ig + npw] = chr->rhog[1][ig];
            }
        }
    }
    else if (PARAM.inp.nspin == 4 && PARAM.inp.mixing_angle <= 0)
    {
        // normal broyden mixing for {rho, mx, my, mz}
        rhog_in = rhogs_in;
        rhog_out = rhogs_out;
        const int npw = this->rhopw->npw;
        auto screen = std::bind(&Charge_Mixing::Kerker_screen_recip, this, std::placeholders::_1); // use old one
        auto twobeta_mix
            = [this, npw](std::complex<double>* out, const std::complex<double>* in, const std::complex<double>* sres) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
                  for (int i = 0; i < npw; ++i)
                  {
                      out[i] = in[i] + this->mixing_beta * sres[i];
                  }
            // magnetism, mx, my, mz
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
                  for (int i = npw; i < 4 * npw; ++i)
                  {
                      out[i] = in[i] + this->mixing_beta_mag * sres[i];
                  }
              };
        this->mixing->push_data(this->rho_mdata, rhog_in, rhog_out, screen, twobeta_mix, true);
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhog_out);
    }
    else if (PARAM.inp.nspin == 4 && PARAM.inp.mixing_angle > 0)
    {
        // special broyden mixing for {rho, |m|} proposed by J. Phys. Soc. Jpn. 82 (2013) 114706
        // here only consider the case of mixing_angle = 1, which mean only change |m| and keep angle fixed
        // old support see mix_rho_recip()
        if ( PARAM.globalv.double_grid)
        {
            ModuleBase::WARNING_QUIT("Charge_Mixing", "double_grid is not supported for new mixing method yet.");
        }
        // allocate memory for rho_magabs and rho_magabs_save
        const int nrxx = this->rhopw->nrxx;
        double* rho_magabs = new double[nrxx];
        double* rho_magabs_save = new double[nrxx];
        ModuleBase::GlobalFunc::ZEROS(rho_magabs, nrxx);
        ModuleBase::GlobalFunc::ZEROS(rho_magabs_save, nrxx);
        // calculate rho_magabs and rho_magabs_save
        for (int ir = 0; ir < nrxx; ir++)
        {
            // |m| for rho
            rho_magabs[ir] = std::sqrt(chr->rho[1][ir] * chr->rho[1][ir] + chr->rho[2][ir] * chr->rho[2][ir] + chr->rho[3][ir] * chr->rho[3][ir]);
            // |m| for rho_save
            rho_magabs_save[ir] = std::sqrt(chr->rho_save[1][ir] * chr->rho_save[1][ir] + chr->rho_save[2][ir] * chr->rho_save[2][ir] + chr->rho_save[3][ir] * chr->rho_save[3][ir]);
        }
        // allocate memory for rhog_magabs and rhog_magabs_save
        const int npw = this->rhopw->npw;
        std::complex<double>* rhog_magabs = new std::complex<double>[npw * 2];
        std::complex<double>* rhog_magabs_save = new std::complex<double>[npw * 2];
        ModuleBase::GlobalFunc::ZEROS(rhog_magabs, npw * 2);
        ModuleBase::GlobalFunc::ZEROS(rhog_magabs_save, npw * 2);
        // calculate rhog_magabs and rhog_magabs_save
        for (int ig = 0; ig < npw; ig++)
        {
            rhog_magabs[ig] = chr->rhog[0][ig]; // rho
            rhog_magabs_save[ig] = chr->rhog_save[0][ig]; // rho_save
        }
        // FT to get rhog_magabs and rhog_magabs_save
        this->rhopw->real2recip(rho_magabs, rhog_magabs + this->rhopw->npw);
        this->rhopw->real2recip(rho_magabs_save, rhog_magabs_save + this->rhopw->npw);
        //
        rhog_in = rhog_magabs_save;
        rhog_out = rhog_magabs;
        auto screen = std::bind(&Charge_Mixing::Kerker_screen_recip, this, std::placeholders::_1); // use old one
        auto twobeta_mix
            = [this, npw](std::complex<double>* out, const std::complex<double>* in, const std::complex<double>* sres) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
                  for (int i = 0; i < npw; ++i)
                  {
                      out[i] = in[i] + this->mixing_beta * sres[i];
                  }
            // magnetism, |m|
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
                  for (int i = npw; i < 2 * npw; ++i)
                  {
                      out[i] = in[i] + this->mixing_beta_mag * sres[i];
                  }
              };
        this->mixing->push_data(this->rho_mdata, rhog_in, rhog_out, screen, twobeta_mix, true);
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhog_out);
        // get new |m| in real space using FT
        this->rhopw->recip2real(rhog_magabs + this->rhopw->npw, rho_magabs);
        // use new |m| and angle to update {mx, my, mz}
        for (int ig = 0; ig < npw; ig++)
        {
            chr->rhog[0][ig] = rhog_magabs[ig]; // rhog
            double norm = std::sqrt(chr->rho[1][ig] * chr->rho[1][ig] + chr->rho[2][ig] * chr->rho[2][ig] + chr->rho[3][ig] * chr->rho[3][ig]);
            if (std::abs(norm) < 1e-10) { continue;
}
            double rescale_tmp = rho_magabs[npw + ig] / norm; 
            chr->rho[1][ig] *= rescale_tmp;
            chr->rho[2][ig] *= rescale_tmp;
            chr->rho[3][ig] *= rescale_tmp;
        }
        // delete
        delete[] rhog_magabs;
        delete[] rhog_magabs_save;
        delete[] rho_magabs;
        delete[] rho_magabs_save;
    }

    if ( PARAM.globalv.double_grid)
    {
        // plain mixing for high_frequencies
        const int ndimhf = (this->rhodpw->npw - this->rhopw->npw) * PARAM.inp.nspin;
        this->mixing_highf->plain_mix(rhoghf_out, rhoghf_in, rhoghf_out, ndimhf, nullptr);

        // combine smooth part and high_frequency part
        combine_data(chr->rhog[0], rhogs_out, rhoghf_out);
        clean_data(rhogs_in, rhoghf_in);
    }

    // rhog to rho
    if (PARAM.inp.nspin == 4 && PARAM.inp.mixing_angle > 0)
    {
        // only tranfer rhog[0]
        // do not support double_grid, use rhopw directly
        chr->rhopw->recip2real(chr->rhog[0], chr->rho[0]);
    }
    else
    {
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            // use rhodpw for double_grid
            // rhodpw is the same as rhopw for ! PARAM.globalv.double_grid
            this->rhodpw->recip2real(chr->rhog[is], chr->rho[is]);
        }
    }

    // For kinetic energy density
    if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
    {
        std::vector<std::complex<double>> kin_g(PARAM.inp.nspin * rhodpw->npw);
        std::vector<std::complex<double>> kin_g_save(PARAM.inp.nspin * rhodpw->npw);
        // FFT to get kin_g and kin_g_save
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            rhodpw->real2recip(chr->kin_r[is], &kin_g[is * rhodpw->npw]);
            rhodpw->real2recip(chr->kin_r_save[is], &kin_g_save[is * rhodpw->npw]);
        }
        // for smooth part, for ! PARAM.globalv.double_grid only have this part
        std::complex<double>*taugs_in = kin_g_save.data(), *taugs_out = kin_g.data();
        // for high frequency part
        std::complex<double>*taughf_in = nullptr, *taughf_out = nullptr;
        if ( PARAM.globalv.double_grid)
        {
            // divide into smooth part and high_frequency part
            divide_data(kin_g_save.data(), taugs_in, taughf_in);
            divide_data(kin_g.data(), taugs_out, taughf_out);
        }

        // Note: there is no kerker modification for tau because I'm not sure
        // if we should have it. If necessary we can try it in the future.
        this->mixing->push_data(this->tau_mdata, taugs_in, taugs_out, nullptr, false);

        this->mixing->mix_data(this->tau_mdata, taugs_out);

        if ( PARAM.globalv.double_grid)
        {
            // simple mixing for high_frequencies
            const int ndimhf = (this->rhodpw->npw - this->rhopw->npw) * PARAM.inp.nspin;
            this->mixing_highf->plain_mix(taughf_out, taughf_in, taughf_out, ndimhf, nullptr);

            // combine smooth part and high_frequency part
            combine_data(kin_g.data(), taugs_out, taughf_out);
            clean_data(taugs_in, taughf_in);
        }

        // kin_g to kin_r
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            rhodpw->recip2real(&kin_g[is * rhodpw->npw], chr->kin_r[is]);
        }
    }

#ifdef USE_PAW
    if(PARAM.inp.use_paw)
    {
        double *nhat_out, *nhat_in;
        nhat_in = chr->nhat_save[0];
        nhat_out = chr->nhat[0];
        // Note: there is no kerker modification for tau because I'm not sure
        // if we should have it. If necessary we can try it in the future.
        this->mixing->push_data(this->nhat_mdata, nhat_in, nhat_out, nullptr, false);

        this->mixing->mix_data(this->nhat_mdata, nhat_out);
    }
#endif

    return;
}

void Charge_Mixing::mix_rho_real(Charge* chr)
{
    double* rhor_in;
    double* rhor_out;
    if (PARAM.inp.nspin == 1)
    {
        rhor_in = chr->rho_save[0];
        rhor_out = chr->rho[0];
        auto screen = std::bind(&Charge_Mixing::Kerker_screen_real, this, std::placeholders::_1);
        this->mixing->push_data(this->rho_mdata, rhor_in, rhor_out, screen, true);    
        auto inner_product
            = std::bind(&Charge_Mixing::inner_product_real, this, std::placeholders::_1, std::placeholders::_2);
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhor_out);
    }
    else if (PARAM.inp.nspin == 2)
    {
        // magnetic density
        double *rho_mag = nullptr;
        double *rho_mag_save = nullptr; 
        const int nrxx = this->rhopw->nrxx;
        // allocate rho_mag[is*nnrx] and rho_mag_save[is*nnrx]
        rho_mag = new double[nrxx * PARAM.inp.nspin];
        rho_mag_save = new double[nrxx * PARAM.inp.nspin];
        ModuleBase::GlobalFunc::ZEROS(rho_mag, nrxx * PARAM.inp.nspin);
        ModuleBase::GlobalFunc::ZEROS(rho_mag_save, nrxx * PARAM.inp.nspin);
        // get rho_mag[is*nnrx] and rho_mag_save[is*nnrx]
        for (int ir = 0; ir < nrxx; ir++)
        {
            rho_mag[ir] = chr->rho[0][ir] + chr->rho[1][ir];
            rho_mag_save[ir] = chr->rho_save[0][ir] + chr->rho_save[1][ir];
        }
        for (int ir = 0; ir < nrxx; ir++)
        {
            rho_mag[ir + nrxx] = chr->rho[0][ir] - chr->rho[1][ir];
            rho_mag_save[ir + nrxx] = chr->rho_save[0][ir] - chr->rho_save[1][ir];
        }
        //
        rhor_in = rho_mag_save;
        rhor_out = rho_mag;
        auto screen = std::bind(&Charge_Mixing::Kerker_screen_real, this, std::placeholders::_1);
        auto twobeta_mix
            = [this, nrxx](double* out, const double* in, const double* sres) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
            for (int i = 0; i < nrxx; ++i)
            {
                out[i] = in[i] + this->mixing_beta * sres[i];
            }
            // magnetism
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
            for (int i = nrxx; i < 2 * nrxx; ++i)
            {
                out[i] = in[i] + this->mixing_beta_mag * sres[i];
            }
        };
        this->mixing->push_data(this->rho_mdata, rhor_in, rhor_out, screen, twobeta_mix, true);
        auto inner_product
            = std::bind(&Charge_Mixing::inner_product_real, this, std::placeholders::_1, std::placeholders::_2);
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhor_out);
        // get new rho[is][nrxx] from rho_mag[is*nrxx]
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(chr->rho[is], nrxx);
            //ModuleBase::GlobalFunc::ZEROS(rho_save[is], nrxx);
        }
        for (int ir = 0; ir < nrxx; ir++)
        {
            chr->rho[0][ir] = 0.5 * (rho_mag[ir] + rho_mag[ir+nrxx]);
            chr->rho[1][ir] = 0.5 * (rho_mag[ir] - rho_mag[ir+nrxx]);
        }
        // delete
        delete[] rho_mag;
        delete[] rho_mag_save;
    }
    else if (PARAM.inp.nspin == 4 && PARAM.inp.mixing_angle <= 0)
    {
        // normal broyden mixing for {rho, mx, my, mz}
        rhor_in = chr->rho_save[0];
        rhor_out = chr->rho[0];
        const int nrxx = this->rhopw->nrxx;
        auto screen = std::bind(&Charge_Mixing::Kerker_screen_real, this, std::placeholders::_1);
        auto twobeta_mix
            = [this, nrxx](double* out, const double* in, const double* sres) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
            for (int i = 0; i < nrxx; ++i)
            {
                out[i] = in[i] + this->mixing_beta * sres[i];
            }
            // magnetism, mx, my, mz
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
            for (int i = nrxx; i < 4 * nrxx; ++i)
            {
                out[i] = in[i] + this->mixing_beta_mag * sres[i];
            }
        };
        this->mixing->push_data(this->rho_mdata, rhor_in, rhor_out, screen, twobeta_mix, true);
        auto inner_product
            = std::bind(&Charge_Mixing::inner_product_real, this, std::placeholders::_1, std::placeholders::_2);
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhor_out);
    }
    else if (PARAM.inp.nspin == 4 && PARAM.inp.mixing_angle > 0)
    {
        // special broyden mixing for {rho, |m|} proposed by J. Phys. Soc. Jpn. 82 (2013) 114706
        // here only consider the case of mixing_angle = 1, which mean only change |m| and keep angle fixed
        const int nrxx = this->rhopw->nrxx;
        // allocate memory for rho_magabs and rho_magabs_save
        double* rho_magabs = new double[nrxx * 2];
        double* rho_magabs_save = new double[nrxx * 2];
        ModuleBase::GlobalFunc::ZEROS(rho_magabs, nrxx * 2);
        ModuleBase::GlobalFunc::ZEROS(rho_magabs_save, nrxx * 2);
        // calculate rho_magabs and rho_magabs_save
        for (int ir = 0; ir < nrxx; ir++)
        {
            rho_magabs[ir] = chr->rho[0][ir]; // rho
            rho_magabs_save[ir] = chr->rho_save[0][ir]; // rho_save
            // |m| for rho
            rho_magabs[nrxx + ir] = std::sqrt(chr->rho[1][ir] * chr->rho[1][ir] + chr->rho[2][ir] * chr->rho[2][ir] + chr->rho[3][ir] * chr->rho[3][ir]);
            // |m| for rho_save
            rho_magabs_save[nrxx + ir] = std::sqrt(chr->rho_save[1][ir] * chr->rho_save[1][ir] + chr->rho_save[2][ir] * chr->rho_save[2][ir] + chr->rho_save[3][ir] * chr->rho_save[3][ir]);
        }
        rhor_in = rho_magabs_save;
        rhor_out = rho_magabs;
        auto screen = std::bind(&Charge_Mixing::Kerker_screen_real, this, std::placeholders::_1);
        auto twobeta_mix
            = [this, nrxx](double* out, const double* in, const double* sres) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
            for (int i = 0; i < nrxx; ++i)
            {
                out[i] = in[i] + this->mixing_beta * sres[i];
            }
            // magnetism, |m|
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
            for (int i = nrxx; i < 2 * nrxx; ++i)
            {
                out[i] = in[i] + this->mixing_beta_mag * sres[i];
            }
        };
        this->mixing->push_data(this->rho_mdata, rhor_in, rhor_out, screen, twobeta_mix, true);
        auto inner_product
            = std::bind(&Charge_Mixing::inner_product_real, this, std::placeholders::_1, std::placeholders::_2);
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhor_out);
        // use new |m| and angle to update {mx, my, mz}
        for (int ir = 0; ir < nrxx; ir++)
        {
            chr->rho[0][ir] = rho_magabs[ir]; // rho
            double norm = std::sqrt(chr->rho[1][ir] * chr->rho[1][ir] + chr->rho[2][ir] * chr->rho[2][ir] + chr->rho[3][ir] * chr->rho[3][ir]);
            if (norm < 1e-10) { continue;
}
            double rescale_tmp = rho_magabs[nrxx + ir] / norm; 
            chr->rho[1][ir] *= rescale_tmp;
            chr->rho[2][ir] *= rescale_tmp;
            chr->rho[3][ir] *= rescale_tmp;
        }
        // delete
        delete[] rho_magabs;
        delete[] rho_magabs_save;
    }
    
    double *taur_out, *taur_in;
    if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
    {
        taur_in = chr->kin_r_save[0];
        taur_out = chr->kin_r[0];
        // Note: there is no kerker modification for tau because I'm not sure
        // if we should have it. If necessary we can try it in the future.
        this->mixing->push_data(this->tau_mdata, taur_in, taur_out, nullptr, false);

        this->mixing->mix_data(this->tau_mdata, taur_out);
    }

}

void Charge_Mixing::mix_rho(Charge* chr)
{
    ModuleBase::TITLE("Charge_Mixing", "mix_rho");
    ModuleBase::timer::tick("Charge", "mix_rho");

    // the charge before mixing.
    const int nrxx = chr->rhopw->nrxx;
    std::vector<double> rho123(PARAM.inp.nspin * nrxx);
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        if (is == 0 || is == 3 || !PARAM.globalv.domag_z)
        {
            double* rho123_is = rho123.data() + is * nrxx;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
            for(int ir = 0 ; ir < nrxx ; ++ir)
            {
                rho123_is[ir] = chr->rho[is][ir];
            }
        }
    }
    std::vector<double> kin_r123;
    if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
    {
        kin_r123.resize(PARAM.inp.nspin * nrxx);
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            double* kin_r123_is = kin_r123.data() + is * nrxx;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
            for(int ir = 0 ; ir < nrxx ; ++ir)
            {
                kin_r123_is[ir] = chr->kin_r[is][ir];
            }
        }
    }
#ifdef USE_PAW
    std::vector<double> nhat_r123;
    if(PARAM.inp.use_paw)
    {
        nhat_r123.resize(PARAM.inp.nspin * nrxx);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
        for(int ir = 0 ; ir < nrxx ; ++ir)
        {
            for(int is = 0; is < PARAM.inp.nspin; ++is)
            {
                nhat_r123[ir+is*nrxx] = chr->nhat[0][ir];
            }
        }
    }        
#endif
    // --------------------Mixing Body--------------------
    if (PARAM.inp.scf_thr_type == 1)
    {
        mix_rho_recip(chr);
    }
    else if (PARAM.inp.scf_thr_type == 2)
    {
        mix_rho_real(chr);
    }
    // ---------------------------------------------------

    // mohan add 2012-06-05
    // rho_save is the charge before mixing
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        if (is == 0 || is == 3 || !PARAM.globalv.domag_z)
        {
            double* rho123_is = rho123.data() + is * nrxx;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
            for(int ir = 0 ; ir < nrxx ; ++ir)
            {
                chr->rho_save[is][ir] = rho123_is[ir];
            }
        }
    }

    if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
    {
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            double* kin_r123_is = kin_r123.data() + is * nrxx;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
            for(int ir = 0 ; ir < nrxx ; ++ir)
            {
                chr->kin_r_save[is][ir] = kin_r123_is[ir];
            }
        }
    }

#ifdef USE_PAW
    if(PARAM.inp.use_paw)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
        for(int ir = 0 ; ir < nrxx ; ++ir)
        {
            for(int is = 0; is < PARAM.inp.nspin; ++is)
            {
                chr->nhat_save[is][ir] = nhat_r123[ir+is*nrxx];
            }
        }
    }
#endif

    if (new_e_iteration) {
        new_e_iteration = false;
}
    ModuleBase::timer::tick("Charge", "mix_rho");
    return;
}