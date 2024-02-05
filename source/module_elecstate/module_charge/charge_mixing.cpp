#include "charge_mixing.h"

#include "module_base/element_elec_config.h"
#include "module_base/inverse_matrix.h"
#include "module_base/module_mixing/broyden_mixing.h"
#include "module_base/module_mixing/plain_mixing.h"
#include "module_base/module_mixing/pulay_mixing.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

Charge_Mixing::Charge_Mixing()
{
}

Charge_Mixing::~Charge_Mixing()
{
    delete this->mixing;
    delete this->mixing_highf;
}

void Charge_Mixing::set_mixing(const std::string& mixing_mode_in,
                               const double& mixing_beta_in,
                               const int& mixing_ndim_in,
                               const double& mixing_gg0_in,
                               const bool& mixing_tau_in,
                               const double& mixing_beta_mag_in,
                               const double& mixing_gg0_mag_in,
                               const double& mixing_gg0_min_in,
                               const double& mixing_angle_in,
                               const bool& mixing_dmr_in)
{
    // get private mixing parameters
    this->mixing_mode = mixing_mode_in;
    this->mixing_beta = mixing_beta_in;
    this->mixing_beta_mag = mixing_beta_mag_in;
    this->mixing_ndim = mixing_ndim_in;
    this->mixing_gg0 = mixing_gg0_in;
    this->mixing_tau = mixing_tau_in;
    this->mixing_gg0_mag = mixing_gg0_mag_in;
    this->mixing_gg0_min = mixing_gg0_min_in;
    this->mixing_angle = mixing_angle_in;
    this->mixing_dmr = mixing_dmr_in;

    // check the paramters
    if (this->mixing_beta > 1.0 || this->mixing_beta < 0.0)
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing", "You'd better set mixing_beta to [0.0, 1.0]!");
    }
    if (GlobalV::NSPIN >= 2 && this->mixing_beta_mag < 0.0)
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing", "You'd better set mixing_beta_mag >= 0.0!");
    }

    if (!(this->mixing_mode == "plain" || this->mixing_mode == "broyden" || this->mixing_mode == "pulay"))
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing", "This Mixing mode is not implemended yet,coming soon.");
    }

    // print into running.log
    GlobalV::ofs_running<<"\n----------- Double Check Mixing Parameters Begin ------------"<<std::endl;
    GlobalV::ofs_running<<"mixing_type: "<< this->mixing_mode <<std::endl;
    GlobalV::ofs_running<<"mixing_beta: "<< this->mixing_beta <<std::endl;
    GlobalV::ofs_running<<"mixing_gg0: "<< this->mixing_gg0 <<std::endl;
    GlobalV::ofs_running<<"mixing_gg0_min: "<< GlobalV::MIXING_GG0_MIN <<std::endl;
    if (GlobalV::NSPIN==2 || GlobalV::NSPIN==4)
    {
        GlobalV::ofs_running<<"mixing_beta_mag: "<< this->mixing_beta_mag <<std::endl;
        GlobalV::ofs_running<<"mixing_gg0_mag: "<< GlobalV::MIXING_GG0_MAG <<std::endl;
    }
    if (GlobalV::MIXING_ANGLE > 0)
    {
        GlobalV::ofs_running<<"mixing_angle: "<< GlobalV::MIXING_ANGLE <<std::endl;
    }
    GlobalV::ofs_running<<"mixing_ndim: "<< this->mixing_ndim <<std::endl;
    GlobalV::ofs_running<<"----------- Double Check Mixing Parameters End ------------"<<std::endl;

    return;
}

void Charge_Mixing::init_mixing()
{
    // this init should be called at the 1-st iteration of each scf loop

    ModuleBase::TITLE("Charge_Mixing", "init_mixing");
    ModuleBase::timer::tick("Charge_Mixing", "init_mixing");

    // (re)construct mixing object
    if (this->mixing_mode == "broyden")
    {
        delete this->mixing;
        this->mixing = new Base_Mixing::Broyden_Mixing(this->mixing_ndim, this->mixing_beta);
    }
    else if (this->mixing_mode == "plain")
    {
        delete this->mixing;
        this->mixing = new Base_Mixing::Plain_Mixing(this->mixing_beta);
    }
    else if (this->mixing_mode == "pulay")
    {
        delete this->mixing;
        this->mixing = new Base_Mixing::Pulay_Mixing(this->mixing_ndim, this->mixing_beta);
    }
    else
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing", "This Mixing mode is not implemended yet,coming soon.");
    }

    if (GlobalV::double_grid)
    {
        // ONLY smooth part of charge density is mixed by specific mixing method
        // The high_frequency part is mixed by plain mixing method.
        delete this->mixing_highf;
        this->mixing_highf = new Base_Mixing::Plain_Mixing(this->mixing_beta);
    }

    // allocate memory for mixing data, if exists, free it first and then allocate new memory
    // initailize rho_mdata
    if (GlobalV::SCF_THR_TYPE == 1)
    {  
        if (GlobalV::NSPIN == 4 && GlobalV::MIXING_ANGLE > 0 )
        {
            this->mixing->init_mixing_data(this->rho_mdata,
                                        this->rhopw->npw * 2,
                                        sizeof(std::complex<double>));
        }
        else
        {
            this->mixing->init_mixing_data(this->rho_mdata,
                                        this->rhopw->npw * GlobalV::NSPIN,
                                        sizeof(std::complex<double>));
        }
    }
    else
    {
        if (GlobalV::NSPIN == 4 && GlobalV::MIXING_ANGLE > 0 )
        {
            this->mixing->init_mixing_data(this->rho_mdata, this->rhopw->nrxx * 2, sizeof(double));
        }
        else
        {
            this->mixing->init_mixing_data(this->rho_mdata, this->rhopw->nrxx * GlobalV::NSPIN, sizeof(double));
        }
    }
    
    // initailize tau_mdata
    if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
    {
        if (GlobalV::SCF_THR_TYPE == 1)
        {
            this->mixing->init_mixing_data(this->tau_mdata,
                                           this->rhopw->npw * GlobalV::NSPIN,
                                           sizeof(std::complex<double>));
        }
        else
        {
            this->mixing->init_mixing_data(this->tau_mdata, this->rhopw->nrxx * GlobalV::NSPIN, sizeof(double));
        }
    }

    // initailize nhat_mdata
#ifdef USE_PAW
    if(GlobalV::use_paw) this->mixing->init_mixing_data(this->nhat_mdata, this->rhopw->nrxx * GlobalV::NSPIN, sizeof(double));
#endif

    ModuleBase::timer::tick("Charge_Mixing", "init_mixing");

    return;
}

void Charge_Mixing::allocate_mixing_dmr(int nnr)
{
    // Note that: we cannot allocate memory for dmr_mdata in set_mixing.
    // since the size of dmr_mdata is given by the size of HContainer.nnr, which is calculated in DensityMatrix::init_DMR().
    // and DensityMatrix::init_DMR() is called in beforescf(). While set_mixing() is called in ESolver_KS::Init().
    ModuleBase::TITLE("Charge_Mixing", "allocate_mixing_dmr");
    ModuleBase::timer::tick("Charge_Mixing", "allocate_mixing_dmr");
    //
    const int dmr_nspin = (GlobalV::NSPIN == 2) ? 2 : 1;
    // allocate memory for dmr_mdata
    if (GlobalV::SCF_THR_TYPE == 1)
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing", "This Mixing of Density Matrix is not supported for PW basis yet");
    }
    else if (GlobalV::SCF_THR_TYPE == 2)
    {
        this->mixing->init_mixing_data(this->dmr_mdata, nnr * dmr_nspin, sizeof(double));
    }

    this->dmr_mdata.reset();
    ModuleBase::timer::tick("Charge_Mixing", "allocate_mixing_dmr");

    return;
}

void Charge_Mixing::set_rhopw(ModulePW::PW_Basis* rhopw_in, ModulePW::PW_Basis* rhodpw_in)
{
    this->rhopw = rhopw_in;
    this->rhodpw = rhodpw_in;
}

double Charge_Mixing::get_drho(Charge* chr, const double nelec)
{
    ModuleBase::TITLE("Charge_Mixing", "get_drho");
    ModuleBase::timer::tick("Charge_Mixing", "get_drho");
    double drho = 0.0;

    if (GlobalV::SCF_THR_TYPE == 1)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            ModuleBase::GlobalFunc::NOTE("Perform FFT on rho(r) to obtain rho(G).");
            chr->rhopw->real2recip(chr->rho[is], chr->rhog[is]);

            ModuleBase::GlobalFunc::NOTE("Perform FFT on rho_save(r) to obtain rho_save(G).");
            chr->rhopw->real2recip(chr->rho_save[is], chr->rhog_save[is]);
        }

        ModuleBase::GlobalFunc::NOTE("Calculate the charge difference between rho(G) and rho_save(G)");
        std::vector<std::complex<double>> drhog(GlobalV::NSPIN * this->rhopw->npw);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 512)
#endif
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            for (int ig = 0; ig < this->rhopw->npw; ig++)
            {
                drhog[is * rhopw->npw + ig] = chr->rhog[is][ig] - chr->rhog_save[is][ig];
            }
        }

        ModuleBase::GlobalFunc::NOTE("Calculate the norm of the Residual std::vector: < R[rho] | R[rho_save] >");
        drho = this->inner_product_recip(drhog.data(), drhog.data());
    }
    else
    {
        // Note: Maybe it is wrong.
        //       The inner_product_real function (L1-norm) is different from that (L2-norm) in mixing.
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            if (is != 0 && is != 3 && GlobalV::DOMAG_Z)
            {
                continue;
            }
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : drho)
#endif
            for (int ir = 0; ir < this->rhopw->nrxx; ir++)
            {
                drho += std::abs(chr->rho[is][ir] - chr->rho_save[is][ir]);
            }
        }
#ifdef __MPI
        Parallel_Reduce::reduce_pool(drho);
#endif
        assert(nelec != 0);
        assert(GlobalC::ucell.omega > 0);
        assert(this->rhopw->nxyz > 0);
        drho *= GlobalC::ucell.omega / static_cast<double>(this->rhopw->nxyz);
        drho /= nelec;
    }

    ModuleBase::timer::tick("Charge_Mixing", "get_drho");
    return drho;
}

void Charge_Mixing::mix_rho_recip_new(Charge* chr)
{
    std::complex<double>* rhog_in = nullptr;
    std::complex<double>* rhog_out = nullptr;
    // for smooth part
    std::complex<double>* rhogs_in = chr->rhog_save[0];
    std::complex<double>* rhogs_out = chr->rhog[0];
    // for high_frequency part
    std::complex<double>* rhoghf_in = nullptr;
    std::complex<double>* rhoghf_out = nullptr;

    if (GlobalV::double_grid)
    {
        // divide into smooth part and high_frequency part
        divide_data(chr->rhog_save[0], rhogs_in, rhoghf_in);
        divide_data(chr->rhog[0], rhogs_out, rhoghf_out);
    }

    //  inner_product_recip_hartree is a hartree-like sum, unit is Ry
    auto inner_product
        = std::bind(&Charge_Mixing::inner_product_recip_hartree, this, std::placeholders::_1, std::placeholders::_2);

    // DIIS Mixing Only for smooth part, while high_frequency part is mixed by plain mixing method.
    if (GlobalV::NSPIN == 1)
    {
        rhog_in = rhogs_in;
        rhog_out = rhogs_out;
        auto screen = std::bind(&Charge_Mixing::Kerker_screen_recip_new, this, std::placeholders::_1);
        this->mixing->push_data(this->rho_mdata, rhog_in, rhog_out, screen, true);
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhog_out);
    }
    else if (GlobalV::NSPIN == 2)
    {
        // magnetic density
        std::complex<double> *rhog_mag = nullptr;
        std::complex<double> *rhog_mag_save = nullptr;
        const int npw = this->rhopw->npw;
        // allocate rhog_mag[is*ngmc] and rhog_mag_save[is*ngmc]
        rhog_mag = new std::complex<double>[npw * GlobalV::NSPIN];
        rhog_mag_save = new std::complex<double>[npw * GlobalV::NSPIN];
        ModuleBase::GlobalFunc::ZEROS(rhog_mag, npw * GlobalV::NSPIN);
        ModuleBase::GlobalFunc::ZEROS(rhog_mag_save, npw * GlobalV::NSPIN);
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
        auto screen = std::bind(&Charge_Mixing::Kerker_screen_recip_new, this, std::placeholders::_1);
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
        for (int is = 0; is < GlobalV::NSPIN; is++)
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
        if (GlobalV::double_grid)
        {
            for (int ig = 0; ig < npw; ig++)
            {
                rhogs_out[ig] = chr->rhog[0][ig];
                rhogs_out[ig + npw] = chr->rhog[1][ig];
            }
        }
    }
    else if (GlobalV::NSPIN == 4 && GlobalV::MIXING_ANGLE <= 0)
    {
        // normal broyden mixing for {rho, mx, my, mz}
        rhog_in = rhogs_in;
        rhog_out = rhogs_out;
        const int npw = this->rhopw->npw;
        auto screen = std::bind(&Charge_Mixing::Kerker_screen_recip_new, this, std::placeholders::_1); // use old one
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
    else if (GlobalV::NSPIN == 4 && GlobalV::MIXING_ANGLE > 0)
    {
        // special broyden mixing for {rho, |m|} proposed by J. Phys. Soc. Jpn. 82 (2013) 114706
        // here only consider the case of mixing_angle = 1, which mean only change |m| and keep angle fixed
        // old support see mix_rho_recip()
        if (GlobalV::double_grid)
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
        auto screen = std::bind(&Charge_Mixing::Kerker_screen_recip_new, this, std::placeholders::_1); // use old one
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
            if (abs(norm) < 1e-10) continue;
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

    if (GlobalV::double_grid)
    {
        // plain mixing for high_frequencies
        const int ndimhf = (this->rhodpw->npw - this->rhopw->npw) * GlobalV::NSPIN;
        this->mixing_highf->plain_mix(rhoghf_out, rhoghf_in, rhoghf_out, ndimhf, nullptr);

        // combine smooth part and high_frequency part
        combine_data(chr->rhog[0], rhogs_out, rhoghf_out);
        clean_data(rhogs_in, rhoghf_in);
    }

    // rhog to rho
    if (GlobalV::NSPIN == 4 && GlobalV::MIXING_ANGLE > 0)
    {
        // only tranfer rhog[0]
        // do not support double_grid, use rhopw directly
        chr->rhopw->recip2real(chr->rhog[0], chr->rho[0]);
    }
    else
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            // use rhodpw for double_grid
            // rhodpw is the same as rhopw for !GlobalV::double_grid
            this->rhodpw->recip2real(chr->rhog[is], chr->rho[is]);
        }
    }

    // For kinetic energy density
    if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
    {
        std::vector<std::complex<double>> kin_g(GlobalV::NSPIN * rhodpw->npw);
        std::vector<std::complex<double>> kin_g_save(GlobalV::NSPIN * rhodpw->npw);
        // FFT to get kin_g and kin_g_save
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            rhodpw->real2recip(chr->kin_r[is], &kin_g[is * rhodpw->npw]);
            rhodpw->real2recip(chr->kin_r_save[is], &kin_g_save[is * rhodpw->npw]);
        }
        // for smooth part, for !GlobalV::double_grid only have this part
        std::complex<double>*taugs_in = kin_g_save.data(), *taugs_out = kin_g.data();
        // for high frequency part
        std::complex<double>*taughf_in = nullptr, *taughf_out = nullptr;
        if (GlobalV::double_grid)
        {
            // divide into smooth part and high_frequency part
            divide_data(kin_g_save.data(), taugs_in, taughf_in);
            divide_data(kin_g.data(), taugs_out, taughf_out);
        }

        // Note: there is no kerker modification for tau because I'm not sure
        // if we should have it. If necessary we can try it in the future.
        this->mixing->push_data(this->tau_mdata, taugs_in, taugs_out, nullptr, false);

        this->mixing->mix_data(this->tau_mdata, taugs_out);

        if (GlobalV::double_grid)
        {
            // simple mixing for high_frequencies
            const int ndimhf = (this->rhodpw->npw - this->rhopw->npw) * GlobalV::NSPIN;
            this->mixing_highf->plain_mix(taughf_out, taughf_in, taughf_out, ndimhf, nullptr);

            // combine smooth part and high_frequency part
            combine_data(kin_g.data(), taugs_out, taughf_out);
            clean_data(taugs_in, taughf_in);
        }

        // kin_g to kin_r
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            rhodpw->recip2real(&kin_g[is * rhodpw->npw], chr->kin_r[is]);
        }
    }

#ifdef USE_PAW
    if(GlobalV::use_paw)
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
    if (GlobalV::NSPIN == 1)
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
    else if (GlobalV::NSPIN == 2)
    {
        // magnetic density
        double *rho_mag = nullptr;
        double *rho_mag_save = nullptr; 
        const int nrxx = this->rhopw->nrxx;
        // allocate rho_mag[is*nnrx] and rho_mag_save[is*nnrx]
        rho_mag = new double[nrxx * GlobalV::NSPIN];
        rho_mag_save = new double[nrxx * GlobalV::NSPIN];
        ModuleBase::GlobalFunc::ZEROS(rho_mag, nrxx * GlobalV::NSPIN);
        ModuleBase::GlobalFunc::ZEROS(rho_mag_save, nrxx * GlobalV::NSPIN);
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
        for (int is = 0; is < GlobalV::NSPIN; is++)
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
    else if (GlobalV::NSPIN == 4 && GlobalV::MIXING_ANGLE <= 0)
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
    else if (GlobalV::NSPIN == 4 && GlobalV::MIXING_ANGLE > 0)
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
            if (norm < 1e-10) continue;
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

void Charge_Mixing::mix_dmr(elecstate::DensityMatrix<double, double>* DM)
{
    // Notice that DensityMatrix object is a Template class
    ModuleBase::TITLE("Charge_Mixing", "mix_dmr");
    ModuleBase::timer::tick("Charge_Mixing", "mix_dmr");
    //
    std::vector<hamilt::HContainer<double>*> dmr = DM->get_DMR_vector();
    std::vector<std::vector<double>> dmr_save = DM->get_DMR_save();
    //
    //const int dmr_nspin = (GlobalV::NSPIN == 2) ? 2 : 1;
    double* dmr_in;
    double* dmr_out;
    if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 4)
    {
        dmr_in = dmr_save[0].data();
        dmr_out = dmr[0]->get_wrapper();
        this->mixing->push_data(this->dmr_mdata, dmr_in, dmr_out, nullptr, false);    
        this->mixing->mix_data(this->dmr_mdata, dmr_out);
    }
    else if (GlobalV::NSPIN == 2)
    {
        // magnetic density matrix
        double* dmr_mag = nullptr;
        double* dmr_mag_save = nullptr; 
        const int nnr = dmr[0]->get_nnr();
        // allocate dmr_mag[is*nnrx] and dmr_mag_save[is*nnrx]
        dmr_mag = new double[nnr * GlobalV::NSPIN];
        dmr_mag_save = new double[nnr * GlobalV::NSPIN];
        ModuleBase::GlobalFunc::ZEROS(dmr_mag, nnr * GlobalV::NSPIN);
        ModuleBase::GlobalFunc::ZEROS(dmr_mag_save, nnr * GlobalV::NSPIN);
        double* dmr_up;
        double* dmr_down;
        // tranfer dmr into dmr_mag
        dmr_up = dmr[0]->get_wrapper();
        dmr_down = dmr[1]->get_wrapper();
        for (int ir = 0; ir < nnr; ir++)
        {
            dmr_mag[ir] = dmr_up[ir] + dmr_down[ir];
            dmr_mag[ir + nnr] = dmr_up[ir] - dmr_down[ir];
        }
        // tranfer dmr_save into dmr_mag_save
        dmr_up = dmr_save[0].data();
        dmr_down = dmr_save[1].data();
        for (int ir = 0; ir < nnr; ir++)
        {
            dmr_mag_save[ir] = dmr_up[ir] + dmr_down[ir];
            dmr_mag_save[ir + nnr] = dmr_up[ir] - dmr_down[ir];
        }
        //
        dmr_in = dmr_mag_save;
        dmr_out = dmr_mag;
        // no kerker in mixing_dmr
        //auto screen = std::bind(&Charge_Mixing::Kerker_screen_real, this, std::placeholders::_1);
        auto twobeta_mix
            = [this, nnr](double* out, const double* in, const double* sres) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
            for (int i = 0; i < nnr; ++i)
            {
                out[i] = in[i] + this->mixing_beta * sres[i];
            }
            // magnetism
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
            for (int i = nnr; i < 2 * nnr; ++i)
            {
                out[i] = in[i] + this->mixing_beta_mag * sres[i];
            }
        };
        this->mixing->push_data(this->dmr_mdata, dmr_in, dmr_out, nullptr, twobeta_mix, false);
        //auto inner_product
        //    = std::bind(&Charge_Mixing::inner_product_real, this, std::placeholders::_1, std::placeholders::_2);
        //this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->dmr_mdata, dmr_out);
        // get new dmr from dmr_mag
        dmr_up = dmr[0]->get_wrapper();
        dmr_down = dmr[1]->get_wrapper();
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(dmr_up, nnr);
            ModuleBase::GlobalFunc::ZEROS(dmr_down, nnr);
        }
        for (int ir = 0; ir < nnr; ir++)
        {
            dmr_up[ir] = 0.5 * (dmr_mag[ir] + dmr_mag[ir+nnr]);
            dmr_down[ir] = 0.5 * (dmr_mag[ir] - dmr_mag[ir+nnr]);
        }
        // delete
        delete[] dmr_mag;
        delete[] dmr_mag_save;
    }

    ModuleBase::timer::tick("Charge_Mixing", "mix_dmr");

    return;
}

void Charge_Mixing::mix_dmr(elecstate::DensityMatrix<std::complex<double>, double>* DM)
{
    // Notice that DensityMatrix object is a Template class
    ModuleBase::TITLE("Charge_Mixing", "mix_dmr");
    ModuleBase::timer::tick("Charge_Mixing", "mix_dmr");
    //
    std::vector<hamilt::HContainer<double>*> dmr = DM->get_DMR_vector();
    std::vector<std::vector<double>> dmr_save = DM->get_DMR_save();
    //
    //const int dmr_nspin = (GlobalV::NSPIN == 2) ? 2 : 1;
    double* dmr_in;
    double* dmr_out;
    if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 4)
    {
        dmr_in = dmr_save[0].data();
        dmr_out = dmr[0]->get_wrapper();
        this->mixing->push_data(this->dmr_mdata, dmr_in, dmr_out, nullptr, false);    
        this->mixing->mix_data(this->dmr_mdata, dmr_out);
    }
    else if (GlobalV::NSPIN == 2)
    {
        // magnetic density matrix
        double* dmr_mag = nullptr;
        double* dmr_mag_save = nullptr; 
        const int nnr = dmr[0]->get_nnr();
        // allocate dmr_mag[is*nnrx] and dmr_mag_save[is*nnrx]
        dmr_mag = new double[nnr * GlobalV::NSPIN];
        dmr_mag_save = new double[nnr * GlobalV::NSPIN];
        ModuleBase::GlobalFunc::ZEROS(dmr_mag, nnr * GlobalV::NSPIN);
        ModuleBase::GlobalFunc::ZEROS(dmr_mag_save, nnr * GlobalV::NSPIN);
        double* dmr_up;
        double* dmr_down;
        // tranfer dmr into dmr_mag
        dmr_up = dmr[0]->get_wrapper();
        dmr_down = dmr[1]->get_wrapper();
        for (int ir = 0; ir < nnr; ir++)
        {
            dmr_mag[ir] = dmr_up[ir] + dmr_down[ir];
            dmr_mag[ir + nnr] = dmr_up[ir] - dmr_down[ir];
        }
        // tranfer dmr_save into dmr_mag_save
        dmr_up = dmr_save[0].data();
        dmr_down = dmr_save[1].data();
        for (int ir = 0; ir < nnr; ir++)
        {
            dmr_mag_save[ir] = dmr_up[ir] + dmr_down[ir];
            dmr_mag_save[ir + nnr] = dmr_up[ir] - dmr_down[ir];
        }
        //
        dmr_in = dmr_mag_save;
        dmr_out = dmr_mag;
        // no kerker in mixing_dmr
        //auto screen = std::bind(&Charge_Mixing::Kerker_screen_real, this, std::placeholders::_1);
        auto twobeta_mix
            = [this, nnr](double* out, const double* in, const double* sres) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
            for (int i = 0; i < nnr; ++i)
            {
                out[i] = in[i] + this->mixing_beta * sres[i];
            }
            // magnetism
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
            for (int i = nnr; i < 2 * nnr; ++i)
            {
                out[i] = in[i] + this->mixing_beta_mag * sres[i];
            }
        };
        this->mixing->push_data(this->dmr_mdata, dmr_in, dmr_out, nullptr, twobeta_mix, false);
        //auto inner_product
        //    = std::bind(&Charge_Mixing::inner_product_real, this, std::placeholders::_1, std::placeholders::_2);
        //this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->dmr_mdata, dmr_out);
        // get new dmr from dmr_mag
        dmr_up = dmr[0]->get_wrapper();
        dmr_down = dmr[1]->get_wrapper();
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(dmr_up, nnr);
            ModuleBase::GlobalFunc::ZEROS(dmr_down, nnr);
        }
        for (int ir = 0; ir < nnr; ir++)
        {
            dmr_up[ir] = 0.5 * (dmr_mag[ir] + dmr_mag[ir+nnr]);
            dmr_down[ir] = 0.5 * (dmr_mag[ir] - dmr_mag[ir+nnr]);
        }
        // delete
        delete[] dmr_mag;
        delete[] dmr_mag_save;
    }

    ModuleBase::timer::tick("Charge_Mixing", "mix_dmr");

    return;
}

void Charge_Mixing::mix_reset()
{
    this->mixing->reset();
    this->rho_mdata.reset();
    // initailize tau_mdata
    if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
    {
        this->tau_mdata.reset();
    }
    // reset for paw
#ifdef USE_PAW
    this->nhat_mdata.reset();
#endif
}

void Charge_Mixing::mix_rho(Charge* chr)
{
    ModuleBase::TITLE("Charge_Mixing", "mix_rho");
    ModuleBase::timer::tick("Charge", "mix_rho");

    // the charge before mixing.
    const int nrxx = chr->rhopw->nrxx;
    std::vector<double> rho123(GlobalV::NSPIN * nrxx);
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        if (is == 0 || is == 3 || !GlobalV::DOMAG_Z)
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
        kin_r123.resize(GlobalV::NSPIN * nrxx);
        for (int is = 0; is < GlobalV::NSPIN; ++is)
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
    if(GlobalV::use_paw)
    {
        nhat_r123.resize(GlobalV::NSPIN * nrxx);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
        for(int ir = 0 ; ir < nrxx ; ++ir)
        {
            for(int is = 0; is < GlobalV::NSPIN; ++is)
            {
                nhat_r123[ir+is*nrxx] = chr->nhat[0][ir];
            }
        }
    }        
#endif
    // --------------------Mixing Body--------------------
    if (GlobalV::SCF_THR_TYPE == 1)
    {
        mix_rho_recip_new(chr);
    }
    else if (GlobalV::SCF_THR_TYPE == 2)
    {
        mix_rho_real(chr);
    }
    // ---------------------------------------------------

    // mohan add 2012-06-05
    // rho_save is the charge before mixing
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        if (is == 0 || is == 3 || !GlobalV::DOMAG_Z)
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
        for (int is = 0; is < GlobalV::NSPIN; ++is)
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
    if(GlobalV::use_paw)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
        for(int ir = 0 ; ir < nrxx ; ++ir)
        {
            for(int is = 0; is < GlobalV::NSPIN; ++is)
            {
                chr->nhat_save[is][ir] = nhat_r123[ir+is*nrxx];
            }
        }
    }
#endif

    if (new_e_iteration)
        new_e_iteration = false;
    ModuleBase::timer::tick("Charge", "mix_rho");
    return;
}

void Charge_Mixing::Kerker_screen_recip(std::complex<double>* drhog)
{
    if (this->mixing_gg0 <= 0.0 || this->mixing_beta <= 0.1)
        return;
    const double fac = this->mixing_gg0;
    const double gg0 = std::pow(fac * 0.529177 / GlobalC::ucell.tpiba, 2);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 512)
#endif
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        for (int ig = 0; ig < this->rhopw->npw; ++ig)
        {
            double gg = this->rhopw->gg[ig];
            double filter_g = std::max(gg / (gg + gg0), GlobalV::MIXING_GG0_MIN / this->mixing_beta);
            drhog[is * this->rhopw->npw + ig] *= filter_g;
        }
    }
    return;
}

void Charge_Mixing::Kerker_screen_recip_new(std::complex<double>* drhog)
{
    if (this->mixing_gg0 <= 0.0 || this->mixing_beta <= 0.1)
        return;
    double fac, gg0, amin;

    // consider a resize for mixing_angle
    int resize_tmp = 1;
    if (GlobalV::NSPIN == 4 && GlobalV::MIXING_ANGLE > 0) resize_tmp = 2;

    // implement Kerker for density and magnetization separately
    for (int is = 0; is < GlobalV::NSPIN / resize_tmp; ++is)
    {
        // new mixing method only support nspin=2 not nspin=4
        if (is >= 1)
        {
            if (GlobalV::MIXING_GG0_MAG <= 0.0001 || GlobalV::MIXING_BETA_MAG <= 0.1)
            {
#ifdef __DEBUG
                assert(is == 1); // make sure break works
#endif
                double is_mag = GlobalV::NSPIN - 1;
                //for (int ig = 0; ig < this->rhopw->npw * is_mag; ig++)
                //{
                //    drhog[is * this->rhopw->npw + ig] *= 1;
                //}
                break;
            }
            fac = GlobalV::MIXING_GG0_MAG;
            amin = GlobalV::MIXING_BETA_MAG;
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
            double filter_g = std::max(gg / (gg + gg0), GlobalV::MIXING_GG0_MIN / amin);
            drhog[is * this->rhopw->npw + ig] *= filter_g;
        }
    }
    return;
}

void Charge_Mixing::Kerker_screen_real(double* drhor)
{
    if (this->mixing_gg0 <= 0.0001 || this->mixing_beta <= 0.1)
        return;
    // consider a resize for mixing_angle
    int resize_tmp = 1;
    if (GlobalV::NSPIN == 4 && GlobalV::MIXING_ANGLE > 0) resize_tmp = 2;
    //
    std::vector<std::complex<double>> drhog(this->rhopw->npw * GlobalV::NSPIN / resize_tmp);
    std::vector<double> drhor_filter(this->rhopw->nrxx * GlobalV::NSPIN / resize_tmp);
    for (int is = 0; is < GlobalV::NSPIN / resize_tmp; ++is)
    {
        // Note after this process some G which is higher than Gmax will be filtered.
        // Thus we cannot use Kerker_screen_recip(drhog.data()) directly after it.
        this->rhopw->real2recip(drhor + is * this->rhopw->nrxx, drhog.data() + is * this->rhopw->npw);
    }
    // implement Kerker for density and magnetization separately
    double fac, gg0, amin;
    for (int is = 0; is < GlobalV::NSPIN / resize_tmp; is++)
    {

        if (is >= 1)
        {
            if (GlobalV::MIXING_GG0_MAG <= 0.0001 || GlobalV::MIXING_BETA_MAG <= 0.1)
            {
#ifdef __DEBUG
                assert(is == 1); // make sure break works
#endif
                double is_mag = GlobalV::NSPIN - 1;
                if (GlobalV::NSPIN == 4 && GlobalV::MIXING_ANGLE > 0) is_mag = 1;
                for (int ig = 0; ig < this->rhopw->npw * is_mag; ig++)
                {
                    drhog[is * this->rhopw->npw + ig] = 0;
                }
                break;
            }
            fac = GlobalV::MIXING_GG0_MAG;
            amin = GlobalV::MIXING_BETA_MAG;
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
            double filter_g = std::max(gg / (gg + gg0), GlobalV::MIXING_GG0_MIN / amin);
            drhog[is * this->rhopw->npw + ig] *= (1 - filter_g);
        }
    }
    // inverse FT
    for (int is = 0; is < GlobalV::NSPIN / resize_tmp; ++is)
    {
        this->rhopw->recip2real(drhog.data() + is * this->rhopw->npw, drhor_filter.data() + is * this->rhopw->nrxx);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
    for (int ir = 0; ir < this->rhopw->nrxx * GlobalV::NSPIN / resize_tmp; ir++)
    {
        drhor[ir] -= drhor_filter[ir];
    }
}

double Charge_Mixing::inner_product_recip(std::complex<double>* rho1, std::complex<double>* rho2)
{
    std::complex<double>** rho1_2d = new std::complex<double>*[GlobalV::NSPIN];
    std::complex<double>** rho2_2d = new std::complex<double>*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        rho1_2d[is] = rho1 + is * this->rhopw->npw;
        rho2_2d[is] = rho2 + is * this->rhopw->npw;
    }
    double result = this->rhog_dot_product(rho1_2d, rho2_2d);
    delete[] rho1_2d;
    delete[] rho2_2d;
    return result;
}

// a simple inner product, now is not used anywhere. For test only.
double Charge_Mixing::inner_product_recip_simple(std::complex<double>* rho1, std::complex<double>* rho2)
{
    ModuleBase::TITLE("Charge_Mixing", "inner_product_recip_simple");
    ModuleBase::timer::tick("Charge_Mixing", "inner_product_recip_simple");

    double rnorm = 0.0;
    // consider a resize for mixing_angle
    int resize_tmp = 1;
    if (GlobalV::NSPIN == 4 && GlobalV::MIXING_ANGLE > 0) resize_tmp = 2;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : rnorm)
#endif
    for (int ig = 0; ig < this->rhopw->npw * GlobalV::NSPIN / resize_tmp; ++ig)
    {
        rnorm += (conj(rho1[ig]) * rho2[ig]).real();
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(rnorm);
#endif

    ModuleBase::timer::tick("Charge_Mixing", "inner_product_recip_simple");

    return rnorm;
}

// a Hartree-like inner product
double Charge_Mixing::inner_product_recip_hartree(std::complex<double>* rhog1, std::complex<double>* rhog2)
{
    ModuleBase::TITLE("Charge_Mixing", "inner_product_recip_hartree");
    ModuleBase::timer::tick("Charge_Mixing", "inner_product_recip_hartree");

    static const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / GlobalC::ucell.tpiba2;
    static const double fac2 = ModuleBase::e2 * ModuleBase::FOUR_PI / (ModuleBase::TWO_PI * ModuleBase::TWO_PI);

    double sum = 0.0;
    const int npw = this->rhopw->npw;

    // a lambda function for summing the charge density
    auto part_of_rho = [&]()
    {
        double sum = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
        for (int ig = 0; ig < this->rhopw->npw; ++ig)
        {
            if (this->rhopw->gg[ig] < 1e-8)
                continue;
            sum += (conj(rhog1[ig]) * rhog2[ig]).real() / this->rhopw->gg[ig];
        }
        sum *= fac;
        return sum;
    };
    
    if (GlobalV::NSPIN==1)
    {
        sum += part_of_rho();
    }
    else if (GlobalV::NSPIN==2)
    {
        // charge density part
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
        for (int ig = 0; ig < this->rhopw->npw; ++ig)
        {
            if (this->rhopw->gg[ig] < 1e-8)
                continue;
            sum += (conj(rhog1[ig]) * (rhog2[ig])).real() / this->rhopw->gg[ig];
        }
        sum *= fac;

        if (GlobalV::GAMMA_ONLY_PW)
        {
            sum *= 2.0;
        }

        // (2) Second part of density error.
        // including |G|=0 term.
        double sum2 = 0.0;

        sum2 += fac2 * (conj(rhog1[0 + this->rhopw->npw]) * rhog2[0 + this->rhopw->npw]).real();

        double mag = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : mag)
#endif
        for (int ig = 0; ig < this->rhopw->npw; ig++)
        {
            mag += (conj(rhog1[ig + this->rhopw->npw]) * rhog2[ig + this->rhopw->npw]).real();
        }
        mag *= fac2;

        if (GlobalV::GAMMA_ONLY_PW)
        {
            mag *= 2.0;
        }

        sum2 += mag;
        sum += sum2;
    }
    else if (GlobalV::NSPIN==4)
    {
        if (!GlobalV::DOMAG && !GlobalV::DOMAG_Z)
        {
            sum += part_of_rho();
        }
        else if (this->mixing_angle <= 0)
        {
            // sum for tradtional mixing
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
            for (int ig = 0; ig < this->rhopw->npw; ig++)
            {
                if (ig == this->rhopw->ig_gge0)
                    continue;
                sum += (conj(rhog1[ig]) * rhog2[ig]).real() / this->rhopw->gg[ig];
            }
            sum *= fac;
            const int ig0 = this->rhopw->ig_gge0;
            if (ig0 > 0)
            {
                sum += fac2
                       * ((conj(rhog1[ig0 + npw]) * rhog2[ig0 + npw]).real() + (conj(rhog1[ig0 + 2*npw]) * rhog2[ig0 + 2*npw]).real()
                          + (conj(rhog1[ig0 + 3*npw]) * rhog2[ig0 + 3*npw]).real());
            }
            double fac3 = fac2;
            if (GlobalV::GAMMA_ONLY_PW)
            {
                fac3 *= 2.0;
            }
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
            for (int ig = 0; ig < this->rhopw->npw; ig++)
            {
                if (ig == ig0)
                    continue;
                sum += fac3
                       * ((conj(rhog1[ig + npw]) * rhog2[ig + npw]).real() + (conj(rhog1[ig + 2*npw]) * rhog2[ig + 2*npw]).real()
                          + (conj(rhog1[ig + 3*npw]) * rhog2[ig + 3*npw]).real());
            }
        }
        else if (this->mixing_angle > 0)
        {
            // sum for angle mixing
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
            for (int ig = 0; ig < this->rhopw->npw; ig++)
            {
                if (ig == this->rhopw->ig_gge0)
                    continue;
                sum += (conj(rhog1[ig]) * rhog2[ig]).real() / this->rhopw->gg[ig];
            }
            sum *= fac;
            const int ig0 = this->rhopw->ig_gge0;
            if (ig0 > 0)
            {
                sum += fac2
                       * ((conj(rhog1[ig0 + this->rhopw->npw]) * rhog2[ig0 + this->rhopw->npw]).real());
            }
            double fac3 = fac2;
            if (GlobalV::GAMMA_ONLY_PW)
            {
                fac3 *= 2.0;
            }
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
            for (int ig = 0; ig < this->rhopw->npw; ig++)
            {
                if (ig == ig0)
                    continue;
                sum += fac3
                       * ((conj(rhog1[ig + this->rhopw->npw]) * rhog2[ig + this->rhopw->npw]).real());
            }
        }
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(sum);
#endif

    ModuleBase::timer::tick("Charge_Mixing", "inner_product_recip_hartree");

    sum *= GlobalC::ucell.omega * 0.5;

    return sum;
}

double Charge_Mixing::inner_product_real(double* rho1, double* rho2)
{
    double rnorm = 0.0;
    // consider a resize for mixing_angle
    int resize_tmp = 1;
    if (GlobalV::NSPIN == 4 && GlobalV::MIXING_ANGLE > 0) resize_tmp = 2;

#ifdef _OPENMP
#pragma omp parallel for reduction(+ : rnorm)
#endif
    for (int ir = 0; ir < this->rhopw->nrxx * GlobalV::NSPIN / resize_tmp; ++ir)
    {
        rnorm += rho1[ir] * rho2[ir];
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(rnorm);
#endif
    return rnorm;
}

double Charge_Mixing::rhog_dot_product(const std::complex<double>* const* const rhog1,
                                       const std::complex<double>* const* const rhog2) const
{
    ModuleBase::TITLE("Charge_Mixing", "rhog_dot_product");
    ModuleBase::timer::tick("Charge_Mixing", "rhog_dot_product");
    static const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / GlobalC::ucell.tpiba2;
    static const double fac2 = ModuleBase::e2 * ModuleBase::FOUR_PI / (ModuleBase::TWO_PI * ModuleBase::TWO_PI);

    double sum = 0.0;

    auto part_of_noncolin = [&]() // Peize Lin change goto to function at 2020.01.31
    {
        double sum = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
        for (int ig = 0; ig < this->rhopw->npw; ++ig)
        {
            if (this->rhopw->gg[ig] < 1e-8)
                continue;
            sum += (conj(rhog1[0][ig]) * rhog2[0][ig]).real() / this->rhopw->gg[ig];
        }
        sum *= fac;
        return sum;
    };

    switch (GlobalV::NSPIN)
    {
    case 1:
        sum += part_of_noncolin();
        break;

    case 2: {
        // (1) First part of density error.
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
        for (int ig = 0; ig < this->rhopw->npw; ++ig)
        {
            if (this->rhopw->gg[ig] < 1e-8)
                continue;
            sum += (conj(rhog1[0][ig] + rhog1[1][ig]) * (rhog2[0][ig] + rhog2[1][ig])).real() / this->rhopw->gg[ig];
        }
        sum *= fac;

        if (GlobalV::GAMMA_ONLY_PW)
        {
            sum *= 2.0;
        }

        // (2) Second part of density error.
        // including |G|=0 term.
        double sum2 = 0.0;

        sum2 += fac2 * (conj(rhog1[0][0] - rhog1[1][0]) * (rhog2[0][0] - rhog2[1][0])).real();

        double mag = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : mag)
#endif
        for (int ig = 0; ig < this->rhopw->npw; ig++)
        {
            mag += (conj(rhog1[0][ig] - rhog1[1][ig]) * (rhog2[0][ig] - rhog2[1][ig])).real();
        }
        mag *= fac2;

        // if(GlobalV::GAMMA_ONLY_PW);
        if (GlobalV::GAMMA_ONLY_PW) // Peize Lin delete ; 2020.01.31
        {
            mag *= 2.0;
        }

        // std::cout << " sum=" << sum << " mag=" << mag << std::endl;
        sum2 += mag;
        sum += sum2;
        break;
    }
    case 4:
        // non-collinear spin, added by zhengdy
        if (!GlobalV::DOMAG && !GlobalV::DOMAG_Z)
            sum += part_of_noncolin();
        else
        {
            // another part with magnetization
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
            for (int ig = 0; ig < this->rhopw->npw; ig++)
            {
                if (ig == this->rhopw->ig_gge0)
                    continue;
                sum += (conj(rhog1[0][ig]) * rhog2[0][ig]).real() / this->rhopw->gg[ig];
            }
            sum *= fac;
            const int ig0 = this->rhopw->ig_gge0;
            if (ig0 > 0)
            {
                sum += fac2
                       * ((conj(rhog1[1][ig0]) * rhog2[1][ig0]).real() + (conj(rhog1[2][ig0]) * rhog2[2][ig0]).real()
                          + (conj(rhog1[3][ig0]) * rhog2[3][ig0]).real());
            }
            double fac3 = fac2;
            if (GlobalV::GAMMA_ONLY_PW)
            {
                fac3 *= 2.0;
            }
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
            for (int ig = 0; ig < this->rhopw->npw; ig++)
            {
                if (ig == ig0)
                    continue;
                sum += fac3
                       * ((conj(rhog1[1][ig]) * rhog2[1][ig]).real() + (conj(rhog1[2][ig]) * rhog2[2][ig]).real()
                          + (conj(rhog1[3][ig]) * rhog2[3][ig]).real());
            }
        }
        break;
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(sum);
#endif
    ModuleBase::timer::tick("Charge_Mixing", "rhog_dot_product");

    sum *= GlobalC::ucell.omega * 0.5;

    return sum;
}

void Charge_Mixing::divide_data(std::complex<double>* data_d,
                                std::complex<double>*& data_s,
                                std::complex<double>*& data_hf)
{
    ModuleBase::TITLE("Charge_Mixing", "divide_data");
    if (GlobalV::NSPIN == 1)
    {
        data_s = data_d;
        data_hf = data_d + this->rhopw->npw;
    }
    else
    {
        const int ndimd = this->rhodpw->npw;
        const int ndims = this->rhopw->npw;
        const int ndimhf = ndimd - ndims;
        data_s = new std::complex<double>[GlobalV::NSPIN * ndims];
        data_hf = nullptr;
        if (ndimhf > 0)
        {
            data_hf = new std::complex<double>[GlobalV::NSPIN * ndimhf];
        }
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            std::memcpy(data_s + is * ndims, data_d + is * ndimd, ndims * sizeof(std::complex<double>));
            std::memcpy(data_hf + is * ndimhf, data_d + is * ndimd + ndims, ndimhf * sizeof(std::complex<double>));
        }
    }
}
void Charge_Mixing::combine_data(std::complex<double>* data_d,
                                 std::complex<double>*& data_s,
                                 std::complex<double>*& data_hf)
{
    ModuleBase::TITLE("Charge_Mixing", "combine_data");
    if (GlobalV::NSPIN == 1)
    {
        data_s = nullptr;
        data_hf = nullptr;
        return;
    }
    else
    {
        const int ndimd = this->rhodpw->npw;
        const int ndims = this->rhopw->npw;
        const int ndimhf = ndimd - ndims;
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            std::memcpy(data_d + is * ndimd, data_s + is * ndims, ndims * sizeof(std::complex<double>));
            std::memcpy(data_d + is * ndimd + ndims, data_hf + is * ndimhf, ndimhf * sizeof(std::complex<double>));
        }
        delete[] data_s;
        delete[] data_hf;
        data_s = nullptr;
        data_hf = nullptr;
    }
}

void Charge_Mixing::clean_data(std::complex<double>*& data_s, std::complex<double>*& data_hf)
{
    ModuleBase::TITLE("Charge_Mixing", "clean_data");
    if (GlobalV::NSPIN == 1)
    {
        data_s = nullptr;
        data_hf = nullptr;
        return;
    }
    else
    {
        delete[] data_s;
        delete[] data_hf;
        data_s = nullptr;
        data_hf = nullptr;
    }
}