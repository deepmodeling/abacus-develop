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
}

void Charge_Mixing::set_mixing(const std::string& mixing_mode_in,
                               const double& mixing_beta_in,
                               const int& mixing_ndim_in,
                               const double& mixing_gg0_in,
                               const bool& mixing_tau_in)
{
    this->mixing_mode = mixing_mode_in;
    this->mixing_beta = mixing_beta_in;
    this->mixing_ndim = mixing_ndim_in;
    this->mixing_gg0 = mixing_gg0_in; // mohan add 2014-09-27
    this->mixing_tau = mixing_tau_in;

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
        this->mixing = new Base_Mixing::Pulay_Mixing(this->mixing_ndim, this->mixing_beta);
    }
    else
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing", "This Mixing mode is not implemended yet,coming soon.");
    }

    if (GlobalV::SCF_THR_TYPE == 1)
    {
        this->mixing->init_mixing_data(this->rho_mdata,
                                       this->rhopw->npw * GlobalV::NSPIN,
                                       sizeof(std::complex<double>));
    }
    else
    {
        this->mixing->init_mixing_data(this->rho_mdata, this->rhopw->nrxx * GlobalV::NSPIN, sizeof(double));
    }
    // Note: we can not init tau_mdata here temporarily, since set_xc_type() is after it.
    // this->mixing->init_mixing_data(this->tau_mdata, this->rhopw->nrxx * GlobalV::NSPIN, sizeof(double));
    return;
}

void Charge_Mixing::set_rhopw(ModulePW::PW_Basis* rhopw_in)
{
    this->rhopw = rhopw_in;
}

void Charge_Mixing::need_auto_set()
{
    this->autoset = true;
}

void Charge_Mixing::auto_set(const double& bandgap_in, const UnitCell& ucell_)
{
    // auto set parameters once
    if (!this->autoset)
    {
        return;
    }
    else
    {
        this->autoset = false;
    }
    GlobalV::ofs_running << "--------------AUTO-SET---------------" << std::endl;
    // 0.2 for metal and 0.7 for others
    if (bandgap_in * ModuleBase::Ry_to_eV < 1.0)
    {
        this->mixing->mixing_beta = this->mixing_beta = 0.2;
    }
    else
    {
        this->mixing->mixing_beta = this->mixing_beta = 0.7;
    }
    GlobalV::ofs_running << "      Autoset mixing_beta to " << this->mixing_beta << std::endl;

    bool has_trans_metal = false;
    // find elements of cell
    for (int it = 0; it < ucell_.ntype; it++)
    {
        if (ModuleBase::IsTransMetal.find(ucell_.atoms[it].ncpp.psd) != ModuleBase::IsTransMetal.end())
        {
            if (ModuleBase::IsTransMetal.at(ucell_.atoms[it].ncpp.psd))
            {
                has_trans_metal = true;
            }
        }
    }
    // auto set kerker mixing for trans metal system
    if (has_trans_metal)
    {
        this->mixing_gg0 = 1.5;
    }
    else
    {
        this->mixing_gg0 = 0.0;
    }
    GlobalV::ofs_running << "      Autoset mixing_gg0 to " << this->mixing_gg0 << std::endl;
    GlobalV::ofs_running << "-------------------------------------" << std::endl;
    // auto set for inhomogeneous system
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

void Charge_Mixing::mix_rho_recip(Charge* chr)
{
    // electronic density
    // rhog and rhog_save are calculated in get_drho() function

    // ONLY smooth part of charge density is mixed by specific mixing method
    // The high_frequency part is mixed by plain mixing method.
    // NOTE: chr->rhopw is dense, while this->rhopw is smooth
    std::vector<std::complex<double>> rhog_in(GlobalV::NSPIN * this->rhopw->npw);
    std::vector<std::complex<double>> rhog_out(GlobalV::NSPIN * this->rhopw->npw);
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        std::memcpy(&rhog_in[is * rhopw->npw], chr->rhog_save[is], rhopw->npw * sizeof(std::complex<double>));
        std::memcpy(&rhog_out[is * rhopw->npw], chr->rhog[is], rhopw->npw * sizeof(std::complex<double>));
    }

    auto screen = std::bind(&Charge_Mixing::Kerker_screen_recip, this, std::placeholders::_1);
    this->mixing->push_data(this->rho_mdata, rhog_in.data(), rhog_out.data(), screen, true);

    auto inner_product
        = std::bind(&Charge_Mixing::inner_product_recip, this, std::placeholders::_1, std::placeholders::_2);
    this->mixing->cal_coef(this->rho_mdata, inner_product);

    this->mixing->mix_data(this->rho_mdata, rhog_out.data());

    // simple mixing for high_frequencies
    this->high_freq_mix(chr->rhog[0], chr->rhog_save[0], GlobalV::NSPIN * chr->rhopw->npw);

    // combine smooth part and high_frequency part
    // rhog to rho
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        std::memcpy(chr->rhog[is], &rhog_out[is * rhopw->npw], rhopw->npw * sizeof(std::complex<double>));
        chr->rhopw->recip2real(chr->rhog[is], chr->rho[is]);
    }

    // For kinetic energy density
    if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
    {
        std::vector<std::complex<double>> kin_g(GlobalV::NSPIN * chr->rhopw->npw);
        std::vector<std::complex<double>> kin_g_save(GlobalV::NSPIN * chr->rhopw->npw);
        std::vector<std::complex<double>> taug_in(GlobalV::NSPIN * this->rhopw->npw);
        std::vector<std::complex<double>> taug_out(GlobalV::NSPIN * this->rhopw->npw);
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            chr->rhopw->real2recip(chr->kin_r[is], &kin_g[is * chr->rhopw->npw]);
            chr->rhopw->real2recip(chr->kin_r_save[is], &kin_g_save[is * chr->rhopw->npw]);
            // similar to rhog and rhog_save
            std::memcpy(&taug_in[is * rhopw->npw],
                        &kin_g_save[is * chr->rhopw->npw],
                        rhopw->npw * sizeof(std::complex<double>));
            std::memcpy(&taug_out[is * rhopw->npw],
                        &kin_g[is * chr->rhopw->npw],
                        rhopw->npw * sizeof(std::complex<double>));
        }
        // Note: there is no kerker modification for tau because I'm not sure
        // if we should have it. If necessary we can try it in the future.
        this->mixing->push_data(this->tau_mdata, taug_in.data(), taug_out.data(), nullptr, false);

        this->mixing->mix_data(this->tau_mdata, taug_out.data());

        // simple mixing for high_frequencies
        this->high_freq_mix(kin_g.data(), kin_g_save.data(), GlobalV::NSPIN * chr->rhopw->npw);

        // combine smooth part and high_frequency part
        // kin_g to kin_r
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            std::memcpy(&kin_g[is * chr->rhopw->npw],
                        &taug_out[is * rhopw->npw],
                        rhopw->npw * sizeof(std::complex<double>));
            chr->rhopw->recip2real(&kin_g[is * chr->rhopw->npw], chr->kin_r[is]);
        }
    }

    return;
}

void Charge_Mixing::mix_rho_real(Charge* chr)
{
    // electronic density
    double* rhor_in = chr->rho_save[0];
    double* rhor_out = chr->rho[0];

    auto screen = std::bind(&Charge_Mixing::Kerker_screen_real, this, std::placeholders::_1);
    this->mixing->push_data(this->rho_mdata, rhor_in, rhor_out, screen, true);

    auto inner_product
        = std::bind(&Charge_Mixing::inner_product_real, this, std::placeholders::_1, std::placeholders::_2);
    this->mixing->cal_coef(this->rho_mdata, inner_product);

    this->mixing->mix_data(this->rho_mdata, rhor_out);

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

void Charge_Mixing::mix_reset()
{
    this->mixing->reset();
    this->rho_mdata.reset();
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
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
        for(int ir = 0 ; ir < GlobalV::NSPIN * nrxx ; ++ir)
        {
            kin_r123[ir] = chr->kin_r[0][ir];
        }
    }

    // --------------------Mixing Body--------------------
    if (GlobalV::SCF_THR_TYPE == 1)
    {
        mix_rho_recip(chr);
    }
    else
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
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
        for(int ir = 0 ; ir < GlobalV::NSPIN * nrxx ; ++ir)
        {
            chr->kin_r_save[0][ir] = kin_r123[ir];
        }
    }

    if (new_e_iteration)
        new_e_iteration = false;
    ModuleBase::timer::tick("Charge", "mix_rho");
    return;
}

void Charge_Mixing::Kerker_screen_recip(std::complex<double>* drhog)
{
    if (this->mixing_gg0 <= 0.0)
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
            double filter_g = std::max(gg / (gg + gg0), 0.1);
            drhog[is * this->rhopw->npw + ig] *= filter_g;
        }
    }
    return;
}

void Charge_Mixing::Kerker_screen_real(double* drhor)
{
    if (this->mixing_gg0 <= 0.0)
        return;
    std::vector<std::complex<double>> drhog(this->rhopw->npw * GlobalV::NSPIN);
    std::vector<double> drhor_filter(this->rhopw->nrxx * GlobalV::NSPIN);
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        // Note after this process some G which is higher than Gmax will be filtered.
        // Thus we cannot use Kerker_screen_recip(drhog.data()) directly after it.
        this->rhopw->real2recip(drhor + is * this->rhopw->nrxx, drhog.data() + is * this->rhopw->npw);
    }
    // Kerker_screen_recip(drhog.data()); Note that we can not use it.
    // However, we should rewrite it:
    const double fac = this->mixing_gg0;
    const double gg0 = std::pow(fac * 0.529177 / GlobalC::ucell.tpiba, 2);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 512)
#endif
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        for (int ig = 0; ig < this->rhopw->npw; ig++)
        {
            double gg = this->rhopw->gg[ig];
            double filter_g = std::max(gg / (gg + gg0), 0.1);
            drhog[is * this->rhopw->npw + ig] *= (1 - filter_g);
        }
    }

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        this->rhopw->recip2real(drhog.data() + is * this->rhopw->npw, drhor_filter.data() + is * this->rhopw->nrxx);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
    for (int ir = 0; ir < this->rhopw->nrxx * GlobalV::NSPIN; ir++)
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

double Charge_Mixing::inner_product_real(double* rho1, double* rho2)
{
    double rnorm = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : rnorm)
#endif
    for (int ir = 0; ir < this->rhopw->nrxx * GlobalV::NSPIN; ++ir)
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

void Charge_Mixing::high_freq_mix(std::complex<double>* data,
                                  const std::complex<double>* data_save,
                                  const int& number) const
{
    ModuleBase::TITLE("Charge_Mixing", "high_freq_mix");
    ModuleBase::timer::tick("Charge_Mixing", "high_freq_mix");

    if (GlobalV::double_grid)
    {
        for (int ig = 0; ig < number; ig++)
        {
            data[ig] += mixing_beta * (data_save[ig] - data[ig]);
        }
    }

    ModuleBase::timer::tick("Charge_Mixing", "high_freq_mix");
}