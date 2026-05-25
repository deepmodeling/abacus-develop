#include "./kedf_extwt.h"

#include "source_io/module_parameter/parameter.h"
#include <iostream>

#include "source_base/global_variable.h"
#include "source_base/parallel_reduce.h"
#include "source_base/tool_quit.h"

namespace
{
constexpr double kDenomEps = 1e-14;
}

/**
 * @brief Set the parameters of ext-WT KEDF, and initialize kernel
 *
 * @param dV the volume of one grid point in real space, omega/nxyz
 * @param alpha
 * @param beta
 * @param nelec the number of electron
 * @param tf_weight
 * @param vw_weight
 * @param of_extwt_kappa
 * @param pw_rho pw_basis
 */
void KEDF_ExtWT::set_para(double dV,
                       double alpha,
                       double beta,
                       double nelec,
                       double tf_weight,
                       double vw_weight,
                       double of_extwt_kappa,
                       ModulePW::PW_Basis* pw_rho)
{
    this->dV_ = dV;
    // this->weightWT = weightWT;
    this->alpha_ = alpha;
    this->beta_ = beta;
    this->kappa_ = of_extwt_kappa;
    // std::cout << "kappa: " << this->kappa_ << std::endl;

    const double rho0_den = pw_rho->nxyz * dV;
    if (std::abs(rho0_den) < kDenomEps)
    {
        ModuleBase::WARNING_QUIT("KEDF_ExtWT", "Zero denominator in rho0 initialization (nxyz * dV)");
    }
    this->rho0_ = nelec / rho0_den;

    this->kf_ = std::pow(3. * std::pow(ModuleBase::PI, 2) * this->rho0_, 1. / 3.);
    this->tkf_ = 2. * this->kf_;

    const double wt_coef_den
        = 9. * this->alpha_ * this->beta_ * std::pow(this->rho0_, this->alpha_ + this->beta_ - 5. / 3.);
    if (std::abs(wt_coef_den) < kDenomEps)
    {
        ModuleBase::WARNING_QUIT("KEDF_ExtWT", "Zero denominator in wt_coef calculation (set_para)");
    }
    this->wt_coef_ = 5. / wt_coef_den;

    delete[] this->kernel_;
    this->kernel_ = new double[pw_rho->npw];
    this->fill_kernel(tf_weight, vw_weight, pw_rho);

    delete[] this->dkernel_deta_;
    this->dkernel_deta_ = new double[pw_rho->npw];

    this->update_dkernel_deta(vw_weight, pw_rho);
}

void KEDF_ExtWT::update_rho0(const double* const* prho,
                            ModulePW::PW_Basis* pw_rho)
{
    // Only for spin unpolarized case, need to be updated for spin polarized case
    int nspin = 1;

    this->sum_rho_kappa_ = 0.;
    this->sum_rho_kappa_plus_one_ = 0.;
    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            this->sum_rho_kappa_ += std::pow(prho[is][ir], this->kappa_);
            this->sum_rho_kappa_plus_one_ += std::pow(prho[is][ir], this->kappa_ + 1.);
        }
    }
    Parallel_Reduce::reduce_all(this->sum_rho_kappa_);
    Parallel_Reduce::reduce_all(this->sum_rho_kappa_plus_one_);
    this->sum_rho_kappa_ *= this->dV_;
    this->sum_rho_kappa_plus_one_ *= this->dV_;

    if (std::abs(this->sum_rho_kappa_) < kDenomEps)
    {
        ModuleBase::WARNING_QUIT("KEDF_ExtWT", "Zero denominator in rho0 update (sum_rho_kappa)");
    }
    this->rho0_ = this->sum_rho_kappa_plus_one_ / this->sum_rho_kappa_;
    // std::cout << "rho0: " << this->rho0_ << std::endl;
}

void KEDF_ExtWT::cal_kernel(double tf_weight,
                double vw_weight,
                double rho0,
                ModulePW::PW_Basis* pw_rho)
{
    this->rho0_ = rho0;

    if (this->rho0_ <= 0.0)
    {
        ModuleBase::WARNING_QUIT("KEDF_ExtWT", "rho0 must be positive in cal_kernel");
    }

    this->kf_ = std::pow(3. * std::pow(ModuleBase::PI, 2) * this->rho0_, 1. / 3.);
    this->tkf_ = 2. * this->kf_;

    const double wt_coef_den
        = 9. * this->alpha_ * this->beta_ * std::pow(this->rho0_, this->alpha_ + this->beta_ - 5. / 3.);
    if (std::abs(wt_coef_den) < kDenomEps)
    {
        ModuleBase::WARNING_QUIT("KEDF_ExtWT", "Zero denominator in wt_coef calculation (cal_kernel)");
    }
    this->wt_coef_ = 5. / wt_coef_den;

    this->fill_kernel(tf_weight, vw_weight, pw_rho);
}

/**
 * @brief Fill the dkernel_deta (this->dkernel_deta_)
 *
 * @param alpha
 * @param beta
 * @param vw_weight
 * @param pw_rho pw_basis
 */
void KEDF_ExtWT::update_dkernel_deta(const double &vw_weight, ModulePW::PW_Basis* pw_rho)
{
    if (std::abs(this->rho0_) < kDenomEps)
    {
        ModuleBase::WARNING_QUIT("KEDF_ExtWT", "Zero denominator in dkernel update (rho0)");
    }
    double eta = 0.;
    for (int ig = 0; ig < pw_rho->npw; ++ig)
    {
        eta = std::sqrt(pw_rho->gg[ig]) * pw_rho->tpiba / this->tkf_;
        this->dkernel_deta_[ig] = - ((this->alpha_ + this->beta_ - 5./3.) * this->kernel_[ig] + 1./3. * eta * this->diff_linhard(eta, vw_weight))
                                     * this->wt_coef_ / this->rho0_;
    }
}

/**
 * @brief Get the energy of ext-WT KEDF
 *
 * @param prho charge density
 * @param pw_rho pw basis
 * @return the energy of ext-WT KEDF
 */
double KEDF_ExtWT::get_energy(const double* const* prho, ModulePW::PW_Basis* pw_rho)
{
    double** kernelRhoBeta = new double*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is) {
        kernelRhoBeta[is] = new double[pw_rho->nrxx];
}
    this->multi_kernel(prho, this->kernel_, kernelRhoBeta, this->beta_, pw_rho);

    double energy = 0.; // in Ry
    if (PARAM.inp.nspin == 1)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            energy += std::pow(prho[0][ir], this->alpha_) * kernelRhoBeta[0][ir];
        }
        energy *= this->dV_ * this->c_tf_;
    }
    else if (PARAM.inp.nspin == 2)
    {
        // for (int is = 0; is < PARAM.inp.nspin; ++is)
        // {
        //     for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        //     {
        //         energy += 2 * pphi[is][ir] * LapPhi[is][ir];
        //     }
        // }
        // energy *= 0.5 * this->dV_ * 0.5;
    }
    this->extwt_energy = energy;
    Parallel_Reduce::reduce_all(this->extwt_energy);

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        delete[] kernelRhoBeta[is];
    }
    delete[] kernelRhoBeta;

    return this->extwt_energy;
}

/**
 * @brief Get the energy density of ext-WT KEDF
 *
 * @param prho charge density
 * @param is the index of spin
 * @param ir the index of real space grid
 * @param pw_rho pw basis
 * @return the energy density of ext-WT KEDF
 */
double KEDF_ExtWT::get_energy_density(const double* const* prho, int is, int ir, ModulePW::PW_Basis* pw_rho)
{
    double** kernelRhoBeta = new double*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is) {
        kernelRhoBeta[is] = new double[pw_rho->nrxx];
}
    this->multi_kernel(prho, this->kernel_, kernelRhoBeta, this->beta_, pw_rho);

    double result = this->c_tf_ * std::pow(prho[is][ir], this->alpha_) * kernelRhoBeta[is][ir];

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        delete[] kernelRhoBeta[is];
    }
    delete[] kernelRhoBeta;
    return result;
}

/**
 * @brief Get the kinetic energy of ext-WT KEDF, and add it onto rtau_extwt
 * 
 * @param prho charge density
 * @param pw_rho pw basis
 * @param rtau_extwt rtau_extwt => rtau_extwt + tau_extwt
 */
void KEDF_ExtWT::tau_extwt(const double* const* prho, ModulePW::PW_Basis* pw_rho, double* rtau_extwt)
{
    double** kernelRhoBeta = new double*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is) {
        kernelRhoBeta[is] = new double[pw_rho->nrxx];
}
    this->multi_kernel(prho, this->kernel_, kernelRhoBeta, this->beta_, pw_rho);

    if (PARAM.inp.nspin == 1)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            rtau_extwt[ir] += std::pow(prho[0][ir], this->alpha_) * kernelRhoBeta[0][ir] * this->c_tf_;
        }
    }
    else if (PARAM.inp.nspin == 2)
    {
        // Waiting for update
    }

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        delete[] kernelRhoBeta[is];
    }
    delete[] kernelRhoBeta;
}

/**
 * @brief Get the potential of ext-WT KEDF, and add it into rpotential,
 * and the ext-WT energy will be calculated and stored in this->extwt_energy
 * r)\rho^{\alpha}(r') dr'}] \f]
 *
 * @param prho charge density
 * @param pw_rho pw basis
 * @param rpotential rpotential => rpotential * 2 * phi + V_{extWT} * 2 * phi
 */
void KEDF_ExtWT::extwt_potential(const double* const* prho, ModulePW::PW_Basis* pw_rho, ModuleBase::matrix& rpotential)
{
    ModuleBase::timer::start("KEDF_ExtWT", "extwt_potential");

    // 1. WT potential
    double** kernelRhoBeta = new double*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is) {
        kernelRhoBeta[is] = new double[pw_rho->nrxx];
    }
    this->multi_kernel(prho, this->kernel_, kernelRhoBeta, this->beta_, pw_rho);

    double** kernelRhoAlpha = new double*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is) {
        kernelRhoAlpha[is] = new double[pw_rho->nrxx];
    }
    this->multi_kernel(prho, this->kernel_, kernelRhoAlpha, this->alpha_, pw_rho);

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            rpotential(is, ir) += this->c_tf_
                                  * (this->alpha_ * std::pow(prho[is][ir], this->alpha_ - 1.) * kernelRhoBeta[is][ir]
                                     + this->beta_ * std::pow(prho[is][ir], this->beta_ - 1.) * kernelRhoAlpha[is][ir])
                                     * 2. * std::sqrt(prho[is][ir]);
        }
    }

    // 2. rho0 part
    // 1) prepare variables
    std::vector<std::vector<double>> rho_alpha(PARAM.inp.nspin, std::vector<double>(pw_rho->nrxx, 0.));
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            rho_alpha[is][ir] = std::pow(prho[is][ir], this->alpha_);
        }
    }

    // 2) calculate the constants
    double coef = 0.;

    double** dkernelRhoBeta = new double*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is) {
        dkernelRhoBeta[is] = new double[pw_rho->nrxx];
    }
    this->multi_kernel(prho, this->dkernel_deta_, dkernelRhoBeta, this->beta_, pw_rho);
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            coef += rho_alpha[is][ir] * dkernelRhoBeta[is][ir];
        }
    }
    Parallel_Reduce::reduce_all(coef);
    coef *= this->dV_ * this->c_tf_;

    // 3) calculate the total potential
    if (std::abs(this->sum_rho_kappa_) < kDenomEps)
    {
        ModuleBase::WARNING_QUIT("KEDF_ExtWT", "Zero denominator in extwt_potential (sum_rho_kappa)");
    }
    int count = 0;
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            double second_term = std::pow(prho[is][ir], this->kappa_ - 1. + 0.5); // 0.5 is for sqrt(rho)
            if (second_term > 100 && count == 0)
            {
                count++;
                std::cout << "Warning: second_term is too large: " << second_term << ", coef = " << coef << std::endl;
            }
            if (second_term >= 10000)
            {
                second_term = 10000;
            }
            rpotential(is, ir) += coef / this->sum_rho_kappa_ * ((this->kappa_ + 1) * std::pow(prho[is][ir], this->kappa_ + 0.5) - this->kappa_ * second_term * this->rho0_) * 2.; // 2 is for 2 * phi
        }
    }

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        delete[] dkernelRhoBeta[is];
    }
    delete[] dkernelRhoBeta;

    // calculate energy
    double energy = 0.; // in Ry
    if (PARAM.inp.nspin == 1)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            energy += std::pow(prho[0][ir], this->alpha_) * kernelRhoBeta[0][ir];
        }
        energy *= this->dV_ * this->c_tf_;
    }
    else if (PARAM.inp.nspin == 2)
    {
        // for (int is = 0; is < PARAM.inp.nspin; ++is)
        // {
        //     for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        //     {
        //         energy += 2 * pphi[is][ir] * LapPhi[is][ir];
        //     }
        // }
        // energy *= 0.5 * this->dV_ * 0.5;
    }
    this->extwt_energy = energy;
    Parallel_Reduce::reduce_all(this->extwt_energy);

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        delete[] kernelRhoBeta[is];
        delete[] kernelRhoAlpha[is];
    }
    delete[] kernelRhoBeta;
    delete[] kernelRhoAlpha;
    ModuleBase::timer::end("KEDF_ExtWT", "extwt_potential");
}

/**
 * @brief Get the stress of ext-WT KEDF, and store it into this->stress
 *
 * @param prho charge density
 * @param pw_rho pw basis
 * @param vw_weight the weight of vW KEDF
 */
void KEDF_ExtWT::get_stress(const double* const* prho, ModulePW::PW_Basis* pw_rho, double vw_weight)
{
    ModuleBase::WARNING_QUIT("KEDF_ExtWT", "ext-WT stress is not implemented yet!");
}

/**
 * @brief Calculate the WT kernel according to Lindhard response function
 *
 * @param eta k / (2 * kF)
 * @param tf_weight
 * @param vw_weight
 * @return W(eta)
 */
double KEDF_ExtWT::extwt_kernel(double eta, double tf_weight, double vw_weight)
{
    if (eta < 0.)
    {
        return 0.;
    }
    // limit for small eta
    else if (eta < 1e-10)
    {
        return 1. - tf_weight + eta * eta * (1. / 3. - 3. * vw_weight);
    }
    // around the singularity
    else if (std::abs(eta - 1.) < 1e-10)
    {
        return 2. - tf_weight - 3. * vw_weight + 20. * (eta - 1);
    }
    // Taylor expansion for high eta
    else if (eta > 3.65)
    {
        double eta2 = eta * eta;
        double invEta2 = 1. / eta2;
        double LindG = 3. * (1. - vw_weight) * eta2
                        -tf_weight-0.6
                        + invEta2 * (-0.13714285714285712
                        + invEta2 * (-6.39999999999999875E-2
                        + invEta2 * (-3.77825602968460128E-2
                        + invEta2 * (-2.51824061652633074E-2
                        + invEta2 * (-1.80879839616166146E-2
                        + invEta2 * (-1.36715733124818332E-2
                        + invEta2 * (-1.07236045520990083E-2
                        + invEta2 * (-8.65192783339199453E-3
                        + invEta2 * (-7.1372762502456763E-3
                        + invEta2 * (-5.9945117538835746E-3
                        + invEta2 * (-5.10997527675418131E-3
                        + invEta2 * (-4.41060829979912465E-3
                        + invEta2 * (-3.84763737842981233E-3
                        + invEta2 * (-3.38745061493813488E-3
                        + invEta2 * (-3.00624946457977689E-3)))))))))))))));
        return LindG;
    }
    else
    {
        return 1. / (0.5 + 0.25 * (1. - eta * eta) / eta * std::log((1 + eta) / std::abs(1 - eta)))
               - 3. * vw_weight * eta * eta - tf_weight;
    }
}

/**
 * @brief The derivative of the WT kernel
 *
 * @param eta k / (2 * kF)
 * @param vw_weight
 * @return d W(eta)/d eta
 */
double KEDF_ExtWT::diff_linhard(double eta, double vw_weight)
{
    if (eta < 0.)
    {
        return 0.;
    }
    else if (eta < 1e-10)
    {
        return 2. * eta * (1. / 3. - 3. * vw_weight);
    }
    else if (std::abs(eta - 1.) < 1e-10)
    {
        return 40.;
    }
    else
    {
        double eta2 = eta * eta;
        return ((eta2 + 1.) * 0.25 / eta2 * std::log(std::abs((1. + eta) / (1. - eta))) - 0.5 / eta)
                   / std::pow((0.5 + 0.25 * (1. - eta2) * std::log((1. + eta) / std::abs(1. - eta)) / eta), 2)
               - 6. * eta * vw_weight;
    }
}

/**
 * @brief Calculate \int{W(r-r')rho^{exponent}(r') dr'}
 *
 * @param [in] prho charge density
 * @param [out] rkernel_rho \int{W(r-r')rho^{exponent}(r') dr'}
 * @param [in] exponent the exponent of rho
 * @param [in] pw_rho pw_basis
 */
void KEDF_ExtWT::multi_kernel(const double* const* prho, const double* kernel, double** rkernel_rho, double exponent, ModulePW::PW_Basis* pw_rho)
{
    std::complex<double>** recipkernelRho = new std::complex<double>*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        recipkernelRho[is] = new std::complex<double>[pw_rho->npw];
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            rkernel_rho[is][ir] = std::pow(prho[is][ir], exponent);
        }
        pw_rho->real2recip(rkernel_rho[is], recipkernelRho[is]);
        for (int ip = 0; ip < pw_rho->npw; ++ip)
        {
            recipkernelRho[is][ip] *= kernel[ip];
            // recipkernelRho[is][ip] *= this->kernel_[ip];
        }
        pw_rho->recip2real(recipkernelRho[is], rkernel_rho[is]);
    }

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        delete[] recipkernelRho[is];
    }
    delete[] recipkernelRho;
}

/**
 * @brief Fill the kernel (this->kernel_)
 *
 * @param tf_weight
 * @param vw_weight
 * @param pw_rho pw_basis
 */
void KEDF_ExtWT::fill_kernel(double tf_weight, double vw_weight, ModulePW::PW_Basis* pw_rho)
{
    double eta = 0.;
    for (int ig = 0; ig < pw_rho->npw; ++ig)
    {
        eta = std::sqrt(pw_rho->gg[ig]) * pw_rho->tpiba / this->tkf_;
        this->kernel_[ig] = this->extwt_kernel(eta, tf_weight, vw_weight) * this->wt_coef_;
    }
}