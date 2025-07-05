#include "./kedf_xwm.h"

#include "module_parameter/parameter.h"
#include "source_base/parallel_reduce.h"
#include "source_base/tool_quit.h"

/**
 * @brief Set the parameters of XWM KEDF, and initialize kernel
 *
 * @param dV the volume of one grid point in real space, omega/nxyz
 * @param rho_ref
 * @param kappa
 * @param nelec the number of electron
 * @param tf_weight
 * @param vw_weight
 * @param pw_rho pw_basis
 */
void KEDF_XWM::set_para(double dV,
                double rho_ref,
                double kappa,
                double nelec,
                double tf_weight,
                double vw_weight,
                ModulePW::PW_Basis* pw_rho)
{
    this->dV_ = dV;
    this->kappa_ = kappa;
    this->kappa_5_6 = this->kappa_ + 5. / 6.;
    this->kappa_11_6 = this->kappa_ + 11. / 6.;
    this->kappa_1_6 = this->kappa_ - 1. / 6.;

    this->rho0_ = 1. / (pw_rho->nxyz * dV) * nelec;
    if (rho_ref != 0)
    {
        this->rho_ref_ = rho_ref;
    }
    else
    {
        this->rho_ref_ = this->rho0_;
    }
    this->c_kernel /= std::pow(this->rho0_, 2. * this->kappa_);

    this->kf_ = std::pow(3. * std::pow(ModuleBase::PI, 2) * this->rho_ref_, 1. / 3.);
    this->tkf_ = 2. * this->kf_;

    double temp1 = 1. / (this->kappa_ + 5. / 6.);
    double temp2 = 1. / (this->kappa_ + 11. / 6.);
    this->c_0 = 0.5 * temp1 * temp1;
    this->c_1 = temp1 * temp2;
    this->c_2 = - temp1 * temp1 * this->rho_ref_;

    this->kernel1_.resize(pw_rho->npw, 0.);
    this->kernel2_.resize(pw_rho->npw, 0.);

    this->fill_kernel(tf_weight, vw_weight, pw_rho);
}

/**
 * @brief Get the energy of XWM KEDF
 * 
 * @param prho charge density
 * @param pw_rho pw basis
 * @return the energy of XWM KEDF
 */
double KEDF_XWM::get_energy(const double* const* prho, ModulePW::PW_Basis* pw_rho)
{
    double** w1Rho5_6 = new double*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        w1Rho5_6[is] = new double[pw_rho->nrxx];
    }
    this->multi_kernel(prho, this->kernel1_.data(), w1Rho5_6, this->kappa_5_6, pw_rho);

    double** w2Rho5_6 = new double*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        w2Rho5_6[is] = new double[pw_rho->nrxx];
    }
    this->multi_kernel(prho, this->kernel2_.data(), w2Rho5_6, this->kappa_5_6, pw_rho);

    double energy = 0.; // in Ry
    if (PARAM.inp.nspin == 1)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            energy += std::pow(prho[0][ir], this->kappa_5_6) * w1Rho5_6[0][ir]
                    + std::pow(prho[0][ir], this->kappa_11_6) * w2Rho5_6[0][ir];
        }
        energy += this->dV_;
    }
    else if (PARAM.inp.nspin == 2)
    {
        // TODO: spin polarized
    }
    this->xwm_energy = energy;
    Parallel_Reduce::reduce_all(this->xwm_energy);

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        delete[] w1Rho5_6[is];
        delete[] w2Rho5_6[is];
    }
    delete[] w1Rho5_6;
    delete[] w2Rho5_6;

    return energy;
}

/**
 * @brief Get the energy density of XWM KEDF
 * 
 * @param prho charge density
 * @param is spin index
 * @param ir grid index
 * @param pw_rho pw basis
 * @return the energy density of XWM KEDF
 */
double KEDF_XWM::get_energy_density(const double* const* prho, int is, int ir, ModulePW::PW_Basis* pw_rho)
{
    double** w1Rho5_6 = new double*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        w1Rho5_6[is] = new double[pw_rho->nrxx];
    }
    this->multi_kernel(prho, this->kernel1_.data(), w1Rho5_6, this->kappa_5_6, pw_rho);

    double** w2Rho5_6 = new double*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        w2Rho5_6[is] = new double[pw_rho->nrxx];
    }
    this->multi_kernel(prho, this->kernel2_.data(), w2Rho5_6, this->kappa_5_6, pw_rho);

    double result = std::pow(prho[is][ir], this->kappa_5_6) * w1Rho5_6[is][ir]
                  + std::pow(prho[is][ir], this->kappa_11_6) * w2Rho5_6[is][ir];
    
    result *= this->dV_;

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        delete[] w1Rho5_6[is];
        delete[] w2Rho5_6[is];
    }
    delete[] w1Rho5_6;
    delete[] w2Rho5_6;

    return result;
}

/**
 * @brief Get the kinetic energy of XWM KEDF, and add it onto rtau_xwm
 * 
 * @param prho charge density
 * @param pw_rho pw basis
 * @param rtau_xwm rtau_xwm => rtau_xwm + tau_xwm
 */
void KEDF_XWM::tau_xwm(const double* const* prho, ModulePW::PW_Basis* pw_rho, double* rtau_xwm)
{
    double** w1Rho5_6 = new double*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        w1Rho5_6[is] = new double[pw_rho->nrxx];
    }
    this->multi_kernel(prho, this->kernel1_.data(), w1Rho5_6, this->kappa_5_6, pw_rho);

    double** w2Rho5_6 = new double*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        w2Rho5_6[is] = new double[pw_rho->nrxx];
    }
    this->multi_kernel(prho, this->kernel2_.data(), w2Rho5_6, this->kappa_5_6, pw_rho);

    if (PARAM.inp.nspin == 1)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            rtau_xwm[ir] += std::pow(prho[0][ir], this->kappa_5_6) * w1Rho5_6[0][ir]
                          + std::pow(prho[0][ir], this->kappa_11_6) * w2Rho5_6[0][ir];
        }
    }
    else if (PARAM.inp.nspin == 2)
    {
        // TODO: spin polarized
    }

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        delete[] w1Rho5_6[is];
        delete[] w2Rho5_6[is];
    }
    delete[] w1Rho5_6;
    delete[] w2Rho5_6;
}

/**
 * @brief Get the potential of XWM KEDF, and add it into rpotential,
 * and the XWM energy is stored in this->xwm_energy
 * 
 * @param prho charge density
 * @param pw_rho pw basis
 * @param rpotential rpotential => rpotential + V_{XWM}
 */
void KEDF_XWM::xwm_potential(const double* const* prho, ModulePW::PW_Basis* pw_rho, ModuleBase::matrix& rpotential)
{
    ModuleBase::timer::tick("KEDF_XWM", "xwm_potential");
    double** w1Rho5_6 = new double*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        w1Rho5_6[is] = new double[pw_rho->nrxx];
    }
    this->multi_kernel(prho, this->kernel1_.data(), w1Rho5_6, this->kappa_5_6, pw_rho);

    double** w2Rho11_6 = new double*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        w2Rho11_6[is] = new double[pw_rho->nrxx];
    }
    this->multi_kernel(prho, this->kernel2_.data(), w2Rho11_6, this->kappa_11_6, pw_rho);

    double** w2Rho5_6 = new double*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        w2Rho5_6[is] = new double[pw_rho->nrxx];
    }
    this->multi_kernel(prho, this->kernel2_.data(), w2Rho5_6, this->kappa_5_6, pw_rho);

    double energy = 0.;
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            double rho1_6 = std::pow(prho[is][ir], this->kappa_1_6);
            double rho5_6 = std::pow(prho[is][ir], this->kappa_5_6);
            double rho11_6 = std::pow(prho[is][ir], this->kappa_11_6);

            rpotential(is, ir) += 2. * this->kappa_5_6 * rho1_6 * w1Rho5_6[is][ir]
                                + this->kappa_5_6 * rho1_6 * w2Rho11_6[is][ir]
                                + this->kappa_11_6 * rho5_6 * w2Rho5_6[is][ir];

            energy += rho5_6 * w1Rho5_6[is][ir] + rho11_6 * w2Rho5_6[is][ir]; // FOR SPIN 1 !!!
        }
    }
    energy *= this->dV_;
    this->xwm_energy = energy;
    Parallel_Reduce::reduce_all(this->xwm_energy);

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        delete[] w1Rho5_6[is];
        delete[] w2Rho11_6[is];
        delete[] w2Rho5_6[is];
    }
    delete[] w1Rho5_6;
    delete[] w2Rho11_6;
    delete[] w2Rho5_6;
    ModuleBase::timer::tick("KEDF_XWM", "xwm_potential");
}

/**
 * @brief Get the stress of XWM KEDF, and store it into this->stress
 * 
 * @param prho charge density
 * @param pw_rho pw basis
 * @param vw_weight the weight of vW KEDF
 */
void KEDF_XWM::get_stress(const double* const* prho, ModulePW::PW_Basis* pw_rho, double vw_weight)
{
    std::cout << "XWM stress is not implemented yet!" << std::endl;
}

/**
 * @brief Calculate \int{W(r-r')rho^{exponent}(r') dr'}
 *
 * @param [in] prho charge density
 * @param [in] kernel W(r-r')
 * @param [out] rkernel_rho \int{W(r-r')rho^{exponent}(r') dr'}
 * @param [in] exponent the exponent of rho
 * @param [in] pw_rho pw_basis
 */
void KEDF_XWM::multi_kernel(const double* const* prho, const double* kernel, double** rkernel_rho, double exponent, ModulePW::PW_Basis* pw_rho)
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
void KEDF_XWM::fill_kernel(double tf_weight, double vw_weight, ModulePW::PW_Basis* pw_rho)
{
    double eta = 0.;
    double coef = - 1./6. / this->rho_ref_;
    for (int ig = 0; ig < pw_rho->npw; ++ig)
    {
        eta = sqrt(pw_rho->gg[ig]) * pw_rho->tpiba / this->tkf_;

        // calculate the lindhard response function and its derivative
        double lindhard = 0.;
        double diff_lindhard = 0.;
        if (eta < 0.)
        {
            lindhard = 0.;
            diff_lindhard = 0.;
        }
        // limit for small eta
        else if (eta < 1e-10)
        {
            lindhard = 1. - tf_weight + eta * eta * (1. / 3. - 3. * vw_weight);
            diff_lindhard = 2. * eta * (1. / 3. - 3. * vw_weight);
        }
        // around the singularity
        else if (std::abs(eta - 1.) < 1e-10)
        {
            lindhard = 2. - tf_weight - 3. * vw_weight + 40. * (eta - 1);
            diff_lindhard = 40.;
        }
        // Taylor expansion for high eta
        else if (eta > 3.65)
        {
            double eta2 = eta * eta;
            double invEta2 = 1. / eta2;
            lindhard = 3. * (1. - vw_weight) * eta2
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
            diff_lindhard = ((eta2 + 1.) * 0.25 / eta2 * std::log(std::abs((1. + eta) / (1. - eta))) - 0.5 / eta)
                            / std::pow((0.5 + 0.25 * (1. - eta2) * std::log((1. + eta) / std::abs(1. - eta)) / eta), 2)
                            - 6. * eta * vw_weight;
        }
        else
        {
            double eta2 = eta * eta;
            double log_term = std::log(std::abs((1. + eta) / (1. - eta)));
            double term1 = 1. / (0.5 + 0.25 * (1. - eta2) / eta * log_term);
            double term2 = 0.25 * (1. + eta2) / eta2 * log_term - 0.5 / eta;

            lindhard = term1 - 3. * vw_weight * eta2 - tf_weight;
            diff_lindhard = term1 * term1 * term2 - 6. * eta * vw_weight;
        }

        diff_lindhard = coef * eta * diff_lindhard;
        this->kernel1_[ig] = (this->c_0 * lindhard + this->c_2 * diff_lindhard) * this->c_kernel;
        this->kernel2_[ig] = this->c_1 * diff_lindhard * this->c_kernel;
    }
}