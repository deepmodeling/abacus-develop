#include "./kedf_wt.h"

#include <iostream>

#include "module_base/global_variable.h"
#include "module_base/parallel_reduce.h"
#include "module_base/tool_quit.h"

/**
 * @brief Set the parameters of WT KEDF, and initialize kernel
 *
 * @param dV the volume of one grid point in real space, omega/nxyz
 * @param alpha
 * @param beta
 * @param nelec the number of electron
 * @param tf_weight
 * @param vw_weight
 * @param of_wt_rho0
 * @param of_hold_rho0
 * @param read_kernel
 * @param kernel_file
 * @param pw_rho pw_basis
 */
void KEDF_WT::set_para(double dV,
                       double alpha,
                       double beta,
                       double nelec,
                       double tf_weight,
                       double vw_weight,
                       double of_wt_rho0,
                       bool of_hold_rho0,
                       bool read_kernel,
                       std::string kernel_file,
                       ModulePW::PW_Basis* pw_rho)
{
    this->dV_ = dV;
    // this->weightWT = weightWT;
    this->alpha_ = alpha;
    this->beta_ = beta;
    this->hold_rho0_ = of_hold_rho0;

    if (of_wt_rho0 != 0)
    {
        this->rho0_ = of_wt_rho0;
    }
    else
    {
        this->rho0_ = 1. / (pw_rho->nxyz * dV) * nelec;
    }

    this->kf_ = std::pow(3. * std::pow(ModuleBase::PI, 2) * this->rho0_, 1. / 3.);
    this->tkf_ = 2. * this->kf_;

    this->wt_coef_
        = 5. / (9. * this->alpha_ * this->beta_ * std::pow(this->rho0_, this->alpha_ + this->beta_ - 5. / 3.));

    delete[] this->kernel_;
    this->kernel_ = new double[pw_rho->npw];

    if (read_kernel)
        this->read_kernel(kernel_file, pw_rho);
    else
        this->fill_kernel(tf_weight, vw_weight, pw_rho);
}

/**
 * @brief Get the energy of WT KEDF
 * \f[ E_{WT} = c_{TF} * \int{\rho^\alpha * W(r - r') * \rho^\beta drdr'} \f]
 *
 * @param prho charge density
 * @param pw_rho pw basis
 * @return the energy of WT KEDF
 */
double KEDF_WT::get_energy(const double* const* prho, ModulePW::PW_Basis* pw_rho)
{
    double** kernelRhoBeta = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
        kernelRhoBeta[is] = new double[pw_rho->nrxx];
    this->multi_kernel(prho, kernelRhoBeta, this->beta_, pw_rho);

    double energy = 0.; // in Ry
    if (GlobalV::NSPIN == 1)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            energy += std::pow(prho[0][ir], this->alpha_) * kernelRhoBeta[0][ir];
        }
        energy *= this->dV_ * this->c_tf_;
    }
    else if (GlobalV::NSPIN == 2)
    {
        // for (int is = 0; is < GlobalV::NSPIN; ++is)
        // {
        //     for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        //     {
        //         energy += 2 * pphi[is][ir] * LapPhi[is][ir];
        //     }
        // }
        // energy *= 0.5 * this->dV_ * 0.5;
    }
    this->wt_energy = energy;
    Parallel_Reduce::reduce_all(this->wt_energy);

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] kernelRhoBeta[is];
    }
    delete[] kernelRhoBeta;

    return energy;
}

/**
 * @brief Get the energy density of LKT KEDF
 * \f[ \tau_{WT} = c_{TF} * \rho^\alpha * \int{W(r - r') * \rho^\beta dr'} \f]
 *
 * @param prho charge density
 * @param is the index of spin
 * @param ir the index of real space grid
 * @param pw_rho pw basis
 * @return the energy density of WT KEDF
 */
double KEDF_WT::get_energy_density(const double* const* prho, int is, int ir, ModulePW::PW_Basis* pw_rho)
{
    double** kernelRhoBeta = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
        kernelRhoBeta[is] = new double[pw_rho->nrxx];
    this->multi_kernel(prho, kernelRhoBeta, this->beta_, pw_rho);

    double result = this->c_tf_ * std::pow(prho[is][ir], this->alpha_) * kernelRhoBeta[is][ir];

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] kernelRhoBeta[is];
    }
    delete[] kernelRhoBeta;
    return result;
}

/**
 * @brief Get the potential of WT KEDF, and add it into rpotential,
 * and the WT energy will be calculated and stored in this->wt_energy
 * \f[ V_{WT} = c_{TF} * [\alpha \rho^{\alpha-1} \int{W(r - r')\rho^{\beta}(r') dr'} + \beta \rho^{\beta-1} \int{W(r' -
 * r)\rho^{\alpha}(r') dr'}] \f]
 *
 * @param prho charge density
 * @param pw_rho pw basis
 * @param rpotential rpotential => rpotential + V_{WT}
 */
void KEDF_WT::wt_potential(const double* const* prho, ModulePW::PW_Basis* pw_rho, ModuleBase::matrix& rpotential)
{
    ModuleBase::timer::tick("KEDF_WT", "wt_potential");

    double** kernelRhoBeta = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
        kernelRhoBeta[is] = new double[pw_rho->nrxx];
    this->multi_kernel(prho, kernelRhoBeta, this->beta_, pw_rho);

    double** kernelRhoAlpha = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
        kernelRhoAlpha[is] = new double[pw_rho->nrxx];
    this->multi_kernel(prho, kernelRhoAlpha, this->alpha_, pw_rho);

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            rpotential(is, ir) += this->c_tf_
                                  * (this->alpha_ * std::pow(prho[is][ir], this->alpha_ - 1.) * kernelRhoBeta[is][ir]
                                     + this->beta_ * std::pow(prho[is][ir], this->beta_ - 1.) * kernelRhoAlpha[is][ir]);
        }
    }

    // calculate energy
    double energy = 0.; // in Ry
    if (GlobalV::NSPIN == 1)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            energy += std::pow(prho[0][ir], this->alpha_) * kernelRhoBeta[0][ir];
        }
        energy *= this->dV_ * this->c_tf_;
    }
    else if (GlobalV::NSPIN == 2)
    {
        // for (int is = 0; is < GlobalV::NSPIN; ++is)
        // {
        //     for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        //     {
        //         energy += 2 * pphi[is][ir] * LapPhi[is][ir];
        //     }
        // }
        // energy *= 0.5 * this->dV_ * 0.5;
    }
    this->wt_energy = energy;
    Parallel_Reduce::reduce_all(this->wt_energy);

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] kernelRhoBeta[is];
        delete[] kernelRhoAlpha[is];
    }
    delete[] kernelRhoBeta;
    delete[] kernelRhoAlpha;
    ModuleBase::timer::tick("KEDF_WT", "wt_potential");
}

/**
 * @brief Get the stress of WT KEDF, and store it into this->stress
 *
 * @param prho charge density
 * @param pw_rho pw basis
 * @param vw_weight the weight of vW KEDF
 */
void KEDF_WT::get_stress(const double* const* prho, ModulePW::PW_Basis* pw_rho, double vw_weight)
{
    double coef = 0.;
    double mult = 0.;
    if (this->hold_rho0_)
    {
        coef = 0.;
        mult = -1. + this->alpha_ + this->beta_;
    }
    else
    {
        coef = 1. / 3.;
        mult = 2. / 3.;
    }

    std::complex<double>** recipRhoAlpha = new std::complex<double>*[GlobalV::NSPIN];
    std::complex<double>** recipRhoBeta = new std::complex<double>*[GlobalV::NSPIN];
    double* tempRho = new double[pw_rho->nrxx];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        recipRhoAlpha[is] = new std::complex<double>[pw_rho->npw];
        recipRhoBeta[is] = new std::complex<double>[pw_rho->npw];

        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            tempRho[ir] = std::pow(prho[is][ir], this->alpha_);
        }
        pw_rho->real2recip(tempRho, recipRhoAlpha[is]);

        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            tempRho[ir] = std::pow(prho[is][ir], this->beta_);
        }
        pw_rho->real2recip(tempRho, recipRhoBeta[is]);
    }

    double eta = 0.;
    double diff = 0.;
    this->stress.zero_out();
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        for (int ip = 0; ip < pw_rho->npw; ++ip)
        {
            eta = sqrt(pw_rho->gg[ip]) * pw_rho->tpiba / this->tkf_;
            diff = this->diff_linhard(eta, vw_weight);
            diff *= eta * (recipRhoAlpha[is][ip] * std::conj(recipRhoBeta[is][ip])).real();
            // cout << "diff    " << diff << endl;
            for (int a = 0; a < 3; ++a)
            {
                for (int b = a; b < 3; ++b)
                {
                    if (pw_rho->gg[ip] == 0.)
                    {
                        continue;
                    }
                    else
                    {
                        this->stress(a, b) += -diff * pw_rho->gcar[ip][a] * pw_rho->gcar[ip][b] / pw_rho->gg[ip];
                        if (a == b)
                            this->stress(a, b) += diff * coef;
                    }
                }
            }
        }
    }

    for (int a = 0; a < 3; ++a)
    {
        for (int b = a; b < 3; ++b)
        {
            Parallel_Reduce::reduce_all(this->stress(a, b));

            if (GlobalV::GAMMA_ONLY_PW)
            {
                this->stress(a, b) *= -std::pow(ModuleBase::PI, 2)
                                      / (this->alpha_ * this->beta_ * this->kf_
                                         * std::pow(this->rho0_, this->alpha_ + this->beta_ - 2.))
                                      * 2.; // multiply by 2 to convert Hartree to Ry
            }
            else
            {
                this->stress(a, b) *= -std::pow(ModuleBase::PI, 2)
                                      / (2. * this->alpha_ * this->beta_ * this->kf_
                                         * std::pow(this->rho0_, this->alpha_ + this->beta_ - 2.))
                                      * 2.; // multiply by 2 to convert Hartree to Ry
            }
        }
    }

    for (int a = 0; a < 3; ++a)
    {
        this->stress(a, a) += mult * this->wt_energy / pw_rho->omega;
        for (int b = 0; b < a; ++b)
        {
            this->stress(a, b) = this->stress(b, a);
        }
    }

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] recipRhoAlpha[is];
        delete[] recipRhoBeta[is];
    }
    delete[] recipRhoAlpha;
    delete[] recipRhoBeta;
    delete[] tempRho;
}

/**
 * @brief Calculate the WT kernel according to Lindhard response function
 *
 * @param eta k / (2 * kF)
 * @param tf_weight
 * @param vw_weight
 * @return W(eta)
 */
double KEDF_WT::wt_kernel(double eta, double tf_weight, double vw_weight)
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
        return 1. / (0.5 + 0.25 * (1. - eta * eta) / eta * log((1 + eta) / std::abs(1 - eta)))
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
double KEDF_WT::diff_linhard(double eta, double vw_weight)
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
        return ((eta2 + 1.) * 0.25 / eta2 * log(std::abs((1. + eta) / (1. - eta))) - 0.5 / eta)
                   / std::pow((0.5 + 0.25 * (1. - eta2) * log((1. + eta) / std::abs(1. - eta)) / eta), 2)
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
void KEDF_WT::multi_kernel(const double* const* prho, double** rkernel_rho, double exponent, ModulePW::PW_Basis* pw_rho)
{
    std::complex<double>** recipkernelRho = new std::complex<double>*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        recipkernelRho[is] = new std::complex<double>[pw_rho->npw];
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            rkernel_rho[is][ir] = std::pow(prho[is][ir], exponent);
        }
        pw_rho->real2recip(rkernel_rho[is], recipkernelRho[is]);
        for (int ip = 0; ip < pw_rho->npw; ++ip)
        {
            recipkernelRho[is][ip] *= this->kernel_[ip];
        }
        pw_rho->recip2real(recipkernelRho[is], rkernel_rho[is]);
    }

    for (int is = 0; is < GlobalV::NSPIN; ++is)
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
void KEDF_WT::fill_kernel(double tf_weight, double vw_weight, ModulePW::PW_Basis* pw_rho)
{
    double eta = 0.;
    for (int ig = 0; ig < pw_rho->npw; ++ig)
    {
        eta = sqrt(pw_rho->gg[ig]) * pw_rho->tpiba / this->tkf_;
        this->kernel_[ig] = this->wt_kernel(eta, tf_weight, vw_weight) * this->wt_coef_;
    }
}

//
/**
 * @brief Read the kernel from file
 *
 * @param file_name the name of the kernel file
 * @param pw_rho pw basis
 */
void KEDF_WT::read_kernel(std::string file_name, ModulePW::PW_Basis* pw_rho)
{
    std::ifstream ifs(file_name.c_str(), std::ios::in);

    GlobalV::ofs_running << "Read WT kernel from " << file_name << std::endl;
    if (!ifs)
        ModuleBase::WARNING_QUIT("kedf_wt.cpp", "The kernel file of WT KEDF not found");

    int kineType = 0;
    double kF_in = 0.;
    double rho0_in = 0.;
    int nq_in = 0;
    double maxEta_in = 0.;

    ifs >> kineType;
    ifs >> kF_in;
    ifs >> rho0_in;
    ifs >> nq_in;

    double* eta_in = new double[nq_in];
    double* w0_in = new double[nq_in];

    for (int iq = 0; iq < nq_in; ++iq)
    {
        ifs >> eta_in[iq] >> w0_in[iq];
    }

    maxEta_in = eta_in[nq_in - 1];

    double eta = 0.;
    double maxEta = 0.;
    int ind1 = 0;
    int ind2 = 0;
    int ind_mid = 0;
    double fac1 = 0.;
    double fac2 = 0.;
    for (int ig = 0; ig < pw_rho->npw; ++ig)
    {
        eta = sqrt(pw_rho->gg[ig]) * pw_rho->tpiba / this->tkf_;
        maxEta = std::max(eta, maxEta);

        if (eta <= eta_in[0])
            this->kernel_[ig] = w0_in[0];
        else if (eta > maxEta_in)
            this->kernel_[ig] = w0_in[nq_in - 1];
        else
        {
            ind1 = 1;
            ind2 = nq_in;
            while (ind1 < ind2 - 1)
            {
                ind_mid = (ind1 + ind2) / 2;
                if (eta > eta_in[ind_mid])
                {
                    ind1 = ind_mid;
                }
                else
                {
                    ind2 = ind_mid;
                }
            }
            fac1 = (eta_in[ind2] - eta) / (eta_in[ind2] - eta_in[ind1]);
            fac2 = (eta - eta_in[ind1]) / (eta_in[ind2] - eta_in[ind1]);
            this->kernel_[ig] = fac1 * w0_in[ind1] + fac2 * w0_in[ind2];
            this->kernel_[ig] *= this->wt_coef_;
        }
    }

    if (maxEta > maxEta_in)
        ModuleBase::WARNING("kedf_wt.cpp", "Please increase the maximal eta value in KEDF kernel file");

    delete[] eta_in;
    delete[] w0_in;
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "FILL WT KERNEL");
}