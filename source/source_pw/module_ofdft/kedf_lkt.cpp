#include "./kedf_lkt.h"

#include "source_io/module_parameter/parameter.h"
#include <iostream>

#include "source_base/parallel_reduce.h"

void KEDF_LKT::set_para(double dV, double lkt_a)
{
    this->dV_ = dV;
    this->lkt_a_ = lkt_a;
}

// Allocate persistent working buffers (called once in init)
void KEDF_LKT::init_buffers(const int nrxx)
{
    free_buffers(); // safety: free any existing buffers first
    as_ = new double[nrxx];
    nabla_term_ = new double[nrxx];
    for (int i = 0; i < 3; ++i)
    {
        nabla_rho_[i] = new double[nrxx];
        div_input_[i] = new double[nrxx];
    }
    buffer_nrxx_ = nrxx;
}

// Free working buffers
void KEDF_LKT::free_buffers()
{
    delete[] as_;
    delete[] nabla_term_;
    for (int i = 0; i < 3; ++i)
    {
        delete[] nabla_rho_[i];
        delete[] div_input_[i];
        nabla_rho_[i] = nullptr;
        div_input_[i] = nullptr;
    }
    as_ = nullptr;
    nabla_term_ = nullptr;
    buffer_nrxx_ = 0;
}

/**
 * @brief Get the energy of LKT KEDF
 * \f[ E_{LKT} = \int{tau_{TF}/\cosh(a * s)}, s = c_s * |\nabla \rho|/\rho^{4/3} \f]
 *
 * @param prho charge density
 * @param pw_rho pw basis
 * @return the energy of LKT KEDF
 */
double KEDF_LKT::get_energy(const double* const* prho, ModulePW::PW_Basis* pw_rho)
{
    double energy = 0.; // in Ry

    if (PARAM.inp.nspin == 1)
    {
        this->nabla(prho[0], pw_rho, nabla_rho_);
        this->get_as(prho[0], nabla_rho_, pw_rho->nrxx, as_);

#ifdef _OPENMP
#pragma omp parallel for reduction(+:energy) schedule(static)
#endif
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            energy += std::pow(prho[0][ir], 5. / 3.) / std::cosh(as_[ir]);
        }
        energy *= this->dV_ * this->c_tf_;
    }
    else if (PARAM.inp.nspin == 2)
    {
        // Waiting for update
    }
    this->lkt_energy = energy;
    Parallel_Reduce::reduce_all(this->lkt_energy);

    return this->lkt_energy;
}

/**
 * @brief Get the energy density of LKT KEDF
 * \f[ \tau_{LKT} = tau_{TF}/\cosh(a * s) \f]
 *
 * @param prho charge density
 * @param is the index of spin
 * @param ir the index of real space grid
 * @param pw_rho pw basis
 * @return the energy density of LKT KEDF
 */
double KEDF_LKT::get_energy_density(const double* const* prho, int is, int ir, ModulePW::PW_Basis* pw_rho)
{
    // Use the pre-allocated buffers for the full gradient computation
    this->nabla(prho[is], pw_rho, nabla_rho_);
    this->get_as(prho[is], nabla_rho_, pw_rho->nrxx, as_);
    
    double energy_den = this->c_tf_ * std::pow(prho[is][ir], 5. / 3.) / std::cosh(as_[ir]);

    return energy_den;
}

void KEDF_LKT::tau_lkt(const double* const* prho, ModulePW::PW_Basis* pw_rho, double* rtau_lkt)
{
    if (PARAM.inp.nspin == 1)
    {
        this->nabla(prho[0], pw_rho, nabla_rho_);
        this->get_as(prho[0], nabla_rho_, pw_rho->nrxx, as_);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            double coshas = std::cosh(as_[ir]);
            rtau_lkt[ir] += std::pow(prho[0][ir], 5. / 3.) / coshas * this->c_tf_;
        }
    }
    else if (PARAM.inp.nspin == 2)
    {
        // Waiting for update
    }
}

/**
 * @brief Get the potential of LKT KEDF, and add it into rpotential,
 * and the LKT energy will be calculated and stored in this->lkt_energy
 * \f[ V_{LKT} =5/3 *\tau_{TF}/\rho * [1/\cosh(as)+5/4 * as * \tanh(as)/\cosh(as)] \f]
 * \f[ +\nabla\cdot(\tau_{TF} * a*\tanh(as)/\cosh(as) * s/|\nabla\rho|^2 * \nabla\rho). \f]
 *
 * Uses pre-allocated buffers and OpenMP parallelization.
 *
 * @param prho charge density
 * @param pw_rho pw basis
 * @param rpotential rpotential => rpotential + V_{LKT}
 */
void KEDF_LKT::lkt_potential(const double* const* prho, ModulePW::PW_Basis* pw_rho, ModuleBase::matrix& rpotential)
{
    ModuleBase::TITLE("KEDF_LKT", "lkt_potential");
    ModuleBase::timer::start("KEDF_LKT", "lkt_potential");
    this->lkt_energy = 0.;

    if (PARAM.inp.nspin == 1)
    {
        this->nabla(prho[0], pw_rho, nabla_rho_);
        this->get_as(prho[0], nabla_rho_, pw_rho->nrxx, as_);

        // OpenMP-parallelized grid loop: prepare div_input_ from nabla_rho_,
        // compute the first term of potential in the same pass
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            double coshas = std::cosh(as_[ir]);
            double tanhas = std::tanh(as_[ir]);

            // add the first term of potential
            rpotential(0, ir) += 5.0 / 3.0 * this->c_tf_ * std::pow(prho[0][ir], 2. / 3.) / coshas
                                 * (1. + 4.0 / 5.0 * as_[ir] * tanhas);

            // prepare div_input_ for divergence()
            for (int i = 0; i < 3; ++i)
            {
                if (as_[ir] == 0)
                {
                    div_input_[i][ir] = 0;
                }
                else
                {
                    div_input_[i][ir] = nabla_rho_[i][ir] * tanhas / coshas / as_[ir] / prho[0][ir]
                                        * this->c_tf_ * std::pow(this->s_coef_ * this->lkt_a_, 2);
                }
            }
        }

        // Energy accumulation (separate loop for correct reduction)
        double energy_local = 0.;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:energy_local) schedule(static)
#endif
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            energy_local += std::pow(prho[0][ir], 5. / 3.) / std::cosh(as_[ir]);
        }

        this->divergence(div_input_, pw_rho, nabla_term_);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            rpotential(0, ir) += nabla_term_[ir];
        }

        this->lkt_energy = energy_local * this->c_tf_ * this->dV_;
        Parallel_Reduce::reduce_all(this->lkt_energy);
    }
    else if (PARAM.inp.nspin == 2)
    {
        // Waiting for update
    }

    ModuleBase::timer::end("KEDF_LKT", "lkt_potential");
}

/**
 * @brief Get the stress of LKT KEDF, and store it into this->stress
 *
 * Uses pre-allocated buffers and OpenMP parallelization.
 *
 * @param prho charge density
 * @param pw_rho pw basis
 */
void KEDF_LKT::get_stress(const double* const* prho, ModulePW::PW_Basis* pw_rho)
{
    if (PARAM.inp.nspin == 1)
    {
        this->nabla(prho[0], pw_rho, nabla_rho_);
        this->get_as(prho[0], nabla_rho_, pw_rho->nrxx, as_);

        // Pre-compute energy for stress diagonal
        double lkt_energy_local = 0.;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:lkt_energy_local) schedule(static)
#endif
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            lkt_energy_local += std::pow(prho[0][ir], 5. / 3.) / std::cosh(as_[ir]);
        }
        lkt_energy_local *= this->c_tf_ * this->dV_;
        Parallel_Reduce::reduce_all(lkt_energy_local);

        for (int alpha = 0; alpha < 3; ++alpha)
        {
            for (int beta = alpha; beta < 3; ++beta)
            {
                this->stress(alpha, beta) = 0;

                if (alpha == beta)
                {
                    this->stress(alpha, beta) = 2.0 / 3.0 / pw_rho->omega * lkt_energy_local;
                }

                double integral_term = 0.;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:integral_term) schedule(static)
#endif
                for (int ir = 0; ir < pw_rho->nrxx; ++ir)
                {
                    double coef = std::tanh(as_[ir]) / std::cosh(as_[ir]);
                    if (as_[ir] != 0.)
                    {
                        integral_term += -nabla_rho_[alpha][ir] * nabla_rho_[beta][ir] / as_[ir] / prho[0][ir]
                                         * std::pow(this->s_coef_ * this->lkt_a_, 2) * coef;
                    }
                    if (alpha == beta)
                    {
                        integral_term += 1.0 / 3.0 * as_[ir] * std::pow(prho[0][ir], 5.0 / 3.0) * coef;
                    }
                }
                Parallel_Reduce::reduce_all(integral_term);
                integral_term *= this->c_tf_ * this->dV_ / pw_rho->omega;

                this->stress(alpha, beta) += integral_term;
            }
        }
        for (int alpha = 1; alpha < 3; ++alpha)
        {
            for (int beta = 0; beta < alpha; ++beta)
            {
                this->stress(alpha, beta) = this->stress(beta, alpha);
            }
        }
    }
    else if (PARAM.inp.nspin == 2)
    {
        // Waiting for update
    }
}

/**
 * @brief Caculate routput = nabla(pinput)
 *
 * @param [in] pinput
 * @param [in] pw_rho pw basis
 * @param [out] routput
 */
void KEDF_LKT::nabla(const double* pinput, ModulePW::PW_Basis* pw_rho, double** routput)
{
    std::complex<double>* recip_data = new std::complex<double>[pw_rho->npw];
    std::complex<double>* recip_nabla = new std::complex<double>[pw_rho->npw];
    pw_rho->real2recip(pinput, recip_data);

    std::complex<double> img{0.0, 1.0};
    for (int j = 0; j < 3; ++j)
    {
        for (int ip = 0; ip < pw_rho->npw; ++ip)
        {
            recip_nabla[ip] = img * pw_rho->gcar[ip][j] * recip_data[ip] * pw_rho->tpiba;
        }
        pw_rho->recip2real(recip_nabla, routput[j]);
    }

    delete[] recip_data;
    delete[] recip_nabla;
}

/**
 * @brief Caculate routput = nabla dot pinput
 *
 * @param [in] pinput
 * @param [in] pw_rho pw basis
 * @param [out] routput
 */
void KEDF_LKT::divergence(const double* const* pinput, ModulePW::PW_Basis* pw_rho, double* routput)
{
    std::complex<double>* recip_container = new std::complex<double>[pw_rho->npw];
    std::complex<double> img{0.0, 1.0};
    ModuleBase::GlobalFunc::ZEROS(routput, pw_rho->nrxx);
    for (int i = 0; i < 3; ++i)
    {
        pw_rho->real2recip(pinput[i], recip_container);
        for (int ip = 0; ip < pw_rho->npw; ++ip)
        {
            recip_container[ip] = img * pw_rho->gcar[ip][i] * pw_rho->tpiba * recip_container[ip];
        }
        pw_rho->recip2real(recip_container, routput, true);
    }

    delete[] recip_container;
}

/**
 * @brief Caculate as = lkt_a * s, s = c_s * |\nabla \rho|/\rho^{4/3}
 *
 * OpenMP parallelized real-space loop.
 *
 * @param [in] prho charge density
 * @param [in] pnabla_rho nabla rho
 * @param [in] nrxx the number of real space grid
 * @param [out] as lkt_a * s
 */
void KEDF_LKT::get_as(const double* prho, const double* const* pnabla_rho, const int nrxx, double* as)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ir = 0; ir < nrxx; ++ir)
    {
        as[ir] = std::sqrt(std::pow(pnabla_rho[0][ir], 2) + std::pow(pnabla_rho[1][ir], 2)
                           + std::pow(pnabla_rho[2][ir], 2))
                 / std::pow(prho[ir], 4.0 / 3.0) * this->s_coef_ * this->lkt_a_;
    }
}
