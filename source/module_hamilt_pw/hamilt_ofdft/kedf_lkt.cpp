#include "./kedf_lkt.h"

#include <iostream>

#include "module_base/parallel_reduce.h"

void KEDF_LKT::set_para(int nx, double dV, double lkt_a)
{
    this->nx = nx;
    this->dV = dV;
    this->lkt_a = lkt_a;
}

//
// ELKT = \int{tau_TF/cosh(a * s)}
//
double KEDF_LKT::get_energy(const double *const *prho, ModulePW::PW_Basis *pw_rho)
{
    double energy = 0.;                // in Ry
    double *as = new double[this->nx]; // a*s
    double **nabla_rho = new double *[3];
    for (int i = 0; i < 3; ++i)
    {
        nabla_rho[i] = new double[this->nx];
    }

    if (GlobalV::NSPIN == 1)
    {
        this->nabla(prho[0], pw_rho, nabla_rho);
        this->get_as(prho[0], nabla_rho, as);
        for (int ir = 0; ir < this->nx; ++ir)
        {
            energy += pow(prho[0][ir], 5. / 3.) / std::cosh(as[ir]);
        }
        energy *= this->dV * this->cTF;
    }
    else if (GlobalV::NSPIN == 2)
    {
        // Waiting for update
    }
    this->LKTenergy = energy;
    Parallel_Reduce::reduce_all(this->LKTenergy);

    delete[] as;
    for (int i = 0; i < 3; ++i)
        delete[] nabla_rho[i];
    delete[] nabla_rho;

    return energy;
}

double KEDF_LKT::get_energy_density(const double *const *prho, int is, int ir, ModulePW::PW_Basis *pw_rho)
{
    double energy_den = 0.;            // in Ry
    double *as = new double[this->nx]; // a*s
    double **nabla_rho = new double *[3];
    for (int i = 0; i < 3; ++i)
    {
        nabla_rho[i] = new double[this->nx];
    }

    this->nabla(prho[is], pw_rho, nabla_rho);
    this->get_as(prho[is], nabla_rho, as);
    energy_den = this->cTF * pow(prho[is][ir], 5. / 3.) / std::cosh(as[ir]);

    delete[] as;
    for (int i = 0; i < 3; ++i)
        delete[] nabla_rho[i];
    delete[] nabla_rho;

    return energy_den;
}

//
// V_LKT =5/3 *tau_TF/rho * [1/cosh(as)+5/4 * as * tanh(as)/cosh(as)]
// +\nabla\cdot(tau_TF * a*tanh(as)/cosh(as) * s/|\nabla\rho|^2 * \nabla\rho).
//
void KEDF_LKT::lkt_potential(const double *const *prho, ModulePW::PW_Basis *pw_rho, ModuleBase::matrix &rpotential)
{
    ModuleBase::timer::tick("KEDF_LKT", "LKT_potential");
    this->LKTenergy = 0.;
    double *as = new double[this->nx]; // a*s
    double **nabla_rho = new double *[3];
    for (int i = 0; i < 3; ++i)
    {
        nabla_rho[i] = new double[this->nx];
    }
    double *nabla_term = new double[this->nx];

    if (GlobalV::NSPIN == 1)
    {
        this->nabla(prho[0], pw_rho, nabla_rho);
        this->get_as(prho[0], nabla_rho, as);

        for (int ir = 0; ir < this->nx; ++ir)
        {
            double coshas = std::cosh(as[ir]);
            double tanhas = std::tanh(as[ir]);

            this->LKTenergy += std::pow(prho[0][ir], 5. / 3.) / coshas;
            // add the first term
            rpotential(0, ir)
                += 5.0 / 3.0 * this->cTF * std::pow(prho[0][ir], 2. / 3.) / coshas * (1. + 4.0 / 5.0 * as[ir] * tanhas);
            // get the nabla_term
            for (int i = 0; i < 3; ++i)
            {
                if (as[ir] == 0)
                {
                    nabla_rho[i][ir] = 0;
                }
                else
                {
                    nabla_rho[i][ir] = nabla_rho[i][ir] * tanhas / coshas / as[ir] / prho[0][ir] * this->cTF
                                       * std::pow(this->s_coef * this->lkt_a, 2);
                }
            }
        }

        this->divergence(nabla_rho, pw_rho, nabla_term);

        for (int ir = 0; ir < this->nx; ++ir)
        {
            rpotential(0, ir) += nabla_term[ir];
        }

        this->LKTenergy *= this->cTF * this->dV;
        Parallel_Reduce::reduce_all(this->LKTenergy);
    }
    else if (GlobalV::NSPIN == 2)
    {
        // Waiting for update
    }

    delete[] as;
    for (int i = 0; i < 3; ++i)
        delete[] nabla_rho[i];
    delete[] nabla_rho;
    delete[] nabla_term;

    ModuleBase::timer::tick("KEDF_LKT", "LKT_potential");
}

void KEDF_LKT::get_stress(const double cell_vol, const double *const *prho, ModulePW::PW_Basis *pw_rho)
{
    double *as = new double[this->nx]; // a*s
    double **nabla_rho = new double *[3];
    for (int i = 0; i < 3; ++i)
    {
        nabla_rho[i] = new double[this->nx];
    }
    double *nabla_term = new double[this->nx];

    if (GlobalV::NSPIN == 1)
    {
        this->nabla(prho[0], pw_rho, nabla_rho);
        this->get_as(prho[0], nabla_rho, as);

        for (int alpha = 0; alpha < 3; ++alpha)
        {
            for (int beta = alpha; beta < 3; ++beta)
            {
                this->stress(alpha, beta) = 0;

                if (alpha == beta)
                {
                    this->stress(alpha, beta) = 2.0 / 3.0 / cell_vol * this->LKTenergy;
                }

                double integral_term = 0.;
                for (int ir = 0; ir < this->nx; ++ir)
                {
                    double coef = std::tanh(as[ir]) / std::cosh(as[ir]);
                    if (as[ir] != 0.)
                    {
                        integral_term += -nabla_rho[alpha][ir] * nabla_rho[beta][ir] / as[ir] / prho[0][ir]
                                         * std::pow(this->s_coef * this->lkt_a, 2) * coef;
                    }
                    if (alpha == beta)
                    {
                        integral_term += 1.0 / 3.0 * as[ir] * std::pow(prho[0][ir], 5.0 / 3.0) * coef;
                    }
                }
                Parallel_Reduce::reduce_all(integral_term);
                integral_term *= this->cTF * this->dV / cell_vol;

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
    else if (GlobalV::NSPIN == 2)
    {
        // Waiting for update
    }

    delete[] as;
    for (int i = 0; i < 3; ++i)
        delete[] nabla_rho[i];
    delete[] nabla_rho;
    delete[] nabla_term;
}

// output = nabla input
void KEDF_LKT::nabla(const double *pinput, ModulePW::PW_Basis *pw_rho, double **routput)
{
    std::complex<double> *recip_data = new std::complex<double>[pw_rho->npw];
    std::complex<double> *recip_nabla = new std::complex<double>[pw_rho->npw];
    pw_rho->real2recip(pinput, recip_data);

    std::complex<double> img {0.0, 1.0};
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

// output = nabla dot input
void KEDF_LKT::divergence(const double *const *pinput, ModulePW::PW_Basis *pw_rho, double *routput)
{
    std::complex<double> *recip_container = new std::complex<double>[pw_rho->npw];
    std::complex<double> img {0.0, 1.0};
    ModuleBase::GlobalFunc::ZEROS(routput, this->nx);
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

// lkt_a * s, s = c_s * |nabla rho|/rho^{4/3}
void KEDF_LKT::get_as(const double *prho, const double *const *pnabla_rho, double *as)
{
    for (int ir = 0; ir < this->nx; ++ir)
    {
        as[ir] = std::sqrt(std::pow(pnabla_rho[0][ir], 2) + std::pow(pnabla_rho[1][ir], 2)
                           + std::pow(pnabla_rho[2][ir], 2))
                 / std::pow(prho[ir], 4.0 / 3.0) * this->s_coef * this->lkt_a;
    }
}
