#include "./kedf_wt.h"
#include <iostream>
#include "module_base/parallel_reduce.h"
#include "module_base/tool_quit.h"
#include "module_base/global_variable.h"

void KEDF_WT::set_para(int nx, double dV, double alpha, double beta, double nelec, double tf_weight, double vw_weight, bool read_kernel, std::string kernel_file, ModulePW::PW_Basis *pw_rho)
{
    this->nx = nx;
    this->dV = dV;
    // this->weightWT = weightWT;
    this->alpha = alpha;
    this->beta = beta;

    if (GlobalV::of_wt_rho0 != 0)
    {
        this->rho0 = GlobalV::of_wt_rho0;
    }
    else
    {
        this->rho0 = 1./(pw_rho->nxyz * dV) * nelec;
    }

    this->kF = std::pow(3. * std::pow(ModuleBase::PI, 2) * this->rho0, 1./3.);
    this->tkF = 2. * this->kF;

    this->WTcoef = 5./(9. * this->alpha * this->beta * std::pow(this->rho0, this->alpha + this->beta - 5./3.));

    delete[] this->kernel;
    this->kernel = new double[pw_rho->npw];

    if (read_kernel)
        this->readKernel(kernel_file, pw_rho);
    else
        this->fillKernel(tf_weight, vw_weight, pw_rho);
}

//
// Ewt = Ctf * \int{rho^alpha * W(r - r') * rho^beta drdr'} + T_vW + T_TF, T_vW and T_TF will be added in ESolver_OF
//
double KEDF_WT::get_energy(const double * const * prho, ModulePW::PW_Basis *pw_rho)
{
    double **kernelRhoBeta = new double* [GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is) kernelRhoBeta[is] = new double[pw_rho->nrxx];
    this->multiKernel(prho, kernelRhoBeta, this->beta, pw_rho);

    double energy = 0.; // in Ry
    if (GlobalV::NSPIN == 1)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            energy += std::pow(prho[0][ir], this->alpha) * kernelRhoBeta[0][ir];
        }
        energy *= this->dV * this->cTF;
    }
    else if (GlobalV::NSPIN == 2)
    {
        // for (int is = 0; is < GlobalV::NSPIN; ++is)
        // {
        //     for (int ir = 0; ir < this->nx; ++ir)
        //     {
        //         energy += 2 * pphi[is][ir] * LapPhi[is][ir];
        //     }
        // }
        // energy *= 0.5 * this->dV * 0.5;
    }
    this->WTenergy = energy;
    Parallel_Reduce::reduce_all(this->WTenergy);

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] kernelRhoBeta[is];
    }
    delete[] kernelRhoBeta;

    return energy;
}

double KEDF_WT::get_energy_density(const double * const *prho, int is, int ir, ModulePW::PW_Basis *pw_rho)
{
    double **kernelRhoBeta = new double* [GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is) kernelRhoBeta[is] = new double[pw_rho->nrxx];
    this->multiKernel(prho, kernelRhoBeta, this->beta, pw_rho);

    double result = this->cTF * std::pow(prho[is][ir], this->alpha) * kernelRhoBeta[is][ir];

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] kernelRhoBeta[is];
    }
    delete[] kernelRhoBeta;
    return result;
}


//
// Vwt = cTF * [alpha rho^{alpha-1} \int{W(r - r')rho^{beta}(r') dr'} + beta rho^{beta-1} \int{W(r' - r)rho^{alpha}(r') dr'}]
//
void KEDF_WT::WT_potential(const double * const *prho, ModulePW::PW_Basis *pw_rho, ModuleBase::matrix &rpotential)
{
    ModuleBase::timer::tick("KEDF_WT", "wt_potential");

    double **kernelRhoBeta = new double* [GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is) kernelRhoBeta[is] = new double[pw_rho->nrxx];
    this->multiKernel(prho, kernelRhoBeta, this->beta, pw_rho);

    double **kernelRhoAlpha = new double* [GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is) kernelRhoAlpha[is] = new double[pw_rho->nrxx];
    this->multiKernel(prho, kernelRhoAlpha, this->alpha, pw_rho);

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            rpotential(is, ir) += this->cTF *
                                    (this->alpha * std::pow(prho[is][ir], this->alpha-1.) * kernelRhoBeta[is][ir]
                                    + this->beta * std::pow(prho[is][ir], this->beta-1.) * kernelRhoAlpha[is][ir]);
        }
    }

    // calculate energy
    double energy = 0.; // in Ry
    if (GlobalV::NSPIN == 1)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            energy += std::pow(prho[0][ir], this->alpha) * kernelRhoBeta[0][ir];
        }
        energy *= this->dV * this->cTF;
    }
    else if (GlobalV::NSPIN == 2)
    {
        // for (int is = 0; is < GlobalV::NSPIN; ++is)
        // {
        //     for (int ir = 0; ir < this->nx; ++ir)
        //     {
        //         energy += 2 * pphi[is][ir] * LapPhi[is][ir];
        //     }
        // }
        // energy *= 0.5 * this->dV * 0.5;
    }
    this->WTenergy = energy;
    Parallel_Reduce::reduce_all(this->WTenergy);

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] kernelRhoBeta[is];
        delete[] kernelRhoAlpha[is];
    }
    delete[] kernelRhoBeta;
    delete[] kernelRhoAlpha;
    ModuleBase::timer::tick("KEDF_WT", "wt_potential");
}

void KEDF_WT::get_stress(double cellVol, const double * const * prho, ModulePW::PW_Basis *pw_rho, double vw_weight)
{
    double coef = 0.;
    double mult = 0.;
    if (GlobalV::of_hold_rho0)
    {
        coef = 0.;
        mult = -1. + this->alpha + this->beta;
    }
    else
    {
        coef = 1./3.;
        mult = 2./3.;
    }

    std::complex<double> **recipRhoAlpha = new std::complex<double> *[GlobalV::NSPIN];
    std::complex<double> **recipRhoBeta = new std::complex<double> *[GlobalV::NSPIN];
    double *tempRho = new double[this->nx];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        recipRhoAlpha[is] = new std::complex<double>[pw_rho->npw];
        recipRhoBeta[is] = new std::complex<double>[pw_rho->npw];

        for (int ir = 0; ir < this->nx; ++ir)
        {
            tempRho[ir] = std::pow(prho[is][ir], this->alpha);
        }
        pw_rho->real2recip(tempRho, recipRhoAlpha[is]);

        for (int ir = 0; ir < this->nx; ++ir)
        {
            tempRho[ir] = std::pow(prho[is][ir], this->beta);
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
            eta = sqrt(pw_rho->gg[ip]) * pw_rho->tpiba / this->tkF;
            diff = this->diffLinhard(eta, vw_weight);
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
                        this->stress(a,b) += - diff * pw_rho->gcar[ip][a] * pw_rho->gcar[ip][b] / pw_rho->gg[ip];
                        if (a == b) this->stress(a,b) += diff * coef;
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
                this->stress(a,b) *= - std::pow(ModuleBase::PI,2) / (this->alpha * this->beta * this->kF * std::pow(this->rho0, this->alpha + this->beta - 2.)) * 2.; // multiply by 2 to convert Hartree to Ry
            }
            else
            {
                this->stress(a,b) *= - std::pow(ModuleBase::PI,2) / (2. * this->alpha * this->beta * this->kF * std::pow(this->rho0, this->alpha + this->beta - 2.)) * 2.; // multiply by 2 to convert Hartree to Ry
            }
        }
    }

    for (int a = 0; a < 3; ++a)
    {
        this->stress(a,a) += mult * this->WTenergy / cellVol;
        for (int b = 0; b < a; ++b)
        {
            this->stress(a,b) = this->stress(b,a);
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

// Calculate WT kernel according to Lindhard response function
double KEDF_WT::WTkernel(double eta, double tf_weight, double vw_weight)
{
    if (eta < 0.)
    {
        return 0.;
    }
    // limit for small eta
    else if (eta < 1e-10)
    {
        return 1. - tf_weight + eta * eta * (1./3. - 3. * vw_weight);
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
        return 1. / (0.5 + 0.25 * (1. - eta * eta) / eta * log((1 + eta)/std::abs(1 - eta)))
                 - 3. * vw_weight * eta * eta - tf_weight;
    }
}

double KEDF_WT::diffLinhard(double eta, double vw_weight)
{
    if (eta < 0.)
    {
        return 0.;
    }
    else if (eta < 1e-10)
    {
        return 2. * eta * (1./3. - 3. * vw_weight);
    }
    else if (std::abs(eta - 1.) < 1e-10)
    {
        return 40.;
    }
    else
    {
        double eta2 = eta * eta;
        return ((eta2 + 1.) * 0.25 / eta2 * log(std::abs((1. + eta)/
             (1.-eta))) - 0.5/eta) / std::pow((0.5 + 0.25 * (1. - eta2)
             * log((1. + eta) / std::abs(1. - eta))/eta) , 2) - 6. * eta * vw_weight;
    }
}

// Calculate \int{W(r-r')rho^{exponent}(r') dr'}
void KEDF_WT::multiKernel(const double * const * prho, double **rkernelRho, double exponent, ModulePW::PW_Basis *pw_rho)
{
    std::complex<double> **recipkernelRho = new std::complex<double> *[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        recipkernelRho[is] = new std::complex<double>[pw_rho->npw];
        for (int ir = 0; ir < this->nx; ++ir)
        {
            rkernelRho[is][ir] = std::pow(prho[is][ir], exponent);
        }
        pw_rho->real2recip(rkernelRho[is], recipkernelRho[is]);
        for (int ip = 0; ip < pw_rho->npw; ++ip)
        {
            recipkernelRho[is][ip] *= this->kernel[ip];
        }
        pw_rho->recip2real(recipkernelRho[is], rkernelRho[is]);
    }

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] recipkernelRho[is];
    }
    delete[] recipkernelRho;
}

void KEDF_WT::fillKernel(double tf_weight, double vw_weight, ModulePW::PW_Basis *pw_rho)
{
    double eta = 0.;
    for (int ig = 0; ig < pw_rho->npw; ++ig)
    {
        eta = sqrt(pw_rho->gg[ig]) * pw_rho->tpiba / this->tkF;
        this->kernel[ig] = this->WTkernel(eta, tf_weight, vw_weight) * this->WTcoef;
    }
}

// Read kernel from file
void KEDF_WT::readKernel(std::string fileName, ModulePW::PW_Basis *pw_rho)
{
    std::ifstream ifs(fileName.c_str(), std::ios::in);

    GlobalV::ofs_running << "Read WT kernel from " << fileName << std::endl;
    if (!ifs) ModuleBase::WARNING_QUIT("kedf_wt.cpp", "The kernel file of WT KEDF not found");

    int kineType = 0;
    double kF_in = 0.;
    double rho0_in = 0.;
    int nq_in = 0;
    double maxEta_in = 0.;

    ifs >> kineType;
    ifs >> kF_in;
    ifs >> rho0_in;
    ifs >> nq_in;

    double *eta_in = new double[nq_in];
    double *w0_in = new double[nq_in];

    for (int iq = 0; iq < nq_in; ++iq)
    {
        ifs >> eta_in[iq] >> w0_in[iq];
    }

    maxEta_in = eta_in[nq_in-1];

    double eta = 0.;
    double maxEta = 0.;
    int ind1 = 0;
    int ind2 = 0;
    int ind_mid = 0;
    double fac1 = 0.;
    double fac2 = 0.;
    for (int ig = 0; ig < pw_rho->npw; ++ig)
    {
        eta = sqrt(pw_rho->gg[ig]) * pw_rho->tpiba / this->tkF;
        maxEta = std::max(eta, maxEta);

        if (eta <= eta_in[0])
            this->kernel[ig] = w0_in[0];
        else if (eta > maxEta_in)
            this->kernel[ig] = w0_in[nq_in-1];
        else
        {
            ind1 = 1;
            ind2 = nq_in;
            while (ind1 < ind2 - 1)
            {
                ind_mid = (ind1 + ind2)/2;
                if (eta > eta_in[ind_mid])
                {
                    ind1 = ind_mid;
                }
                else
                {
                    ind2 = ind_mid;
                }
            }
            fac1 = (eta_in[ind2] - eta)/(eta_in[ind2] - eta_in[ind1]);
            fac2 = (eta - eta_in[ind1])/(eta_in[ind2] - eta_in[ind1]);
            this->kernel[ig] = fac1 * w0_in[ind1] + fac2 * w0_in[ind2];
            this->kernel[ig] *= this->WTcoef;
        }
    }

    if (maxEta > maxEta_in) ModuleBase::WARNING("kedf_wt.cpp", "Please increase the maximal eta value in KEDF kernel file");

    delete[] eta_in;
    delete[] w0_in;
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "FILL WT KERNEL");
}