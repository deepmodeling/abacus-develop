#include "cal_mlkedf_descriptors.h"

namespace ModuleIO
{

void Cal_MLKEDF_Descriptors::set_para(
    const int &nx,
    const double &nelec, 
    const double &tf_weight, 
    const double &vw_weight,
    const double &chi_p,
    const double &chi_q,
    const std::vector<double> &chi_xi,
    const std::vector<double> &chi_pnl,
    const std::vector<double> &chi_qnl,
    const int &nkernel,
    const std::vector<int> &kernel_type,
    const std::vector<double> &kernel_scaling,
    const std::vector<double> &yukawa_alpha,
    const std::vector<std::string> &kernel_file,
    const double &omega,
    ModulePW::PW_Basis *pw_rho
)
{
    this->nx = nx;
    this->nkernel = nkernel;
    this->chi_p = chi_p;
    this->chi_q = chi_q;
    this->chi_xi = chi_xi;
    this->chi_pnl = chi_pnl;
    this->chi_qnl = chi_qnl;

    this->kernel_type = kernel_type;
    this->kernel_scaling = kernel_scaling;
    this->yukawa_alpha = yukawa_alpha;
    this->kernel_file = kernel_file;
    std::cout << "nkernel = " << nkernel << std::endl;

    if (PARAM.inp.of_wt_rho0 != 0)
    {
        this->rho0 = PARAM.inp.of_wt_rho0;
    }
    else
    {
        this->rho0 = 1./omega * nelec;
    }

    this->kF = std::pow(3. * std::pow(ModuleBase::PI, 2) * this->rho0, 1./3.);
    this->tkF = 2. * this->kF;

    this->kernel = std::vector<std::vector<double>>(this->nkernel);
    for (int ik = 0; ik < this->nkernel; ++ik)
    {
        // delete[] this->kernel[ik];
        this->kernel[ik] = std::vector<double>(pw_rho->npw, 0.0);
        if (this->kernel_type[ik] == 3 || this->kernel_type[ik] == 4) // 3 for TKK Al, and 4 for TKK Si
        {
            this->read_kernel(this->kernel_file[ik], this->kernel_scaling[ik], pw_rho, this->kernel[ik].data());
        }
        else
        {
            double eta = 0.;
            for (int ip = 0; ip < pw_rho->npw; ++ip)
            {
                eta = sqrt(pw_rho->gg[ip]) * pw_rho->tpiba / this->tkF * this->kernel_scaling[ik];
                if (this->kernel_type[ik] == 1)
                {
                    this->kernel[ik][ip] = this->MLkernel(eta, tf_weight, vw_weight);
                }
                else if (this->kernel_type[ik] == 2)
                {
                    this->kernel[ik][ip] = this->MLkernel_yukawa(eta, this->yukawa_alpha[ik]);
                }
            }
        }
    }
}

double Cal_MLKEDF_Descriptors::MLkernel(double eta, double tf_weight, double vw_weight)
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
        return 1. / (0.5 + 0.25 * (1. - eta * eta) / eta * std::log((1 + eta)/std::abs(1 - eta)))
                 - 3. * vw_weight * eta * eta - tf_weight;
    }
}

double Cal_MLKEDF_Descriptors::MLkernel_yukawa(double eta, double alpha)
{
    return (eta == 0 && alpha == 0) ? 0. : M_PI / (eta * eta + alpha * alpha / 4.);
}

// Read kernel from file
void Cal_MLKEDF_Descriptors::read_kernel(const std::string &fileName, const double& scaling, ModulePW::PW_Basis *pw_rho, double* kernel_)
{
    std::ifstream ifs(fileName.c_str(), std::ios::in);

    GlobalV::ofs_running << "Read WT kernel from " << fileName << std::endl;
    if (!ifs) ModuleBase::WARNING_QUIT("cal_mlkedf_descriptors.cpp", "The kernel file not found");

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
        eta = eta * scaling;
        maxEta = std::max(eta, maxEta);

        if (eta <= eta_in[0])
            kernel_[ig] = w0_in[0];
        else if (eta > maxEta_in)
            kernel_[ig] = w0_in[nq_in-1];
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
            kernel_[ig] = fac1 * w0_in[ind1] + fac2 * w0_in[ind2];
        }
    }

    if (maxEta > maxEta_in) ModuleBase::WARNING("cal_mlkedf_descriptors.cpp", "Please increase the maximal eta value in KEDF kernel file");

    delete[] eta_in;
    delete[] w0_in;
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "FILL WT KERNEL");
}

void Cal_MLKEDF_Descriptors::multiKernel(const int ikernel, double *pinput, ModulePW::PW_Basis *pw_rho, double *routput)
{
    std::complex<double> *recipOutput = new std::complex<double>[pw_rho->npw];

    pw_rho->real2recip(pinput, recipOutput);
    for (int ip = 0; ip < pw_rho->npw; ++ip)
    {
        recipOutput[ip] *= this->kernel[ikernel][ip];
    }
    pw_rho->recip2real(recipOutput, routput);

    delete[] recipOutput;
}

void Cal_MLKEDF_Descriptors::Laplacian(double * pinput, ModulePW::PW_Basis *pw_rho, double * routput)
{
    std::complex<double> *recipContainer = new std::complex<double>[pw_rho->npw];

    pw_rho->real2recip(pinput, recipContainer);
    for (int ip = 0; ip < pw_rho->npw; ++ip)
    {
        recipContainer[ip] *= - pw_rho->gg[ip] * pw_rho->tpiba2;
    }
    pw_rho->recip2real(recipContainer, routput);

    delete[] recipContainer;
}

void Cal_MLKEDF_Descriptors::divergence(double ** pinput, ModulePW::PW_Basis *pw_rho, double * routput)
{
    std::complex<double> *recipContainer = new std::complex<double>[pw_rho->npw];
    std::complex<double> img(0.0, 1.0);
    ModuleBase::GlobalFunc::ZEROS(routput, this->nx);
    for (int i = 0; i < 3; ++i)
    {
        pw_rho->real2recip(pinput[i], recipContainer);
        for (int ip = 0; ip < pw_rho->npw; ++ip)
        {
            recipContainer[ip] = img * pw_rho->gcar[ip][i] * pw_rho->tpiba * recipContainer[ip];
        }
        pw_rho->recip2real(recipContainer, routput, true);
    }

    delete[] recipContainer;
}

void Cal_MLKEDF_Descriptors::tanh(std::vector<double> &pinput, std::vector<double> &routput, double chi)
{
    for (int i = 0; i < this->nx; ++i)
    {
        routput[i] = std::tanh(pinput[i] * chi);
    }
}

double Cal_MLKEDF_Descriptors::dtanh(double tanhx, double chi)
{
    return (1. - tanhx * tanhx) * chi;
}

void Cal_MLKEDF_Descriptors::getGamma(const double * const *prho, std::vector<double> &rgamma)
{
    for(int ir = 0; ir < this->nx; ++ir)
    {
        rgamma[ir] = std::pow(prho[0][ir]/this->rho0, 1./3.);
    }
}

void Cal_MLKEDF_Descriptors::getP(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<std::vector<double>> &pnablaRho, std::vector<double> &rp)
{
    for(int ir = 0; ir < this->nx; ++ir)
    {
        rp[ir] = 0.;
        for (int j = 0; j < 3; ++j)
        {
            rp[ir] += std::pow(pnablaRho[j][ir], 2);
        }
        rp[ir] *= this->pqcoef / std::pow(prho[0][ir], 8.0/3.0);
    }
}

void Cal_MLKEDF_Descriptors::getQ(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rq)
{
    // get Laplacian rho
    std::complex<double> *recipRho = new std::complex<double>[pw_rho->npw];
    pw_rho->real2recip(prho[0], recipRho);
    for (int ip = 0; ip < pw_rho->npw; ++ip)
    {
        recipRho[ip] *= - pw_rho->gg[ip] * pw_rho->tpiba2;
    }
    pw_rho->recip2real(recipRho, rq.data());

    for (int ir = 0; ir < this->nx; ++ir)
    {
        rq[ir] *= this->pqcoef / std::pow(prho[0][ir], 5.0/3.0);
    }

    delete[] recipRho;
}

void Cal_MLKEDF_Descriptors::getGammanl(const int ikernel, std::vector<double> &pgamma, ModulePW::PW_Basis *pw_rho, std::vector<double> &rgammanl)
{
    this->multiKernel(ikernel, pgamma.data(), pw_rho, rgammanl.data());
}

void Cal_MLKEDF_Descriptors::getPnl(const int ikernel, std::vector<double> &pp, ModulePW::PW_Basis *pw_rho, std::vector<double> &rpnl)
{
    this->multiKernel(ikernel, pp.data(), pw_rho, rpnl.data());
}

void Cal_MLKEDF_Descriptors::getQnl(const int ikernel, std::vector<double> &pq, ModulePW::PW_Basis *pw_rho, std::vector<double> &rqnl)
{
    this->multiKernel(ikernel, pq.data(), pw_rho, rqnl.data());
}

// xi = gammanl/gamma
void Cal_MLKEDF_Descriptors::getXi(std::vector<double> &pgamma, std::vector<double> &pgammanl, std::vector<double> &rxi)
{
    for (int ir = 0; ir < this->nx; ++ir)
    {
        if (pgamma[ir] == 0)
        {
            std::cout << "WARNING: gamma=0" << std::endl;
            rxi[ir] = 0.;
        }
        else
        {
            rxi[ir] = pgammanl[ir]/pgamma[ir];
        }
    }
}

// tanhxi = tanh(gammanl/gamma)
void Cal_MLKEDF_Descriptors::getTanhXi(const int ikernel, std::vector<double> &pgamma, std::vector<double> &pgammanl, std::vector<double> &rtanhxi)
{
    for (int ir = 0; ir < this->nx; ++ir)
    {
        if (pgamma[ir] == 0)
        {
            std::cout << "WARNING: gamma=0" << std::endl;
            rtanhxi[ir] = 0.;
        }
        else
        {
            rtanhxi[ir] = std::tanh(pgammanl[ir]/pgamma[ir] * this->chi_xi[ikernel]);
        }
    }
}

// tanh(p)
void Cal_MLKEDF_Descriptors::getTanhP(std::vector<double> &pp, std::vector<double> &rtanhp)
{
    this->tanh(pp, rtanhp, this->chi_p);
}

// tanh(q)
void Cal_MLKEDF_Descriptors::getTanhQ(std::vector<double> &pq, std::vector<double> &rtanhq)
{
    this->tanh(pq, rtanhq, this->chi_q);
}

// tanh(pnl)
void Cal_MLKEDF_Descriptors::getTanh_Pnl(const int ikernel, std::vector<double> &ppnl, std::vector<double> &rtanh_pnl)
{
    this->tanh(ppnl, rtanh_pnl, this->chi_pnl[ikernel]);
}

// tanh(qnl)
void Cal_MLKEDF_Descriptors::getTanh_Qnl(const int ikernel, std::vector<double> &pqnl, std::vector<double> &rtanh_qnl)
{
    this->tanh(pqnl, rtanh_qnl, this->chi_qnl[ikernel]);
}

// tanh(p)_nl
void Cal_MLKEDF_Descriptors::getTanhP_nl(const int ikernel, std::vector<double> &ptanhp, ModulePW::PW_Basis *pw_rho, std::vector<double> &rtanhp_nl)
{
    this->multiKernel(ikernel, ptanhp.data(), pw_rho, rtanhp_nl.data());
}

// tanh(q)_nl
void Cal_MLKEDF_Descriptors::getTanhQ_nl(const int ikernel, std::vector<double> &ptanhq, ModulePW::PW_Basis *pw_rho, std::vector<double> &rtanhq_nl)
{
    this->multiKernel(ikernel, ptanhq.data(), pw_rho, rtanhq_nl.data());
}

// (tanhxi)_nl
void Cal_MLKEDF_Descriptors::getTanhXi_nl(const int ikernel, std::vector<double> &ptanhxi, ModulePW::PW_Basis *pw_rho, std::vector<double> &rtanhxi_nl)
{
    this->multiKernel(ikernel, ptanhxi.data(), pw_rho, rtanhxi_nl.data());
}

void Cal_MLKEDF_Descriptors::getF_KS(
    psi::Psi<std::complex<double>> *psi,
    elecstate::ElecState *pelec,
    ModulePW::PW_Basis_K *pw_psi,
    ModulePW::PW_Basis *pw_rho,
    UnitCell& ucell,
    const std::vector<std::vector<double>> &nablaRho,
    std::vector<double> &rF,
    std::vector<double> &rpauli
)
{
    double *pauliED = new double[this->nx]; // Pauli Energy Density
    ModuleBase::GlobalFunc::ZEROS(pauliED, this->nx);

    double *pauliPot = new double[this->nx];
    ModuleBase::GlobalFunc::ZEROS(pauliPot, this->nx);

    std::complex<double> *wfcr = new std::complex<double>[this->nx];
    ModuleBase::GlobalFunc::ZEROS(wfcr, this->nx);

    double epsilonM = pelec->ekb(0,0);
    assert(PARAM.inp.nspin == 1);

    base_device::DEVICE_CPU* ctx;

    // calculate positive definite kinetic energy density
    for (int ik = 0; ik < psi->get_nk(); ++ik)
    {
        psi->fix_k(ik);
        int ikk = psi->get_current_k();
        assert(ikk == ik);
        int npw = psi->get_current_nbas();
        int nbands = psi->get_nbands();
        for (int ibnd = 0; ibnd < nbands; ++ibnd)
        {
            if (pelec->wg(ik, ibnd) < ModuleBase::threshold_wg) {
                continue;
            }

            pw_psi->recip_to_real(ctx, &psi->operator()(ibnd,0), wfcr, ik);
            const double w1 = pelec->wg(ik, ibnd) / ucell.omega;
            
            // output one wf, to check KS equation
            if (ik == 0 && ibnd == 0)
            {
                std::vector<double> wf_real = std::vector<double>(this->nx);
                std::vector<double> wf_imag = std::vector<double>(this->nx);
                for (int ir = 0; ir < this->nx; ++ir)
                {
                    wf_real[ir] = wfcr[ir].real();
                    wf_imag[ir] = wfcr[ir].imag();
                }
                const long unsigned cshape[] = {(long unsigned) this->nx}; // shape of container and containernl
            }

            if (w1 != 0.0)
            {
                // Find the energy of HOMO
                if (pelec->ekb(ik,ibnd) > epsilonM)
                {
                    epsilonM = pelec->ekb(ik,ibnd);
                }
                // The last term of Pauli potential
                for (int ir = 0; ir < pelec->charge->nrxx; ir++)
                {
                    pauliPot[ir] -= w1 * pelec->ekb(ik,ibnd) * norm(wfcr[ir]);
                }
            }

            for (int j = 0; j < 3; ++j)
            {
                ModuleBase::GlobalFunc::ZEROS(wfcr, pelec->charge->nrxx);
                for (int ig = 0; ig < npw; ig++)
                {
                    double fact
                        = pw_psi->getgpluskcar(ik, ig)[j] * ucell.tpiba;
                    wfcr[ig] = psi->operator()(ibnd, ig) * std::complex<double>(0.0, fact);
                }

                pw_psi->recip2real(wfcr, wfcr, ik);
                
                for (int ir = 0; ir < this->nx; ++ir)
                {
                    pauliED[ir] += w1 * norm(wfcr[ir]); // actually, here should be w1/2 * norm(wfcr[ir]), but we multiply 2 to convert Ha to Ry.
                }
            }
        }
    }

    for (int j = 0; j < 3; ++j)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            pauliED[ir] -= nablaRho[j][ir] * nablaRho[j][ir] / (8. * pelec->charge->rho[0][ir]) * 2.; // convert Ha to Ry.
        }
    }

    for (int ir = 0; ir < this->nx; ++ir)
    {
        rF[ir] = pauliED[ir] / (this->cTF * std::pow(pelec->charge->rho[0][ir], 5./3.));
        rpauli[ir] = (pauliED[ir] + pauliPot[ir])/pelec->charge->rho[0][ir] + epsilonM;
    }
}

void Cal_MLKEDF_Descriptors::getNablaRho(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<std::vector<double>> &rnablaRho)
{
    std::complex<double> *recipRho = new std::complex<double>[pw_rho->npw];
    std::complex<double> *recipNablaRho = new std::complex<double>[pw_rho->npw];
    pw_rho->real2recip(prho[0], recipRho);
    
    std::complex<double> img(0.0, 1.0);
    for (int j = 0; j < 3; ++j)
    {
        for (int ip = 0; ip < pw_rho->npw; ++ip)
        {
            recipNablaRho[ip] = img * pw_rho->gcar[ip][j] * recipRho[ip] * pw_rho->tpiba;
        }
        pw_rho->recip2real(recipNablaRho, rnablaRho[j].data());
    }

    delete[] recipRho;
    delete[] recipNablaRho;
}

}