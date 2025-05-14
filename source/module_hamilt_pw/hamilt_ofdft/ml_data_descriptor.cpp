#ifdef __MLKEDF

#include "ml_data.h"

void ML_data::getGamma(const double * const *prho, std::vector<double> &rgamma)
{
    for(int ir = 0; ir < this->nx; ++ir)
    {
        rgamma[ir] = std::pow(prho[0][ir]/this->rho0, 1./3.);
    }
}

void ML_data::getP(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<std::vector<double>> &pnablaRho, std::vector<double> &rp)
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

void ML_data::getQ(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rq)
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

void ML_data::getGammanl(const int ikernel, std::vector<double> &pgamma, ModulePW::PW_Basis *pw_rho, std::vector<double> &rgammanl)
{
    this->multiKernel(ikernel, pgamma.data(), pw_rho, rgammanl.data());
}

void ML_data::getPnl(const int ikernel, std::vector<double> &pp, ModulePW::PW_Basis *pw_rho, std::vector<double> &rpnl)
{
    this->multiKernel(ikernel, pp.data(), pw_rho, rpnl.data());
}

void ML_data::getQnl(const int ikernel, std::vector<double> &pq, ModulePW::PW_Basis *pw_rho, std::vector<double> &rqnl)
{
    this->multiKernel(ikernel, pq.data(), pw_rho, rqnl.data());
}

// xi = gammanl/gamma
void ML_data::getXi(std::vector<double> &pgamma, std::vector<double> &pgammanl, std::vector<double> &rxi)
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
void ML_data::getTanhXi(const int ikernel, std::vector<double> &pgamma, std::vector<double> &pgammanl, std::vector<double> &rtanhxi)
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
void ML_data::getTanhP(std::vector<double> &pp, std::vector<double> &rtanhp)
{
    this->tanh(pp, rtanhp, this->chi_p);
}

// tanh(q)
void ML_data::getTanhQ(std::vector<double> &pq, std::vector<double> &rtanhq)
{
    this->tanh(pq, rtanhq, this->chi_q);
}

// tanh(pnl)
void ML_data::getTanh_Pnl(const int ikernel, std::vector<double> &ppnl, std::vector<double> &rtanh_pnl)
{
    this->tanh(ppnl, rtanh_pnl, this->chi_pnl[ikernel]);
}

// tanh(qnl)
void ML_data::getTanh_Qnl(const int ikernel, std::vector<double> &pqnl, std::vector<double> &rtanh_qnl)
{
    this->tanh(pqnl, rtanh_qnl, this->chi_qnl[ikernel]);
}

// tanh(p)_nl
void ML_data::getTanhP_nl(const int ikernel, std::vector<double> &ptanhp, ModulePW::PW_Basis *pw_rho, std::vector<double> &rtanhp_nl)
{
    this->multiKernel(ikernel, ptanhp.data(), pw_rho, rtanhp_nl.data());
}

// tanh(q)_nl
void ML_data::getTanhQ_nl(const int ikernel, std::vector<double> &ptanhq, ModulePW::PW_Basis *pw_rho, std::vector<double> &rtanhq_nl)
{
    this->multiKernel(ikernel, ptanhq.data(), pw_rho, rtanhq_nl.data());
}

// (tanhxi)_nl
void ML_data::getTanhXi_nl(const int ikernel, std::vector<double> &ptanhxi, ModulePW::PW_Basis *pw_rho, std::vector<double> &rtanhxi_nl)
{
    this->multiKernel(ikernel, ptanhxi.data(), pw_rho, rtanhxi_nl.data());
}

void ML_data::getPauli_WT(KEDF_WT &wt, KEDF_TF &tf, const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rpauli)
{
    ModuleBase::matrix potential(1, this->nx, true);

    tf.tf_potential(prho, potential);
    wt.wt_potential(prho, pw_rho, potential);

    for (int ir = 0; ir < this->nx; ++ir){
        rpauli[ir] = potential(0, ir);
    }
}

void ML_data::getF_WT(KEDF_WT &wt, KEDF_TF &tf, const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rF)
{
    double wtden = 0.;
    double tfden = 0.;
    for (int ir = 0; ir < this->nx; ++ir)
    {
        wtden = wt.get_energy_density(prho, 0, ir, pw_rho);
        tfden = tf.get_energy_density(prho, 0, ir);
        rF[ir] = 1. + wtden/tfden;
        // if (wtden < 0) std::cout << wtden/tfden << std::endl;
    }
}

void ML_data::getF_KS1(
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
                    wfcr[ig] = psi->operator()(ibnd, ig) * complex<double>(0.0, fact);
                }

                pw_psi->recip2real(wfcr, wfcr, ik);
                
                for (int ir = 0; ir < this->nx; ++ir)
                {
                    pauliED[ir] += w1 * norm(wfcr[ir]); // actually, here should be w1/2 * norm(wfcr[ir]), but we multiply 2 to convert Ha to Ry.
                }
            }
        }
    }

    std::cout << "(1) epsilon max = " << epsilonM << std::endl;
    // calculate the positive definite vW energy density
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

void ML_data::getF_KS2(
    psi::Psi<std::complex<double>> *psi,
    elecstate::ElecState *pelec,
    ModulePW::PW_Basis_K *pw_psi,
    ModulePW::PW_Basis *pw_rho,
    UnitCell& ucell,
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

    std::complex<double> *wfcg = nullptr;
    std::complex<double> *LapWfcr = new std::complex<double>[this->nx];
    ModuleBase::GlobalFunc::ZEROS(LapWfcr, this->nx);

    double epsilonM = pelec->ekb(0,0);
    assert(PARAM.inp.nspin == 1);

    base_device::DEVICE_CPU* ctx;

    // calculate kinetic energy density
    for (int ik = 0; ik < psi->get_nk(); ++ik)
    {
        psi->fix_k(ik);
        int ikk = psi->get_current_k();
        assert(ikk == ik);
        int npw = psi->get_current_nbas();
        int nbands = psi->get_nbands();
        delete[] wfcg;
        wfcg = new std::complex<double>[npw];
        for (int ibnd = 0; ibnd < nbands; ++ibnd)
        {
            if (pelec->wg(ik, ibnd) < ModuleBase::threshold_wg) {
                continue;
            }

            pw_psi->recip_to_real(ctx, &psi->operator()(ibnd,0), wfcr, ik);
            const double w1 = pelec->wg(ik, ibnd) / ucell.omega;

            if (pelec->ekb(ik,ibnd) > epsilonM)
            {
                epsilonM = pelec->ekb(ik,ibnd);
            }
            for (int ir = 0; ir < pelec->charge->nrxx; ir++)
            {
                pauliPot[ir] -= w1 * pelec->ekb(ik,ibnd) * norm(wfcr[ir]);
            }

            ModuleBase::GlobalFunc::ZEROS(wfcg, npw);
            for (int ig = 0; ig < npw; ig++)
            {
                double fact = pw_psi->getgk2(ik, ig) * ucell.tpiba2;
                wfcg[ig] = - psi->operator()(ibnd, ig) * fact;
            }

            pw_psi->recip2real(wfcg, LapWfcr, ik);
            
            for (int ir = 0; ir < this->nx; ++ir)
            {
                pauliED[ir] += - w1 * (conj(wfcr[ir]) * LapWfcr[ir]).real(); // actually, here should be w1/2 * norm(wfcr[ir]), but we multiply 2 to convert Ha to Ry.
            }
        }
    }

    std::cout << "(2) epsilon max = " << epsilonM << std::endl;
    // calculate the positive definite vW energy density
    double *phi = new double[this->nx];
    double *LapPhi = new double[this->nx];
    for (int ir = 0; ir < this->nx; ++ir)
    {
        phi[ir] = sqrt(pelec->charge->rho[0][ir]);
    }
    this->Laplacian(phi, pw_rho, LapPhi);

    for (int ir = 0; ir < this->nx; ++ir)
    {
        pauliED[ir] -= - phi[ir] * LapPhi[ir]; // convert Ha to Ry.
    }

    for (int ir = 0; ir < this->nx; ++ir)
    {
        rF[ir] = pauliED[ir] / (this->cTF * std::pow(pelec->charge->rho[0][ir], 5./3.));
        rpauli[ir] = (pauliED[ir] + pauliPot[ir])/pelec->charge->rho[0][ir] + epsilonM;
    }
}

void ML_data::getNablaRho(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<std::vector<double>> &rnablaRho)
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
#endif