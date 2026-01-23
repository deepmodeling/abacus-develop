#include "ml_base.h"

#ifdef __MLALGO

double ML_Base::potGammaTerm(int ir)
{
    return (this->ml_gamma) ? 1./3. * gamma[ir] * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["gamma"][0]] : 0.;
}
double ML_Base::potPTerm1(int ir)
{
    return (this->ml_p) ? - 8./3. * p[ir] * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["p"][0]] : 0.;
}
double ML_Base::potQTerm1(int ir)
{
    return (this->ml_q) ? - 5./3. * q[ir] * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["q"][0]] : 0.;
}
double ML_Base::potXiTerm1(int ir)
{
    double result = 0.;
    for (int ik = 0; ik < this->descriptor2kernel["xi"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["xi"][ik];
        int d2i = this->descriptor2index["xi"][ik];
        result += -1./3. * xi[d2k][ir] * this->gradient_cpu_ptr[ir * this->ninput + d2i];
    }
    return result;
}
double ML_Base::potTanhxiTerm1(int ir)
{
    double result = 0.;
    for (int ik = 0; ik < this->descriptor2kernel["tanhxi"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["tanhxi"][ik];
        int d2i = this->descriptor2index["tanhxi"][ik];
        result += -1./3. * xi[d2k][ir] * this->cal_tool->dtanh(this->tanhxi[d2k][ir], this->chi_xi[d2k])
                                    * this->gradient_cpu_ptr[ir * this->ninput + d2i];
    }
    return result;
}
double ML_Base::potTanhpTerm1(int ir)
{
    return (this->ml_tanhp) ? - 8./3. * p[ir] * this->cal_tool->dtanh(this->tanhp[ir], this->chi_p)
                                 * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["tanhp"][0]] : 0.;
}
double ML_Base::potTanhqTerm1(int ir)
{
    return (this->ml_tanhq) ? - 5./3. * q[ir] * this->cal_tool->dtanh(this->tanhq[ir], this->chi_q)
                                 * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["tanhq"][0]] : 0.;
}

// Implementations of nl terms using energy_prefactor/exponent logic
void ML_Base::potGammanlTerm(const double * const *prho, const std::vector<double> &tau_lda, const ModulePW::PW_Basis *pw_rho, std::vector<double> &rGammanlTerm)
{
    double *dFdgammanl = new double[this->nx];
    for (int ik = 0; ik < this->descriptor2kernel["gammanl"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["gammanl"][ik];
        int d2i = this->descriptor2index["gammanl"][ik];
        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdgammanl[ir] = tau_lda[ir] * this->gradient_cpu_ptr[ir * this->ninput + d2i];
        }
        this->cal_tool->multiKernel(d2k, dFdgammanl, pw_rho, dFdgammanl);
        for (int ir = 0; ir < this->nx; ++ir)
        {
            rGammanlTerm[ir] += 1./3. * this->gamma[ir] / prho[0][ir] * dFdgammanl[ir];
        }
    }
    delete[] dFdgammanl;
}

void ML_Base::potXinlTerm(const double * const *prho, const std::vector<double> &tau_lda, const ModulePW::PW_Basis *pw_rho, std::vector<double> &rXinlTerm)
{
    double *dFdxi = new double[this->nx];
    for (int ik = 0; ik < this->descriptor2kernel["xi"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["xi"][ik];
        int d2i = this->descriptor2index["xi"][ik];
        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdxi[ir] = tau_lda[ir] / std::pow(prho[0][ir], 1./3.) * this->gradient_cpu_ptr[ir * this->ninput + d2i];
        }
        this->cal_tool->multiKernel(d2k, dFdxi, pw_rho, dFdxi);
        for (int ir = 0; ir < this->nx; ++ir)
        {
            rXinlTerm[ir] += 1./3. * std::pow(prho[0][ir], -2./3.) * dFdxi[ir];
        }
    }
    delete[] dFdxi;
}

void ML_Base::potTanhxinlTerm(const double * const *prho, const std::vector<double> &tau_lda, const ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhxinlTerm)
{
    double *dFdtanhxi = new double[this->nx];
    for (int ik = 0; ik < this->descriptor2kernel["tanhxi"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["tanhxi"][ik];
        int d2i = this->descriptor2index["tanhxi"][ik];
        for (int ir = 0; ir < this->nx; ++ir)
        {
             dFdtanhxi[ir] = tau_lda[ir] / std::pow(prho[0][ir], 1./3.) * this->cal_tool->dtanh(this->tanhxi[d2k][ir], this->chi_xi[d2k])
                        * this->gradient_cpu_ptr[ir * this->ninput + d2i];
        }
        this->cal_tool->multiKernel(d2k, dFdtanhxi, pw_rho, dFdtanhxi);
        for (int ir = 0; ir < this->nx; ++ir)
        {
            rTanhxinlTerm[ir] += 1./3. * std::pow(prho[0][ir], -2./3.) * dFdtanhxi[ir];
        }
    }
    delete[] dFdtanhxi;
}

void ML_Base::potTanhxi_nlTerm(const double * const *prho, const std::vector<double> &tau_lda, const ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhxi_nlTerm)
{
    double *dFdtanhxi_nl = new double[this->nx];
    double *dFdtanhxi_nl_nl = new double[this->nx];
    std::vector<double> result(this->nx, 0.);
    for (int ik = 0; ik < this->descriptor2kernel["tanhxi_nl"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["tanhxi_nl"][ik];
        int d2i = this->descriptor2index["tanhxi_nl"][ik];
        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdtanhxi_nl[ir] = tau_lda[ir] * this->gradient_cpu_ptr[ir * this->ninput + d2i];
        }
        this->cal_tool->multiKernel(d2k, dFdtanhxi_nl, pw_rho, dFdtanhxi_nl);
        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdtanhxi_nl[ir] *= this->cal_tool->dtanh(this->tanhxi[d2k][ir], this->chi_xi[d2k]) / std::pow(prho[0][ir], 1./3.);
        }
        this->cal_tool->multiKernel(d2k, dFdtanhxi_nl, pw_rho, dFdtanhxi_nl_nl);
        for (int ir = 0; ir < this->nx; ++ir)
        {
            result[ir] += dFdtanhxi_nl_nl[ir] - dFdtanhxi_nl[ir] * this->xi[d2k][ir];
        }
    }
    for (int ir = 0; ir < this->nx; ++ir)
    {
        rTanhxi_nlTerm[ir] += result[ir] * 1./3. * std::pow(prho[0][ir], -2./3.);
    }
    delete[] dFdtanhxi_nl;
    delete[] dFdtanhxi_nl_nl;
}

// get contribution of p and pnl
void ML_Base::potPPnlTerm(const double * const *prho, const std::vector<double> &tau_lda, const ModulePW::PW_Basis *pw_rho, std::vector<double> &rPPnlTerm)
{
    double *dFdpnl = new double[this->nx];
    std::vector<double> dFdpnl_tot(this->nx, 0.);
    std::vector<double> result(this->nx, 0.);
    for (int ik = 0; ik < this->descriptor2index["pnl"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["pnl"][ik];
        int d2i = this->descriptor2index["pnl"][ik];

        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdpnl[ir] = tau_lda[ir] * this->gradient_cpu_ptr[ir * this->ninput + d2i];
        }
        this->cal_tool->multiKernel(d2k, dFdpnl, pw_rho, dFdpnl);
        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdpnl_tot[ir] += dFdpnl[ir];
        }
    }
    delete[] dFdpnl;

    double ** tempP = new double*[3];
    for (int i = 0; i < 3; ++i)
    {
        tempP[i] = new double[this->nx];
        for (int ir = 0; ir < this->nx; ++ir)
        {
            tempP[i][ir] = (this->ml_p)? - this->pqcoef * 2. * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["p"][0]] * this->nablaRho[i][ir] * tau_lda[ir] / std::pow(prho[0][ir], 8./3.): 0.;
            if (this->ml_pnl)
            {
                tempP[i][ir] += - this->pqcoef * 2. * this->nablaRho[i][ir] / std::pow(prho[0][ir], 8./3.) * dFdpnl_tot[ir];
            }
        }
    }
    this->cal_tool->divergence(tempP, pw_rho, result.data());

    if (this->ml_pnl)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            result[ir] += -8./3. * this->p[ir]/prho[0][ir] * dFdpnl_tot[ir];
        }
    }

    for (int ir = 0; ir < this->nx; ++ir)
    {
        rPPnlTerm[ir] += result[ir];
    }

    for (int i = 0; i < 3; ++i)
    { 
        delete[] tempP[i];
    }
    delete[] tempP;
}


void ML_Base::potQQnlTerm(const double * const *prho, const std::vector<double> &tau_lda, const ModulePW::PW_Basis *pw_rho, std::vector<double> &rQQnlTerm)
{
    double *dFdqnl = new double[this->nx];
    std::vector<double> dFdqnl_tot(this->nx, 0.);
    std::vector<double> result(this->nx, 0.);
    for (int ik = 0; ik < this->descriptor2index["qnl"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["qnl"][ik];
        int d2i = this->descriptor2index["qnl"][ik];

        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdqnl[ir] = tau_lda[ir] * this->gradient_cpu_ptr[ir * this->ninput + d2i];
        }
        this->cal_tool->multiKernel(d2k, dFdqnl, pw_rho, dFdqnl);
        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdqnl_tot[ir] += dFdqnl[ir];
        }
    }
    delete[] dFdqnl;

    double * tempQ = new double[this->nx];
    for (int ir = 0; ir < this->nx; ++ir)
    {
        tempQ[ir] = (this->ml_q)? this->pqcoef * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["q"][0]] * tau_lda[ir] / std::pow(prho[0][ir], 5./3.) : 0.;
        if (this->ml_qnl)
        {
            tempQ[ir] += this->pqcoef / std::pow(prho[0][ir], 5./3.) * dFdqnl_tot[ir];
        }
    }
    this->cal_tool->Laplacian(tempQ, pw_rho, result.data());
    if (this->ml_qnl)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            result[ir] += -5./3. * this->q[ir] / prho[0][ir] * dFdqnl_tot[ir];
        }
    }

    for (int ir = 0; ir < this->nx; ++ir)
    {
        rQQnlTerm[ir] += result[ir];
    }
    delete[] tempQ;
}


void ML_Base::potTanhpTanh_pnlTerm(const double * const *prho, const std::vector<double> &tau_lda, const ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhpTanh_pnlTerm)
{
    // Note we assume that tanhp_nl and tanh_pnl will NOT be used together.
    double *dFdpnl = new double[this->nx];
    std::vector<double> dFdpnl_tot(this->nx, 0.);
    std::vector<double> result(this->nx, 0.);
    for (int ik = 0; ik < this->descriptor2kernel["tanh_pnl"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["tanh_pnl"][ik];
        int d2i = this->descriptor2index["tanh_pnl"][ik];

        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdpnl[ir] = tau_lda[ir] * this->cal_tool->dtanh(this->tanh_pnl[d2k][ir], this->chi_pnl[d2k])
                         * this->gradient_cpu_ptr[ir * this->ninput + d2i];
        }
        this->cal_tool->multiKernel(d2k, dFdpnl, pw_rho, dFdpnl);
        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdpnl_tot[ir] += dFdpnl[ir];
        }
    }
    delete[] dFdpnl;

    double ** tempP = new double*[3];
    for (int i = 0; i < 3; ++i)
    {
        tempP[i] = new double[this->nx];
        for (int ir = 0; ir < this->nx; ++ir)
        {
            tempP[i][ir] = (this->ml_tanhp)? - this->pqcoef * 2. * this->cal_tool->dtanh(this->tanhp[ir], this->chi_p)
                           * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["tanhp"][0]] * this->nablaRho[i][ir] * tau_lda[ir] / std::pow(prho[0][ir], 8./3.) : 0.;
            if (this->ml_tanh_pnl)
            {
                tempP[i][ir] += - this->pqcoef * 2. * this->nablaRho[i][ir] / std::pow(prho[0][ir], 8./3.) * dFdpnl_tot[ir];
            }
        }
    }
    this->cal_tool->divergence(tempP, pw_rho, result.data());

    if (this->ml_tanh_pnl)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            result[ir] += -8./3. * this->p[ir]/prho[0][ir] * dFdpnl_tot[ir];
        }
    }

    for (int ir = 0; ir < this->nx; ++ir)
    {
        rTanhpTanh_pnlTerm[ir] += result[ir];
    }
    for (int i = 0; i < 3; ++i) 
    { 
        delete[] tempP[i];
    }
    delete[] tempP;
}

void ML_Base::potTanhqTanh_qnlTerm(const double * const *prho, const std::vector<double> &tau_lda, const ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhqTanh_qnlTerm)
{
    // Note we assume that tanhq_nl and tanh_qnl will NOT be used together.
    double *dFdqnl = new double[this->nx];
    std::vector<double> dFdqnl_tot(this->nx, 0.);
    std::vector<double> result(this->nx, 0.);
    for (int ik = 0; ik < this->descriptor2kernel["tanh_qnl"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["tanh_qnl"][ik];
        int d2i = this->descriptor2index["tanh_qnl"][ik];

        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdqnl[ir] = tau_lda[ir] * this->cal_tool->dtanh(this->tanh_qnl[d2k][ir], this->chi_qnl[d2k])
                         * this->gradient_cpu_ptr[ir * this->ninput + d2i];
        }
        this->cal_tool->multiKernel(d2k, dFdqnl, pw_rho, dFdqnl);
        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdqnl_tot[ir] += dFdqnl[ir];
        }
    }
    delete[] dFdqnl;

    double * tempQ = new double[this->nx];
    for (int ir = 0; ir < this->nx; ++ir)
    {
        tempQ[ir] = (this->ml_tanhq)? this->pqcoef * this->cal_tool->dtanh(this->tanhq[ir], this->chi_q)
                    * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["tanhq"][0]] * tau_lda[ir] / std::pow(prho[0][ir], 5./3.) : 0.;
        if (this->ml_tanh_qnl)
        {
            tempQ[ir] += this->pqcoef / std::pow(prho[0][ir], 5./3.) * dFdqnl_tot[ir];
        }
    }
    this->cal_tool->Laplacian(tempQ, pw_rho, result.data());
    if (this->ml_tanh_qnl)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            result[ir] += -5./3. * this->q[ir] / prho[0][ir] * dFdqnl_tot[ir];
        }
    }

    for (int ir = 0; ir < this->nx; ++ir)
    {
        rTanhqTanh_qnlTerm[ir] += result[ir];
    }
    delete[] tempQ;
}

// Note we assume that tanhp_nl and tanh_pnl will NOT be used together.
void ML_Base::potTanhpTanhp_nlTerm(const double * const *prho, const std::vector<double> &tau_lda, const ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhpTanhp_nlTerm)
{
    double *dFdpnl = new double[this->nx];
    std::vector<double> dFdpnl_tot(this->nx, 0.);
    std::vector<double> result(this->nx, 0.);
    for (int ik = 0; ik < this->descriptor2kernel["tanhp_nl"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["tanhp_nl"][ik];
        int d2i = this->descriptor2index["tanhp_nl"][ik];

        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdpnl[ir] = tau_lda[ir]
                         * this->gradient_cpu_ptr[ir * this->ninput + d2i];
        }
        this->cal_tool->multiKernel(d2k, dFdpnl, pw_rho, dFdpnl);
        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdpnl_tot[ir] += this->cal_tool->dtanh(this->tanhp[ir], this->chi_p) * dFdpnl[ir];
        }
    }
    delete[] dFdpnl;

    double ** tempP = new double*[3];
    for (int i = 0; i < 3; ++i)
    {
        tempP[i] = new double[this->nx];
        for (int ir = 0; ir < this->nx; ++ir)
        {
            tempP[i][ir] = (this->ml_tanhp)? - this->pqcoef * 2. * this->cal_tool->dtanh(this->tanhp[ir], this->chi_p)
                           * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["tanhp"][0]] * this->nablaRho[i][ir] * tau_lda[ir] / std::pow(prho[0][ir], 8./3.) : 0.;
            if (this->ml_tanhp_nl)
            {
                tempP[i][ir] += - this->pqcoef * 2. * this->nablaRho[i][ir] / std::pow(prho[0][ir], 8./3.) * dFdpnl_tot[ir];
            }
        }
    }
    this->cal_tool->divergence(tempP, pw_rho, result.data());

    if (this->ml_tanhp_nl)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            result[ir] += -8./3. * this->p[ir]/prho[0][ir] * dFdpnl_tot[ir];
        }
    }

    for (int ir = 0; ir < this->nx; ++ir)
    {
        rTanhpTanhp_nlTerm[ir] += result[ir];
    }
    for (int i = 0; i < 3; ++i) 
    { 
        delete[] tempP[i];
    }
    delete[] tempP;
}

void ML_Base::potTanhqTanhq_nlTerm(const double * const *prho, const std::vector<double> &tau_lda, const ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhqTanhq_nlTerm)
{
    double *dFdqnl = new double[this->nx];
    std::vector<double> dFdqnl_tot(this->nx, 0.);
    std::vector<double> result(this->nx, 0.);
    for (int ik = 0; ik < this->descriptor2kernel["tanhq_nl"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["tanhq_nl"][ik];
        int d2i = this->descriptor2index["tanhq_nl"][ik];

        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdqnl[ir] = tau_lda[ir]
                         * this->gradient_cpu_ptr[ir * this->ninput + d2i];
        }
        this->cal_tool->multiKernel(d2k, dFdqnl, pw_rho, dFdqnl);
        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdqnl_tot[ir] += dFdqnl[ir];
        }
    }
    delete[] dFdqnl;

    double * tempQ = new double[this->nx];
    for (int ir = 0; ir < this->nx; ++ir)
    {
        tempQ[ir] = (this->ml_tanhq)? this->pqcoef * this->cal_tool->dtanh(this->tanhq[ir], this->chi_q)
                    * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["tanhq"][0]] * tau_lda[ir] / std::pow(prho[0][ir], 5./3.) : 0.;
        if (this->ml_tanhq_nl)
        {
            dFdqnl_tot[ir] *= this->cal_tool->dtanh(this->tanhq[ir], this->chi_q);
            tempQ[ir] += this->pqcoef / std::pow(prho[0][ir], 5./3.) * dFdqnl_tot[ir];
        }
    }
    this->cal_tool->Laplacian(tempQ, pw_rho, result.data());

    if (this->ml_tanhq_nl)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            result[ir] += -5./3. * this->q[ir] / prho[0][ir] * dFdqnl_tot[ir];
        }
    }

    for (int ir = 0; ir < this->nx; ++ir)
    {
        rTanhqTanhq_nlTerm[ir] += result[ir];
    }
    delete[] tempQ;
}
#endif