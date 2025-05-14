#ifdef __MLKEDF

#include "kedf_ml.h"

#include "module_base/parallel_reduce.h"
#include "module_base/global_function.h"

/**
 * @brief Calculate the Pauli potential of ML KEDF
 * 
 * @param prho charge density
 * @param pw_rho PW_Basis
 * @param rpotential rpotential => rpotential + V_{ML}
 */
void KEDF_ML::get_potential_(const double * const * prho, ModulePW::PW_Basis *pw_rho, ModuleBase::matrix &rpotential)
{
    // get potential
    ModuleBase::timer::tick("KEDF_ML", "Pauli Potential");

    std::vector<double> pauli_potential(this->nx, 0.);

    if (this->ml_gammanl){
        this->potGammanlTerm(prho, pw_rho, pauli_potential);
    }
    if (this->ml_xi){
        this->potXinlTerm(prho, pw_rho, pauli_potential);
    }
    if (this->ml_tanhxi){
        this->potTanhxinlTerm(prho, pw_rho, pauli_potential);
    }
    if (this->ml_tanhxi_nl){
        this->potTanhxi_nlTerm(prho, pw_rho, pauli_potential);
    }
    if (this->ml_p || this->ml_pnl){
        this->potPPnlTerm(prho, pw_rho, pauli_potential);
    }
    if (this->ml_q || this->ml_qnl){
        this->potQQnlTerm(prho, pw_rho, pauli_potential);
    }
    if (this->ml_tanh_pnl){
        this->potTanhpTanh_pnlTerm(prho, pw_rho, pauli_potential);
    }
    if (this->ml_tanh_qnl){
        this->potTanhqTanh_qnlTerm(prho, pw_rho, pauli_potential);
    }
    if ((this->ml_tanhp || this->ml_tanhp_nl) && !this->ml_tanh_pnl){
        this->potTanhpTanhp_nlTerm(prho, pw_rho, pauli_potential);
    }
    if ((this->ml_tanhq || this->ml_tanhq_nl) && !this->ml_tanh_qnl){
        this->potTanhqTanhq_nlTerm(prho, pw_rho, pauli_potential);
    }

    for (int ir = 0; ir < this->nx; ++ir)
    {
        pauli_potential[ir] += this->cTF * std::pow(prho[0][ir], 5./3.) / prho[0][ir] *
                      (5./3. * this->enhancement_cpu_ptr[ir] + this->potGammaTerm(ir) + this->potPTerm1(ir) + this->potQTerm1(ir)
                      + this->potXiTerm1(ir) + this->potTanhxiTerm1(ir) + this->potTanhpTerm1(ir) + this->potTanhqTerm1(ir));
        rpotential(0, ir) += pauli_potential[ir];
    }
    ModuleBase::timer::tick("KEDF_ML", "Pauli Potential");
}


double KEDF_ML::potGammaTerm(int ir)
{
    return (this->ml_gamma) ? 1./3. * gamma[ir] * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["gamma"][0]] : 0.;
}

double KEDF_ML::potPTerm1(int ir)
{
    return (this->ml_p) ? - 8./3. * p[ir] * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["p"][0]] : 0.;
}

double KEDF_ML::potQTerm1(int ir)
{
    return (this->ml_q) ? - 5./3. * q[ir] * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["q"][0]] : 0.;
}

double KEDF_ML::potXiTerm1(int ir)
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

double KEDF_ML::potTanhxiTerm1(int ir)
{
    double result = 0.;
    for (int ik = 0; ik < this->descriptor2kernel["tanhxi"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["tanhxi"][ik];
        int d2i = this->descriptor2index["tanhxi"][ik];
        result += -1./3. * xi[d2k][ir] * this->ml_data->dtanh(this->tanhxi[d2k][ir], this->chi_xi[d2k])
                                    * this->gradient_cpu_ptr[ir * this->ninput + d2i];
    }
    return result;
}

double KEDF_ML::potTanhpTerm1(int ir)
{
    return (this->ml_tanhp) ? - 8./3. * p[ir] * this->ml_data->dtanh(this->tanhp[ir], this->chi_p)
                                 * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["tanhp"][0]] : 0.;
}

double KEDF_ML::potTanhqTerm1(int ir)
{
    return (this->ml_tanhq) ? - 5./3. * q[ir] * this->ml_data->dtanh(this->tanhq[ir], this->chi_q)
                                 * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["tanhq"][0]] : 0.;
}

void KEDF_ML::potGammanlTerm(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rGammanlTerm)
{
    double *dFdgammanl = new double[this->nx];
    for (int ik = 0; ik < this->descriptor2kernel["gammanl"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["gammanl"][ik];
        int d2i = this->descriptor2index["gammanl"][ik];
        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdgammanl[ir] = this->cTF * std::pow(prho[0][ir], 5./3.) * this->gradient_cpu_ptr[ir * this->ninput + d2i];
        }
        this->ml_data->multiKernel(d2k, dFdgammanl, pw_rho, dFdgammanl);
        for (int ir = 0; ir < this->nx; ++ir)
        {
            rGammanlTerm[ir] += 1./3. * this->gamma[ir] / prho[0][ir] * dFdgammanl[ir];
        }
    }
    delete[] dFdgammanl;
}

void KEDF_ML::potXinlTerm(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rXinlTerm)
{
    double *dFdxi = new double[this->nx];
    for (int ik = 0; ik < this->descriptor2kernel["xi"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["xi"][ik];
        int d2i = this->descriptor2index["xi"][ik];
        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdxi[ir] = this->cTF * std::pow(prho[0][ir], 4./3.) * this->gradient_cpu_ptr[ir * this->ninput + d2i];
        }
        this->ml_data->multiKernel(d2k, dFdxi, pw_rho, dFdxi);
        for (int ir = 0; ir < this->nx; ++ir)
        {
            rXinlTerm[ir] += 1./3. * std::pow(prho[0][ir], -2./3.) * dFdxi[ir];
        }
    }
    delete[] dFdxi;
}

void KEDF_ML::potTanhxinlTerm(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhxinlTerm)
{
    double *dFdtanhxi = new double[this->nx];
    for (int ik = 0; ik < this->descriptor2kernel["tanhxi"].size(); ++ik)
    {
        int d2k = this->descriptor2kernel["tanhxi"][ik];
        int d2i = this->descriptor2index["tanhxi"][ik];
        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdtanhxi[ir] = this->cTF * std::pow(prho[0][ir], 4./3.) * this->ml_data->dtanh(this->tanhxi[d2k][ir], this->chi_xi[d2k])
                        * this->gradient_cpu_ptr[ir * this->ninput + d2i];
        }
        this->ml_data->multiKernel(d2k, dFdtanhxi, pw_rho, dFdtanhxi);
        for (int ir = 0; ir < this->nx; ++ir)
        {
            rTanhxinlTerm[ir] += 1./3. * std::pow(prho[0][ir], -2./3.) * dFdtanhxi[ir];
        }
    }
    delete[] dFdtanhxi;
}

void KEDF_ML::potTanhxi_nlTerm(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhxi_nlTerm)
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
            dFdtanhxi_nl[ir] = this->cTF * std::pow(prho[0][ir], 5./3.) * this->gradient_cpu_ptr[ir * this->ninput + d2i];
        }
        this->ml_data->multiKernel(d2k, dFdtanhxi_nl, pw_rho, dFdtanhxi_nl);
        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdtanhxi_nl[ir] *= this->ml_data->dtanh(this->tanhxi[d2k][ir], this->chi_xi[d2k]) / std::pow(prho[0][ir], 1./3.);
        }
        this->ml_data->multiKernel(d2k, dFdtanhxi_nl, pw_rho, dFdtanhxi_nl_nl);
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
void KEDF_ML::potPPnlTerm(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rPPnlTerm)
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
            dFdpnl[ir] = this->cTF * std::pow(prho[0][ir], 5./3.) * this->gradient_cpu_ptr[ir * this->ninput + d2i];
        }
        this->ml_data->multiKernel(d2k, dFdpnl, pw_rho, dFdpnl);
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
            tempP[i][ir] = (this->ml_p)? - 3./20. * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["p"][0]] * nablaRho[i][ir] / prho[0][ir] * /*Ha2Ry*/ 2. : 0.;
            if (this->ml_pnl)
            {
                tempP[i][ir] += - this->pqcoef * 2. * this->nablaRho[i][ir] / std::pow(prho[0][ir], 8./3.) * dFdpnl_tot[ir];
            }
        }
    }
    this->ml_data->divergence(tempP, pw_rho, result.data());

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

void KEDF_ML::potQQnlTerm(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rQQnlTerm)
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
            dFdqnl[ir] = this->cTF * std::pow(prho[0][ir], 5./3.) * this->gradient_cpu_ptr[ir * this->ninput + d2i];
        }
        this->ml_data->multiKernel(d2k, dFdqnl, pw_rho, dFdqnl);
        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdqnl_tot[ir] += dFdqnl[ir];
        }
    }
    delete[] dFdqnl;

    double * tempQ = new double[this->nx];
    for (int ir = 0; ir < this->nx; ++ir)
    {
        tempQ[ir] = (this->ml_q)? 3./40. * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["q"][0]] * /*Ha2Ry*/ 2. : 0.;
        if (this->ml_qnl)
        {
            tempQ[ir] += this->pqcoef / std::pow(prho[0][ir], 5./3.) * dFdqnl_tot[ir];
        }
    }
    this->ml_data->Laplacian(tempQ, pw_rho, result.data());

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

void KEDF_ML::potTanhpTanh_pnlTerm(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhpTanh_pnlTerm)
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
            dFdpnl[ir] = this->cTF * std::pow(prho[0][ir], 5./3.) * this->ml_data->dtanh(this->tanh_pnl[d2k][ir], this->chi_pnl[d2k])
                         * this->gradient_cpu_ptr[ir * this->ninput + d2i];
        }
        this->ml_data->multiKernel(d2k, dFdpnl, pw_rho, dFdpnl);
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
            tempP[i][ir] = (this->ml_tanhp)? - 3./20. * this->ml_data->dtanh(this->tanhp[ir], this->chi_p)
                           * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["tanhp"][0]] * nablaRho[i][ir] / prho[0][ir] * /*Ha2Ry*/ 2. : 0.;
            if (this->ml_tanh_pnl)
            {
                tempP[i][ir] += - this->pqcoef * 2. * this->nablaRho[i][ir] / std::pow(prho[0][ir], 8./3.) * dFdpnl_tot[ir];
            }
        }
    }
    this->ml_data->divergence(tempP, pw_rho, result.data());

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

void KEDF_ML::potTanhqTanh_qnlTerm(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhqTanh_qnlTerm)
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
            dFdqnl[ir] = this->cTF * std::pow(prho[0][ir], 5./3.) * this->ml_data->dtanh(this->tanh_qnl[d2k][ir], this->chi_qnl[d2k])
                         * this->gradient_cpu_ptr[ir * this->ninput + d2i];
        }
        this->ml_data->multiKernel(d2k, dFdqnl, pw_rho, dFdqnl);
        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdqnl_tot[ir] += dFdqnl[ir];
        }
    }
    delete[] dFdqnl;

    double * tempQ = new double[this->nx];
    for (int ir = 0; ir < this->nx; ++ir)
    {
        tempQ[ir] = (this->ml_tanhq)? 3./40. * this->ml_data->dtanh(this->tanhq[ir], this->chi_q)
                    * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["tanhq"][0]] * /*Ha2Ry*/ 2. : 0.;
        if (this->ml_tanh_qnl)
        {
            tempQ[ir] += this->pqcoef / std::pow(prho[0][ir], 5./3.) * dFdqnl_tot[ir];
        }
    }
    this->ml_data->Laplacian(tempQ, pw_rho, result.data());

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
void KEDF_ML::potTanhpTanhp_nlTerm(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhpTanhp_nlTerm)
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
            dFdpnl[ir] = this->cTF * std::pow(prho[0][ir], 5./3.)
                         * this->gradient_cpu_ptr[ir * this->ninput + d2i];
        }
        this->ml_data->multiKernel(d2k, dFdpnl, pw_rho, dFdpnl);
        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdpnl_tot[ir] += this->ml_data->dtanh(this->tanhp[ir], this->chi_p) * dFdpnl[ir];
        }
    }
    delete[] dFdpnl;

    double ** tempP = new double*[3];
    for (int i = 0; i < 3; ++i)
    {
        tempP[i] = new double[this->nx];
        for (int ir = 0; ir < this->nx; ++ir)
        {
            tempP[i][ir] = (this->ml_tanhp)? - 3./20. * this->ml_data->dtanh(this->tanhp[ir], this->chi_p)
                           * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["tanhp"][0]] * nablaRho[i][ir] / prho[0][ir] * /*Ha2Ry*/ 2. : 0.;
            if (this->ml_tanhp_nl)
            {
                tempP[i][ir] += - this->pqcoef * 2. * this->nablaRho[i][ir] / std::pow(prho[0][ir], 8./3.) * dFdpnl_tot[ir];
            }
        }
    }
    this->ml_data->divergence(tempP, pw_rho, result.data());

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

void KEDF_ML::potTanhqTanhq_nlTerm(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhqTanhq_nlTerm)
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
            dFdqnl[ir] = this->cTF * std::pow(prho[0][ir], 5./3.)
                         * this->gradient_cpu_ptr[ir * this->ninput + d2i];
        }
        this->ml_data->multiKernel(d2k, dFdqnl, pw_rho, dFdqnl);
        for (int ir = 0; ir < this->nx; ++ir)
        {
            dFdqnl_tot[ir] += dFdqnl[ir];
        }
    }
    delete[] dFdqnl;

    double * tempQ = new double[this->nx];
    for (int ir = 0; ir < this->nx; ++ir)
    {
        tempQ[ir] = (this->ml_tanhq)? 3./40. * this->ml_data->dtanh(this->tanhq[ir], this->chi_q)
                    * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["tanhq"][0]] * /*Ha2Ry*/ 2. : 0.;
        if (this->ml_tanhq_nl)
        {
            dFdqnl_tot[ir] *= this->ml_data->dtanh(this->tanhq[ir], this->chi_q);
            tempQ[ir] += this->pqcoef / std::pow(prho[0][ir], 5./3.) * dFdqnl_tot[ir];
        }
    }
    this->ml_data->Laplacian(tempQ, pw_rho, result.data());

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