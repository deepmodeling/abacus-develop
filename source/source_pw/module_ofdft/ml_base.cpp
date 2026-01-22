#include "ml_base.h"
#include "source_base/parallel_reduce.h"
#include "source_base/global_function.h"
#include "source_pw/module_pwdft/global.h"
#include "npy.hpp"

#ifdef __MLALGO

ML_Base::ML_Base(){}

ML_Base::~ML_Base()
{
    if (this->cal_tool) delete this->cal_tool;
}

void ML_Base::set_device(std::string device_inpt)
{
    if (device_inpt == "cpu")
    {
        std::cout << "------------------- Running NN on CPU -------------------" << std::endl;
        this->device_type = torch::kCPU;
    }
    else if (device_inpt == "gpu")
    {
        if (torch::cuda::cudnn_is_available())
        {
            std::cout << "------------------- Running NN on GPU -------------------" << std::endl;
            this->device_type = torch::kCUDA;
        }
        else
        {
            std::cout << "--------------- Warning: GPU is unaviable ---------------" << std::endl;
            std::cout << "------------------- Running NN on CPU -------------------" << std::endl;
            this->device_type = torch::kCPU;
        }
    }
    this->device = torch::Device(this->device_type);
}

void ML_Base::updateInput(const double * const * prho, const ModulePW::PW_Basis *pw_rho)
{
    ModuleBase::timer::tick("ML_Base", "updateInput");
    if (this->gene_data_label["gamma"][0])
    {   
        this->cal_tool->getGamma(prho, this->gamma);
    }
    if (this->gene_data_label["p"][0])
    {
        this->cal_tool->getNablaRho(prho, pw_rho, this->nablaRho);
        this->cal_tool->getP(prho, pw_rho, this->nablaRho, this->p);
    }
    if (this->gene_data_label["q"][0])
    {
        this->cal_tool->getQ(prho, pw_rho, this->q);
    }
    if (this->gene_data_label["tanhp"][0])
    {
        this->cal_tool->getTanhP(this->p, this->tanhp);
    }
    if (this->gene_data_label["tanhq"][0])
    {
        this->cal_tool->getTanhQ(this->q, this->tanhq);
    }

    for (int ik = 0; ik < nkernel; ++ik)
    {
        if (this->gene_data_label["gammanl"][ik]){
            this->cal_tool->getGammanl(ik, this->gamma, pw_rho, this->gammanl[ik]);
        }
        if (this->gene_data_label["pnl"][ik]){
            this->cal_tool->getPnl(ik, this->p, pw_rho, this->pnl[ik]);
        }
        if (this->gene_data_label["qnl"][ik]){
            this->cal_tool->getQnl(ik, this->q, pw_rho, this->qnl[ik]);
        }
        if (this->gene_data_label["xi"][ik]){
            this->cal_tool->getXi(this->gamma, this->gammanl[ik], this->xi[ik]);
        }
        if (this->gene_data_label["tanhxi"][ik]){
            this->cal_tool->getTanhXi(ik, this->gamma, this->gammanl[ik], this->tanhxi[ik]);
        }
        if (this->gene_data_label["tanhxi_nl"][ik]){
            this->cal_tool->getTanhXi_nl(ik, this->tanhxi[ik], pw_rho, this->tanhxi_nl[ik]);
        }
        if (this->gene_data_label["tanh_pnl"][ik]){
            this->cal_tool->getTanh_Pnl(ik, this->pnl[ik], this->tanh_pnl[ik]);
        }
        if (this->gene_data_label["tanh_qnl"][ik]){
            this->cal_tool->getTanh_Qnl(ik, this->qnl[ik], this->tanh_qnl[ik]);
        }
        if (this->gene_data_label["tanhp_nl"][ik]){
            this->cal_tool->getTanhP_nl(ik, this->tanhp, pw_rho, this->tanhp_nl[ik]);
        }
        if (this->gene_data_label["tanhq_nl"][ik]){
            this->cal_tool->getTanhQ_nl(ik, this->tanhq, pw_rho, this->tanhq_nl[ik]);
        }
    }
    ModuleBase::timer::tick("ML_Base", "updateInput");
}

void ML_Base::NN_forward(const double * const * prho, const ModulePW::PW_Basis *pw_rho, bool cal_grad)
{
    ModuleBase::timer::tick("ML_Base", "Forward");

    this->nn->zero_grad();
    this->nn->inputs.requires_grad_(false);
    this->nn->set_data(this, this->descriptor_type, this->kernel_index, this->nn->inputs);
    this->nn->inputs.requires_grad_(true);

    this->nn->F = this->nn->forward(this->nn->inputs);    
    if (this->nn->inputs.grad().numel()) 
    {
        this->nn->inputs.grad().zero_(); 
    }

    if (PARAM.inp.of_ml_feg != 3)
    {
        this->nn->F = torch::softplus(this->nn->F);
    }
    if (PARAM.inp.of_ml_feg == 1)
    {
        this->nn->F = this->nn->F - this->feg_net_F + 1.;
    }
    else if (PARAM.inp.of_ml_feg == 3)
    {
        this->nn->F = torch::softplus(this->nn->F - this->feg_net_F + this->feg3_correct);
    }
    ModuleBase::timer::tick("ML_Base", "Forward");

    if (cal_grad)
    {
        ModuleBase::timer::tick("ML_Base", "Backward");
        this->nn->F.backward(torch::ones({this->nx, 1}, this->device_type));
        ModuleBase::timer::tick("ML_Base", "Backward");
    }
}

torch::Tensor ML_Base::get_data(std::string parameter, const int ikernel){

    if (parameter == "gamma") return torch::tensor(this->gamma, this->device_type);
    if (parameter == "p") return torch::tensor(this->p, this->device_type);
    if (parameter == "q") return torch::tensor(this->q, this->device_type);
    if (parameter == "tanhp") return torch::tensor(this->tanhp, this->device_type);
    if (parameter == "tanhq") return torch::tensor(this->tanhq, this->device_type);
    if (parameter == "gammanl") return torch::tensor(this->gammanl[ikernel], this->device_type);
    if (parameter == "pnl") return torch::tensor(this->pnl[ikernel], this->device_type);
    if (parameter == "qnl") return torch::tensor(this->qnl[ikernel], this->device_type);
    if (parameter == "xi") return torch::tensor(this->xi[ikernel], this->device_type);
    if (parameter == "tanhxi") return torch::tensor(this->tanhxi[ikernel], this->device_type);
    if (parameter == "tanhxi_nl") return torch::tensor(this->tanhxi_nl[ikernel], this->device_type);
    if (parameter == "tanh_pnl") return torch::tensor(this->tanh_pnl[ikernel], this->device_type);
    if (parameter == "tanh_qnl") return torch::tensor(this->tanh_qnl[ikernel], this->device_type);
    if (parameter == "tanhp_nl") return torch::tensor(this->tanhp_nl[ikernel], this->device_type);
    if (parameter == "tanhq_nl") return torch::tensor(this->tanhq_nl[ikernel], this->device_type);
    return torch::zeros({});
}

void ML_Base::get_potential_(const double * const * prho, const ModulePW::PW_Basis *pw_rho, ModuleBase::matrix &rpotential)
{
    ModuleBase::timer::tick("ML_Base", "Pauli Potential");

    std::vector<double> pauli_potential(this->nx, 0.);
    std::vector<double> tau_lda(this->nx, 0.); // Dummy or calculated inside
    for (int ir = 0; ir < this->nx; ++ir)
    {
        tau_lda[ir] = this->energy_prefactor * std::pow(prho[0][ir], this->energy_exponent);
    }

    if (this->ml_gammanl) this->potGammanlTerm(prho, tau_lda, pw_rho, pauli_potential);
    if (this->ml_xi) this->potXinlTerm(prho, tau_lda, pw_rho, pauli_potential);
    if (this->ml_tanhxi) this->potTanhxinlTerm(prho, tau_lda, pw_rho, pauli_potential);
    if (this->ml_tanhxi_nl) this->potTanhxi_nlTerm(prho, tau_lda, pw_rho, pauli_potential);
    if (this->ml_p || this->ml_pnl) this->potPPnlTerm(prho, tau_lda, pw_rho, pauli_potential);
    if (this->ml_q || this->ml_qnl) this->potQQnlTerm(prho, tau_lda, pw_rho, pauli_potential);
    if (this->ml_tanh_pnl) this->potTanhpTanh_pnlTerm(prho, tau_lda, pw_rho, pauli_potential);
    if (this->ml_tanh_qnl) this->potTanhqTanh_qnlTerm(prho, tau_lda, pw_rho, pauli_potential);
    if ((this->ml_tanhp || this->ml_tanhp_nl) && !this->ml_tanh_pnl) this->potTanhpTanhp_nlTerm(prho, tau_lda, pw_rho, pauli_potential);
    if ((this->ml_tanhq || this->ml_tanhq_nl) && !this->ml_tanh_qnl) this->potTanhqTanhq_nlTerm(prho, tau_lda, pw_rho, pauli_potential);

    // Common output logic, using virtual local term functions
    for (int ir = 0; ir < this->nx; ++ir)
    {
        double factor = tau_lda[ir] / prho[0][ir];       
        pauli_potential[ir] += factor *
                      (this->energy_exponent * this->enhancement_cpu_ptr[ir] + this->potGammaTerm(ir) + this->potPTerm1(ir) + this->potQTerm1(ir)
                      + this->potXiTerm1(ir) + this->potTanhxiTerm1(ir) + this->potTanhpTerm1(ir) + this->potTanhqTerm1(ir));
        rpotential(0, ir) += pauli_potential[ir];
    }
    ModuleBase::timer::tick("ML_Base", "Pauli Potential");
}

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

// IO tools
void ML_Base::loadVector(std::string filename, std::vector<double> &data)
{
    npy::npy_data<double> d = npy::read_npy<double>(filename);
    data = d.data;
}

void ML_Base::dumpVector(std::string filename, const std::vector<double> &data)
{
    npy::npy_data_ptr<double> d;
    d.data_ptr = data.data();
    d.shape = {(long unsigned) this->cal_tool->nx};
    d.fortran_order = false;
    npy::write_npy(filename, d);
}

void ML_Base::dumpTensor(std::string filename, const torch::Tensor &data)
{
    std::cout << "Dumping " << filename << std::endl;
    torch::Tensor data_cpu = data.to(this->device_CPU).contiguous();
    std::vector<double> v(data_cpu.data_ptr<double>(), data_cpu.data_ptr<double>() + data_cpu.numel());
    this->dumpVector(filename, v);
}

void ML_Base::dumpMatrix(std::string filename, const ModuleBase::matrix &data)
{
    std::cout << "Dumping " << filename << std::endl;
    std::vector<double> v(data.c, data.c + this->nx);
    this->dumpVector(filename, v);
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
        // tempQ[ir] = (this->ml_q)? 3./40. * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["q"][0]] * /*Ha2Ry*/ 2. : 0.;
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
            // dFdpnl[ir] = this->cDirac * std::pow(prho[0][ir], 5./3.) * this->cal_tool->dtanh(this->tanh_pnl[d2k][ir], this->chi_pnl[d2k])
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
            // tempP[i][ir] = (this->ml_tanhp)? - 3./20. * this->cal_tool->dtanh(this->tanhp[ir], this->chi_p)
                        //    * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["tanhp"][0]] * this->nablaRho[i][ir] / prho[0][ir] * /*Ha2Ry*/ 2. : 0.;
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
            // dFdqnl[ir] = this->cDirac * std::pow(prho[0][ir], 5./3.) * this->cal_tool->dtanh(this->tanh_qnl[d2k][ir], this->chi_qnl[d2k])
                        //  * this->gradient_cpu_ptr[ir * this->ninput + d2i];
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
        // tempQ[ir] = (this->ml_tanhq)? 3./40. * this->cal_tool->dtanh(this->tanhq[ir], this->chi_q)
                    // * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["tanhq"][0]] * /*Ha2Ry*/ 2. : 0.;
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
            // dFdpnl[ir] = this->cDirac * std::pow(prho[0][ir], 5./3.)
            //              * this->gradient_cpu_ptr[ir * this->ninput + d2i];
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
            // tempP[i][ir] = (this->ml_tanhp)? - 3./20. * this->cal_tool->dtanh(this->tanhp[ir], this->chi_p)
                        //    * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["tanhp"][0]] * this->nablaRho[i][ir] / prho[0][ir] * /*Ha2Ry*/ 2. : 0.;
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
            // dFdqnl[ir] = this->cDirac * std::pow(prho[0][ir], 5./3.)
                        //  * this->gradient_cpu_ptr[ir * this->ninput + d2i];
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
        // tempQ[ir] = (this->ml_tanhq)? 3./40. * this->cal_tool->dtanh(this->tanhq[ir], this->chi_q)
                    // * this->gradient_cpu_ptr[ir * this->ninput + this->descriptor2index["tanhq"][0]] * /*Ha2Ry*/ 2. : 0.;
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
