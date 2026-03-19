#include "ml_base.h"
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

torch::Tensor ML_Base::get_data(std::string parameter, const int ikernel) const {

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

#endif
