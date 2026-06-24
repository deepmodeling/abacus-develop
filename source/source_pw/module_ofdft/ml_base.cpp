#include "ml_base.h"
#include "npy.hpp"

#ifdef __MLALGO

ML_Base::ML_Base(){}

ML_Base::~ML_Base()
{
    if (this->cal_tool) delete this->cal_tool;
}

void ML_Base::set_device(const std::string& device_inpt, std::ostream& ofs_running)
{
    ModuleBase::TITLE("ML_Base", "set_device");
    if (device_inpt == "cpu")
    {
        ofs_running << "------------------- Running Neural Network on CPU -------------------" << std::endl;
        this->device_type = torch::kCPU;
    }
    else if (device_inpt == "gpu")
    {
        if (torch::cuda::cudnn_is_available())
        {
            ofs_running << "------------------- Running Neural Network on GPU -------------------" << std::endl;
            this->device_type = torch::kCUDA;
        }
        else
        {
	    std::cout << "--------------- Warning: GPU is unavailable ---------------" << std::endl;

            ofs_running << "--------------- Warning: GPU is unavailable ---------------" << std::endl;
            ofs_running << "------------------- Running Neural Network on CPU -------------------" << std::endl;
            this->device_type = torch::kCPU;
        }
    }
    this->device = torch::Device(this->device_type);
}

void ML_Base::update_input(const double * const * prho, const ModulePW::PW_Basis *pw_rho)
{
    ModuleBase::TITLE("ML_Base", "update_input");
    ModuleBase::timer::start("ML_Base", "update_input");
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
    ModuleBase::timer::end("ML_Base", "update_input");
}

void ML_Base::nn_forward(const double * const * prho, const ModulePW::PW_Basis *pw_rho, bool cal_grad)
{
    ModuleBase::TITLE("ML_Base", "nn_forward");
    ModuleBase::timer::start("ML_Base", "nn_forward");

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
    ModuleBase::timer::end("ML_Base", "nn_forward");

    if (cal_grad)
    {
        ModuleBase::timer::start("ML_Base", "backward");
        this->nn->F.backward(torch::ones({this->nx, 1}, this->device_type));
        ModuleBase::timer::end("ML_Base", "backward");
    }
}

torch::Tensor ML_Base::get_data(std::string parameter, const int ikernel) const
{
    if (parameter == "gamma")
    {
        return torch::tensor(this->gamma, this->device_type);
    }
    else if (parameter == "p")
    {
        return torch::tensor(this->p, this->device_type);
    }
    else if (parameter == "q")
    {
        return torch::tensor(this->q, this->device_type);
    }
    else if (parameter == "tanhp")
    {
        return torch::tensor(this->tanhp, this->device_type);
    }
    else if (parameter == "tanhq")
    {
        return torch::tensor(this->tanhq, this->device_type);
    }
    else if (parameter == "gammanl")
    {
        return torch::tensor(this->gammanl[ikernel], this->device_type);
    }
    else if (parameter == "pnl")
    {
        return torch::tensor(this->pnl[ikernel], this->device_type);
    }
    else if (parameter == "qnl")
    {
        return torch::tensor(this->qnl[ikernel], this->device_type);
    }
    else if (parameter == "xi")
    {
        return torch::tensor(this->xi[ikernel], this->device_type);
    }
    else if (parameter == "tanhxi")
    {
        return torch::tensor(this->tanhxi[ikernel], this->device_type);
    }
    else if (parameter == "tanhxi_nl")
    {
        return torch::tensor(this->tanhxi_nl[ikernel], this->device_type);
    }
    else if (parameter == "tanh_pnl")
    {
        return torch::tensor(this->tanh_pnl[ikernel], this->device_type);
    }
    else if (parameter == "tanh_qnl")
    {
        return torch::tensor(this->tanh_qnl[ikernel], this->device_type);
    }
    else if (parameter == "tanhp_nl")
    {
        return torch::tensor(this->tanhp_nl[ikernel], this->device_type);
    }
    else if (parameter == "tanhq_nl")
    {
        return torch::tensor(this->tanhq_nl[ikernel], this->device_type);
    }
    else
    {
        return torch::zeros({});
    }
}

void ML_Base::get_potential_(const double * const * prho, const ModulePW::PW_Basis *pw_rho, ModuleBase::matrix &rpotential)
{
    ModuleBase::TITLE("ML_Base", "get_potential_");
    ModuleBase::timer::start("ML_Base", "pauli_potential");

    std::vector<double> pauli_potential(this->nx, 0.);
    std::vector<double> tau_lda(this->nx, 0.); // Dummy or calculated inside
    for (int ir = 0; ir < this->nx; ++ir)
    {
        tau_lda[ir] = this->energy_prefactor * std::pow(prho[0][ir], this->energy_exponent);
    }

    if (this->ml_gammanl)
    {
        this->pot_gammanl_term(prho, tau_lda, pw_rho, pauli_potential);
    }
    if (this->ml_xi)
    {
        this->pot_xi_nl_term(prho, tau_lda, pw_rho, pauli_potential);
    }
    if (this->ml_tanhxi)
    {
        this->pot_tanhxi_nl_term(prho, tau_lda, pw_rho, pauli_potential);
    }
    if (this->ml_tanhxi_nl)
    {
        this->pot_tanhxi_nl_nl_term(prho, tau_lda, pw_rho, pauli_potential);
    }
    if (this->ml_p || this->ml_pnl)
    {
        this->pot_p_pnl_term(prho, tau_lda, pw_rho, pauli_potential);
    }
    if (this->ml_q || this->ml_qnl)
    {
        this->pot_q_qnl_term(prho, tau_lda, pw_rho, pauli_potential);
    }
    if (this->ml_tanh_pnl)
    {
        this->pot_tanhp_tanh_pnl_term(prho, tau_lda, pw_rho, pauli_potential);
    }
    if (this->ml_tanh_qnl)
    {
        this->pot_tanhq_tanh_qnl_term(prho, tau_lda, pw_rho, pauli_potential);
    }
    if ((this->ml_tanhp || this->ml_tanhp_nl) && !this->ml_tanh_pnl)
    {
        this->pot_tanhp_tanhp_nl_term(prho, tau_lda, pw_rho, pauli_potential);
    }
    if ((this->ml_tanhq || this->ml_tanhq_nl) && !this->ml_tanh_qnl)
    {
        this->pot_tanhq_tanhq_nl_term(prho, tau_lda, pw_rho, pauli_potential);
    }

    for (int ir = 0; ir < this->nx; ++ir)
    {
        double factor = tau_lda[ir] / prho[0][ir];       

        pauli_potential[ir] += factor *
                      (this->energy_exponent * this->enhancement_cpu_ptr[ir] 
		       + this->pot_gamma_term(ir) + this->pot_p_term_1(ir) + this->pot_q_term_1(ir)
                      + this->pot_xi_term_1(ir) + this->pot_tanhxi_term_1(ir) 
		      + this->pot_tanhp_term_1(ir) + this->pot_tanhq_term_1(ir));

        rpotential(0, ir) += pauli_potential[ir];
    }

    ModuleBase::timer::end("ML_Base", "pauli_potential");
}

// IO tools
void ML_Base::load_vector(std::string filename, std::vector<double> &data)
{
    npy::npy_data<double> d = npy::read_npy<double>(filename);
    data = d.data;
}

void ML_Base::dump_vector(std::string filename, const std::vector<double> &data)
{
    npy::npy_data_ptr<double> d;
    d.data_ptr = data.data();
    d.shape = {(long unsigned) this->cal_tool->nx};
    d.fortran_order = false;
    npy::write_npy(filename, d);
}

void ML_Base::dump_tensor(std::string filename, const torch::Tensor &data)
{
    std::cout << "Dumping " << filename << std::endl;
    torch::Tensor data_cpu = data.to(this->device_CPU).contiguous();
    std::vector<double> v(data_cpu.data_ptr<double>(), data_cpu.data_ptr<double>() + data_cpu.numel());
    this->dump_vector(filename, v);
}

void ML_Base::dump_matrix(std::string filename, const ModuleBase::matrix &data)
{
    std::cout << "Dumping " << filename << std::endl;
    std::vector<double> v(data.c, data.c + this->nx);
    this->dump_vector(filename, v);
}

#endif
