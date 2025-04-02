#ifdef __MLKEDF

#include "kedf_ml.h"

#include "module_base/parallel_reduce.h"
#include "module_base/global_function.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

void KEDF_ML::set_para(
    const int nx, 
    const double dV, 
    const double nelec, 
    const double tf_weight, 
    const double vw_weight, 
    const double chi_p,
    const double chi_q,
    const std::vector<double> &chi_xi,
    const std::vector<double> &chi_pnl,
    const std::vector<double> &chi_qnl,
    const int &nkernel,
    const std::vector<int> &kernel_type,
    const std::vector<double> &kernel_scaling,
    const std::vector<double> &yukawa_alpha,
    const std::vector<std::string> &kernel_file,
    const bool &of_ml_gamma,
    const bool &of_ml_p,
    const bool &of_ml_q,
    const bool &of_ml_tanhp,
    const bool &of_ml_tanhq,
    const std::vector<int> &of_ml_gammanl,
    const std::vector<int> &of_ml_pnl,
    const std::vector<int> &of_ml_qnl,
    const std::vector<int> &of_ml_xi,
    const std::vector<int> &of_ml_tanhxi,
    const std::vector<int> &of_ml_tanhxi_nl,
    const std::vector<int> &of_ml_tanh_pnl,
    const std::vector<int> &of_ml_tanh_qnl,
    const std::vector<int> &of_ml_tanhp_nl,
    const std::vector<int> &of_ml_tanhq_nl,
    const std::string device_inpt,
    ModulePW::PW_Basis *pw_rho
)
{
    torch::set_default_dtype(caffe2::TypeMeta::fromScalarType(torch::kDouble));
    auto output = torch::get_default_dtype();
    std::cout << "Default type: " << output << std::endl;

    this->set_device(device_inpt);

    this->nx = nx;
    this->nx_tot = nx;
    this->dV = dV;
    this->nkernel = nkernel;

    this->init_data(
        nkernel,
        of_ml_gamma,
        of_ml_p,
        of_ml_q,
        of_ml_tanhp,
        of_ml_tanhq,
        of_ml_gammanl,
        of_ml_pnl,
        of_ml_qnl,
        of_ml_xi,
        of_ml_tanhxi,
        of_ml_tanhxi_nl,
        of_ml_tanh_pnl,
        of_ml_tanh_qnl,
        of_ml_tanhp_nl,
        of_ml_tanhq_nl);

    std::cout << "ninput = " << ninput << std::endl;

    if (PARAM.inp.of_kinetic == "ml")
    {
        int nnode = 100;
        int nlayer = 3;
        this->nn = std::make_shared<NN_OFImpl>(this->nx, 0, this->ninput, nnode, nlayer, this->device);
        torch::load(this->nn, "net.pt", this->device_type);
        std::cout << "load net done" << std::endl;
        if (PARAM.inp.of_ml_feg != 0)
        {
            torch::Tensor feg_inpt = torch::zeros(this->ninput, this->device_type);
            for (int i = 0; i < this->ninput; ++i)
            {
                if (this->descriptor_type[i] == "gamma") feg_inpt[i] = 1.;
            }

            if (PARAM.inp.of_ml_feg == 1) 
                this->feg_net_F = torch::softplus(this->nn->forward(feg_inpt)).to(this->device_CPU).contiguous().data_ptr<double>()[0];
            else
            {
                this->feg_net_F = this->nn->forward(feg_inpt).to(this->device_CPU).contiguous().data_ptr<double>()[0];
            }

            std::cout << "feg_net_F = " << this->feg_net_F << std::endl;
        }
    } 
    
    if (PARAM.inp.of_kinetic == "ml" || PARAM.inp.of_ml_gene_data == 1)
    {
        this->ml_data = new ML_data;

        this->chi_p = chi_p;
        this->chi_q = chi_q;
        this->chi_xi = chi_xi;
        this->chi_pnl = chi_pnl;
        this->chi_qnl = chi_qnl;

        this->ml_data->set_para(nx, nelec, tf_weight, vw_weight, chi_p, chi_q,
                                chi_xi, chi_pnl, chi_qnl, nkernel, kernel_type, kernel_scaling, yukawa_alpha, kernel_file, this->dV * pw_rho->nxyz, pw_rho);
    }
}

/**
 * @brief Get the energy of ML KEDF
 * \f[ E_{ML} = c_{TF} * \int{F(\rho) \rho^{5/3} dr} \f]
 * 
 * @param prho charge density
 * @param pw_rho PW_Basis
 * @return the energy of ML KEDF
 */
double KEDF_ML::get_energy(const double * const * prho, ModulePW::PW_Basis *pw_rho)
{
    this->updateInput(prho, pw_rho);

    this->NN_forward(prho, pw_rho, false);
    
    torch::Tensor enhancement_cpu_tensor = this->nn->F.to(this->device_CPU).contiguous();
    this->enhancement_cpu_ptr = enhancement_cpu_tensor.data_ptr<double>();

    double energy = 0.;
    for (int ir = 0; ir < this->nx; ++ir)
    {
        energy += enhancement_cpu_ptr[ir] * std::pow(prho[0][ir], 5./3.);
    }
    std::cout << "energy" << energy << std::endl;
    energy *= this->dV * this->cTF;
    this->ml_energy = energy;
    Parallel_Reduce::reduce_all(this->ml_energy);
    return this->ml_energy;
}

/**
 * @brief Get the potential of ML KEDF, and add it into rpotential
 * 
 * @param prho charge density
 * @param pw_rho PW_Basis
 * @param rpotential rpotential => rpotential + V_{ML}
 */
void KEDF_ML::ml_potential(const double * const * prho, ModulePW::PW_Basis *pw_rho, ModuleBase::matrix &rpotential)
{
    this->updateInput(prho, pw_rho);

    this->NN_forward(prho, pw_rho, true);
    
    torch::Tensor enhancement_cpu_tensor = this->nn->F.to(this->device_CPU).contiguous();
    this->enhancement_cpu_ptr = enhancement_cpu_tensor.data_ptr<double>();
    torch::Tensor gradient_cpu_tensor = this->nn->inputs.grad().to(this->device_CPU).contiguous();
    this->gradient_cpu_ptr = gradient_cpu_tensor.data_ptr<double>();

    this->get_potential_(prho, pw_rho, rpotential);

    // get energy
    ModuleBase::timer::tick("KEDF_ML", "Pauli Energy");
    double energy = 0.;
    for (int ir = 0; ir < this->nx; ++ir)
    {
        energy += enhancement_cpu_ptr[ir] * std::pow(prho[0][ir], 5./3.);
    }
    energy *= this->dV * this->cTF;
    this->ml_energy = energy;
    Parallel_Reduce::reduce_all(this->ml_energy);
    ModuleBase::timer::tick("KEDF_ML", "Pauli Energy");
}

/**
 * @brief Generate training data for ML KEDF
 * 
 * @param prho charge density
 * @param wt KEDF_WT
 * @param tf KEDF_TF
 * @param pw_rho PW_Basis
 * @param veff effective potential
 */
void KEDF_ML::generateTrainData(const double * const *prho, KEDF_WT &wt, KEDF_TF &tf,  ModulePW::PW_Basis *pw_rho, const double *veff)
{
    this->ml_data->generateTrainData_WT(prho, wt, tf, pw_rho, veff);
    if (PARAM.inp.of_kinetic == "ml")
    {
        this->updateInput(prho, pw_rho);

        this->NN_forward(prho, pw_rho, true);
        
        torch::Tensor enhancement_cpu_tensor = this->nn->F.to(this->device_CPU).contiguous();
        this->enhancement_cpu_ptr = enhancement_cpu_tensor.data_ptr<double>();
        torch::Tensor gradient_cpu_tensor = this->nn->inputs.grad().to(this->device_CPU).contiguous();
        this->gradient_cpu_ptr = gradient_cpu_tensor.data_ptr<double>();

        torch::Tensor enhancement = this->nn->F.reshape({this->nx});
        ModuleBase::matrix potential(1, this->nx);

        this->get_potential_(prho, pw_rho, potential);

        std::cout << "dumpdump\n";
        this->dumpTensor(enhancement, "enhancement.npy");
        this->dumpMatrix(potential, "potential.npy");
    }
}

/**
 * @brief For test
 * 
 * @param prho charge density
 * @param pw_rho PW_Basis
 */
void KEDF_ML::localTest(const double * const *pprho, ModulePW::PW_Basis *pw_rho)
{
    // for test =====================
    std::vector<long unsigned int> cshape = {(long unsigned) this->nx};
    bool fortran_order = false;

    std::vector<double> temp_prho(this->nx);
    this->ml_data->loadVector("dir_of_input_rho", temp_prho);
    double ** prho = new double *[1];
    prho[0] = new double[this->nx];
    for (int ir = 0; ir < this->nx; ++ir) prho[0][ir] = temp_prho[ir];
    for (int ir = 0; ir < this->nx; ++ir) 
    {
        if (prho[0][ir] == 0.){
            std::cout << "WARNING: rho = 0" << std::endl;
        }
    };
    std::cout << "Load rho done" << std::endl;
    // ==============================

    this->updateInput(prho, pw_rho);
    std::cout << "update done" << std::endl;

    this->NN_forward(prho, pw_rho, true);
    
    torch::Tensor enhancement_cpu_tensor = this->nn->F.to(this->device_CPU).contiguous();
    this->enhancement_cpu_ptr = enhancement_cpu_tensor.data_ptr<double>();
    torch::Tensor gradient_cpu_tensor = this->nn->inputs.grad().to(this->device_CPU).contiguous();
    this->gradient_cpu_ptr = gradient_cpu_tensor.data_ptr<double>();

    std::cout << "enhancement done" << std::endl;

    torch::Tensor enhancement = this->nn->F.reshape({this->nx});
    ModuleBase::matrix potential(1, this->nx);

    this->get_potential_(prho, pw_rho, potential);
    std::cout << "potential done" << std::endl;

    this->dumpTensor(enhancement, "enhancement-abacus.npy");
    this->dumpMatrix(potential, "potential-abacus.npy");
    exit(0);
}

/**
 * @brief Set the device for ML KEDF
 * 
 * @param device_inpt "cpu" or "gpu"
 */
void KEDF_ML::set_device(std::string device_inpt)
{
    if (device_inpt == "cpu")
    {
        std::cout << "---------- Running NN on CPU ----------" << std::endl;
        this->device_type = torch::kCPU;
    }
    else if (device_inpt == "gpu")
    {
        if (torch::cuda::cudnn_is_available())
        {
            std::cout << "---------- Running NN on GPU ----------" << std::endl;
            this->device_type = torch::kCUDA;
        }
        else
        {
            std::cout << "------ Warning: GPU is unaviable ------" << std::endl;
            std::cout << "---------- Running NN on CPU ----------" << std::endl;
            this->device_type = torch::kCPU;
        }
    }
    this->device = torch::Device(this->device_type);
}

/**
 * @brief Interface to Neural Network forward
 * 
 * @param prho charge density
 * @param pw_rho PW_Basis
 * @param cal_grad whether to calculate the gradient
 */
void KEDF_ML::NN_forward(const double * const * prho, ModulePW::PW_Basis *pw_rho, bool cal_grad)
{
    ModuleBase::timer::tick("KEDF_ML", "Forward");

    this->nn->zero_grad();
    this->nn->inputs.requires_grad_(false);
    this->nn->set_data(this, this->descriptor_type, this->kernel_index, this->nn->inputs);
    this->nn->inputs.requires_grad_(true);

    this->nn->F = this->nn->forward(this->nn->inputs);    
    if (this->nn->inputs.grad().numel()) 
    {
        this->nn->inputs.grad().zero_(); // In the first step, inputs.grad() returns an undefined Tensor, so that numel() = 0.
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
    ModuleBase::timer::tick("KEDF_ML", "Forward");

    if (cal_grad)
    {
        ModuleBase::timer::tick("KEDF_ML", "Backward");
        this->nn->F.backward(torch::ones({this->nx, 1}, this->device_type));
        ModuleBase::timer::tick("KEDF_ML", "Backward");
    }
}

/**
 * @brief Dump the torch::Tensor into .npy file
 * 
 * @param data torch::Tensor
 * @param filename file name
 */
void KEDF_ML::dumpTensor(const torch::Tensor &data, std::string filename)
{
    std::cout << "Dumping " << filename << std::endl;
    torch::Tensor data_cpu = data.to(this->device_CPU).contiguous();
    std::vector<double> v(data_cpu.data_ptr<double>(), data_cpu.data_ptr<double>() + data_cpu.numel());
    // for (int ir = 0; ir < this->nx; ++ir) assert(v[ir] == data[ir].item<double>());
    this->ml_data->dumpVector(filename, v);
}

/**
 * @brief Dump the matrix into .npy file
 * 
 * @param data matrix
 * @param filename file name
 */
void KEDF_ML::dumpMatrix(const ModuleBase::matrix &data, std::string filename)
{
    std::cout << "Dumping " << filename << std::endl;
    std::vector<double> v(data.c, data.c + this->nx);
    // for (int ir = 0; ir < this->nx; ++ir) assert(v[ir] == data[ir].item<double>());
    this->ml_data->dumpVector(filename, v);
}

/**
 * @brief Update the desciptors for ML KEDF
 * 
 * @param prho charge density
 * @param pw_rho PW_Basis
 */
void KEDF_ML::updateInput(const double * const * prho, ModulePW::PW_Basis *pw_rho)
{
    ModuleBase::timer::tick("KEDF_ML", "updateInput");
    // std::cout << "updata_input" << std::endl;
    if (this->gene_data_label["gamma"][0])
    {   
        this->ml_data->getGamma(prho, this->gamma);
    }
    if (this->gene_data_label["p"][0])
    {
        this->ml_data->getNablaRho(prho, pw_rho, this->nablaRho);
        this->ml_data->getP(prho, pw_rho, this->nablaRho, this->p);
    }
    if (this->gene_data_label["q"][0])
    {
        this->ml_data->getQ(prho, pw_rho, this->q);
    }
    if (this->gene_data_label["tanhp"][0])
    {
        this->ml_data->getTanhP(this->p, this->tanhp);
    }
    if (this->gene_data_label["tanhq"][0])
    {
        this->ml_data->getTanhQ(this->q, this->tanhq);
    }

    for (int ik = 0; ik < nkernel; ++ik)
    {
        if (this->gene_data_label["gammanl"][ik]){
            this->ml_data->getGammanl(ik, this->gamma, pw_rho, this->gammanl[ik]);
        }
        if (this->gene_data_label["pnl"][ik]){
            this->ml_data->getPnl(ik, this->p, pw_rho, this->pnl[ik]);
        }
        if (this->gene_data_label["qnl"][ik]){
            this->ml_data->getQnl(ik, this->q, pw_rho, this->qnl[ik]);
        }
        if (this->gene_data_label["xi"][ik]){
            this->ml_data->getXi(this->gamma, this->gammanl[ik], this->xi[ik]);
        }
        if (this->gene_data_label["tanhxi"][ik]){
            this->ml_data->getTanhXi(ik, this->gamma, this->gammanl[ik], this->tanhxi[ik]);
        }
        if (this->gene_data_label["tanhxi_nl"][ik]){
            this->ml_data->getTanhXi_nl(ik, this->tanhxi[ik], pw_rho, this->tanhxi_nl[ik]);
        }
        if (this->gene_data_label["tanh_pnl"][ik]){
            this->ml_data->getTanh_Pnl(ik, this->pnl[ik], this->tanh_pnl[ik]);
        }
        if (this->gene_data_label["tanh_qnl"][ik]){
            this->ml_data->getTanh_Qnl(ik, this->qnl[ik], this->tanh_qnl[ik]);
        }
        if (this->gene_data_label["tanhp_nl"][ik]){
            this->ml_data->getTanhP_nl(ik, this->tanhp, pw_rho, this->tanhp_nl[ik]);
        }
        if (this->gene_data_label["tanhq_nl"][ik]){
            this->ml_data->getTanhQ_nl(ik, this->tanhq, pw_rho, this->tanhq_nl[ik]);
        }
    }
    ModuleBase::timer::tick("KEDF_ML", "updateInput");
}

/**
 * @brief Return the descriptors for ML KEDF
 * 
 * @param parameter "gamma", "p", "q", "tanhp", "tanhq", "gammanl", "pnl", "qnl", "xi", "tanhxi", "tanhxi_nl", "tanh_pnl", "tanh_qnl", "tanhp_nl", "tanhq_nl"
 * @param ikernel kernel index
 */
torch::Tensor KEDF_ML::get_data(std::string parameter, const int ikernel){

    if (parameter == "gamma")
    {
        return torch::tensor(this->gamma, this->device_type);
    }
    if (parameter == "p")
    {
        return torch::tensor(this->p, this->device_type);
    }
    if (parameter == "q")
    {
        return torch::tensor(this->q, this->device_type);
    }
    if (parameter == "tanhp")
    {
        return torch::tensor(this->tanhp, this->device_type);
    }
    if (parameter == "tanhq")
    {
        return torch::tensor(this->tanhq, this->device_type);
    }
    if (parameter == "gammanl")
    {
        return torch::tensor(this->gammanl[ikernel], this->device_type);
    }
    if (parameter == "pnl")
    {
        return torch::tensor(this->pnl[ikernel], this->device_type);
    }
    if (parameter == "qnl")
    {
        return torch::tensor(this->qnl[ikernel], this->device_type);
    }
    if (parameter == "xi")
    {
        return torch::tensor(this->xi[ikernel], this->device_type);
    }
    if (parameter == "tanhxi")
    {
        return torch::tensor(this->tanhxi[ikernel], this->device_type);
    }
    if (parameter == "tanhxi_nl")
    {
        return torch::tensor(this->tanhxi_nl[ikernel], this->device_type);
    }
    if (parameter == "tanh_pnl")
    {
        return torch::tensor(this->tanh_pnl[ikernel], this->device_type);
    }
    if (parameter == "tanh_qnl")
    {
        return torch::tensor(this->tanh_qnl[ikernel], this->device_type);
    }
    if (parameter == "tanhp_nl")
    {
        return torch::tensor(this->tanhp_nl[ikernel], this->device_type);
    }
    if (parameter == "tanhq_nl")
    {
        return torch::tensor(this->tanhq_nl[ikernel], this->device_type);
    }
    return torch::zeros({});
}
#endif