#ifdef __MLALGO

#include "kedf_ml.h"

#include "npy.hpp"
#include "source_base/parallel_reduce.h"
#include "source_base/global_function.h"

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
    ModulePW::PW_Basis *pw_rho,
    std::ostream& ofs_running
)
{
    ModuleBase::TITLE("KEDF_ML", "set_para");
    torch::set_default_dtype(caffe2::TypeMeta::fromScalarType(torch::kDouble));
    auto output = torch::get_default_dtype();
    ofs_running << " Default type: " << output << std::endl;

    this->set_device(device_inpt, ofs_running);

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
        of_ml_tanhq_nl,
        ofs_running);

    ofs_running << " ninput = " << ninput << " (number of descriptors)" << std::endl;
    ofs_running << " nkernel = " << this->nkernel << " (number of kernel functions)" << std::endl;

    if (PARAM.inp.of_kinetic == "ml")
    {
        int nnode = 100;
        int nlayer = 3;
        this->nn = std::make_shared<NN_OFImpl>(this->nx, 0, this->ninput, nnode, nlayer, this->device, ofs_running);
        try
        {
            torch::load(this->nn, "net.pt", this->device_type);
        }
        catch (const std::exception& e)
        {
            ModuleBase::WARNING_QUIT("KEDF_ML::set_para", 
                                    "Failed to load neural network model from net.pt: " + std::string(e.what()));
        }
        ofs_running << " load net done (ML KEDF neural network model loaded successfully)" << std::endl;
        if (PARAM.inp.of_ml_feg != 0)
        {
            torch::Tensor feg_inpt = torch::zeros(this->ninput, this->device_type);
            for (int i = 0; i < this->ninput; ++i)
            {
                if (this->descriptor_type[i] == "gamma") feg_inpt[i] = 1.;
            }

            if (PARAM.inp.of_ml_feg == 1)
            {
                this->feg_net_F = torch::softplus(this->nn->forward(feg_inpt)).to(this->device_CPU).contiguous().data_ptr<double>()[0];
            }
            else
            {
                this->feg_net_F = this->nn->forward(feg_inpt).to(this->device_CPU).contiguous().data_ptr<double>()[0];
            }

            ofs_running << " feg_net_F = " << this->feg_net_F 
		    << " (Pauli energy enhancement factor in free electron gas)" << std::endl << std::endl;
        }
    }
    else
    {
        ofs_running << " ML KEDF not enabled (of_kinetic != \"ml\")" << std::endl;
    }
    
    if (PARAM.inp.of_kinetic == "ml" || PARAM.inp.of_ml_gene_data == 1)
    {
        this->cal_tool = new ModuleIO::Cal_MLKEDF_Descriptors;

        this->chi_p = chi_p;
        this->chi_q = chi_q;
        this->chi_xi = chi_xi;
        this->chi_pnl = chi_pnl;
        this->chi_qnl = chi_qnl;

        this->cal_tool->set_para(nx, nelec, tf_weight, vw_weight, chi_p, chi_q,
                                chi_xi, chi_pnl, chi_qnl, nkernel, kernel_type, 
				kernel_scaling, yukawa_alpha, kernel_file, 
				this->dV * pw_rho->nxyz, pw_rho, ofs_running);
    }
    else
    {
        ofs_running << " ML descriptor calculator not initialized (neither ml kinetic nor gene_data enabled)" << std::endl;
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
    ModuleBase::TITLE("KEDF_ML", "get_energy");
    this->update_input(prho, pw_rho);

    this->nn_forward(prho, pw_rho, false);
    
    torch::Tensor enhancement_cpu_tensor = this->nn->F.to(this->device_CPU).contiguous();
    this->enhancement_cpu_ptr = enhancement_cpu_tensor.data_ptr<double>();

    double energy = 0.;
    for (int ir = 0; ir < this->nx; ++ir)
    {
        energy += enhancement_cpu_ptr[ir] * std::pow(prho[0][ir], this->energy_exponent);
    }
    energy *= this->dV * this->energy_prefactor;
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
    ModuleBase::TITLE("KEDF_ML", "ml_potential");
    ModuleBase::timer::start("KEDF_ML", "ml_potential");

    this->update_input(prho, pw_rho);

    this->nn_forward(prho, pw_rho, true);
    
    torch::Tensor enhancement_cpu_tensor = this->nn->F.to(this->device_CPU).contiguous();

    this->enhancement_cpu_ptr = enhancement_cpu_tensor.data_ptr<double>();

    torch::Tensor gradient_cpu_tensor = this->nn->inputs.grad().to(this->device_CPU).contiguous();

    this->gradient_cpu_ptr = gradient_cpu_tensor.data_ptr<double>();

    this->get_potential_(prho, pw_rho, rpotential);

    // Calculate Pauli energy (ml_energy) from enhancement factor
    // E_pauli = c_TF * ∫ F(ρ) * ρ^(5/3) dr
    double energy = 0.;
    for (int ir = 0; ir < this->nx; ++ir)
    {
        energy += enhancement_cpu_ptr[ir] * std::pow(prho[0][ir], this->energy_exponent);
    }
    energy *= this->dV * this->energy_prefactor;
    this->ml_energy = energy;
    Parallel_Reduce::reduce_all(this->ml_energy);

    ModuleBase::timer::end("KEDF_ML", "ml_potential");
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
void KEDF_ML::gen_training_data(const double * const *prho, ModulePW::PW_Basis *pw_rho, const double *veff)
{
    ModuleBase::TITLE("KEDF_ML", "gen_training_data");
    // this->cal_tool->generateTrainData_WT(prho, wt, tf, pw_rho, veff); // Will be fixed in next pr
    if (PARAM.inp.of_kinetic == "ml")
    {
        this->update_input(prho, pw_rho);

        this->nn_forward(prho, pw_rho, true);
        
        torch::Tensor enhancement_cpu_tensor = this->nn->F.to(this->device_CPU).contiguous();
        this->enhancement_cpu_ptr = enhancement_cpu_tensor.data_ptr<double>();
        torch::Tensor gradient_cpu_tensor = this->nn->inputs.grad().to(this->device_CPU).contiguous();
        this->gradient_cpu_ptr = gradient_cpu_tensor.data_ptr<double>();

        torch::Tensor enhancement = this->nn->F.reshape({this->nx});
        ModuleBase::matrix potential(1, this->nx);

        this->get_potential_(prho, pw_rho, potential);

        this->dump_tensor("enhancement.npy", enhancement);
        this->dump_matrix("potential.npy", potential);
    }
    else
    {
        std::cout << " Warning: gen_training_data skipped (of_kinetic != \"ml\")" << std::endl;
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
    ModuleBase::TITLE("KEDF_ML", "local_test");
    // for test =====================
    std::vector<long unsigned int> cshape = {(long unsigned) this->nx};
    bool fortran_order = false;

    std::vector<double> temp_prho(this->nx);
    this->load_vector("dir_of_input_rho", temp_prho);
    double ** prho = new double *[1];
    prho[0] = new double[this->nx];
    for (int ir = 0; ir < this->nx; ++ir) prho[0][ir] = temp_prho[ir];
    for (int ir = 0; ir < this->nx; ++ir) 
    {
        if (prho[0][ir] == 0.)
        {
            std::cout << "WARNING: rho = 0 at grid point " << ir << std::endl;
        }
        else
        {
            // Normal case: non-zero density
        }
    };
    // ==============================

    this->update_input(prho, pw_rho);

    this->nn_forward(prho, pw_rho, true);
    
    torch::Tensor enhancement_cpu_tensor = this->nn->F.to(this->device_CPU).contiguous();
    this->enhancement_cpu_ptr = enhancement_cpu_tensor.data_ptr<double>();
    torch::Tensor gradient_cpu_tensor = this->nn->inputs.grad().to(this->device_CPU).contiguous();
    this->gradient_cpu_ptr = gradient_cpu_tensor.data_ptr<double>();

    torch::Tensor enhancement = this->nn->F.reshape({this->nx});
    ModuleBase::matrix potential(1, this->nx);

    this->get_potential_(prho, pw_rho, potential);

    this->dump_tensor("enhancement-abacus.npy", enhancement);
    this->dump_matrix("potential-abacus.npy", potential);
    exit(0);
}
#endif
