#ifdef __MLALGO

#include "pot_ml_exx.h"

#include "npy.hpp"
#include "source_base/parallel_reduce.h"
#include "source_base/global_function.h"
#include "source_pw/module_pwdft/global.h"

namespace elecstate
{

ML_EXX::ML_EXX()
{
    this->energy_prefactor = - 3. / 4. * std::pow(3. / M_PI, 1./3.) * 2;
    this->energy_exponent = 4. / 3.;
}

ML_EXX::~ML_EXX(){}

void ML_EXX::set_para(const Input_para& inp, const UnitCell* ucell_in, const ModulePW::PW_Basis* rho_basis_in)
{
    torch::set_default_dtype(caffe2::TypeMeta::fromScalarType(torch::kDouble));
    auto output = torch::get_default_dtype();
    std::cout << "Default type: " << output << std::endl;

    this->set_device(inp.of_ml_device);

    this->nx = rho_basis_in->nrxx;
    this->nx_tot = rho_basis_in->nrxx;
    this->dV = ucell_in->omega / rho_basis_in->nxyz;
    this->nkernel = inp.of_ml_nkernel;

    this->init_data(
        this->nkernel,
        inp.of_ml_gamma,
        inp.of_ml_p,
        inp.of_ml_q,
        inp.of_ml_tanhp,
        inp.of_ml_tanhq,
        inp.of_ml_gammanl,
        inp.of_ml_pnl,
        inp.of_ml_qnl,
        inp.of_ml_xi,
        inp.of_ml_tanhxi,
        inp.of_ml_tanhxi_nl,
        inp.of_ml_tanh_pnl,
        inp.of_ml_tanh_qnl,
        inp.of_ml_tanhp_nl,
        inp.of_ml_tanhq_nl);

    std::cout << "ninput = " << this->ninput << std::endl;

    if (PARAM.inp.ml_exx)
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
    
    if (PARAM.inp.ml_exx || PARAM.inp.of_ml_gene_data == 1)
    {
        this->cal_tool = new ModuleIO::Cal_MLKEDF_Descriptors;

        this->chi_p = inp.of_ml_chi_p;
        this->chi_q = inp.of_ml_chi_q;
        this->chi_xi = inp.of_ml_chi_xi;
        this->chi_pnl = inp.of_ml_chi_pnl;
        this->chi_qnl = inp.of_ml_chi_qnl;

        this->cal_tool->set_para(this->nx, inp.nelec, inp.of_tf_weight, inp.of_vw_weight, this->chi_p, this->chi_q,
                                this->chi_xi, this->chi_pnl, this->chi_qnl, this->nkernel, inp.of_ml_kernel, inp.of_ml_kernel_scaling, inp.of_ml_yukawa_alpha, inp.of_ml_kernel_file, this->dV * rho_basis_in->nxyz, rho_basis_in);
    }
}


/**
 * @brief Get the potential of ML KEDF, and add it into rpotential
 * 
 * @param prho charge density
 * @param pw_rho PW_Basis
 * @param rpotential rpotential => rpotential + V_{ML}
 */
void ML_EXX::ml_potential(const double * const * prho, const ModulePW::PW_Basis *pw_rho, ModuleBase::matrix &rpotential)
{
    double* rho_data = new double[this->nx];
    const double** prho_mod = new const double*[1];
    prho_mod[0] = rho_data;

    for (int ir = 0; ir < this->nx; ++ir)
    {
        rho_data[ir] = std::abs(prho[0][ir]);
    }

    this->updateInput(prho_mod, pw_rho);

    this->NN_forward(prho_mod, pw_rho, true);
    
    torch::Tensor enhancement_cpu_tensor = this->nn->F.to(this->device_CPU).contiguous();
    this->enhancement_cpu_ptr = enhancement_cpu_tensor.data_ptr<double>();
    torch::Tensor gradient_cpu_tensor = this->nn->inputs.grad().to(this->device_CPU).contiguous();
    this->gradient_cpu_ptr = gradient_cpu_tensor.data_ptr<double>();

    this->get_potential_(prho_mod, pw_rho, rpotential);

    // get energy
    ModuleBase::timer::tick("ML_EXX", "Pauli Energy");
    double energy = 0.;
    for (int ir = 0; ir < this->nx; ++ir)
    {
        energy += this->enhancement_cpu_ptr[ir] * std::pow(prho_mod[0][ir], this->energy_exponent);
    }
    energy *= this->dV * this->energy_prefactor;
    this->ml_exx_energy = energy;
    Parallel_Reduce::reduce_pool(this->ml_exx_energy);
    ModuleBase::timer::tick("ML_EXX", "Pauli Energy");

    delete[] rho_data;
    delete[] prho_mod;
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
void ML_EXX::generateTrainData(const double * const *prho, const ModulePW::PW_Basis *pw_rho, const double *veff)
{
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

        this->dumpTensor("enhancement.npy", enhancement);
        this->dumpMatrix("potential.npy", potential);
    }
}

/**
 * @brief For test
 * 
 * @param prho charge density
 * @param pw_rho PW_Basis
 */
void ML_EXX::localTest(const double * const *pprho, const ModulePW::PW_Basis *pw_rho)
{
    // for test =====================
    std::vector<long unsigned int> cshape = {(long unsigned) this->nx};
    bool fortran_order = false;

    std::vector<double> temp_prho(this->nx);
    this->loadVector("path_to_rho_file", temp_prho);
    
    double ** prho = new double *[1];
    prho[0] = new double[this->nx];
    for (int ir = 0; ir < this->nx; ++ir) prho[0][ir] = temp_prho[ir];
    for (int ir = 0; ir < this->nx; ++ir) 
    {
        if (prho[0][ir] == 0.){
            std::cout << "WARNING: rho = 0" << std::endl;
        }
    };
    // ==============================
    this->updateInput(prho, pw_rho);

    this->NN_forward(prho, pw_rho, true);
    
    torch::Tensor enhancement_cpu_tensor = this->nn->F.to(this->device_CPU).contiguous();
    this->enhancement_cpu_ptr = enhancement_cpu_tensor.data_ptr<double>();
    torch::Tensor gradient_cpu_tensor = this->nn->inputs.grad().to(this->device_CPU).contiguous();
    this->gradient_cpu_ptr = gradient_cpu_tensor.data_ptr<double>();

    torch::Tensor enhancement = this->nn->F.reshape({this->nx});
    ModuleBase::matrix potential(1, this->nx);

    this->get_potential_(prho, pw_rho, potential);

    this->dumpTensor("enhancement-abacus.npy", enhancement);
    this->dumpMatrix("potential-abacus.npy", potential);
    exit(0);
}

}   // namespace elecstate
#endif
