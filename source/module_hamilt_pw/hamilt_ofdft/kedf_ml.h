#ifndef KEDF_ML_H
#define KEDF_ML_H

#ifdef __MLKEDF

#include "ml_data.h"

#include <vector>
#include "module_hamilt_pw/hamilt_ofdft/kedf_wt.h"
#include "module_hamilt_pw/hamilt_ofdft/kedf_tf.h"
#include "./ml_tools/nn_of.h"

class KEDF_ML
{
public:
    KEDF_ML()
    {
        // this->stress.create(3,3);
    }
    ~KEDF_ML()
    {
        delete this->ml_data;
    }

    void set_para(
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
        ModulePW::PW_Basis *pw_rho);

    void set_device(std::string device_inpt);

    double get_energy(const double * const * prho, ModulePW::PW_Basis *pw_rho);
    // double get_energy_density(const double * const *prho, int is, int ir, ModulePW::PW_Basis *pw_rho);
    void ml_potential(const double * const * prho, ModulePW::PW_Basis *pw_rho, ModuleBase::matrix &rpotential);
    // void get_stress(double cellVol, const double * const * prho, ModulePW::PW_Basis *pw_rho, double vw_weight);
    // double diffLinhard(double eta, double vw_weight);

    // output all parameters
    void generateTrainData(const double * const *prho, KEDF_WT &wt, KEDF_TF &tf, ModulePW::PW_Basis *pw_rho, const double *veff);
    void localTest(const double * const *prho, ModulePW::PW_Basis *pw_rho);

    // interface to NN
    void NN_forward(const double * const * prho, ModulePW::PW_Basis *pw_rho, bool cal_grad);

    void get_potential_(const double * const * prho, ModulePW::PW_Basis *pw_rho, ModuleBase::matrix &rpotential);

    // potentials
    double potGammaTerm(int ir);
    double potPTerm1(int ir);
    double potQTerm1(int ir);
    double potXiTerm1(int ir);
    double potTanhxiTerm1(int ir);
    double potTanhpTerm1(int ir);
    double potTanhqTerm1(int ir);
    void potGammanlTerm(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rGammanlTerm);
    void potXinlTerm(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rXinlTerm);
    void potTanhxinlTerm(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhxinlTerm);
    void potTanhxi_nlTerm(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhxi_nlTerm); // 2023-03-20 for tanhxi_nl
    void potPPnlTerm(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rPPnlTerm);
    void potQQnlTerm(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rQQnlTerm);
    void potTanhpTanh_pnlTerm(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhpTanh_pnlTerm);
    void potTanhqTanh_qnlTerm(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhqTanh_qnlTerm);
    void potTanhpTanhp_nlTerm(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhpTanhp_nlTerm);
    void potTanhqTanhq_nlTerm(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhqTanhq_nlTerm);
    // tools
    void dumpTensor(const torch::Tensor &data, std::string filename);
    void dumpMatrix(const ModuleBase::matrix &data, std::string filename);
    void updateInput(const double * const * prho, ModulePW::PW_Basis *pw_rho);

    ML_data *ml_data = nullptr;

    int nx = 0; // number of grid points
    int nx_tot = 0; // equal to nx (called by NN)
    double dV = 0.;
    // double weightml = 1.;
    const double cTF = 3.0/10.0 * std::pow(3*std::pow(M_PI, 2.0), 2.0/3.0) * 2; // 10/3*(3*pi^2)^{2/3}, multiply by 2 to convert unit from Hartree to Ry, finally in Ry*Bohr^(-2)
    const double pqcoef = 1.0 / (4.0 * std::pow(3*std::pow(M_PI, 2.0), 2.0/3.0)); // coefficient of p and q
    double ml_energy = 0.;
    // ModuleBase::matrix stress;
    double feg_net_F = 0.;
    double feg3_correct = 0.541324854612918; // ln(e - 1)

    // Descriptors and hyperparameters
    int ninput = 0; // number of descriptors
    std::vector<double> gamma = {};
    std::vector<double> p = {};
    std::vector<double> q = {};
    std::vector<std::vector<double>> gammanl = {};
    std::vector<std::vector<double>> pnl = {};
    std::vector<std::vector<double>> qnl = {};
    std::vector<std::vector<double>> nablaRho = {};
    // new parameters 2023-02-13
    std::vector<double> chi_xi = {1.0};
    double chi_p = 1.;
    double chi_q = 1.;
    std::vector<std::vector<double>> xi = {}; // we assume ONLY ONE of them is used.
    std::vector<std::vector<double>> tanhxi = {};
    std::vector<std::vector<double>> tanhxi_nl= {}; // 2023-03-20
    std::vector<double> tanhp = {};
    std::vector<double> tanhq = {};
    // plan 1
    std::vector<double> chi_pnl = {1.0};
    std::vector<double> chi_qnl = {1.0};
    std::vector<std::vector<double>> tanh_pnl = {};
    std::vector<std::vector<double>> tanh_qnl = {};
    // plan 2
    std::vector<std::vector<double>> tanhp_nl = {};
    std::vector<std::vector<double>> tanhq_nl = {};
    // GPU
    torch::DeviceType device_type = torch::kCPU;
    torch::Device device = torch::Device(torch::kCPU);
    torch::Device device_CPU = torch::Device(torch::kCPU);

    // Nueral Network
    std::shared_ptr<NN_OFImpl> nn;
    double* enhancement_cpu_ptr = nullptr;
    double* gradient_cpu_ptr = nullptr;

    int nkernel = 1; // number of kernels

    // maps
    void init_data(
        const int &nkernel,
        const bool &of_ml_gamma,
        const bool &of_ml_p,
        const bool &of_ml_q,
        const bool &of_ml_tanhp,
        const bool &of_ml_tanhq,
        const std::vector<int> &of_ml_gammanl_,
        const std::vector<int> &of_ml_pnl,
        const std::vector<int> &of_ml_qnl,
        const std::vector<int> &of_ml_xi,
        const std::vector<int> &of_ml_tanhxi,
        const std::vector<int> &of_ml_tanhxi_nl,
        const std::vector<int> &of_ml_tanh_pnl,
        const std::vector<int> &of_ml_tanh_qnl,
        const std::vector<int> &of_ml_tanhp_nl,
        const std::vector<int> &of_ml_tanhq_nl
    );
    
    // Whether to use corresponding descriptors
    bool ml_gamma = false;
    bool ml_p = false;
    bool ml_q = false;
    bool ml_tanhp = false;
    bool ml_tanhq = false;
    bool ml_gammanl = false;
    bool ml_pnl = false;
    bool ml_qnl = false;
    bool ml_xi = false;
    bool ml_tanhxi = false;
    bool ml_tanhxi_nl = false;
    bool ml_tanh_pnl = false;
    bool ml_tanh_qnl = false;
    bool ml_tanhp_nl = false;
    bool ml_tanhq_nl = false;

    std::vector<std::string> descriptor_type = {};                  // the descriptors used
    std::vector<int> kernel_index = {};                             // the index of the kernel used
    std::map<std::string, std::vector<int>> descriptor2kernel = {}; // the map from descriptor to kernel index
    std::map<std::string, std::vector<int>> descriptor2index = {};  // the map from descriptor to index
    std::map<std::string, std::vector<bool>> gene_data_label = {};  // the map from descriptor to gene label

    torch::Tensor get_data(std::string parameter, const int ikernel);   // get the descriptor data for the ikernel-th kernel
};

#endif
#endif