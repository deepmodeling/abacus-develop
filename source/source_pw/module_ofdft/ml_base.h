#ifndef ML_BASE_H
#define ML_BASE_H

#include <vector>
#include <string>
#include <map>
#include <memory>
#include <cmath>

#ifdef __MLALGO
#include "source_pw/module_ofdft/ml_tools/nn_of.h"
#include "source_io/module_ml/cal_mlkedf_descriptors.h"

// The ML_Base class encapsulates common functionality for Machine Learning based
// constructs in OFDFT and EXX.
class ML_Base
{
public:
    ML_Base();
    ~ML_Base();

    // Common Interface
    void set_device(std::string device_inpt);
    
    // Tools
    void loadVector(std::string filename, std::vector<double> &data);
    void dumpVector(std::string filename, const std::vector<double> &data);
    void dumpTensor(std::string filename, const torch::Tensor &data);
    void dumpMatrix(std::string filename, const ModuleBase::matrix &data);

    int nx_tot = 0; // equal to nx (called by NN)
    torch::Tensor get_data(std::string parameter, const int ikernel) const;

protected:
    void updateInput(const double * const * prho, const ModulePW::PW_Basis *pw_rho);
    void NN_forward(const double * const * prho, const ModulePW::PW_Basis *pw_rho, bool cal_grad);
    void get_potential_(const double * const * prho, const ModulePW::PW_Basis *pw_rho, ModuleBase::matrix &rpotential);

    // Potential Terms - these appear identical in both classes or are intended to be shared
    double potGammaTerm(int ir);
    double potPTerm1(int ir);
    double potQTerm1(int ir);
    double potXiTerm1(int ir);
    double potTanhxiTerm1(int ir);
    double potTanhpTerm1(int ir);
    double potTanhqTerm1(int ir);

    // Derived classes should ensure they can work with these signatures.
    // Note: ML_EXX originally passed tau_lda for some of these. 
    // If tau_lda is needed, derived classes can override or we can add it to member variables.
    // For now, keeping signatures compatible with member access.
    void potGammanlTerm(const double * const *prho, const std::vector<double> &tau_lda, const ModulePW::PW_Basis *pw_rho, std::vector<double> &rGammanlTerm);
    void potXinlTerm(const double * const *prho, const std::vector<double> &tau_lda, const ModulePW::PW_Basis *pw_rho, std::vector<double> &rXinlTerm);
    void potTanhxinlTerm(const double * const *prho, const std::vector<double> &tau_lda, const ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhxinlTerm);
    void potTanhxi_nlTerm(const double * const *prho, const std::vector<double> &tau_lda, const ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhxi_nlTerm); 
    void potPPnlTerm(const double * const *prho, const std::vector<double> &tau_lda, const ModulePW::PW_Basis *pw_rho, std::vector<double> &rPPnlTerm);
    void potQQnlTerm(const double * const *prho, const std::vector<double> &tau_lda, const ModulePW::PW_Basis *pw_rho, std::vector<double> &rQQnlTerm);
    void potTanhpTanh_pnlTerm(const double * const *prho, const std::vector<double> &tau_lda, const ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhpTanh_pnlTerm);
    void potTanhqTanh_qnlTerm(const double * const *prho, const std::vector<double> &tau_lda, const ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhqTanh_qnlTerm);
    void potTanhpTanhp_nlTerm(const double * const *prho, const std::vector<double> &tau_lda, const ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhpTanhp_nlTerm);
    void potTanhqTanhq_nlTerm(const double * const *prho, const std::vector<double> &tau_lda, const ModulePW::PW_Basis *pw_rho, std::vector<double> &rTanhqTanhq_nlTerm);

protected: 
    // --- Member Variables (Common) ---

    ModuleIO::Cal_MLKEDF_Descriptors *cal_tool = nullptr;

    int nx = 0; // number of grid points
    double dV = 0.;
    
    // Constants
    double pqcoef = 1.0 / (4.0 * std::pow(3*std::pow(M_PI, 2.0), 2.0/3.0)); // coefficient of p and q
    double feg_net_F = 0.0;
    double feg3_correct = 0.541324854612918; // ln(e - 1)
    double energy_prefactor = 0.0; // cTF for KEDF, cDirac for EXX
    double energy_exponent = 0.0; // 5/3 for KEDF, 4/3 for EXX

    // Descriptors and hyperparameters
    int ninput = 0; // number of descriptors
    std::vector<double> gamma;
    std::vector<double> p;
    std::vector<double> q;
    std::vector<std::vector<double>> gammanl;
    std::vector<std::vector<double>> pnl;
    std::vector<std::vector<double>> qnl;
    std::vector<std::vector<double>> nablaRho;
    
    // Parameters
    std::vector<double> chi_xi;
    double chi_p = 1.0;
    double chi_q = 1.0;
    std::vector<std::vector<double>> xi; 
    std::vector<std::vector<double>> tanhxi;
    std::vector<std::vector<double>> tanhxi_nl; 
    std::vector<double> tanhp;
    std::vector<double> tanhq;
    
    // plan 1
    std::vector<double> chi_pnl;
    std::vector<double> chi_qnl;
    std::vector<std::vector<double>> tanh_pnl;
    std::vector<std::vector<double>> tanh_qnl;
    // plan 2
    std::vector<std::vector<double>> tanhp_nl;
    std::vector<std::vector<double>> tanhq_nl;

    // GPU / Device
    torch::DeviceType device_type = torch::kCPU;
    torch::Device device = torch::Device(torch::kCPU);
    torch::Device device_CPU = torch::Device(torch::kCPU);

    // Neural Network
    std::shared_ptr<NN_OFImpl> nn;
    double* enhancement_cpu_ptr = nullptr;
    double* gradient_cpu_ptr = nullptr;
    int nkernel = 1;

    // Switch flags
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

    // Maps
    std::vector<std::string> descriptor_type;                  
    std::vector<int> kernel_index;                             
    std::map<std::string, std::vector<int>> descriptor2kernel; 
    std::map<std::string, std::vector<int>> descriptor2index;  
    std::map<std::string, std::vector<bool>> gene_data_label;  
};

#endif // __MLALGO
#endif // ML_BASE_H
