#ifndef CAL_MLKEDF_DESCRIPTORS_H
#define CAL_MLKEDF_DESCRIPTORS_H

#include <vector>
#include "source_base/global_function.h"
#include "source_pw/module_pwdft/global.h"
#include "source_io/module_parameter/parameter.h"

namespace ModuleIO
{

/**
 * @brief A class to calculate the descriptors for ML KEDF.
 * Sun, Liang, and Mohan Chen. Physical Review B 109.11 (2024): 115135.
 * Sun, Liang, and Mohan Chen. Electronic Structure 6.4 (2024): 045006.
 * @author sunliang on 2025-06-12
 */
class Cal_MLKEDF_Descriptors
{
public:
    ~Cal_MLKEDF_Descriptors() {}

    void set_para(
        const int &nx,
        const double &nelec, 
        const double &tf_weight, 
        const double &vw_weight,
        const double &chi_p,
        const double &chi_q,
        const std::vector<double> &chi_xi,
        const std::vector<double> &chi_pnl,
        const std::vector<double> &chi_qnl,
        const int &nkernel,
        const std::vector<int> &kernel_type,
        const std::vector<double> &kernel_scaling,
        const std::vector<double> &yukawa_alpha,
        const std::vector<std::string> &kernel_file,
        const double &omega,
        ModulePW::PW_Basis *pw_rho);
    // get input parameters
    void getGamma(const double * const *prho, std::vector<double> &rgamma);
    void getP(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<std::vector<double>> &pnablaRho, std::vector<double> &rp);
    void getQ(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<double> &rq);
    void getGammanl(const int ikernel, std::vector<double> &pgamma, ModulePW::PW_Basis *pw_rho, std::vector<double> &rgammanl);
    void getPnl(const int ikernel, std::vector<double> &pp, ModulePW::PW_Basis *pw_rho, std::vector<double> &rpnl);
    void getQnl(const int ikernel, std::vector<double> &pq, ModulePW::PW_Basis *pw_rho, std::vector<double> &rqnl);
    // new parameters 2023-02-03
    void getXi(std::vector<double> &pgamma, std::vector<double> &pgammanl, std::vector<double> &rxi);
    void getTanhXi(const int ikernel, std::vector<double> &pgamma, std::vector<double> &pgammanl, std::vector<double> &rtanhxi);
    void getTanhP(std::vector<double> &pp, std::vector<double> &rtanhp);
    void getTanhQ(std::vector<double> &pq, std::vector<double> &rtanhq);
    void getTanh_Pnl(const int ikernel, std::vector<double> &ppnl, std::vector<double> &rtanh_pnl);
    void getTanh_Qnl(const int ikernel, std::vector<double> &pqnl, std::vector<double> &rtanh_qnl);
    void getTanhP_nl(const int ikernel, std::vector<double> &ptanhp, ModulePW::PW_Basis *pw_rho, std::vector<double> &rtanhp_nl);
    void getTanhQ_nl(const int ikernel, std::vector<double> &ptanhq, ModulePW::PW_Basis *pw_rho, std::vector<double> &rtanhq_nl);
    // 2023-03-20
    void getTanhXi_nl(const int ikernel, std::vector<double> &ptanhxi, ModulePW::PW_Basis *pw_rho, std::vector<double> &rtanhxi_nl);

    void getF_KS(
        psi::Psi<std::complex<double>> *psi,
        elecstate::ElecState *pelec,
        ModulePW::PW_Basis_K *pw_psi,
        ModulePW::PW_Basis *pw_rho,
        UnitCell& ucell,
        const std::vector<std::vector<double>> &nablaRho,
        std::vector<double> &rF,
        std::vector<double> &rpauli
    );
    // get intermediate variables of V_Pauli
    void getNablaRho(const double * const *prho, ModulePW::PW_Basis *pw_rho, std::vector<std::vector<double>> &rnablaRho);

    // tools
    double MLkernel(double eta, double tf_weight, double vw_weight);
    double MLkernel_yukawa(double eta, double alpha);
    void read_kernel(const std::string &fileName, const double& scaling, ModulePW::PW_Basis *pw_rho, double* kernel_);
    void multiKernel(const int ikernel, double *pinput, ModulePW::PW_Basis *pw_rho, double *routput);
    void Laplacian(double * pinput, ModulePW::PW_Basis *pw_rho, double * routput);
    void divergence(double ** pinput, ModulePW::PW_Basis *pw_rho, double * routput);

    void tanh(std::vector<double> &pinput, std::vector<double> &routput, double chi=1.);
    double dtanh(double tanhx, double chi=1.);

    // new parameters 2023-02-13
    std::vector<double> chi_xi = {1.0};
    double chi_p = 1.;
    double chi_q = 1.;
    std::vector<double> chi_pnl = {1.0};
    std::vector<double> chi_qnl = {1.0};

    int nx = 0;
    double dV = 0.;
    double rho0 = 0.; // average rho
    double kF = 0.; // Fermi vector kF = (3 pi^2 rho0)^(1/3)
    double tkF = 0.; // 2 * kF
    // double weightml = 1.;
    const double cTF = 3.0/10.0 * std::pow(3*std::pow(M_PI, 2.0), 2.0/3.0) * 2; // 10/3*(3*pi^2)^{2/3}, multiply by 2 to convert unit from Hartree to Ry, finally in Ry*Bohr^(-2)
    const double pqcoef = 1.0 / (4.0 * std::pow(3*std::pow(M_PI, 2.0), 2.0/3.0)); // coefficient of p and q
    
    int nkernel = 1;
    std::vector<int> kernel_type = {1};
    std::vector<double> kernel_scaling = {1.0};
    std::vector<double> yukawa_alpha = {1.0};
    std::vector<std::string> kernel_file = {"none"};
    std::vector<std::vector<double>> kernel = {}; // kernel[ikernel][ipw] = kernel value for ikernel and ipw
};

} // namespace ModuleIO

#endif