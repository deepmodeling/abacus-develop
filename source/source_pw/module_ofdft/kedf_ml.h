#ifndef KEDF_ML_H
#define KEDF_ML_H

#ifdef __MLALGO

#include "ml_base.h"

class KEDF_ML : public ML_Base
{
public:
    KEDF_ML()
    {
        this->energy_prefactor = 3. /10. * std::pow(3*std::pow(M_PI, 2.0), 2.0/3.0) * 2; // 10/3*(3*pi^2)^{2/3}, multiply by 2 to convert unit from Hartree to Ry, finally in Ry*Bohr^(-2)
        this->energy_exponent = 5. / 3.;
        // this->stress.create(3,3);
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

    double get_energy(const double * const * prho, ModulePW::PW_Basis *pw_rho);
    // double get_energy_density(const double * const *prho, int is, int ir, ModulePW::PW_Basis *pw_rho);
    void ml_potential(const double * const * prho, ModulePW::PW_Basis *pw_rho, ModuleBase::matrix &rpotential);
    // void get_stress(double cellVol, const double * const * prho, ModulePW::PW_Basis *pw_rho, double vw_weight);

    // output all parameters
    void generateTrainData(const double * const *prho, ModulePW::PW_Basis *pw_rho, const double *veff);
    void localTest(const double * const *prho, ModulePW::PW_Basis *pw_rho);

    double ml_energy = 0.;
    // ModuleBase::matrix stress;
    
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
};

#endif
#endif