#ifndef KEDF_XWM_H
#define KEDF_XWM_H

#include "source_base/matrix.h"
#include "source_base/timer.h"
#include "source_basis/module_pw/pw_basis.h"

/**
 * @brief A class which calculates the kinetic energy, potential, and stress with Xie-Wang-Morales (XWM) KEDF.
 * See Xu Q, Wang Y, Ma Y. Physical Review B, 2019, 100(20): 205132.
 * @author sunliang on 2025-01
 */
class KEDF_XWM
{
  public:
    KEDF_XWM()
    {
        this->stress.create(3, 3);
    }
    ~KEDF_XWM(){}

    void set_para(double dV,
                  double rho_ref,
                  double kappa,
                  double nelec,
                  double tf_weight,
                  double vw_weight,
                  ModulePW::PW_Basis* pw_rho);


    double get_energy(const double* const* prho, ModulePW::PW_Basis* pw_rho);
    double get_energy_density(const double* const* prho, int is, int ir, ModulePW::PW_Basis* pw_rho);
    void tau_xwm(const double* const* prho, ModulePW::PW_Basis* pw_rho, double* rtau_xwm);
    void xwm_potential(const double* const* prho, ModulePW::PW_Basis* pw_rho, ModuleBase::matrix& rpotential);
    void get_stress(const double* const* prho, ModulePW::PW_Basis* pw_rho, double vw_weight);
    double xwm_energy = 0.;
    ModuleBase::matrix stress;

  private:
    void multi_kernel(const double* const* prho, const double* kernel, double** rkernel_rho, double exponent, ModulePW::PW_Basis* pw_rho);
    void fill_kernel(double tf_weight, double vw_weight, ModulePW::PW_Basis* pw_rho);

    double dV_ = 0.;
    double rho0_ = 0.; // average rho
    double rho_ref_ = 0.; // reference rho
    double kf_ = 0.;  // Fermi vector kF = (3 pi^2 rho_star_)^(1/3)
    double tkf_ = 0.; // 2 * kF
    double kappa_ = 0.;
    double c_kernel = std::pow(ModuleBase::PI, 4./3.) / std::pow(3., 1. / 3.) * 2.; // multiply by 2 to convert unit from Hartree to Ry
    double c_0 = 0.; // coef of T0, c_0 = 18 / (6 * kappa + 5) ^ 2
    double c_1 = 0.; // coef of T1, c_1 = [(kappa + 5/6)(kappa + 11/6)] ^ -1
    double c_2 = 0.; // coef of T1, c_2 = - (kappa + 5/6) ^ -2 * rho_ref_
    double kappa_5_6 = 0.; // kappa + 5/6
    double kappa_11_6 = 0.; // kappa + 11/6
    double kappa_1_6 = 0.; // kappa - 1/6
    std::vector<double> kernel1_ = {}; // w1 in PHYSICAL REVIEW B 110, 085113 (2024)
    std::vector<double> kernel2_ = {}; // w2 in PHYSICAL REVIEW B 110, 085113 (2024)
};
#endif