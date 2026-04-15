#ifndef KEDF_EXTWT_H
#define KEDF_EXTWT_H
#include <cmath>
#include <cstdio>

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/matrix.h"
#include "source_base/timer.h"
#include "source_basis/module_pw/pw_basis.h"

/**
 * @brief A class which calculates the kinetic energy, potential, and stress with extended Wang-Teter (ext-WT) KEDF.
 * See Sun L, Chen M. Physical Review B, 2026, 113(16): L161107.
 * @author sunliang on 2026-04
 */
class KEDF_ExtWT
{
  public:
    KEDF_ExtWT()
    {
        this->stress.create(3, 3);
    }
    ~KEDF_ExtWT()
    {
        delete[] this->kernel_;
        delete[] this->dkernel_deta_;
    }

    void set_para(double dV,
                  double alpha,
                  double beta,
                  double nelec,
                  double tf_weight,
                  double vw_weight,
                  double of_extwt_kappa,
                  ModulePW::PW_Basis* pw_rho);

    void update_rho0(const double* const* prho,
                     ModulePW::PW_Basis* pw_rho);
    void cal_kernel(double tf_weight,
                    double vw_weight,
                    double rho0,
                    ModulePW::PW_Basis* pw_rho);

    void update_dkernel_deta(const double &vw_weight, ModulePW::PW_Basis* pw_rho);

    double get_energy(const double* const* prho, ModulePW::PW_Basis* pw_rho);
    double get_energy_density(const double* const* prho, int is, int ir, ModulePW::PW_Basis* pw_rho);
    void tau_extwt(const double* const* prho, ModulePW::PW_Basis* pw_rho, double* rtau_extwt);
    void extwt_potential(const double* const* prho, ModulePW::PW_Basis* pw_rho, ModuleBase::matrix& rpotential);
    void get_stress(const double* const* prho, ModulePW::PW_Basis* pw_rho, double vw_weight);
    double extwt_energy = 0.;
    ModuleBase::matrix stress;

    double rho0_ = 0.; // average rho

  private:
    double extwt_kernel(double eta, double tf_weight, double vw_weight);
    double diff_linhard(double eta, double vw_weight);
    void multi_kernel(const double* const* prho, const double* kernel, double** rkernel_rho, double exponent, ModulePW::PW_Basis* pw_rho);
    void fill_kernel(double tf_weight, double vw_weight, ModulePW::PW_Basis* pw_rho);

    double dV_ = 0.;
    double kf_ = 0.;  // Fermi vector kF = (3 pi^2 rho)^(1/3)
    double tkf_ = 0.; // 2 * kF
    double alpha_ = 5. / 6.;
    double beta_ = 5. / 6.;
    // double weightWT = 1.;
    const double c_tf_
        = 3.0 / 10.0 * std::pow(3 * std::pow(M_PI, 2.0), 2.0 / 3.0)
          * 2; // 10/3*(3*pi^2)^{2/3}, multiply by 2 to convert unit from Hartree to Ry, finally in Ry*Bohr^(-2)
    double wt_coef_ = 0.; // coefficient of WT kernel
    double* kernel_ = nullptr;
    double* dkernel_deta_ = nullptr; // \partial w/ \partial rho0 = coef * ((alpha + beta - 5/3)F(eta) + 1/3 eta F'(eta))
    double sum_rho_kappa_ = 0.;
    double sum_rho_kappa_plus_one_ = 0.;
    double kappa_ = 1.0 / (2.0 * std::pow(4./3., 1./3.) - 1.0);
};
#endif