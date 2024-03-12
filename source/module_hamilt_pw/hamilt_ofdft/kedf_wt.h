#ifndef KEDF_WT_H
#define KEDF_WT_H
#include <math.h>
#include <stdio.h>

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_base/timer.h"
#include "module_basis/module_pw/pw_basis.h"

/**
 * @brief A class which calculates the kinetic energy, potential, and stress with Wang-Teter (WT) KEDF.
 * See Wang L W, Teter M P. Physical Review B, 1992, 45(23): 13196.
 * @author sunliang on 2022-06
 */
class KEDF_WT
{
  public:
    KEDF_WT()
    {
        this->stress.create(3, 3);
    }
    ~KEDF_WT()
    {
        delete[] this->kernel_;
    }

    void set_para(double dV,
                  double alpha,
                  double beta,
                  double nelec,
                  double tf_weight,
                  double vw_weight,
                  double of_wt_rho0,
                  bool of_hold_rho0,
                  bool read_kernel,
                  std::string kernel_file,
                  ModulePW::PW_Basis* pw_rho);

    double get_energy(const double* const* prho, ModulePW::PW_Basis* pw_rho);
    double get_energy_density(const double* const* prho, int is, int ir, ModulePW::PW_Basis* pw_rho);
    void wt_potential(const double* const* prho, ModulePW::PW_Basis* pw_rho, ModuleBase::matrix& rpotential);
    void get_stress(const double* const* prho, ModulePW::PW_Basis* pw_rho, double vw_weight);
    double wt_energy = 0.;
    ModuleBase::matrix stress;

  private:
    double wt_kernel(double eta, double tf_weight, double vw_weight);
    double diff_linhard(double eta, double vw_weight);
    void multi_kernel(const double* const* prho, double** rkernel_rho, double exponent, ModulePW::PW_Basis* pw_rho);
    void read_kernel(std::string file_name, ModulePW::PW_Basis* pw_rho);
    void fill_kernel(double tf_weight, double vw_weight, ModulePW::PW_Basis* pw_rho);

    double dV_ = 0.;
    double rho0_ = 0.; // average rho
    bool hold_rho0_ = false;
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
};
#endif