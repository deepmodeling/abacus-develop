#ifndef KEDF_LKT_H
#define KEDF_LKT_H
#include <math.h>
#include <stdio.h>

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_base/timer.h"
#include "module_basis/module_pw/pw_basis.h"

/**
 * @brief A class which calculates the kinetic energy, potential, and stress with Luo-Karasiev-Trickey (LKT) KEDF.
 * See Luo K, Karasiev V V, Trickey S B. Physical Review B, 2018, 98(4): 041111.
 * @author sunliang on 2023-04-28
 */
class KEDF_LKT
{
  public:
    KEDF_LKT()
    {
        this->stress.create(3, 3);
    }
    ~KEDF_LKT(){}

    void set_para(double dV, double lkt_a);

    double get_energy(const double *const *prho, ModulePW::PW_Basis *pw_rho);
    double get_energy_density(const double *const *prho, int is, int ir, ModulePW::PW_Basis *pw_rho);
    void lkt_potential(const double *const *prho, ModulePW::PW_Basis *pw_rho, ModuleBase::matrix &rpotential);
    void get_stress(const double *const *prho, ModulePW::PW_Basis *pw_rho);

    double lkt_energy = 0.; // LKT energy
    ModuleBase::matrix stress;

  private:
    void nabla(const double *pinput, ModulePW::PW_Basis *pw_rho, double **routput);
    void divergence(const double *const *pinput, ModulePW::PW_Basis *pw_rho, double *routput);
    void get_as(const double *prho, const double *const *pnabla_rho, const int nrxx, double *as);

    double dV_ = 0.; // volume element = V/nxyz
    const double c_tf_
        = 3.0 / 10.0 * std::pow(3 * std::pow(M_PI, 2.0), 2.0 / 3.0)
          * 2; // 10/3*(3*pi^2)^{2/3}, multiply by 2 to convert unit from Hartree to Ry, finally in Ry*Bohr^(-2)
    const double s_coef_
        = 1.0 / (2. * std::pow(3 * std::pow(M_PI, 2.0), 1.0 / 3.0)); // coef of s, s=s_coef * |nabla rho|/rho^{4/3}
    double lkt_a_ = 1.3;
};
#endif