#ifndef KEDF_TF_H
#define KEDF_TF_H
#include <math.h>
#include <stdio.h>

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_base/timer.h"

/**
 * @brief A class which calculates the kinetic energy, potential, and stress with Thomas-Fermi (TF) KEDF.
 * See Fermi E. Rend. Accad. Naz. Lincei, 1927, 6(602-607): 5.
 *     Thomas L H. Cambridge University Press, 1927, 23(5): 542-548.
 * @author sunliang on 2022-05
 */
class KEDF_TF
{
  public:
    KEDF_TF()
    {
        this->stress.create(3, 3);
    }
    ~KEDF_TF(){}

    void set_para(int nx, double dV, double tf_weight);

    double get_energy(const double *const *prho);
    double get_energy_density(const double *const *prho, int is, int ir);
    void tf_potential(const double *const *prho, ModuleBase::matrix &rpotential);
    void get_stress(double cell_vol);

    double tf_energy = 0.; // TF energy
    ModuleBase::matrix stress;

  private:
    int nx_ = 0;            // number of real space points in current core
    double dV_ = 0.;        // volume element = V/nxyz
    double tf_weight_ = 1.; // weight of TF KEDF
    const double c_tf_
        = 3.0 / 10.0 * std::pow(3 * std::pow(M_PI, 2.0), 2.0 / 3.0)
          * 2; // 10/3*(3*pi^2)^{2/3}, multiply by 2 to convert unit from Hartree to Ry, finally in Ry*Bohr^(-2)
};
#endif