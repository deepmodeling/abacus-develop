#ifndef EVOLVE_PHI_H
#define EVOLVE_PHI_H
#include <math.h>
#include <stdio.h>

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/matrix.h"
#include "source_base/timer.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_estate/elecstate.h" // electronic states
#include "source_estate/module_charge/charge.h"

/**
 * @brief TDOFDFT
 * @author liyuanbo on 2025-09
 */
class EVOLVE_PHI
{
  public:
    EVOLVE_PHI()
    {
    }
    ~EVOLVE_PHI()
    {
    }
    void propagate_psi(elecstate::ElecState* pelec, const Charge& chr, UnitCell& ucell, std::complex<double>** pphi_, ModulePW::PW_Basis* pw_rho);

  private:
    const double c_tf_
        = 3.0 / 10.0 * std::pow(3 * std::pow(M_PI, 2.0), 2.0 / 3.0)
          * 2; // 10/3*(3*pi^2)^{2/3}, multiply by 2 to convert unit from Hartree to Ry, finally in Ry*Bohr^(-2)

    void get_Hpsi(elecstate::ElecState* pelec, const Charge& chr, UnitCell& ucell, const std::complex<double>* const* psi_, ModulePW::PW_Basis* pw_rho, std::complex<double>** Hpsi);
    void get_tf_potential(const double* const* prho, ModulePW::PW_Basis* pw_rho, ModuleBase::matrix& rpot);
    void get_vw_potential_phi(const std::complex<double>* const* pphi, ModulePW::PW_Basis* pw_rho, std::complex<double>** Hpsi); // -1/2 \nabla^2 \phi
    void get_CD_potential(const std::complex<double>* const* psi_, ModulePW::PW_Basis* pw_rho, ModuleBase::matrix& rpot);
 
};
#endif