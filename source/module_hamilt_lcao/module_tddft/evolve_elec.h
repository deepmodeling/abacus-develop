#ifndef EVOLVE_ELEC_H
#define EVOLVE_ELEC_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_esolver/esolver_ks_lcao.h"
#include "module_esolver/esolver_ks_lcao_tddft.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_psi/psi.h"

//-----------------------------------------------------------
// mohan add 2021-02-09
// This class is used to evolve the electronic wave functions
// in TDDFT in terms of the multiple k points
// k is the index for the points in the first Brillouin zone
//-----------------------------------------------------------

namespace module_tddft
{
class Evolve_elec
{

    friend class ELEC_scf;
    friend class ModuleESolver::ESolver_KS_LCAO;
    friend class ModuleESolver::ESolver_KS_LCAO_TDDFT;

  public:
    Evolve_elec();
    ~Evolve_elec();

    static double td_force_dt;
    static bool td_vext;
    static std::vector<int> td_vext_dire_case;
    static bool out_dipole;
    static bool out_efield;

    static double td_print_eij; // the threshold to output Eij elements
    static int td_edm;          // 0: new edm method   1: old edm method

  private:
    static void solve_psi(const int& istep,
                          const int nband,
                          const int nlocal,
                          hamilt::Hamilt<double>* phm,
                          Local_Orbital_wfc& lowf,
                          psi::Psi<std::complex<double>>* psi,
                          psi::Psi<std::complex<double>>* psi_laststep,
                          std::complex<double>** Hk_laststep,
                          std::complex<double>** Sk_laststep,
                          ModuleBase::matrix& ekb,
                          int htype,
                          int propagator,
                          const int& nks);
};
} // namespace module_tddft
#endif
