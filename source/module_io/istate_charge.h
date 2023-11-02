#ifndef ISTATE_CHARGE_H
#define ISTATE_CHARGE_H
#include <module_base/complexmatrix.h>
#include <module_base/matrix.h>

#include <stdexcept>
#include <vector>

#include "module_basis/module_pw/pw_basis.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_psi/psi.h"

class IState_Charge
{
  public:
    IState_Charge(psi::Psi<double>* psi_gamma_in, Local_Orbital_Charge& loc_in);
    IState_Charge(psi::Psi<std::complex<double>>* psi_k_in, Local_Orbital_Charge& loc_in)
    {
        throw std::logic_error("IState_Charge for multi-k is not implemented.");
    };

    ~IState_Charge();

    void begin(Gint_Gamma& gg,
               elecstate::ElecState* pelec,
               const ModulePW::PW_Basis* rhopw,
               const ModulePW::PW_Basis_Big* bigpw,
               const bool gamma_only_local,
               const int nbands_istate,
               const int nbands,
               const double nelec,
               const int nspin,
               const std::string& global_out_dir,
               const int my_rank,
               std::ofstream& ofs_warning);

  private:
    int* bands_picked;

#ifdef __MPI
    void idmatrix(const int& ib, elecstate::ElecState* pelec, const int nspin, const double nelec, const int nlocal);
#endif
    psi::Psi<double>* psi_gamma;
    Local_Orbital_Charge* loc;
};
#endif
