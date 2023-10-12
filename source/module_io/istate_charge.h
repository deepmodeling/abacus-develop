#ifndef ISTATE_CHARGE_H
#define ISTATE_CHARGE_H
#include <module_base/complexmatrix.h>
#include <module_base/matrix.h>

#include <vector>

#include "module_basis/module_pw/pw_basis.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_psi/psi.h"
#include <stdexcept>

class IState_Charge
{
public:
    IState_Charge(psi::Psi<double>* psi_gamma_in,
        Local_Orbital_Charge& loc_in);
    IState_Charge(psi::Psi<std::complex<double>>* psi_k_in,
        Local_Orbital_Charge& loc_in)
    {
        throw std::logic_error("IState_Charge for multi-k is not implemented.");
    };

    ~IState_Charge();

    void begin(Gint_Gamma& gg,
               elecstate::ElecState* pelec,
               const ModulePW::PW_Basis* rhopw,
               const ModulePW::PW_Basis_Big* bigpw);

  private:
    int *bands_picked;

#ifdef __MPI
	void idmatrix(const int &ib, elecstate::ElecState* pelec);
#endif
    psi::Psi<double>* psi_gamma;
    Local_Orbital_Charge* loc;

};
#endif
