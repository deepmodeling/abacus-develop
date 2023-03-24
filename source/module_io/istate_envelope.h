#ifndef ISTATE_ENVELOPE_H
#define ISTATE_ENVELOPE_H
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_wfc.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_psi/psi.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_elecstate/elecstate.h"

class IState_Envelope
{
public:
    IState_Envelope(const elecstate::ElecState* pes_in);
    ~IState_Envelope();

    /// for gamma_only
    void begin(const psi::Psi<double>* psid, Local_Orbital_wfc& lowf, Gint_Gamma& gg, int& out_wfc_pw, int& out_wfc_r);
    /// for multi-k
    void begin(const psi::Psi<std::complex<double>>* psi, Local_Orbital_wfc& lowf, Gint_k& gk, int& out_wfc_pw, int& out_wfc_r);

private:
    bool* bands_picked = nullptr;
    const elecstate::ElecState* pes = nullptr;

    void set_pw_wfc(ModulePW::PW_Basis_K* wfc_basis,
        const int& ik, const int& ib,
        const int& nspin,
        const double* const* const rho,
        psi::Psi<std::complex<double>> &wfc_g);

};
#endif
