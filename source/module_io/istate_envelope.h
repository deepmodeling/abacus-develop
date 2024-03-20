#ifndef ISTATE_ENVELOPE_H
#define ISTATE_ENVELOPE_H
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_cell/klist.h"
#include "module_elecstate/elecstate.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_wfc.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_psi/psi.h"
#include <stdexcept>
class IState_Envelope
{
public:
    IState_Envelope(const elecstate::ElecState* pes_in);
    ~IState_Envelope();

    /// for gamma_only
    void begin(const psi::Psi<double>* psid,
               const ModulePW::PW_Basis* rhopw,
               const ModulePW::PW_Basis_K* wfcpw,
               const ModulePW::PW_Basis_Big* bigpw,
               Local_Orbital_wfc& lowf,
               Gint_Gamma& gg,
               int& out_wfc_pw,
               int& out_wfc_r,
        const K_Vectors& kv);


    /// tmp, delete after Gint is refactored.
    void begin(const psi::Psi<double>* psid,
        const ModulePW::PW_Basis* rhopw,
        const ModulePW::PW_Basis_K* wfcpw,
        const ModulePW::PW_Basis_Big* bigpw,
        Local_Orbital_wfc& lowf,
        Gint_k& gg,
        int& out_wfc_pw,
        int& out_wfc_r,
        const K_Vectors& kv) {
        throw std::logic_error("gint_k should use with complex psi.");
    };
    /// for multi-k
    void begin(const psi::Psi<std::complex<double>>* psi,
               const ModulePW::PW_Basis* rhopw,
               const ModulePW::PW_Basis_K* wfcpw,
               const ModulePW::PW_Basis_Big* bigpw,
               Local_Orbital_wfc& lowf,
               Gint_k& gk,
               int& out_wfc_pw,
               int& out_wfc_r,
               const K_Vectors& kv);

    /// tmp, delete after Gint is refactored.
    void begin(const psi::Psi<std::complex<double>>* psi,
        const ModulePW::PW_Basis* rhopw,
        const ModulePW::PW_Basis_K* wfcpw,
        const ModulePW::PW_Basis_Big* bigpw,
        Local_Orbital_wfc& lowf,
        Gint_Gamma& gk,
        int& out_wfc_pw,
        int& out_wfc_r,
        const K_Vectors& kv) {
        throw std::logic_error("gint_gamma should use with real psi.");
    };

  private:
    bool* bands_picked = nullptr;
    const elecstate::ElecState* pes = nullptr;

    void set_pw_wfc(const ModulePW::PW_Basis_K* wfcpw,
                    const int& ik,
                    const int& ib,
                    const int& nspin,
                    const double* const* const rho,
                    psi::Psi<std::complex<double>>& wfc_g);
};
#endif
