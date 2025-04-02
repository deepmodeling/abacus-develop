#ifndef ISTATE_ENVELOPE_H
#define ISTATE_ENVELOPE_H
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_cell/klist.h"
#include "module_elecstate/elecstate.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_psi/psi.h"

#include <stdexcept>
class IState_Envelope
{
  public:
    IState_Envelope(const elecstate::ElecState* pes);
    ~IState_Envelope();

    /// For gamma_only
    void begin(const UnitCell& ucell,
               const psi::Psi<double>* psid,
               const ModulePW::PW_Basis* pw_rhod,
               const ModulePW::PW_Basis_K* pw_wfc,
               const ModulePW::PW_Basis_Big* pw_big,
               const Parallel_Grid& pgrid,
               const Parallel_Orbitals& para_orb,
               Gint_Gamma& gg,
               const int& out_wfc_pw,
               const int& out_wfc_r,
               const K_Vectors& kv,
               const double nelec,
               const int nbands_istate,
               const std::vector<int>& out_wfc_norm,
               const std::vector<int>& out_wfc_re_im,
               const int nbands,
               const int nspin,
               const int nlocal,
               const std::string& global_out_dir);

    /// tmp, delete after Gint is refactored.
    void begin(const UnitCell& ucell,
               const psi::Psi<double>* psid,
               const ModulePW::PW_Basis* pw_rhod,
               const ModulePW::PW_Basis_K* pw_wfc,
               const ModulePW::PW_Basis_Big* pw_big,
               const Parallel_Grid& pgrid,
               const Parallel_Orbitals& para_orb,
               Gint_k& gg,
               const int& out_wfc_pw,
               const int& out_wfc_r,
               const K_Vectors& kv,
               const double nelec,
               const int nbands_istate,
               const std::vector<int>& out_wfc_norm,
               const std::vector<int>& out_wfc_re_im,
               const int nbands,
               const int nspin,
               const int nlocal,
               const std::string& global_out_dir)
    {
        throw std::logic_error("gint_k should use with complex psi.");
    };

    /// For multi-k
    void begin(const UnitCell& ucell,
               const psi::Psi<std::complex<double>>* psi,
               const ModulePW::PW_Basis* pw_rhod,
               const ModulePW::PW_Basis_K* pw_wfc,
               const ModulePW::PW_Basis_Big* pw_big,
               const Parallel_Grid& pgrid,
               const Parallel_Orbitals& para_orb,
               Gint_k& gk,
               const int& out_wfc_pw,
               const int& out_wfc_r,
               const K_Vectors& kv,
               const double nelec,
               const int nbands_istate,
               const std::vector<int>& out_wfc_norm,
               const std::vector<int>& out_wfc_re_im,
               const int nbands,
               const int nspin,
               const int nlocal,
               const std::string& global_out_dir);

    /// tmp, delete after Gint is refactored.
    void begin(const UnitCell& ucell,
               const psi::Psi<std::complex<double>>* psi,
               const ModulePW::PW_Basis* pw_rhod,
               const ModulePW::PW_Basis_K* pw_wfc,
               const ModulePW::PW_Basis_Big* pw_big,
               const Parallel_Grid& pgrid,
               const Parallel_Orbitals& para_orb,
               Gint_Gamma& gk,
               const int& out_wfc_pw,
               const int& out_wfc_r,
               const K_Vectors& kv,
               const double nelec,
               const int nbands_istate,
               const std::vector<int>& out_wfc_norm,
               const std::vector<int>& out_wfc_re_im,
               const int nbands,
               const int nspin,
               const int nlocal,
               const std::string& global_out_dir)
    {
        throw std::logic_error("gint_gamma should use with real psi.");
    };

  private:
    void select_bands(const int nbands_istate,
                      const std::vector<int>& out_wfc_kb,
                      const int nbands,
                      const double nelec,
                      const int mode,
                      const int fermi_band);

    void set_pw_wfc(const ModulePW::PW_Basis_K* pw_wfc,
                    const int& ik,
                    const int& ib,
                    const int& nspin,
                    const double* const* const rho,
                    psi::Psi<std::complex<double>>& wfc_g);

    int globalIndex(int localindex, int nblk, int nprocs, int myproc);

    int localIndex(int globalindex, int nblk, int nprocs, int& myproc);

#ifdef __MPI
    template <typename T>
    int set_wfc_grid(const int naroc[2],
                     const int nb,
                     const int dim0,
                     const int dim1,
                     const int iprow,
                     const int ipcol,
                     const T* in,
                     T** out,
                     const std::vector<int>& trace_lo);
    template <typename T>
    void wfc_2d_to_grid(const T* wfc_2d, const Parallel_Orbitals& pv, T** wfc_grid, const std::vector<int>& trace_lo);
#endif

    std::vector<int> bands_picked_;
    const elecstate::ElecState* pes_ = nullptr;
};
#endif
