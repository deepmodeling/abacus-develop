#ifndef TD_CURRENT_PW_H
#define TD_CURRENT_PW_H

#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_estate/elecstate.h"
#include "source_psi/psi.h"
#include "source_pw/module_pwdft/op_pw_vel.h"
#include "source_pw/module_pwdft/vnl_pw.h"

#include <memory>
#include <string>

namespace ModuleIO
{
/** @brief Reduced current densities in Hartree atomic units. */
struct PWCurrentResult
{
    ModuleBase::Vector3<double> total;
    std::vector<double> per_k;
};

/** @brief Write already reduced current densities; performs no collective operations. */
void write_pw_current(const PWCurrentResult& current,
                      const int istep,
                      const int nk_per_spin,
                      const bool out_current_k,
                      const std::string& out_dir);

/** @brief Persistent one-k-point current workspace for a propagation run. */
template <typename FPTYPE, typename Device>
class CurrentPW
{
  private:
    std::unique_ptr<hamilt::Velocity<FPTYPE, Device>> velocity_;
    ct::Tensor vpsi_;
    ct::Tensor dots_;
    std::vector<std::complex<FPTYPE>> band_current_;

  public:
    /** @brief Compute globally reduced current densities without file output. */
    PWCurrentResult calculate(const UnitCell& ucell,
                              const ModulePW::PW_Basis_K* wfcpw,
                              psi::Psi<std::complex<FPTYPE>, Device>* psi,
                              const elecstate::ElecState* pelec,
                              const K_Vectors& kv,
                              pseudopot_cell_vnl* ppcell,
                              const int gauge,
                              const ModuleBase::Vector3<double>& A_right_ha);
    /** @brief Evaluate and write the current using native wavefunctions. */
    void write(const int istep,
               const UnitCell& ucell,
               const ModulePW::PW_Basis_K* wfcpw,
               psi::Psi<std::complex<FPTYPE>, Device>* psi,
               const elecstate::ElecState* pelec,
               const K_Vectors& kv,
               pseudopot_cell_vnl* ppcell,
               const int gauge,
               const ModuleBase::Vector3<double>& vector_potential,
               const bool out_current_k,
               const std::string& out_dir,
               const int world_rank);
};
} // namespace ModuleIO
#endif
