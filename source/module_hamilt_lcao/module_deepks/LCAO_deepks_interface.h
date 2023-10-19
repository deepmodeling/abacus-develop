#ifndef LCAO_DEEPKS_INTERFACE_H
#define LCAO_DEEPKS_INTERFACE_H

#ifdef __DEEPKS
#include "LCAO_deepks.h"
#include "module_base/complexmatrix.h"
#include "module_base/matrix.h"
#include <memory>

class LCAO_Deepks_Interface
{
  public:
    /// @brief Constructor for LCAO_Deepks_Interface
    /// @param ld_in
    LCAO_Deepks_Interface(std::shared_ptr<LCAO_Deepks> ld_in);
    /// @brief output deepks-related labels, descriptors and energy corrections
    /// @param[in] etot
    /// @param[in] nks
    /// @param[in] nat
    /// @param[in] ekb
    /// @param[in] kvec_d
    /// @param[in] ucell
    /// @param[in] orb
    /// @param[in] GridD
    /// @param[in] ParaV
    /// @param[in] psi
    /// @param[in] psid
    /// @param[in] dm_gamma
    /// @param[in] dm_k
    // for Gamma-only
    void out_deepks_labels(double etot,
                           int nks,
                           int nat,
                           const ModuleBase::matrix& ekb,
                           const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                           const UnitCell& ucell,
                           const LCAO_Orbitals& orb,
                           Grid_Driver& GridD,
                           const Parallel_Orbitals* ParaV,
                           const psi::Psi<double>& psid,
                           const std::vector<std::vector<double>>& dm_gamma);
  // for multi-k
  void out_deepks_labels(double etot,
                           int nks,
                           int nat,
                           const ModuleBase::matrix& ekb,
                           const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                           const UnitCell& ucell,
                           const LCAO_Orbitals& orb,
                           Grid_Driver& GridD,
                           const Parallel_Orbitals* ParaV,
                           const psi::Psi<std::complex<double>>& psi,
                           const std::vector<std::vector<std::complex<double>>>& dm_k);

  private:
    std::shared_ptr<LCAO_Deepks> ld;
};

#endif
#endif