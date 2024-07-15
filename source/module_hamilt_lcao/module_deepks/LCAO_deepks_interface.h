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
    /// @param[in] nlocal
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
    /// @param[in] deepks_v_delta
    // for Gamma-only
    void out_deepks_labels(const double& etot,
                           const int& nks,
                           const int& nat,
                           const int& nlocal,
                           const ModuleBase::matrix& ekb,
                           const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                           const UnitCell& ucell,
                           const LCAO_Orbitals& orb,
                           Grid_Driver& GridD,
                           const Parallel_Orbitals* ParaV,
                           const psi::Psi<double>& psid,
                           const elecstate::DensityMatrix<double, double>* dm,
                           const int& deepks_v_delta);
  // for multi-k
  void out_deepks_labels(const double& etot,
                           const int& nks,
                           const int& nat,
                           const int& nlocal,
                           const ModuleBase::matrix& ekb,
                           const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                           const UnitCell& ucell,
                           const LCAO_Orbitals& orb,
                           Grid_Driver& GridD,
                           const Parallel_Orbitals* ParaV,
                           const psi::Psi<std::complex<double>>& psi,
                           const elecstate::DensityMatrix<std::complex<double>, double>* dm,
                           const int& deepks_v_delta);

  private:
    std::shared_ptr<LCAO_Deepks> ld;
};

#endif
#endif