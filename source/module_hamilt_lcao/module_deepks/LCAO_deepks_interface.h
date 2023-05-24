#ifndef LCAO_DEEPKS_INTERFACE_H
#define LCAO_DEEPKS_INTERFACE_H

#ifdef __DEEPKS
#include "LCAO_deepks.h"
#include "module_base/complexmatrix.h"
#include "module_base/matrix.h"

class LCAO_Deepks_Interface
{
  public:
    LCAO_Deepks_Interface(LCAO_Deepks* ld_in);
    ~LCAO_Deepks_Interface();
    LCAO_Deepks* ld = nullptr;
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
                           const psi::Psi<double>& psid,
                           const std::vector<ModuleBase::matrix>& dm_gamma,
                           const std::vector<ModuleBase::ComplexMatrix>& dm_k);
};

#endif
#endif