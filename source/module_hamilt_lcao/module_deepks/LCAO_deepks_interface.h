#ifndef LCAO_DEEPKS_INTERFACE_H
#define LCAO_DEEPKS_INTERFACE_H

#ifdef __DEEPKS
#include "LCAO_deepks.h"

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
                             const Parallel_Orbitals& ParaV,
                             psi::Psi<std::complex<double>>* psi,
                             psi::Psi<double>* psid);
    //void print_projected_dm();
    //void print_descriptors();
    //void perform_deepks_scf();
};

#endif
#endif