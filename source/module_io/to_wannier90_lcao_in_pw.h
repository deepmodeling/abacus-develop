#ifndef TOWannier90_LCAO_IN_PW_H
#define TOWannier90_LCAO_IN_PW_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "to_wannier90.h"
#include "to_wannier90_pw.h"

#include "module_base/complexmatrix.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/lapack_connector.h"
#include "module_base/matrix.h"
#include "module_base/matrix3.h"
#include "module_cell/klist.h"
#include "module_hamilt_lcao/hamilt_lcaodft/wavefunc_in_pw.h"
#include "module_psi/psi.h"

#include "module_basis/module_ao/ORB_table_phi.h"
#include "module_basis/module_ao/ORB_gaunt_table.h"
#include "module_basis/module_ao/ORB_atomic_lm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_base/vector3.h"
#include "module_base/abfs-vector3_order.h"
#include "module_base/ylm.h"
#include "single_R_io.h"

#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

#ifdef __LCAO
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_wfc.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"


class toWannier90_LCAO_IN_PW : public toWannier90_PW
{
  public:
    toWannier90_LCAO_IN_PW(
      const bool &out_wannier_mmn, 
      const bool &out_wannier_amn, 
      const bool &out_wannier_unk, 
      const bool &out_wannier_eig,
      const bool &out_wannier_wvfn_formatted, 
      const std::string &nnkpfile,
      const std::string &wannier_spin
    );
    ~toWannier90_LCAO_IN_PW();

    void calculate(
      const ModuleBase::matrix& ekb,
      const ModulePW::PW_Basis_K* wfcpw,
      const ModulePW::PW_Basis_Big* bigpw,
      const Structure_Factor& sf,
      const K_Vectors& kv,
      const psi::Psi<std::complex<double>>* psi,
      const Parallel_Orbitals *pv
    );

    void calculate(
      const ModuleBase::matrix& ekb,
      const ModulePW::PW_Basis_K* wfcpw,
      const ModulePW::PW_Basis_Big* bigpw,
      const Structure_Factor& sf,
      const K_Vectors& kv,
      const psi::Psi<double>* psi,
      const Parallel_Orbitals *pv
    )
    {
      throw std::logic_error("The wave function of toWannier90_LCAO_IN_PW is generally a std::complex<double> type.");
    }

  protected:
    const Parallel_Orbitals* ParaV;

    psi::Psi<std::complex<double>>* get_unk_from_lcao(
      const psi::Psi<std::complex<double>>& psi_in, 
      const ModulePW::PW_Basis_K* wfcpw, 
      const Structure_Factor& sf, 
      const K_Vectors& kv
    );

    void produce_local_basis_in_pw(
      const int& ik,
      const ModulePW::PW_Basis_K* wfc_basis,
      const Structure_Factor& sf,
      ModuleBase::ComplexMatrix& psi,
      const ModuleBase::realArray& table_local
    );

    void get_lcao_wfc_global_ik(const int ik, const psi::Psi<std::complex<double>>& psi_in, ModuleBase::ComplexMatrix &lcao_wfc_global);
};
#endif

#endif
