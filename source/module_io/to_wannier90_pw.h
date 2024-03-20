#ifndef TOWannier90_PW_H
#define TOWannier90_PW_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "to_wannier90.h"

#include "module_base/complexmatrix.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/lapack_connector.h"
#include "module_base/matrix.h"
#include "module_base/matrix3.h"
#include "module_cell/klist.h"
#include "module_hamilt_lcao/hamilt_lcaodft/wavefunc_in_pw.h"
#include "module_psi/psi.h"

class toWannier90_PW : public toWannier90
{
  public:
    toWannier90_PW(
      const bool &out_wannier_mmn, 
      const bool &out_wannier_amn, 
      const bool &out_wannier_unk, 
      const bool &out_wannier_eig,
      const bool &out_wannier_wvfn_formatted, 
      const std::string &nnkpfile,
      const std::string &wannier_spin
    );
    ~toWannier90_PW();

    void calculate(
      const ModuleBase::matrix& ekb,
      const ModulePW::PW_Basis_K* wfcpw,
      const ModulePW::PW_Basis_Big* bigpw,
      const K_Vectors& kv,
      const psi::Psi<std::complex<double>>* psi
    );

    void calculate(
      const ModuleBase::matrix& ekb,
      const ModulePW::PW_Basis_K* wfcpw,
      const ModulePW::PW_Basis_Big* bigpw,
      const K_Vectors& kv,
      const psi::Psi<double>* psi
    )
    {
      throw std::logic_error("The wave function of toWannier90_PW is generally a std::complex<double> type.");
    }

    void cal_Amn(const psi::Psi<std::complex<double>>& psi_pw, const ModulePW::PW_Basis_K* wfcpw);
    void cal_Mmn(const psi::Psi<std::complex<double>>& psi_pw, const ModulePW::PW_Basis_K* wfcpw);
    void out_unk(
      const psi::Psi<std::complex<double>>& psi_pw,
      const ModulePW::PW_Basis_K* wfcpw,
      const ModulePW::PW_Basis_Big* bigpw
    );

  protected:
    // Radial section of trial orbitals
    const int mesh_r = 333;
    const double dx = 0.025;
    const double x_min = -6.0;

    void unkdotkb(
      const psi::Psi<std::complex<double>>& psi_pw, 
      const ModulePW::PW_Basis_K* wfcpw, 
      const int& ik, 
      const int& ikb, 
      const ModuleBase::Vector3<double> G, 
      ModuleBase::ComplexMatrix &Mmn
    );

    void gen_radial_function_in_q(std::vector<ModuleBase::matrix> &radial_in_q);

    void integral(
      const int meshr,
      const double *psir,
      const double *r,
      const double *rab,
      const int &l,
      double *table
    );

    void produce_trial_in_pw(
      const psi::Psi<std::complex<double>>& psi_pw,
      const int& ik,
      const ModulePW::PW_Basis_K* wfcpw,
      const std::vector<ModuleBase::matrix> &radial_in_q,
      ModuleBase::ComplexMatrix& trial_orbitals_k
    );

    void get_trial_orbitals_lm_k(
      const int &orbital_L,
      const int &orbital_m,
      const ModuleBase::matrix &ylm,
      const ModuleBase::Vector3<double> *gk,
      const int &npw,
      double *radial_in_q_single,
      std::complex<double> *orbital_in_G_single
    );

    void unkdotW_A(
      const psi::Psi<std::complex<double>>& psi_pw, 
      const ModulePW::PW_Basis_K* wfcpw, 
      const int& ik, 
      const std::vector<ModuleBase::matrix> &radial_in_q, 
      ModuleBase::ComplexMatrix &Amn
    );

};

#endif
