#ifndef DEEPKS_FORCE_H 
#define DEEPKS_FORCE_H 

#ifdef __DEEPKS


#include "module_base/complexmatrix.h"
#include "module_base/intarray.h"
#include "module_base/matrix.h"
#include "module_base/timer.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/module_dm/density_matrix.h"

#include <torch/script.h>
#include <torch/torch.h>
#include <unordered_map>

namespace DeePKS_domain
{
    //------------------------
    // LCAO_deepks_fgamma.cpp
    // LCAO_deepks_fk.cpp
    //------------------------

    // This file contains subroutines for calculating F_delta,
    // which is defind as sum_mu,nu rho_mu,nu d/dX (<chi_mu|alpha>V(D)<alpha|chi_nu>)

    // There are 3 subroutines in this file:
    // 1. cal_f_delta_gamma, which is used for gamma point calculation
    // 2. cal_f_delta_k, which is used for multi-k calculation
    // 3. check_f_delta, which prints F_delta into F_delta.dat for checking

    // for gamma only, pulay and HF terms of force are calculated together
void cal_f_delta_gamma(const std::vector<std::vector<double>>& dm,
                       const UnitCell& ucell,
                       const LCAO_Orbitals& orb,
                       const Grid_Driver& gd,
                       const Parallel_Orbitals& pv,
                       const int lmaxd,
                       std::vector<std::vector<std::unordered_map<int, std::vector<std::vector<double>>>>>& nlm_save,
                       double** gedm,
                       ModuleBase::IntArray* inl_index,
                       ModuleBase::matrix& f_delta,
                       const bool isstress,
                       ModuleBase::matrix& svnl_dalpha);

// for multi-k, pulay and HF terms of force are calculated together

typedef std::tuple<int, int, int, int> key_tuple;

void cal_f_delta_k(
    const std::vector<std::vector<std::complex<double>>>& dm, /**<[in] density matrix*/
    const UnitCell& ucell,
    const LCAO_Orbitals& orb,
    const Grid_Driver& GridD,
    const Parallel_Orbitals& pv,
    const int lmaxd,
    const int nks,
    const std::vector<ModuleBase::Vector3<double>>& kvec_d,
    std::vector<std::map<key_tuple, std::unordered_map<int, std::vector<std::vector<double>>>>>& nlm_save_k,
    double** gedm,
    ModuleBase::IntArray* inl_index,
    ModuleBase::matrix& f_delta,
    const bool isstress,
    ModuleBase::matrix& svnl_dalpha);

void check_f_delta(const int nat, ModuleBase::matrix& f_delta, ModuleBase::matrix& svnl_dalpha);
}

#endif
#endif
