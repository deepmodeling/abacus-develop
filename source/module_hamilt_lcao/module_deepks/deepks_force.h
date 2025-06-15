#ifndef DEEPKS_FORCE_H
#define DEEPKS_FORCE_H

#ifdef __MLALGO

#include "module_base/complexmatrix.h"
#include "module_base/intarray.h"
#include "module_base/matrix.h"
#include "module_base/timer.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/module_dm/density_matrix.h"

namespace DeePKS_domain
{
//------------------------
// deepks_force.cpp
//------------------------

// This file contains subroutines for calculating F_delta,
// which is defind as sum_mu,nu rho_mu,nu d/dX (<chi_mu|alpha>V(D)<alpha|chi_nu>)

// There are 1 subroutine in this file:
// 1. cal_f_delta, which is used for F_delta calculation

template <typename TK>
void cal_f_delta(const hamilt::HContainer<double>* dmr,
                 const UnitCell& ucell,
                 const LCAO_Orbitals& orb,
                 const Grid_Driver& GridD,
                 const Parallel_Orbitals& pv,
                 const int nks,
                 const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                 std::vector<hamilt::HContainer<double>*> phialpha,
                 double** gedm,
                 ModuleBase::IntArray* inl_index,
                 ModuleBase::matrix& f_delta,
                 const bool isstress,
                 ModuleBase::matrix& svnl_dalpha);
} // namespace DeePKS_domain

#endif
#endif
