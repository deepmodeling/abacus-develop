#ifndef DEEPKS_PDM_H
#define DEEPKS_PDM_H

#ifdef __DEEPKS

#include "module_base/complexmatrix.h"
#include "module_base/matrix.h"
#include "module_base/timer.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"

#include <torch/script.h>
#include <torch/torch.h>

namespace DeePKS_domain
{
//-------------------
// deepks_pdm.cpp
//-------------------

// This file contains subroutines for calculating pdm,
// which is defind as sum_mu,nu rho_mu,nu (<chi_mu|alpha><alpha|chi_nu>);
// It also contains subroutines for printing pdm for checking purpose

// There are 3 subroutines in this file:
// 1. read_pdm, which reads pdm from file
// 2. cal_pdm, which is used for calculating pdm
// 3. check_pdm, which prints pdm to descriptor.dat

// read pdm from file, do it only once in whole calculation
void read_pdm(bool read_pdm_file,
              bool is_equiv,
              bool& init_pdm,
              const int inlmax,
              const int lmaxd,
              const int* inl_l,
              const Numerical_Orbital& alpha,
              std::vector<torch::Tensor>& pdm);

// calculate projected density matrix: pdm = sum_i,occ <phi_i|alpha1><alpha2|phi_k>
// 3 cases to skip calculation of pdm:
//   - NSCF calculation of DeePKS, init_chg = file and pdm has been read
//   - SCF calculation of DeePKS with init_chg = file and pdm has been read for restarting SCF
//   - Relax/Cell-Relax/MD calculation, non-first step will use the convergence pdm from the last step as initial pdm
template <typename TK>
void cal_pdm(bool& init_pdm,
             const int inlmax,
             const int lmaxd,
             const int* inl_l,
             const ModuleBase::IntArray* inl_index,
             const elecstate::DensityMatrix<TK, double>* dm,
             const std::vector<hamilt::HContainer<double>*> phialpha,
             const UnitCell& ucell,
             const LCAO_Orbitals& orb,
             const Grid_Driver& GridD,
             const Parallel_Orbitals& pv,
             std::vector<torch::Tensor>& pdm);

void check_pdm(const int inlmax, const int* inl_l, const std::vector<torch::Tensor>& pdm);
} // namespace DeePKS_domain

#endif
#endif
