#ifndef DEEPKS_FPRE_H
#define DEEPKS_FPRE_H

#ifdef __MLALGO

#include "module_base/complexmatrix.h"
#include "module_base/intarray.h"
#include "module_base/matrix.h"
#include "module_base/timer.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"

#include <torch/script.h>
#include <torch/torch.h>

namespace DeePKS_domain
{
//------------------------
// deepks_fpre.cpp
//------------------------

// This file contains 2 subroutines for calculating,
// 1. cal_gdmx, calculating gdmx
// 2. cal_gvx : gvx is used for training with force label, which is gradient of descriptors,
//       calculated by d(des)/dX = d(pdm)/dX * d(des)/d(pdm) = gdmx * gvdm
//       using einsum

// calculate the gradient of pdm with regard to atomic positions
// d/dX D_{Inl,mm'}
template <typename TK>
void cal_gdmx(const int lmaxd,
              const int inlmax,
              const int nks,
              const std::vector<ModuleBase::Vector3<double>>& kvec_d,
              std::vector<hamilt::HContainer<double>*> phialpha,
              const ModuleBase::IntArray* inl_index,
              const hamilt::HContainer<double>* dmr,
              const UnitCell& ucell,
              const LCAO_Orbitals& orb,
              const Parallel_Orbitals& pv,
              const Grid_Driver& GridD,
              torch::Tensor& gdmx);

/// calculates gradient of descriptors w.r.t atomic positions
///----------------------------------------------------
/// m, n: 2*l+1
/// v: eigenvalues of dm , 2*l+1
/// a,b: natom
///  - (a: the center of descriptor orbitals
///  - b: the atoms whose force being calculated)
/// gevdm*gdmx->gvx
///----------------------------------------------------
void cal_gvx(const int nat,
             const int inlmax,
             const int des_per_atom,
             const std::vector<int>& inl2l,
             const std::vector<torch::Tensor>& gevdm,
             const torch::Tensor& gdmx,
             torch::Tensor& gvx,
             const int rank);

} // namespace DeePKS_domain
#endif
#endif