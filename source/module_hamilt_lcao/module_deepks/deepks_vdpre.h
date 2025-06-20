#ifndef DEEPKS_VDPRE_H
#define DEEPKS_VDPRE_H

#ifdef __MLALGO

#include "source_base/complexmatrix.h"
#include "source_base/intarray.h"
#include "source_base/matrix.h"
#include "source_base/timer.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"

#include <torch/script.h>
#include <torch/torch.h>

namespace DeePKS_domain
{
//------------------------
// deepks_vdpre.cpp
//------------------------

// This file contains 3 subroutines for calculating v_delta,
// 1. cal_v_delta_precalc : v_delta_precalc is used for training with v_delta label,
//                         which equals gevdm * v_delta_pdm,
//                         v_delta_pdm = overlap * overlap
// 2. prepare_phialpha : prepare phialpha for outputting npy file
// 3. prepare_gevdm : prepare gevdm for outputting npy file

// for deepks_v_delta = 1
// calculates v_delta_precalc
template <typename TK>
void cal_v_delta_precalc(const int nlocal,
                         const int lmaxd,
                         const int inlmax,
                         const int nat,
                         const int nks,
                         const std::vector<int>& inl2l,
                         const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                         const std::vector<hamilt::HContainer<double>*> phialpha,
                         const std::vector<torch::Tensor> gevdm,
                         const ModuleBase::IntArray* inl_index,
                         const UnitCell& ucell,
                         const LCAO_Orbitals& orb,
                         const Parallel_Orbitals& pv,
                         const Grid_Driver& GridD,
                         torch::Tensor& v_delta_precalc);

// for deepks_v_delta = 2
// prepare phialpha for outputting npy file
template <typename TK>
void prepare_phialpha(const int nlocal,
                      const int lmaxd,
                      const int inlmax,
                      const int nat,
                      const int nks,
                      const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                      const std::vector<hamilt::HContainer<double>*> phialpha,
                      const UnitCell& ucell,
                      const LCAO_Orbitals& orb,
                      const Parallel_Orbitals& pv,
                      const Grid_Driver& GridD,
                      torch::Tensor& phialpha_out);

// prepare gevdm for outputting npy file
void prepare_gevdm(const int nat,
                   const int lmaxd,
                   const int inlmax,
                   const LCAO_Orbitals& orb,
                   const std::vector<torch::Tensor>& gevdm_in,
                   torch::Tensor& gevdm_out);
} // namespace DeePKS_domain
#endif
#endif