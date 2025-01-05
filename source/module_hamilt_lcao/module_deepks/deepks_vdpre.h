#ifndef DEEPKS_VDPRE_H
#define DEEPKS_VDPRE_H

#ifdef __DEEPKS

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
// deepks_vdpre.cpp
//------------------------

// This file contains 6 subroutines for calculating v_delta,
// 1. cal_v_delta_precalc : v_delta_precalc is used for training with v_delta label,
//                         which equals gvdm * v_delta_pdm,
//                         v_delta_pdm = overlap * overlap
// 2. check_v_delta_precalc : check v_delta_precalc
// 3. prepare_phialpha : prepare phialpha for outputting npy file
// 4. check_vdp_phialpha : check phialpha
// 5. prepare_gevdm : prepare gevdm for outputting npy file
// 6. check_vdp_gevdm : check gevdm

// for deepks_v_delta = 1
// calculates v_delta_precalc
template <typename TK>
void cal_v_delta_precalc(const int nlocal,
                         const int lmaxd,
                         const int inlmax,
                         const int nat,
                         const int nks,
                         const int* inl_l,
                         const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                         const std::vector<hamilt::HContainer<double>*> phialpha,
                         const std::vector<torch::Tensor> gevdm,
                         const ModuleBase::IntArray* inl_index,
                         const UnitCell& ucell,
                         const LCAO_Orbitals& orb,
                         const Parallel_Orbitals& pv,
                         const Grid_Driver& GridD,
                         torch::Tensor& v_delta_precalc);

template <typename TK>
void check_v_delta_precalc(const int nat,
                           const int nks,
                           const int nlocal,
                           const int des_per_atom,
                           const torch::Tensor& v_delta_precalc);

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

template <typename TK>
void check_vdp_phialpha(const int nat,
                        const int nks,
                        const int nlocal,
                        const int inlmax,
                        const int lmaxd,
                        const torch::Tensor& phialpha_out);

// prepare gevdm for outputting npy file
void prepare_gevdm(const int nat,
                   const int lmaxd,
                   const int inlmax,
                   const LCAO_Orbitals& orb,
                   const std::vector<torch::Tensor>& gevdm_in,
                   torch::Tensor& gevdm_out);

void check_vdp_gevdm(const int nat, const int lmaxd, const int inlmax, const torch::Tensor& gevdm);
} // namespace DeePKS_domain
#endif
#endif