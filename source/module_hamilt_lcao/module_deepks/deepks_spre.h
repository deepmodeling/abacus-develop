#ifndef DEEPKS_SPRE_H
#define DEEPKS_SPRE_H

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
// deepks_spre.cpp
//------------------------

// This file contains 4 subroutines for calculating,
// 1. cal_gdmepsl, calculating gdmepsl
// 2. check_gdmepsl, which prints gdmepsl to a series of .dat files
// 3. cal_gvepsl : gvepsl is used for training with stress label, which is derivative of
//       descriptors wrt strain tensor, calculated by
//       d(des)/d\epsilon_{ab} = d(pdm)/d\epsilon_{ab} * d(des)/d(pdm) = gdmepsl * gvdm
//       using einsum
// 4. check_gvepsl : prints gvepsl into gvepsl.dat for checking

// calculate the gradient of pdm with regard to atomic virial stress tensor
// d/d\epsilon D_{Inl,mm'}
template <typename TK>
void cal_gdmepsl( // const ModuleBase::matrix& dm,
    const int lmaxd,
    const int inlmax,
    const int nks,
    const std::vector<ModuleBase::Vector3<double>>& kvec_d,
    std::vector<hamilt::HContainer<double>*> phialpha,
    const ModuleBase::IntArray* inl_index,
    const std::vector<std::vector<TK>>& dm,
    const UnitCell& ucell,
    const LCAO_Orbitals& orb,
    const Parallel_Orbitals& pv,
    const Grid_Driver& GridD,
    torch::Tensor& gdmepsl);

void check_gdmepsl(const torch::Tensor& gdmepsl);

void cal_gvepsl(const int nat,
                const int inlmax,
                const int des_per_atom,
                const int* inl_l,
                const std::vector<torch::Tensor>& gevdm,
                const torch::Tensor& gdmepsl,
                torch::Tensor& gvepsl);

void check_gvepsl(const torch::Tensor& gvepsl);
} // namespace DeePKS_domain
#endif
#endif