#ifndef DEEPKS_VDRPRE_H
#define DEEPKS_VDRPRE_H

#ifdef __MLALGO

#include "deepks_param.h"
#include "source_base/complexmatrix.h"
#include "source_base/intarray.h"
#include "source_base/matrix.h"
#include "source_base/timer.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_lcao/module_hcontainer/hcontainer.h"

#include <torch/script.h>
#include <torch/torch.h>

namespace DeePKS_domain
{
//------------------------
// deepks_vdrpre.cpp
//------------------------

// This file contains 1 subroutine for calculating v_delta,
// cal_vdr_precalc : v_delta_r_precalc is used for training with v_delta_r label,
//                         which equals gevdm * v_delta_pdm,
//                         v_delta_pdm = overlap * overlap

// for deepks_v_delta = -1
// calculates v_delta_r_precalc

void prepare_phialpha_iRmat(const int nlocal,
                            const int R_size,
                            const DeePKS_Param& deepks_param,
                            const std::vector<hamilt::HContainer<double>*> phialpha,
                            const UnitCell& ucell,
                            const LCAO_Orbitals& orb,
                            const Grid_Driver& GridD,
                            torch::Tensor& overlap,
                            torch::Tensor& iRmat);

void cal_vdr_precalc(const int nlocal,
                     const int nat,
                     const int nks,
                     const int R_size,
                     const DeePKS_Param& deepks_param,
                     const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                     const std::vector<hamilt::HContainer<double>*> phialpha,
                     const std::vector<torch::Tensor> gevdm,
                     const UnitCell& ucell,
                     const LCAO_Orbitals& orb,
                     const Parallel_Orbitals& pv,
                     const Grid_Driver& GridD,
                     torch::Tensor& vdr_precalc);

int mapping_R(int R);
torch::Tensor mapping_R(const torch::Tensor& R_tensor);

template <typename T>
int get_R_size(const hamilt::HContainer<T>& hcontainer);
} // namespace DeePKS_domain
#endif
#endif