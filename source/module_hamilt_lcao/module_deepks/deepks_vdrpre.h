#ifndef DEEPKS_VDRPRE_H
#define DEEPKS_VDRPRE_H

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
// deepks_vdrpre.cpp
//------------------------

// This file contains 1 subroutine for calculating v_delta,
// cal_vdr_precalc : v_delta_r_precalc is used for training with v_delta_r label,
//                         which equals gevdm * v_delta_pdm,
//                         v_delta_pdm = overlap * overlap

// for deepks_v_delta = -1
// calculates v_delta_r_precalc
void prepare_phialpha_r(const int nlocal,
                        const int lmaxd,
                        const int inlmax,
                        const int nat,
                        const std::vector<hamilt::HContainer<double>*> phialpha,
                        const UnitCell& ucell,
                        const LCAO_Orbitals& orb,
                        const Parallel_Orbitals& pv,
                        const Grid_Driver& GridD,
                        torch::Tensor& phialpha_r_out,
                        torch::Tensor& R_query);

} // namespace DeePKS_domain
#endif
#endif