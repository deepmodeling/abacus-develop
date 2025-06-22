#ifndef DEEPKS_ORBPRE_H
#define DEEPKS_ORBPRE_H

#ifdef __MLALGO

#include "source_base/complexmatrix.h"
#include "source_base/intarray.h"
#include "source_base/matrix.h"
#include "source_base/timer.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"

#include <torch/script.h>
#include <torch/torch.h>

namespace DeePKS_domain
{
//------------------------
// deepks_orbpre.cpp
//------------------------

// This file contains one subroutine for calculating orbital_precalc,
// which is defind as gevdm * dm_hl * overlap * overlap

template <typename TK, typename TH>
void cal_orbital_precalc(const std::vector<TH>& dm_hl,
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
                         torch::Tensor& orbital_precalc);
} // namespace DeePKS_domain
#endif
#endif
