#ifndef DEEPKS_PHIALPHA_H
#define DEEPKS_PHIALPHA_H

#ifdef __DEEPKS

#include "module_base/complexmatrix.h"
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
// This file contains 3 subroutines:
// 1. allocate_phialpha, which allocates memory for phialpha
// 2. build_phialpha, which calculates the overlap
// between atomic basis and projector alpha : <phi_mu|alpha>
// which will be used in calculating pdm, gdmx, H_V_delta, F_delta;
// 3. check_phialpha, which prints the results into .dat files
// for checking

// calculates <chi|alpha>
void allocate_phialpha(const bool& cal_deri,
                       const UnitCell& ucell,
                       const LCAO_Orbitals& orb,
                       const Grid_Driver& GridD,
                       const Parallel_Orbitals* pv,
                       std::vector<hamilt::HContainer<double>*>& phialpha);

void build_phialpha(const bool& cal_deri /**< [in] 0 for 2-center intergration, 1 for its derivation*/,
                    const UnitCell& ucell,
                    const LCAO_Orbitals& orb,
                    const Grid_Driver& GridD,
                    const Parallel_Orbitals* pv,
                    const TwoCenterIntegrator& overlap_orb_alpha,
                    std::vector<hamilt::HContainer<double>*>& phialpha);

void check_phialpha(const bool& cal_deri /**< [in] 0 for 2-center intergration, 1 for its derivation*/,
                    const UnitCell& ucell,
                    const LCAO_Orbitals& orb,
                    const Grid_Driver& GridD,
                    const Parallel_Orbitals* pv,
                    std::vector<hamilt::HContainer<double>*>& phialpha,
                    const int rank);
} // namespace DeePKS_domain

#endif
#endif
