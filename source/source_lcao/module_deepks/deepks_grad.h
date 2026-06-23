#ifndef DEEPKS_GRAD_H
#define DEEPKS_GRAD_H

#ifdef __MLALGO

#include "deepks_param.h"
#include "source_base/timer.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_lcao/module_hcontainer/hcontainer.h"

#include <torch/torch.h>

namespace DeePKS_domain
{
//------------------------
// deepks_grad.cpp
//------------------------
// Building block for the descriptor-gradient label g* = (B^T B)^{-1} B^T H.
// B: vdr_precalc, g: gradient (\partial E_delta / \partial descriptor), H: Hamiltonian(HR)
//
// B_{(mu,nu,R),(I,nl,k)} = sum_{m,m'} v^I_{nl,k,m} * Theta^I_{nl,mu,nu}(R)_{mm'} * v^I_{nl,k,m'}
//   where Theta^I_{nl,mu,nu}(R)_{mm'} = sum_{dR1-dR2=R} phi[mu,m,dR1] * phi[nu,m',dR2]
//
// cal_phialpha_hamilt_proj  (Step 1)
//   dot_phialpha_hamilt[inl, m, m']
//     = <phialpha^I_{nl,m} | H | phialpha^I_{nl,m'}>
//     = sum_{mu,nu,R} phi[mu,m,dR1] * H[mu,nu](dR2-dR1) * phi[nu,m',dR2]
//   Shape: (inlmax, nm_max, nm_max).
//
// cal_vdrpre_square: full B^T B matrix for the descriptor-gradient label.
//     (B^T B)_{a1,a2} = sum_{R,mu,nu} B_{(mu,nu,R),a1} * B_{(mu,nu,R),a2}
//                     = flat(vdr_precalc)^T @ flat(vdr_precalc)
//   Shape: (nat * des_per_atom, nat * des_per_atom).
//
// Requires gevdm (deepks_scf must be active).
void cal_vdrpre_square(const int nlocal,
                       const int nat,
                       const int R_size,
                       const DeePKS_Param& deepks_param,
                       const std::vector<hamilt::HContainer<double>*>& phialpha,
                       const std::vector<torch::Tensor>& gevdm,
                       const UnitCell& ucell,
                       const LCAO_Orbitals& orb,
                       const Parallel_Orbitals& pv,
                       const Grid_Driver& GridD,
                       torch::Tensor& vdrpre_square);

template <typename TR>
void cal_phialpha_hamilt_proj(const int nlocal,
                               const int nat,
                               const DeePKS_Param& deepks_param,
                               const std::vector<hamilt::HContainer<double>*>& phialpha,
                               const hamilt::HContainer<TR>& hR,
                               const UnitCell& ucell,
                               const LCAO_Orbitals& orb,
                               const Parallel_Orbitals& pv,
                               const Grid_Driver& GridD,
                               torch::Tensor& dot_phialpha_hamilt);


} // namespace DeePKS_domain

#endif
#endif
