#pragma once

#include <cuda_runtime.h>

namespace ModuleGint
{

__global__ void set_phi_kernel(
    const int nwmax,
    const int mgrids_num,
    const int nrmax,
    const double dr_uniform,
    const double* __restrict__ ylmcoef,
    const int* __restrict__ ucell_atom_nwl,
    const bool* __restrict__ atom_iw2_new,
    const int* __restrict__ atom_iw2_ylm,
    const int* __restrict__ atom_nw,
    const int* __restrict__ iat2it,
    const double* __restrict__ rcut,
    const double* __restrict__ psi_u,
    const double* __restrict__ dpsi_u,
    const double3* __restrict__ mgrids_pos,
    const int* __restrict__ atoms_iat,
    const double3* __restrict__ atoms_bgrids_rcoords,
    const int2* __restrict__ atoms_num_info,
    const int* __restrict__ atoms_phi_start,
    const int* __restrict__ bgrids_phi_len,
    double* __restrict__ phi);

__global__ void set_phi_dphi_kernel(
    const int nwmax,
    const int mgrids_num,
    const int nrmax,
    const double dr_uniform,
    const double* __restrict__ ylmcoef,
    const int* __restrict__ ucell_atom_nwl,
    const bool* __restrict__ atom_iw2_new,
    const int* __restrict__ atom_iw2_ylm,
    const int* __restrict__ atom_iw2_l,
    const int* __restrict__ atom_nw,
    const int* __restrict__ iat2it,
    const double* __restrict__ rcut,
    const double* __restrict__ psi_u,
    const double* __restrict__ dpsi_u,
    const double3* __restrict__ mgrids_pos,
    const int* __restrict__ atoms_iat,
    const double3* __restrict__ atoms_bgrids_rcoords,
    const int2* __restrict__ atoms_num_info,
    const int* __restrict__ atoms_phi_start,
    const int* __restrict__ bgrids_phi_len,
    double* __restrict__ phi,
    double* __restrict__ dphi_x,
    double* __restrict__ dphi_y,
    double* __restrict__ dphi_z);

__global__ void set_ddphi_kernel(
    const int nwmax,
    const int mgrids_num,
    const int nrmax,
    const double dr_uniform,
    const double* __restrict__ ylmcoef,
    const int* __restrict__ ucell_atom_nwl,
    const bool* __restrict__ atom_iw2_new,
    const int* __restrict__ atom_iw2_ylm,
    const int* __restrict__ atom_iw2_l,
    const int* __restrict__ atom_nw,
    const int* __restrict__ iat2it,
    const double* __restrict__ rcut,
    const double* __restrict__ psi_u,
    const double* __restrict__ dpsi_u,
    const double3* __restrict__ mgrids_pos,
    const int* __restrict__ atoms_iat,
    const double3* __restrict__ atoms_bgrids_rcoords,
    const int2* __restrict__ atoms_num_info,
    const int* __restrict__ atoms_phi_start,
    const int* __restrict__ bgrids_phi_len,
    double* __restrict__ ddphi_xx,
    double* __restrict__ ddphi_xy,
    double* __restrict__ ddphi_xz,
    double* __restrict__ ddphi_yy,
    double* __restrict__ ddphi_yz,
    double* __restrict__ ddphi_zz);

__global__ void phi_mul_vldr3_kernel(
    const double* __restrict__ vl,
    const double dr3,
    const double* __restrict__ phi,
    const int mgrids_per_bgrid,
    const int* __restrict__ mgrids_local_idx,
    const int* __restrict__ bgrids_phi_len,
    const int* __restrict__ bgrids_phi_start,
    double* __restrict__ result);

// rho(ir) = \sum_{iwt} \phi_i(ir,iwt) * \phi_j^*(ir,iwt)
// each block calculate the dot product of phi_i and phi_j of a meshgrid
__global__ void phi_dot_phi_kernel(
    const double* __restrict__ phi_i,           // phi_i(ir,iwt)
    const double* __restrict__ phi_j,           // phi_j(ir,iwt)
    const int mgrids_per_bgrid,                 // the number of mgrids of each biggrid
    const int* __restrict__ mgrids_local_idx,   // the idx of mgrid in local cell
    const int* __restrict__ bgrids_phi_len,     // the length of phi on a mgrid of a biggrid
    const int* __restrict__ bgrids_phi_start,   // the start idx in phi of each biggrid
    double* __restrict__ rho);                  // rho(ir)

__global__ void phi_dot_dphi_kernel(
    const double* __restrict__ phi,
    const double* __restrict__ dphi_x,
    const double* __restrict__ dphi_y,
    const double* __restrict__ dphi_z,
    const int mgrids_per_bgrid,
    const int* __restrict__ bgrids_phi_len,
    const int2* __restrict__ atoms_num_info,
    const int* __restrict__ atoms_phi_start,
    const int* __restrict__ atoms_iat,
    const int* __restrict__ iat2it,
    const int* __restrict__ atom_nw,
    double* force);

__global__ void phi_dot_dphi_r_kernel(
    const double* __restrict__ phi,
    const double* __restrict__ dphi_x,
    const double* __restrict__ dphi_y,
    const double* __restrict__ dphi_z,
    const int mgrids_per_bgrid,
    const int* __restrict__ bgrids_phi_len,
    const int2* __restrict__ atoms_num_info,
    const int* __restrict__ atoms_phi_start,
    const int* __restrict__ atoms_iat,
    const double3* __restrict__ atoms_bgrids_rcoords,
    const double3* __restrict__ mgrids_pos,
    const int* __restrict__ iat2it,
    const int* __restrict__ atom_nw,
    double* __restrict__ svl);
    
}
