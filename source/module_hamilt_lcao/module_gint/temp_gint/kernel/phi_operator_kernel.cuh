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

__global__ void phi_mul_vldr3_kernel(
    const double* __restrict__ vl,
    const double dr3,
    const double* __restrict__ phi,
    const int mgrids_per_bgrid,
    const int* __restrict__ mgrids_local_idx,
    const int* __restrict__ bgrids_phi_len,
    const int* __restrict__ bgrids_phi_start,
    double* __restrict__ result);

}