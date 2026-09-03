#pragma once

#include <cuda_runtime.h>

#include "gint_helper.cuh"
#include "sph.cuh"
#include "source_base/module_device/kernel_compat.h"

// The template kernels below are defined in this header (not in the .cu) on
// purpose: in whole-program compilation mode (-rdc=false, the default), nvcc
// gives the host-side stubs of __global__ function templates internal linkage
// (forced since CUDA 13 via -static-global-template-stub=true), so every
// <<<...>>> launch of a template kernel must see its definition in the same
// translation unit. Non-template kernels keep external linkage and stay in
// phi_operator_kernel.cu.

namespace ModuleGint
{

// Templated version: internal computation in double, output cast to Real
template<typename Real>
__global__ void set_phi_kernel(
    const int nwmax,
    const int mgrids_num,
    const int nrmax,
    const double dr_uniform,
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
    const double3* __restrict__ atom_rcoords,
    const int2* __restrict__ atoms_num_info,
    const int* __restrict__ atom_phi_start,
    const int* __restrict__ bgrid_phi_len,
    Real* __restrict__ phi)
{
    const int bgrid_id = blockIdx.y;
    const int mgrid_id = blockIdx.x;
    const int atoms_num = atoms_num_info[bgrid_id].x;
    const int pre_atoms_num = atoms_num_info[bgrid_id].y;
    const double3 mgrid_pos = mgrids_pos[mgrid_id];

    for (int atom_id = threadIdx.x; atom_id < atoms_num; atom_id += blockDim.x)
    {
        const int atom_type = iat2it[atoms_iat[atom_id + pre_atoms_num]];
        const double3 rcoord = atom_rcoords[atom_id + pre_atoms_num];       // rcoord is the ralative coordinate of an atom and a biggrid
        const double3 coord = make_double3(mgrid_pos.x-rcoord.x,                    // coord is the relative coordinate of an atom and a meshgrid
                                           mgrid_pos.y-rcoord.y,
                                           mgrid_pos.z-rcoord.z);
        // Preserve the existing near-origin behavior. Only the exact
        // atomic grid point follows the CPU direct-recurrence semantics.
        const bool exact_origin
            = (coord.x == 0.0 && coord.y == 0.0 && coord.z == 0.0);

        double dist = norm3d(coord.x, coord.y, coord.z);
        if (dist < rcut[atom_type])
        {
            if (dist < 1.0E-9)
            { dist += 1.0E-9; }
            // since nwl is less or equal than 5, the size of ylma is (5+1)^2
            double ylma[36];
            const int nwl = ucell_atom_nwl[atom_type];
            if (exact_origin)
            {
                ModuleBase::sph_harm_direct(
                    nwl, 0.0, 0.0, 0.0, ylma);
            }
            else
            {
                sph_harm(
                    nwl,
                    coord.x/dist,
                    coord.y/dist,
                    coord.z/dist,
                    ylma);
            }

            const double pos = dist / dr_uniform;
            const int ip = static_cast<int>(pos);
            const double dx = pos - ip;
            const double dx2 = dx * dx;
            const double dx3 = dx2 * dx;

            const double c3 = 3.0 * dx2 - 2.0 * dx3;
            const double c1 = 1.0 - c3;
            const double c2 = (dx - 2.0 * dx2 + dx3) * dr_uniform;
            const double c4 = (dx3 - dx2) * dr_uniform;

            double psi = 0;
            const int it_nw = atom_type * nwmax;
            int iw_nr = it_nw * nrmax + ip;
            int phi_idx = atom_phi_start[atom_id + pre_atoms_num] +
                          bgrid_phi_len[bgrid_id] * mgrid_id;

            for (int iw = 0; iw < atom_nw[atom_type]; iw++, iw_nr += nrmax)
            {
                if (atom_iw2_new[it_nw + iw])
                {
                    psi = c1 * psi_u[iw_nr] + c2 * dpsi_u[iw_nr]
                          + c3 * psi_u[iw_nr + 1] + c4 * dpsi_u[iw_nr + 1];
                }
                phi[phi_idx + iw] = static_cast<Real>(psi * ylma[atom_iw2_ylm[it_nw + iw]]);
            }
        }
        else
        {
            int phi_idx = atom_phi_start[atom_id + pre_atoms_num] +
                          bgrid_phi_len[bgrid_id] * mgrid_id;
            for (int iw = 0; iw < atom_nw[atom_type]; iw++)
            {
                phi[phi_idx + iw] = Real(0.0);
            }
        }
    }
}

// WantPhi == false: skip phi[] writes entirely (callers like gint_tau pass nullptr).
template<bool WantPhi>
__global__ void set_phi_dphi_kernel(
    const int nwmax,
    const int mgrids_num,
    const int nrmax,
    const double dr_uniform,
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
    const double3* __restrict__ atom_rcoords,
    const int2* __restrict__ atoms_num_info,
    const int* __restrict__ atom_phi_start,
    const int* __restrict__ bgrid_phi_len,
    double* __restrict__ phi,
    double* __restrict__ dphi_x,
    double* __restrict__ dphi_y,
    double* __restrict__ dphi_z)
{
    const int bgrid_id = blockIdx.y;
    const int mgrid_id = blockIdx.x;
    const int atoms_num = atoms_num_info[bgrid_id].x;
    const int pre_atoms_num = atoms_num_info[bgrid_id].y;
    const double3 mgrid_pos = mgrids_pos[mgrid_id];

    for (int atom_id = threadIdx.x; atom_id < atoms_num; atom_id += blockDim.x)
    {
        const int atom_type = iat2it[atoms_iat[atom_id + pre_atoms_num]];
        const double3 rcoord = atom_rcoords[atom_id + pre_atoms_num];
        const double3 coord = make_double3(mgrid_pos.x-rcoord.x,
                                           mgrid_pos.y-rcoord.y,
                                           mgrid_pos.z-rcoord.z);
        double dist = norm3d(coord.x, coord.y, coord.z);
        if (dist < rcut[atom_type])
        {
            if (dist < 1.0E-9)
            { dist += 1.0E-9; }
            // since nwl is less or equal than 5, the size of rly is (5+1)^2
            // size of grly = 36 * 3
            double rly[36];
            double grly[36 * 3];
            const int nwl = ucell_atom_nwl[atom_type];
            grad_rl_sph_harm(nwl, coord.x, coord.y, coord.z, rly, grly);

            // interpolation
            const double inv_dist = 1.0 / dist;  // hoisted: re-used by every iw below
            const double pos = dist / dr_uniform;
            const int ip = static_cast<int>(pos);
            const double x0 = pos - ip;
            const double x1 = 1.0 - x0;
            const double x2 = 2.0 - x0;
            const double x3 = 3.0 - x0;
            const double x12 = x1 * x2 / 6;
            const double x03 = x0 * x3 / 2;
            double tmp = 0;
            double dtmp = 0;
            const int it_nw = atom_type * nwmax;
            int iw_nr = it_nw * nrmax + ip;
            int phi_idx = atom_phi_start[atom_id + pre_atoms_num] +
                          bgrid_phi_len[bgrid_id] * mgrid_id;
            for (int iw = 0; iw < atom_nw[atom_type]; iw++, iw_nr += nrmax)
            {
                if (atom_iw2_new[it_nw + iw])
                {
                    tmp = x12 * (psi_u[iw_nr] * x3 + psi_u[iw_nr + 3] * x0)
                        + x03 * (psi_u[iw_nr + 1] * x2 - psi_u[iw_nr + 2] * x1);
                    dtmp = x12 * (dpsi_u[iw_nr] * x3 + dpsi_u[iw_nr + 3] * x0)
                         + x03 * (dpsi_u[iw_nr + 1] * x2 - dpsi_u[iw_nr + 2] * x1);
                }
                const int iw_l = atom_iw2_l[it_nw + iw];
                const int idx_ylm = atom_iw2_ylm [it_nw + iw];
                const double rl = ::pow_int(dist, iw_l);
                const double inv_rl = 1.0 / rl;
                const double tmprl = tmp * inv_rl;

                if (WantPhi)
                {
                    phi[phi_idx + iw] = tmprl * rly[idx_ylm];
                }
                // derivative of wave functions with respect to atom positions.
                // (dtmp - tmp*iw_l/dist) / rl * rly / dist  ==  (dtmp*inv_dist - tmp*iw_l*inv_dist^2) * inv_rl * rly
                const double tmpdphi_rly = (dtmp * inv_dist - tmp * iw_l * inv_dist * inv_dist)
                                           * inv_rl * rly[idx_ylm];

                dphi_x[phi_idx + iw] =  tmpdphi_rly * coord.x + tmprl * grly[idx_ylm * 3 + 0];
                dphi_y[phi_idx + iw] =  tmpdphi_rly * coord.y + tmprl * grly[idx_ylm * 3 + 1];
                dphi_z[phi_idx + iw] =  tmpdphi_rly * coord.z + tmprl * grly[idx_ylm * 3 + 2];
            }
        }
        else
        {
            int phi_idx = atom_phi_start[atom_id + pre_atoms_num] +
                          bgrid_phi_len[bgrid_id] * mgrid_id;
            for (int iw = 0; iw < atom_nw[atom_type]; iw++)
            {
                if (WantPhi)
                {
                    phi[phi_idx + iw] = 0.0;
                }
                dphi_x[phi_idx + iw] = 0.0;
                dphi_y[phi_idx + iw] = 0.0;
                dphi_z[phi_idx + iw] = 0.0;
            }
        }
    }
}

__global__ void set_ddphi_kernel(
    const int nwmax,
    const int mgrids_num,
    const int nrmax,
    const double dr_uniform,
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
    const double3* __restrict__ atom_rcoords,
    const int2* __restrict__ atoms_num_info,
    const int* __restrict__ atom_phi_start,
    const int* __restrict__ bgrid_phi_len,
    double* __restrict__ ddphi_xx,
    double* __restrict__ ddphi_xy,
    double* __restrict__ ddphi_xz,
    double* __restrict__ ddphi_yy,
    double* __restrict__ ddphi_yz,
    double* __restrict__ ddphi_zz);

template<typename Real>
__global__ void phi_mul_vldr3_kernel(
    const Real* __restrict__ vl,
    const Real dr3,
    const Real* __restrict__ phi,
    const int mgrids_per_bgrid,
    const int* __restrict__ mgrid_lidx,
    const int* __restrict__ bgrid_phi_len,
    const int* __restrict__ bgrid_phi_start,
    Real* __restrict__ result)
{
    const int bgrid_id = blockIdx.y;
    const int mgrid_id = blockIdx.x;
    const int phi_len = bgrid_phi_len[bgrid_id];
    const int phi_start = bgrid_phi_start[bgrid_id] + mgrid_id * phi_len;
    const int batch_mgrid_id = bgrid_id * mgrids_per_bgrid + mgrid_id;
    const Real vldr3 =  vl[mgrid_lidx[batch_mgrid_id]] * dr3;
    for(int i = threadIdx.x; i < phi_len; i += blockDim.x)
    {
        result[phi_start + i] = phi[phi_start + i] * vldr3;
    }
}

// rho(ir) = \sum_{iwt} \phi_i(ir,iwt) * \phi_j^*(ir,iwt)
// each block calculate the dot product of phi_i and phi_j of a meshgrid.
// Inputs phi_i and phi_j can have different element types: in the rho path
// phi_i is fp32 (Real) while phi_j (phi_dm) is fp64; in the tau path both are
// fp64. The per-block reduction and atomicAdd to rho run in fp64 regardless.
template<typename Tin_a, typename Tin_b>
__global__ void phi_dot_phi_kernel(
    const Tin_a* __restrict__ phi_i,          // phi_i(ir,iwt)
    const Tin_b* __restrict__ phi_j,          // phi_j(ir,iwt)
    const int mgrids_per_bgrid,                 // the number of mgrids of each biggrid
    const int* __restrict__ mgrid_lidx,   // the idx of mgrid in local cell
    const int* __restrict__ bgrid_phi_len,     // the length of phi on a mgrid of a biggrid
    const int* __restrict__ bgrid_phi_start,   // the start idx in phi of each biggrid
    double* __restrict__ rho)                // rho(ir)
{
    __shared__ double s_data[32];  // the length of s_data equals the max warp num of a block
    const int bgrid_id = blockIdx.y;
    const int mgrid_id = blockIdx.x;
    const int phi_len = bgrid_phi_len[bgrid_id];
    const int phi_start = bgrid_phi_start[bgrid_id] + mgrid_id * phi_len;
    const Tin_a* phi_i_mgrid = phi_i + phi_start;
    const Tin_b* phi_j_mgrid = phi_j + phi_start;
    const int batch_mgrid_id = bgrid_id * mgrids_per_bgrid + mgrid_id;
    const int mgrid_local_idx = mgrid_lidx[batch_mgrid_id];
    const int tid = threadIdx.x;
    const int warp_id = tid / 32;
    const int lane_id = tid % 32;
    double tmp_sum = 0.0;

    for (int i = tid; i < phi_len; i += blockDim.x)
    {
        tmp_sum += phi_i_mgrid[i] * phi_j_mgrid[i];
    }

    tmp_sum = warpReduceSum(tmp_sum);

    if (lane_id == 0)
    {
        s_data[warp_id] = tmp_sum;
    }
    __syncthreads();

    tmp_sum = (tid < blockDim.x / 32) ? s_data[tid] : 0.0;
    if(warp_id == 0)
    {
        tmp_sum = warpReduceSum(tmp_sum);
    }

    if(tid == 0)
    {
        atomicAdd(&rho[mgrid_local_idx], tmp_sum);
    }
}

__global__ void phi_dot_dphi_kernel(
    const double* __restrict__ phi,
    const double* __restrict__ dphi_x,
    const double* __restrict__ dphi_y,
    const double* __restrict__ dphi_z,
    const int mgrids_per_bgrid,
    const int* __restrict__ bgrid_phi_len,
    const int2* __restrict__ atoms_num_info,
    const int* __restrict__ atom_phi_start,
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
    const int* __restrict__ bgrid_phi_len,
    const int2* __restrict__ atoms_num_info,
    const int* __restrict__ atom_phi_start,
    const int* __restrict__ atoms_iat,
    const double3* __restrict__ atom_rcoords,
    const double3* __restrict__ mgrids_pos,
    const int* __restrict__ iat2it,
    const int* __restrict__ atom_nw,
    double* __restrict__ svl);

}
