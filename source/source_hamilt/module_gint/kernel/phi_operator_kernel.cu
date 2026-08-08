#include "phi_operator_kernel.cuh"
#include "gint_helper.cuh"
#include "sph.cuh"
#include "source_base/module_device/device.h"
#include "source_base/module_device/kernel_compat.h"

// The template kernels (set_phi_kernel, set_phi_dphi_kernel,
// phi_mul_vldr3_kernel, phi_dot_phi_kernel) are defined in
// phi_operator_kernel.cuh so that their <<<...>>> launches in
// phi_operator_gpu.cu see the definitions in the same translation unit; only
// the non-template kernels live here (see the note in the header).

namespace ModuleGint
{

// The code for `set_ddphi_kernel` is quite difficult to understand.
// To grasp it, you better refer to the CPU function `set_ddphi`
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
    double* __restrict__ ddphi_zz)
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
        double coord[3]{mgrid_pos.x-rcoord.x,
                        mgrid_pos.y-rcoord.y,
                        mgrid_pos.z-rcoord.z};
        double dist = norm3d(coord[0], coord[1], coord[2]);
        if (dist < rcut[atom_type])
        {
            int phi_idx = atom_phi_start[atom_id + pre_atoms_num] +
                          bgrid_phi_len[bgrid_id] * mgrid_id;
            for(int i = 0; i < 6; i++)
            {
                const double eps = (i & 1) ? -0.0001 : 0.0001;
                coord[i/2] += eps;
                double dist = norm3d(coord[0], coord[1], coord[2]);
                if (dist < 1.0E-9)
                { dist += 1.0E-9; }
                // since nwl is less or equal than 5, the size of rly is (5+1)^2
                // size of grly = 36 * 3
                double rly[36];
                double grly[36 * 3];
                const int nwl = ucell_atom_nwl[atom_type];
                grad_rl_sph_harm(nwl, coord[0], coord[1], coord[2], rly, grly);

                // interpolation
                const double inv_dist = 1.0 / dist;  // hoisted: re-used by every iw
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
                    const double rl = pow_int(dist, iw_l);
                    const double inv_rl = 1.0 / rl;
                    const double tmprl = tmp * inv_rl;
                    const double tmpdphi_rly = (dtmp * inv_dist - tmp * iw_l * inv_dist * inv_dist)
                                               * inv_rl * rly[idx_ylm];

                    double dphi[3];
                    dphi[0] = tmpdphi_rly * coord[0] + tmprl * grly[idx_ylm * 3 + 0];
                    dphi[1] = tmpdphi_rly * coord[1] + tmprl * grly[idx_ylm * 3 + 1];
                    dphi[2] = tmpdphi_rly * coord[2] + tmprl * grly[idx_ylm * 3 + 2];

                    if (i == 0)
                    {
                        ddphi_xx[phi_idx + iw] += dphi[0];
                        ddphi_xy[phi_idx + iw] += dphi[1];
                        ddphi_xz[phi_idx + iw] += dphi[2];
                    } else if (i == 1)
                    {
                        ddphi_xx[phi_idx + iw] -= dphi[0];
                        ddphi_xy[phi_idx + iw] -= dphi[1];
                        ddphi_xz[phi_idx + iw] -= dphi[2];
                    } else if (i == 2)
                    {
                        ddphi_xy[phi_idx + iw] += dphi[0];
                        ddphi_yy[phi_idx + iw] += dphi[1];
                        ddphi_yz[phi_idx + iw] += dphi[2];
                    } else if (i == 3)
                    {
                        ddphi_xy[phi_idx + iw] -= dphi[0];
                        ddphi_yy[phi_idx + iw] -= dphi[1];
                        ddphi_yz[phi_idx + iw] -= dphi[2];
                    } else if (i == 4)
                    {
                        ddphi_xz[phi_idx + iw] += dphi[0];
                        ddphi_yz[phi_idx + iw] += dphi[1];
                        ddphi_zz[phi_idx + iw] += dphi[2];
                    } else // i == 5
                    {
                        ddphi_xz[phi_idx + iw] -= dphi[0];
                        ddphi_yz[phi_idx + iw] -= dphi[1];
                        ddphi_zz[phi_idx + iw] -= dphi[2];
                    }
                }
                coord[i/2] -= eps;  // recover coord
            }

            for (int iw = 0; iw < atom_nw[atom_type]; iw++)
            {
                ddphi_xx[phi_idx + iw] /= 0.0002;
                ddphi_xy[phi_idx + iw] /= 0.0004;
                ddphi_xz[phi_idx + iw] /= 0.0004;
                ddphi_yy[phi_idx + iw] /= 0.0002;
                ddphi_yz[phi_idx + iw] /= 0.0004;
                ddphi_zz[phi_idx + iw] /= 0.0002;
            }
        }
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
    double* force)
{
    // NOTE: this kernel assumes blockDim.x == 32 (a single warp). If the launch
    // configuration is ever changed, the reduce below needs a shared-memory stage.
    const int bgrid_id = blockIdx.y;
    const int atoms_num = atoms_num_info[bgrid_id].x;
    const int pre_atoms_num = atoms_num_info[bgrid_id].y;
    const int b_phi_len = bgrid_phi_len[bgrid_id];
    const int tid = threadIdx.x;
    const int lane_id = tid;  // blockDim.x == 32

    for (int atom_id = blockIdx.x; atom_id < atoms_num; atom_id += gridDim.x)
    {
        const int a_phi_start = atom_phi_start[atom_id + pre_atoms_num];
        const int iat = atoms_iat[atom_id + pre_atoms_num];
        const int nw = atom_nw[iat2it[iat]];
        double f[3] = {0.0, 0.0, 0.0};
        for (int mgrid_id = 0; mgrid_id < mgrids_per_bgrid; mgrid_id++)
        {
            const int phi_start = a_phi_start + mgrid_id * b_phi_len;
            for (int iw = tid; iw < nw; iw += blockDim.x)
            {
                int phi_idx = phi_start + iw;
                const double p = phi[phi_idx];
                f[0] += p * dphi_x[phi_idx];
                f[1] += p * dphi_y[phi_idx];
                f[2] += p * dphi_z[phi_idx];
            }
        }

        // single-warp reduce
        f[0] = warpReduceSum(f[0]);
        f[1] = warpReduceSum(f[1]);
        f[2] = warpReduceSum(f[2]);

        if (lane_id == 0)
        {
            atomicAdd(&force[iat * 3 + 0], f[0] * 2);
            atomicAdd(&force[iat * 3 + 1], f[1] * 2);
            atomicAdd(&force[iat * 3 + 2], f[2] * 2);
        }
    }
}

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
    double* __restrict__ svl)
{
    // NOTE: this kernel assumes blockDim.x == 32 (a single warp). If the launch
    // configuration is ever changed, the reduce below needs a shared-memory stage.
    const int tid = threadIdx.x;
    const int bgrid_id = blockIdx.y;
    const int atoms_num = atoms_num_info[bgrid_id].x;
    const int pre_atoms_num = atoms_num_info[bgrid_id].y;
    const int b_phi_len = bgrid_phi_len[bgrid_id];
    const int lane_id = tid;  // blockDim.x == 32

    double stress[6]{0.0};
    for (int mgrid_id = blockIdx.x; mgrid_id < mgrids_per_bgrid; mgrid_id += gridDim.x)
    {
        const double3 mgrid_pos = mgrids_pos[mgrid_id];
        for (int atom_id = 0; atom_id < atoms_num; atom_id++)
        {
            const int phi_start = atom_phi_start[atom_id + pre_atoms_num] + mgrid_id * b_phi_len;
            const int iat = atoms_iat[atom_id + pre_atoms_num];
            const int nw = atom_nw[iat2it[iat]];
            const double3 rcoord = atom_rcoords[atom_id + pre_atoms_num];       // rcoord is the ralative coordinate of an atom and a biggrid
            const double3 coord = make_double3(mgrid_pos.x-rcoord.x,                    // coord is the relative coordinate of an atom and a meshgrid
                                               mgrid_pos.y-rcoord.y,
                                               mgrid_pos.z-rcoord.z);
            for (int iw = tid; iw < nw; iw += blockDim.x)
            {
                int phi_idx = phi_start + iw;
                const double p   = phi[phi_idx];
                const double pdx = p * dphi_x[phi_idx];
                const double pdy = p * dphi_y[phi_idx];
                const double pdz = p * dphi_z[phi_idx];
                stress[0] += pdx * coord.x;
                stress[1] += pdx * coord.y;
                stress[2] += pdx * coord.z;
                stress[3] += pdy * coord.y;
                stress[4] += pdy * coord.z;
                stress[5] += pdz * coord.z;
            }
        }
    }

    // single-warp reduce
    #pragma unroll
    for (int i = 0; i < 6; i++)
    {
        stress[i] = warpReduceSum(stress[i]);
    }
    if (lane_id == 0)
    {
        #pragma unroll
        for (int i = 0; i < 6; i++)
        {
            atomicAdd(&svl[i], stress[i] * 2);
        }
    }
}

}
