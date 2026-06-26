#include "phi_operator_kernel.cuh"
#include "gint_helper.cuh"
#include "sph.cuh"
#include "source_base/module_device/device.h"
#include "source_base/module_device/kernel_compat.h"

namespace ModuleGint
{

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
        double dist = norm3d(coord.x, coord.y, coord.z);
        if (dist < rcut[atom_type])
        {
            if (dist < 1.0E-9)
            { dist += 1.0E-9; }
            // since nwl is less or equal than 5, the size of ylma is (5+1)^2
            double ylma[36];
            const int nwl = ucell_atom_nwl[atom_type];
            sph_harm(nwl, coord.x/dist, coord.y/dist, coord.z/dist, ylma);

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

// Explicit instantiations for set_phi_kernel
template __global__ void set_phi_kernel<double>(
    const int, const int, const int, const double,
    const int*, const bool*, const int*, const int*, const int*,
    const double*, const double*, const double*, const double3*,
    const int*, const double3*, const int2*, const int*, const int*,
    double*);
template __global__ void set_phi_kernel<float>(
    const int, const int, const int, const double,
    const int*, const bool*, const int*, const int*, const int*,
    const double*, const double*, const double*, const double3*,
    const int*, const double3*, const int2*, const int*, const int*,
    float*);

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
                const double rl = pow_int(dist, iw_l);
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

// Explicit instantiations for set_phi_dphi_kernel
template __global__ void set_phi_dphi_kernel<true>(
    const int, const int, const int, const double,
    const int*, const bool*, const int*, const int*, const int*, const int*,
    const double*, const double*, const double*, const double3*,
    const int*, const double3*, const int2*, const int*, const int*,
    double*, double*, double*, double*);
template __global__ void set_phi_dphi_kernel<false>(
    const int, const int, const int, const double,
    const int*, const bool*, const int*, const int*, const int*, const int*,
    const double*, const double*, const double*, const double3*,
    const int*, const double3*, const int2*, const int*, const int*,
    double*, double*, double*, double*);

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

// Explicit instantiations for phi_mul_vldr3_kernel
template __global__ void phi_mul_vldr3_kernel<double>(
    const double*, const double, const double*, const int,
    const int*, const int*, const int*, double*);
template __global__ void phi_mul_vldr3_kernel<float>(
    const float*, const float, const float*, const int,
    const int*, const int*, const int*, float*);

// rho(ir) = \sum_{iwt} \phi_i(ir,iwt) * \phi_j^*(ir,iwt)
// each block calculate the dot product of phi_i and phi_j of a meshgrid.
// The per-thread/warp/block reduction is in double regardless of input types
// so that fp32 inputs are summed without catastrophic precision loss.
template<typename Tin_a, typename Tin_b>
__global__ void phi_dot_phi_kernel(
    const Tin_a* __restrict__ phi_i,
    const Tin_b* __restrict__ phi_j,
    const int mgrids_per_bgrid,
    const int* __restrict__ mgrid_lidx,
    const int* __restrict__ bgrid_phi_len,
    const int* __restrict__ bgrid_phi_start,
    double* __restrict__ rho)
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

// Explicit instantiations for phi_dot_phi_kernel
template __global__ void phi_dot_phi_kernel<double, double>(
    const double*, const double*, const int,
    const int*, const int*, const int*, double*);
template __global__ void phi_dot_phi_kernel<float, double>(
    const float*, const double*, const int,
    const int*, const int*, const int*, double*);

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
