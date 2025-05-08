#include "phi_operator_kernel.cuh"
#include "gint_helper.cuh"
#include "sph.cuh"

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
    double* __restrict__ phi)
{
    const int bgrid_id = blockIdx.y;
    const int mgrid_id = blockIdx.x;
    const int atoms_num = atoms_num_info[bgrid_id].x;
    const int pre_atoms_num = atoms_num_info[bgrid_id].y;
    const double3 mgrid_pos = mgrids_pos[mgrid_id];
    
    for (int atom_id = threadIdx.x; atom_id < atoms_num; atom_id += blockDim.x)
    {
        const int atom_type = iat2it[atoms_iat[atom_id + pre_atoms_num]];
        const double3 rcoord = atoms_bgrids_rcoords[atom_id + pre_atoms_num];
        const double3 coord = make_double3(mgrid_pos.x-rcoord.x,
                                           mgrid_pos.y-rcoord.y,
                                           mgrid_pos.z-rcoord.z);
        double dist = sqrt(coord.x * coord.x + coord.y * coord.y + coord.z * coord.z);
        if (dist < rcut[atom_type])
        {
            if (dist < 1.0E-9)
            { dist += 1.0E-9; }
            // since nwl is less or equal than 5, the size of ylma is (5+1)^2
            double ylma[36];
            const int nwl = ucell_atom_nwl[atom_type];
            sph_harm(nwl, ylmcoef, coord.x/dist, coord.y/dist, coord.z/dist, ylma);

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
            int phi_idx = atoms_phi_start[atom_id + pre_atoms_num] +
                          bgrids_phi_len[bgrid_id] * mgrid_id;

            for (int iw = 0; iw < atom_nw[atom_type]; iw++)
            {
                if (atom_iw2_new[it_nw + iw])
                {
                    psi = c1 * psi_u[iw_nr] + c2 * dpsi_u[iw_nr]
                          + c3 * psi_u[iw_nr + 1] + c4 * dpsi_u[iw_nr + 1];
                }
                phi[phi_idx + iw] = psi * ylma[atom_iw2_ylm[it_nw + iw]];
                iw_nr += nrmax;
            }
        }
        else
        {
            int phi_idx = atoms_phi_start[atom_id + pre_atoms_num] +
                          bgrids_phi_len[bgrid_id] * mgrid_id;
            for (int iw = 0; iw < atom_nw[atom_type]; iw++)
            {
                phi[phi_idx + iw] = 0.0;
            }
        }
    }
}

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
        const double3 rcoord = atoms_bgrids_rcoords[atom_id + pre_atoms_num];
        const double3 coord = make_double3(mgrid_pos.x-rcoord.x,
                                           mgrid_pos.y-rcoord.y,
                                           mgrid_pos.z-rcoord.z);
        double dist = sqrt(coord.x * coord.x + coord.y * coord.y + coord.z * coord.z);
        if (dist < rcut[atom_type])
        {
            if (dist < 1.0E-9)
            { dist += 1.0E-9; }
            // since nwl is less or equal than 5, the size of rly is (5+1)^2
            // size of grly = 36 * 3
            double rly[36];
            double grly[36 * 3];
            const int nwl = ucell_atom_nwl[atom_type];
            grad_rl_sph_harm(nwl, ylmcoef, coord.x, coord.y, coord.z, rly, grly);

            // interpolation
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
            int phi_idx = atoms_phi_start[atom_id + pre_atoms_num] +
                          bgrids_phi_len[bgrid_id] * mgrid_id;
            for (int iw = 0; iw < atom_nw[atom_type]; iw++)
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
                const double tmprl = tmp / rl;

                // if phi == nullptr, it means that we only need dphi.
                if(phi != nullptr)
                {
                    phi[phi_idx + iw] = tmprl * rly[idx_ylm];
                }
                // derivative of wave functions with respect to atom positions.
                const double tmpdphi_rly = (dtmp - tmp * iw_l / dist) / rl * rly[idx_ylm] / dist;

                dphi_x[phi_idx + iw] =  tmpdphi_rly * coord.x + tmprl * grly[idx_ylm * 3 + 0];
                dphi_y[phi_idx + iw] =  tmpdphi_rly * coord.y + tmprl * grly[idx_ylm * 3 + 1];
                dphi_z[phi_idx + iw] =  tmpdphi_rly * coord.z + tmprl * grly[idx_ylm * 3 + 2];
            }
        }
        else
        {
            int phi_idx = atoms_phi_start[atom_id + pre_atoms_num] +
                          bgrids_phi_len[bgrid_id] * mgrid_id;
            for (int iw = 0; iw < atom_nw[atom_type]; iw++)
            {
                if(phi != nullptr)
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

__global__ void phi_mul_vldr3_kernel(
    const double* __restrict__ vl,
    const double dr3,
    const double* __restrict__ phi,
    const int mgrids_per_bgrid,
    const int* __restrict__ mgrids_local_idx,
    const int* __restrict__ bgrids_phi_len,
    const int* __restrict__ bgrids_phi_start,
    double* __restrict__ result)
{
    const int bgrid_id = blockIdx.y;
    const int mgrid_id = blockIdx.x;
    const int phi_len = bgrids_phi_len[bgrid_id];
    const int phi_start = bgrids_phi_start[bgrid_id] + mgrid_id * phi_len;
    const int mgrid_id_in_batch = bgrid_id * mgrids_per_bgrid + mgrid_id;
    const double vldr3 =  vl[mgrids_local_idx[mgrid_id_in_batch]] * dr3;
    for(int i = threadIdx.x; i < phi_len; i += blockDim.x)
    {
        result[phi_start + i] = phi[phi_start + i] * vldr3;
    }
}

// rho(ir) = \sum_{iwt} \phi_i(ir,iwt) * \phi_j^*(ir,iwt)
// each block calculate the dot product of phi_i and phi_j of a meshgrid
__global__ void phi_dot_phi_kernel(
    const double* __restrict__ phi_i,
    const double* __restrict__ phi_j,
    const int mgrids_per_bgrid,
    const int* __restrict__ mgrids_local_idx,
    const int* __restrict__ bgrids_phi_len,
    const int* __restrict__ bgrids_phi_start,
    double* __restrict__ rho)
{
    __shared__ double s_data[32];    // the length of s_data equals the max warp num of a block
    const int bgrid_id = blockIdx.y;
    const int mgrid_id = blockIdx.x;
    const int phi_len = bgrids_phi_len[bgrid_id];
    const int phi_start = bgrids_phi_start[bgrid_id] + mgrid_id * phi_len;
    const double* phi_i_mgrid = phi_i + phi_start;
    const double* phi_j_mgrid = phi_j + phi_start;
    const int mgrid_id_in_batch = bgrid_id * mgrids_per_bgrid + mgrid_id;
    const int mgrid_local_idx = mgrids_local_idx[mgrid_id_in_batch];
    const int tid = threadIdx.x;
    const int warp_id = tid / 32;
    const int lane_id = tid % 32;
    double tmp_sum = 0;

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

    tmp_sum = (tid < blockDim.x / 32) ? s_data[tid] : 0;
    if(warp_id == 0)
    {
        tmp_sum = warpReduceSum(tmp_sum);
    }

    if(tid == 0)
    {
        rho[mgrid_local_idx] = tmp_sum;
    }
}

}