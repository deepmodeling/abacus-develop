#include "phi_operator_kernel.cuh"
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
            double dr[3] = { coord.x / dist, coord.y / dist, coord.z / dist };
            // since nwl is less or equal than 5, the size of ylma is (5+1)^2
            double ylma[36];
            const int nwl = ucell_atom_nwl[atom_type];
            sph_harm(nwl, ylmcoef, dr, ylma);

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
}