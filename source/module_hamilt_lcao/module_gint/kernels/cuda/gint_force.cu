#include "sph.cuh"
#include "interp.cuh"
#include "gint_force.cuh"
#include "cuda_tools.cuh"
#include "module_base/module_device/device.h"
// CUDA kernel to calculate psi and force
namespace GintKernel
{
__inline__ __device__ double warpReduceSum(double val)
{   
    val += __shfl_xor_sync(0xffffffff, val, 16, 32);
    val += __shfl_xor_sync(0xffffffff, val, 8, 32);
    val += __shfl_xor_sync(0xffffffff, val, 4, 32);
    val += __shfl_xor_sync(0xffffffff, val, 2, 32);
    val += __shfl_xor_sync(0xffffffff, val, 1, 32);
    return val;
}


__global__ void get_psi_force(double* ylmcoef,
                              double delta_r,
                              int bxyz,
                              const int nwmax,
                              const int max_atom,
                              const int* const ucell_atom_nwl,
                              const bool* const atom_iw2_new,
                              const int* const atom_iw2_ylm,
                              const int* const atom_iw2_l,
                              const int* const atom_nw,
                              const double* const rcut,
                              const int nr_max,
                              const double* const psi_u,
                              const double* const mcell_pos,
                              const double* const dr_part,
                              const double* const vldr3,
                              const uint8_t* const atoms_type,
                              const int* const atoms_num_info,
                              double* psi,
                              double* dpsi,
                              double* d2psi)
{
    const int bcell_id = blockIdx.x;
    const int num_atoms = atoms_num_info[2 * bcell_id];
    const int pre_atoms = atoms_num_info[2 * bcell_id + 1];
    const int mcell_id = blockIdx.y;
    const double vldr3_value = vldr3[bcell_id*bxyz + mcell_id];
    const double mcell_pos_x = mcell_pos[3 * mcell_id];
    const double mcell_pos_y = mcell_pos[3 * mcell_id + 1];
    const double mcell_pos_z = mcell_pos[3 * mcell_id + 2];

    for(int atom_id = threadIdx.x; atom_id < num_atoms; atom_id += blockDim.x)
    {
        const int dr_start = 3 * (pre_atoms + atom_id);
        const double dr_x = dr_part[dr_start] + mcell_pos_x;
        const double dr_y = dr_part[dr_start + 1] + mcell_pos_y;
        const double dr_z = dr_part[dr_start + 2] + mcell_pos_z;
        double dist = sqrt(dr_x * dr_x + dr_y * dr_y + dr_z * dr_z);
        const int atype = __ldg(atoms_type + pre_atoms + atom_id);
        if(dist < rcut[atype])
        {
            if (dist < 1.0E-9)
            {
                dist += 1.0E-9;
            }
            // dr is different from that in interp_rho and interp_vl
            double dr[3] = {dr_x, dr_y, dr_z};
            double ylma[49];
            double grly[49][3];
            const int nwl = __ldg(ucell_atom_nwl + atype);
            spherical_harmonics_d(dr, dist*dist, grly, nwl, ylma, ylmcoef);
            int psi_idx = ((pre_atoms + atom_id) * bxyz + mcell_id) * nwmax;
            interp_f(dist,
                     delta_r,
                     atype,
                     nwmax,
                     nr_max,
                     atom_nw,
                     atom_iw2_new,
                     psi_u,
                     ylma,
                     atom_iw2_l,
                     atom_iw2_ylm,
                     vldr3_value,
                     dr,
                     grly,
                     psi_idx,
                     psi,
                     dpsi,
                     d2psi);
        }
    }
}


__global__ void dot_product_stress(const double* d2psi,
                                   const double* psi_dm,
                                   const int size,
                                   double* stress)
{
    __shared__ double cache[32 * 6];
    const int tid = threadIdx.x;
    const int stride = blockDim.x * gridDim.x;
    const int warp_id = tid / 32;
    const int lane_id = tid % 32;
    double tmp[6] = {0.0};
    for(int id = threadIdx.x + blockIdx.x * blockDim.x; id < size; id += stride)
    {   
        const double psi_dm_2 = psi_dm[id] * 2;
        const int id_stress = id * 6;
        tmp[0] += d2psi[id_stress] * psi_dm_2;
        tmp[1] += d2psi[id_stress + 1] * psi_dm_2;
        tmp[2] += d2psi[id_stress + 2] * psi_dm_2;
        tmp[3] += d2psi[id_stress + 3] * psi_dm_2;
        tmp[4] += d2psi[id_stress + 4] * psi_dm_2;
        tmp[5] += d2psi[id_stress + 5] * psi_dm_2;
    }

    for(int i = 0; i<6; i++)
    {
        tmp[i] = warpReduceSum(tmp[i]);
    }

    if (lane_id == 0)
    {
        for (int i = 0; i < 6; i++)
        {
            cache[warp_id * 6 + i] = tmp[i];
        }
    }
    __syncthreads();

    for (int i = 0; i < 6; i++)
    {
        tmp[i] = (tid < blockDim.x / 32) ? cache[tid * 6 + i] : 0;
    }

    if(warp_id == 0)
    {
        for (int i = 0; i < 6; i++)
        {
            tmp[i] = warpReduceSum(tmp[i]);
        }
    }

    if (tid == 0)
    {
        for (int i = 0; i < 6; i++)
        {
            atomicAdd(&stress[i], tmp[i]); // Use atomicAdd() instead of atomic_add().
        }
    }
}


__global__ void dot_product_force(const int bxyz,
                                  const int nwmax,
                                  const int *atoms_num_info,
                                  const int *iat_on_nbz,
                                  const double* dpsi,
                                  const double* psi_dm,
                                  double* force)
{
    __shared__ double cache[32 * 3];
    const int tid = threadIdx.x;
    const int bcell_id = blockIdx.x;
    const int warp_id = tid / 32;
    const int lane_id = tid % 32;
    const int vec_size = bxyz * nwmax;
    const int atom_num = atoms_num_info[2 * bcell_id];
    const int pre_atoms = atoms_num_info[2 * bcell_id + 1];

    for(int k = 0; k < atom_num; k++)
    {
        const int atom_id = pre_atoms + k;
        const int offset = atom_id * vec_size;
        const int iat = iat_on_nbz[atom_id];
        double force_iat[3] = {0.0};

        for(int i =tid; i < vec_size; i += blockDim.x)
        {
            int psi_offset = offset + i;
            double psi_dm_2 = psi_dm[psi_offset] * 2;
            force_iat[0] += dpsi[psi_offset * 3] * psi_dm_2;
            force_iat[1] += dpsi[psi_offset * 3 + 1] * psi_dm_2;
            force_iat[2] += dpsi[psi_offset * 3 + 2] * psi_dm_2;
        }

        for (int i = 0; i < 3; i++)
        {
            force_iat[i] = warpReduceSum(force_iat[i]);
        }

        if (lane_id == 0)
        {
            for (int i = 0; i < 3; i++)
            {
                cache[warp_id * 3 + i] = force_iat[i];
            }
        }
        __syncthreads();

        for (int i = 0; i < 3; i++)
        {
            force_iat[i] = (tid < blockDim.x / 32) ? cache[tid * 3 + i] : 0;
        }

        if (warp_id == 0)
        {
            for (int i = 0; i < 3; i++)
            {
                force_iat[i] = warpReduceSum(force_iat[i]);
            }
        }

        if (tid == 0)
        {
            for (int i = 0; i < 3; i++)
            {
                atomicAdd(&force[iat * 3 + i], force_iat[i]);
            }
        }
    }
}

} // namespace GintKernel
