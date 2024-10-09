#include "interp.cuh"
#include "gint_rho.cuh"
#include "sph.cuh"
#include "cuda_tools.cuh"

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


/*
    each block calculates the wavefunction on a meshcell,
    and each thread loops over the atoms on a meshcell.
*/
__global__ void get_psi(const double* const ylmcoef,
                        const double delta_r,
                        const int bxyz,
                        const int nwmax,
                        const int max_atom,
                        const int* const ucell_atom_nwl,
                        const bool* const atom_iw2_new,
                        const int* const atom_iw2_ylm,
                        const int* const atom_nw,
                        const double* const rcut,
                        const int nr_max,
                        const double* const psi_u,
                        const double* const mcell_pos,
                        const double* const dr_part,
                        const uint8_t* const atoms_type,
                        const int* const atoms_num_info,
                        double* psi)
{
    const int bcell_id = blockIdx.x;
    const int num_atoms = atoms_num_info[2 * bcell_id];
    const int pre_atoms = atoms_num_info[2 * bcell_id + 1];
    const int mcell_id = blockIdx.y;
    const double mcell_pos_x = mcell_pos[3 * mcell_id];
    const double mcell_pos_y = mcell_pos[3 * mcell_id + 1];
    const double mcell_pos_z = mcell_pos[3 * mcell_id + 2];

    for(int atom_id = threadIdx.x; atom_id < num_atoms; atom_id += blockDim.x)
    {
        const int aid = pre_atoms + atom_id;
        const double dr_x = dr_part[aid * 3] + mcell_pos_x;
        const double dr_y = dr_part[aid * 3 + 1] + mcell_pos_y;
        const double dr_z = dr_part[aid * 3 + 2] + mcell_pos_z;
        double dist = sqrt(dr_x * dr_x + dr_y * dr_y + dr_z * dr_z);
        const int atype = __ldg(atoms_type + aid);
        if(dist < rcut[atype])
        {
            if (dist < 1.0E-9)
            {
                dist += 1.0E-9;
            }
            double dr[3] = {dr_x / dist, dr_y / dist, dr_z / dist};
            double ylma[49];
            const int nwl = __ldg(ucell_atom_nwl + atype);
            int psi_idx = (pre_atoms * bxyz + mcell_id * num_atoms + atom_id) * nwmax;
            spherical_harmonics(dr, nwl, ylma, ylmcoef);
            interp_rho(dist,
                       delta_r,
                       atype,
                       nwmax,
                       nr_max,
                       atom_nw,
                       atom_iw2_new,
                       psi_u,
                       ylma,
                       atom_iw2_ylm,
                       psi,
                       psi_idx);
        }
    }
}

/*
    Each block calculates the dot product on a meshcell,
    and each thread loops over the wavefunction of atoms on a meshcell.
*/
__global__ void psir_dot(const int bxyz,
                         const int nwmax,
                         const int* atoms_num_info,
                         const double* __restrict__ vec_a_g,
                         const double* __restrict__  vec_b_g,
                         double** results_g)
{
    __shared__ double s_data[32];
    const int tid = threadIdx.x;
    const int bcell_id = blockIdx.x;
    const int mcell_id = blockIdx.y;
    const int vec_size = atoms_num_info[2 * bcell_id] * nwmax;
    const int offset = atoms_num_info[2 * bcell_id + 1] * nwmax * bxyz + mcell_id * vec_size;
    const double* vec_a_mcell = vec_a_g + offset;
    const double* vec_b_mcell = vec_b_g + offset;
    const int warp_id = tid / 32;
    const int lane_id = tid % 32;
    double mySum = 0;

    for (int k = tid; k < vec_size; k += blockDim.x)
    {
        mySum += vec_a_mcell[k] * vec_b_mcell[k];
    }    

    mySum = warpReduceSum(mySum);

    if (lane_id == 0)
    {
        s_data[warp_id] = mySum;
    }
    __syncthreads();

    mySum = (tid < blockDim.x / 32) ? s_data[tid] : 0;
    if (warp_id == 0)
    {
        mySum = warpReduceSum(mySum);
    }

    if (tid == 0) {
        *results_g[bcell_id*bxyz + mcell_id] = mySum;
    }
}
} // namespace GintKernel