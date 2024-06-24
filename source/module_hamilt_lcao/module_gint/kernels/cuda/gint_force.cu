#include "sph.cuh"
#include "interp.cuh"
#include "gint_force.cuh"
#include "cuda_tools.cuh"
#include "cuda_runtime.h"
// CUDA kernel to calculate psi and force
namespace GintKernel
{
__global__ void get_psi_force(double* ylmcoef,
                              double delta_r_g,
                              int bxyz_g,
                              const int nwmax_g,
                              double* __restrict__ psi_input_double,
                              int* __restrict__ psi_input_int,
                              int* __restrict__ atom_num_per_bcell,
                              int* __restrict__ start_idx_per_bcell,
                              const int* __restrict__ ucell_atom_nwl,
                              const bool* __restrict__ atom_iw2_new,
                              const int* __restrict__ atom_iw2_ylm,
                              const int* __restrict__ atom_iw2_l,
                              const int* __restrict__ atom_nw,
                              int nr_max,
                              const double* __restrict__ psi_u,
                              double* psi,
                              double* dpsi,
                              double* d2psi)
{
    const int size = atom_num_per_bcell[blockIdx.x];
    const int bcell_start = start_idx_per_bcell[blockIdx.x];
    const int end_index = bcell_start + size;
    const int start_index = bcell_start + threadIdx.x + blockDim.x * blockIdx.y;
    for (int index = start_index; index < end_index;
         index += blockDim.x * gridDim.y)
    {
        double dr[3];
        const int index_double = index * 5;
        dr[0] = psi_input_double[index_double];
        dr[1] = psi_input_double[index_double + 1];
        dr[2] = psi_input_double[index_double + 2];
        const double distance = psi_input_double[index_double + 3];
        const double vlbr3_value = psi_input_double[index_double + 4];
        double ylma[49]; // Attention!!! At present, we only use L=5 at
                         // most. So (L+1) * (L+1)=36
        double grly[49][3];
        const int index_int = index * 2;
        const int it = psi_input_int[index_int];
        int dist_tmp = psi_input_int[index_int + 1];

        const int nwl = ucell_atom_nwl[it];
        spherical_harmonics_d(dr, distance*distance, grly, nwl, ylma, ylmcoef);
        interpolate_f(distance,
                      delta_r_g,
                      it,
                      nwmax_g,
                      nr_max,
                      atom_nw,
                      atom_iw2_new,
                      psi_u,
                      atom_iw2_l,
                      atom_iw2_ylm,
                      psi,
                      dist_tmp,
                      ylma,
                      vlbr3_value,
                      dpsi,
                      dr,
                      grly,
                      d2psi);
    }
}

__global__ void dot_product_stress(double* d2psi,
                                   double* psir_ylm_dm,
                                   double* stress_dot,
                                   int elements_num)
{
    __shared__ double cache[256][6]; 
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int cacheIndex = threadIdx.x;
    double tmp[6] = {0.0};
    while (tid < elements_num)
    {   
        double psi_dm_2 = psir_ylm_dm[tid] * 2;
        const int tid_stress = tid * 6;
        tmp[0] += d2psi[tid_stress] * psi_dm_2;
        tmp[1] += d2psi[tid_stress + 1] * psi_dm_2;
        tmp[2] += d2psi[tid_stress + 2] * psi_dm_2;
        tmp[3] += d2psi[tid_stress + 3] * psi_dm_2;
        tmp[4] += d2psi[tid_stress + 4] * psi_dm_2;
        tmp[5] += d2psi[tid_stress + 5] * psi_dm_2;
        tid += blockDim.x * gridDim.x;
    }

    for (int i = 0; i < 6; i++)
    {
        cache[cacheIndex][i] = tmp[i];
    }
    __syncthreads();

    int i = blockDim.x / 2;
    while (i != 0)
    {
        if (cacheIndex < i)
        {
            for (int index = 0; index < 6; index++)
            {
                cache[cacheIndex][index] += cache[cacheIndex + i][index];
            }
        }
        __syncthreads();
        i /= 2;
    }

    if (cacheIndex == 0){
        for (int index = 0; index < 6; index++)
        {
            atomicAdd(&stress_dot[index], cache[0][index]); // Use atomicAdd() instead of atomic_add().
            // stress_dot[blockIdx.x + gridDim.x * index] = cache[0][index];
        }
    }
}

__global__ void dot_product_force(double* __restrict__ dpsi,
                                  double* __restrict__ psir_ylm_dm,
                                  double* force_dot,
                                  int* iat,
                                  int nwmax)
{
    extern __shared__ double localsum[];
    int tid = threadIdx.x;
    int bid = blockIdx.x;
    int iat_on_nbz = iat[bid];
    if(iat_on_nbz <= -1)
    {
        return;
    }

    int offset = bid * nwmax;
    localsum[tid * 3] = 0.0;
    localsum[tid * 3 + 1] = 0.0;
    localsum[tid * 3 + 2] = 0.0;

    for (int i = tid; i < nwmax; i += blockDim.x)
    {
        int ls_offset = tid * 3;
        int psi_offset = offset + i;
        int psi_offset_force = psi_offset * 3;
        double psi_dm_2 = psir_ylm_dm[psi_offset] * 2;
        localsum[ls_offset] += dpsi[psi_offset_force] * psi_dm_2;
        localsum[ls_offset + 1] += dpsi[psi_offset_force + 1] * psi_dm_2;
        localsum[ls_offset + 2] += dpsi[psi_offset_force + 2] * psi_dm_2;
    }
    __syncthreads();
    
    for (int i = blockDim.x / 2; i > 0; i >>= 1)
    {
        if (tid < i)
        {
            int ls_offset = tid * 3;
            localsum[ls_offset] += localsum[ls_offset + i * 3];
            localsum[ls_offset + 1] += localsum[ls_offset + i * 3 + 1];
            localsum[ls_offset + 2] += localsum[ls_offset + i * 3 + 2];
        }
        __syncthreads();
    }

    if(tid == 0)
    {
        for (int i = 0; i < 3; i++)
        {
            atomicAdd(&force_dot[iat_on_nbz*3 + i], localsum[i]);
        }
    }
}
} // namespace GintKernel
