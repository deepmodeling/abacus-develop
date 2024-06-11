#include "gint_force.cuh"
#include "interp.cuh"
#include "module_hamilt_lcao/module_gint/gint_force_gpu.h"
#include "cuda_tools.cuh"
#include "gint_force.cuh"
#include "sph.cuh"
#include "cuda_runtime.h"
// CUDA kernel to calculate psi and force
namespace GintKernel
{
__global__ void get_psi_force(double* ylmcoef,
                              double delta_r_g,
                              int bxyz_g,
                              double nwmax_g,
                              double* __restrict__ input_double,
                              int* __restrict__ input_int,
                              int* __restrict__ num_psir,
                              int psi_size_max,
                              const int* __restrict__ ucell_atom_nwl,
                              const bool* __restrict__ atom_iw2_new,
                              const int* __restrict__ atom_iw2_ylm,
                              const int* __restrict__ atom_iw2_l,
                              const int* __restrict__ atom_nw,
                              int nr_max,
                              const double* __restrict__ psi_u,
                              double* psir_r,
                              double* psir_lx,
                              double* psir_ly,
                              double* psir_lz,
                              double* psir_lxx,
                              double* psir_lxy,
                              double* psir_lxz,
                              double* psir_lyy,
                              double* psir_lyz,
                              double* psir_lzz)
{
    int size = num_psir[blockIdx.x];
    int start_index = psi_size_max * blockIdx.x;
    int end_index = start_index + size;
    start_index += threadIdx.x + blockDim.x * blockIdx.y;
    for (int index = start_index; index < end_index;
         index += blockDim.x * gridDim.y)
    {
        double dr[3];
        int index_double = index * 5;
        dr[0] = input_double[index_double];
        dr[1] = input_double[index_double + 1];
        dr[2] = input_double[index_double + 2];
        double distance = input_double[index_double + 3];
        double vlbr3_value = input_double[index_double + 4];
        double ylma[49]; // Attention!!! At present, we only use L=5 at
                         // most. So (L+1) * (L+1)=36
        double grly[49][3];
        int index_int = index * 2;
        int it = input_int[index_int];
        int dist_tmp = input_int[index_int + 1];

        int nwl = ucell_atom_nwl[it];
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
                      psir_r,
                      dist_tmp,
                      ylma,
                      vlbr3_value,
                      psir_lx,
                      dr,
                      grly,
                      psir_ly,
                      psir_lz,
                      psir_lxx,
                      psir_lxy,
                      psir_lxz,
                      psir_lyy,
                      psir_lyz,
                      psir_lzz);
    }
}


/**
 * \brief Compute dot product of stress components and partial derivatives.
 *
 * This CUDA kernel computes the dot product of stress components and partial
 * derivatives based on the input arrays.
 */

__global__ void dot_product_stress(double* psir_lxx,
                                   double* psir_lxy,
                                   double* psir_lxz,
                                   double* psir_lyy,
                                   double* psir_lyz,
                                   double* psir_lzz,
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
        tmp[0] += psir_lxx[tid] * psi_dm_2;
        tmp[1] += psir_lxy[tid] * psi_dm_2;
        tmp[2] += psir_lxz[tid] * psi_dm_2;
        tmp[3] += psir_lyy[tid] * psi_dm_2;
        tmp[4] += psir_lyz[tid] * psi_dm_2;
        tmp[5] += psir_lzz[tid] * psi_dm_2;
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

/**
 * @brief Calculate the dot product force.
 *
 * This function calculates the dot product force based on the provided
 * parameters.
 *
 * @param psir_lx Pointer to the array of psir_lx values.
 * @param psir_ly Pointer to the array of psir_ly values.
 * @param psir_lz Pointer to the array of psir_lz values.
 * @param psir_ylm_dm Pointer to the array of psir_ylm_dm values.
 * @param force_dot Pointer to the array where the calculated force will be
 * stored.
 * @param iat Pointer to the array of iat values.
 * @param nwmax Maximum value for nwmax.
 */

__global__ void dot_product_force(double* __restrict__ psir_lx,
                                  double* __restrict__ psir_ly,
                                  double* __restrict__ psir_lz,
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

    for (int i = tid; i < nwmax; i += blockDim.x)
    {
        int ls_offset = tid * 3;
        int psi_offset = offset + i;
        double psi_dm_2 = psir_ylm_dm[psi_offset] * 2;
        localsum[ls_offset] += psir_lx[psi_offset] * psi_dm_2;
        localsum[ls_offset + 1] += psir_ly[psi_offset] * psi_dm_2;
        localsum[ls_offset + 2] += psir_lz[psi_offset] * psi_dm_2;
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
