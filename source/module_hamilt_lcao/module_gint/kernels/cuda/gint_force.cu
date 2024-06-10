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
                              double* input_double,
                              int* input_int,
                              int* num_psir,
                              int psi_size_max,
                              int* ucell_atom_nwl,
                              bool* atom_iw2_new,
                              int* atom_iw2_ylm,
                              int* atom_iw2_l,
                              int* atom_nw,
                              int nr_max,
                              double* psi_u,
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
        tmp[0] += psir_lxx[tid] * psir_ylm_dm[tid] * 2;
        tmp[1] += psir_lxy[tid] * psir_ylm_dm[tid] * 2;
        tmp[2] += psir_lxz[tid] * psir_ylm_dm[tid] * 2;
        tmp[3] += psir_lyy[tid] * psir_ylm_dm[tid] * 2;
        tmp[4] += psir_lyz[tid] * psir_ylm_dm[tid] * 2;
        tmp[5] += psir_lzz[tid] * psir_ylm_dm[tid] * 2;
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
 * @param max_size Maximum size for arrays.
 * @param elements_num Number of elements to process.
 */

__global__ void dot_product_force(double* psir_lx,
                                  double* psir_ly,
                                  double* psir_lz,
                                  double* psir_ylm_dm,
                                  double* force_dot,
                                  int* iat,
                                  int nwmax,
                                  int elements_num)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < elements_num)
    {
        int iat_on_nbz = iat[tid];
        if (iat_on_nbz <= -1)
        {
            tid += blockDim.x * gridDim.x;
            continue;
        }

        int iat_index = tid * 3;
        int dist = tid * nwmax;
        double tmp[3] = {0.0};

        for (int i = 0; i < nwmax; i++)
        {
            tmp[0] += psir_lx[dist + i] * psir_ylm_dm[dist + i] * 2;
            tmp[1] += psir_ly[dist + i] * psir_ylm_dm[dist + i] * 2;
            tmp[2] += psir_lz[dist + i] * psir_ylm_dm[dist + i] * 2;
        }
        
        for (int i = 0; i < 3; i++)
        {
            atomicAdd(&force_dot[iat_on_nbz*3 + i], tmp[i]);
        }
        tid += blockDim.x * gridDim.x;
    }
}
} // namespace GintKernel
