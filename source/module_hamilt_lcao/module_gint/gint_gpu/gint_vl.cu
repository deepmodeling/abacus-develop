#include "module_hamilt_lcao/module_gint/gint_gpu/gint_vl.cuh"
#include "module_hamilt_lcao/module_gint/gint_gpu/gint_vl.h"
#include "module_hamilt_lcao/module_gint/gint_gpu/vbatch_matrix_multiple/cuda_tools.cuh"
#include "sph.cuh"
namespace lcaoCudaKernel{

__global__ void get_psi_and_vldr3(double *ylmcoef,
                                  double delta_r_g,
                                  double bxyz_g,
                                  double nwmax_g,
                                  double *input_double,
                                  int *input_int,
                                  int *num_psir,
                                  int psi_size_max,
                                  int *ucell_atom_nwl,
                                  bool *atom_iw2_new,
                                  int *atom_iw2_ylm,
                                  int *atom_nw,
                                  int nr_max,
                                  double *psi_u,
                                  double *psir_ylm_left,
                                  double *psir_ylm_right)
{
    int size = num_psir[blockIdx.x];
    int start_index = psi_size_max * blockIdx.x;
    int end_index = start_index + size;
    start_index += threadIdx.x + blockDim.x * blockIdx.y;
    for (int index = start_index; index < end_index; index += blockDim.x * gridDim.y)
    {
        double dr[3];
        int index_double = index * 5;
        dr[0] = input_double[index_double];
        dr[1] = input_double[index_double + 1];
        dr[2] = input_double[index_double + 2];
        double distance = input_double[index_double + 3];
        double vlbr3_value = input_double[index_double + 4];
        double ylma[49];
        int index_int = index * 2;
        int it = input_int[index_int];
        int dist_tmp = input_int[index_int + 1];
        int nwl = ucell_atom_nwl[it];
        spherical_harmonics(dr,distance,nwl,ylma,ylmcoef);
    
        distance /= delta_r_g;

        int ip = (int)(distance);
        double dx = distance - ip;
        double dx2 = dx * dx;
        double dx3 = dx2 * dx;

        double c3 = 3.0 * dx2 - 2.0 * dx3;
        double c1 = 1.0 - c3;
        double c2 = (dx - 2.0 * dx2 + dx3) * delta_r_g;
        double c4 = (dx3 - dx2) * delta_r_g;

        double phi = 0.0;
        int it_nw = it * nwmax_g;
        int iw_nr = (it_nw * nr_max + ip) * 2;
        int it_nw_iw = it_nw;
        for (int iw = 0; iw < atom_nw[it]; ++iw)
        {
            if (atom_iw2_new[it_nw_iw])
            {
                phi = c1 * psi_u[iw_nr] + c2 * psi_u[iw_nr + 1] + c3 * psi_u[iw_nr + 2] + c4 * psi_u[iw_nr + 3];
            }
            double temp = phi * ylma[atom_iw2_ylm[it_nw_iw]];
            psir_ylm_right[dist_tmp] = temp * vlbr3_value;
            psir_ylm_left[dist_tmp] = temp;
            dist_tmp += bxyz_g;
            iw_nr += nr_max;
            iw_nr += nr_max;
            it_nw_iw++;
        }
    }
}

__global__ void psi_multiple(int bxyz_g, const int* m, int* n,
                                double  const * const * global_A_array,
                                double const * const * global_B_array,
                                double ** global_C_array)
{
    int atom_pair_index = blockIdx.x;
    int nw_mul = m[atom_pair_index] * n[atom_pair_index];
    int atom_nw2 = n[atom_pair_index];
    const double * atom_left = global_A_array[atom_pair_index];
    const double * atom_right = global_B_array[atom_pair_index];
    double * output = global_C_array[atom_pair_index];
    #pragma unroll
    for (int iw_index = threadIdx.x; iw_index < nw_mul; iw_index += blockDim.x)
    {
        int iw1 = iw_index / atom_nw2 * bxyz_g;
        int iw2 = iw_index % atom_nw2;
        double v2 = 0.0;
        const double * left = &atom_left[iw1];
        const double * right = &atom_right[iw2 * bxyz_g];
        #pragma unroll
        for (int ib = 0; ib <  bxyz_g; ++ib)
        {
            v2 += left[ib] * right[ib];
        }
        atomicAdd(&(output[iw1 + iw2]), v2);
    }
}

} // namespace lcaoCudaKernel