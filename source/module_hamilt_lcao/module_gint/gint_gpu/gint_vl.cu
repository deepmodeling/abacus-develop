#include "module_hamilt_lcao/module_gint/gint_gpu/gint_vl.cuh"
#include "module_hamilt_lcao/module_gint/gint_gpu/gint_vl.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_lcao/module_gint/gint_gpu/vbatch_matrix_multiple/cuda_tools.cuh"

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
        // begin calculation
        double ylma[49]; // Attention!!! At present, we only use L=5 at most. So (L+1) * (L+1)=36
        int index_int = index * 2;
        int it = input_int[index_int];
        int dist_tmp = input_int[index_int + 1];

        /***************************
        L = 0
        ***************************/
        ylma[0] = ylmcoef[0]; // l=0, m=0
        int nwl = ucell_atom_nwl[it];
        double tmp0;
        if (nwl == 0)
            goto YLM_END;

        /***************************
        L = 1
        ***************************/
        ylma[1] = ylmcoef[1] * dr[2];  // l=1, m=0
        ylma[2] = -ylmcoef[1] * dr[0]; // l=1, m=1
        ylma[3] = -ylmcoef[1] * dr[1]; // l=1, m=-1
        if (nwl == 1)
            goto YLM_END;

        /***************************
        L = 2
        ***************************/
        ylma[4] = ylmcoef[2] * dr[2] * ylma[1] - ylmcoef[3] * ylma[0]; // l=2, m=0
        tmp0 = ylmcoef[4] * dr[2];
        ylma[5] = tmp0 * ylma[2]; // l=2,m=1
        ylma[6] = tmp0 * ylma[3]; // l=2,m=-1

        tmp0 = ylmcoef[4] * dr[0];
        ylma[7] = ylmcoef[5] * ylma[4] - ylmcoef[6] * ylma[0] - tmp0 * ylma[2]; // l=2,m=2
        ylma[8] = -tmp0 * ylma[3];
        if (nwl == 2)
            goto YLM_END;

        /***************************
        L = 3
        ***************************/
        ylma[9] = ylmcoef[7] * dr[2] * ylma[4] - ylmcoef[8] * ylma[1]; // l=3, m=0

        tmp0 = ylmcoef[9] * dr[2];
        ylma[10] = tmp0 * ylma[5] - ylmcoef[10] * ylma[2]; // l=3,m=1
        ylma[11] = tmp0 * ylma[6] - ylmcoef[10] * ylma[3]; // l=3,m=-1

        tmp0 = ylmcoef[11] * dr[2];
        ylma[12] = tmp0 * ylma[7]; // l=3,m=2
        ylma[13] = tmp0 * ylma[8]; // l=3,m=-2

        tmp0 = ylmcoef[14] * dr[0];
        ylma[14] = ylmcoef[12] * ylma[10] - ylmcoef[13] * ylma[2] - tmp0 * ylma[7]; // l=3,m=3
        ylma[15] = ylmcoef[12] * ylma[11] - ylmcoef[13] * ylma[3] - tmp0 * ylma[8]; // l=3,m=-3
        if (nwl == 3)
            goto YLM_END;

        /***************************
        L = 4
        ***************************/
        ylma[16] = ylmcoef[15] * dr[2] * ylma[9] - ylmcoef[16] * ylma[4]; // l=4,m=0

        tmp0 = ylmcoef[17] * dr[2];
        ylma[17] = tmp0 * ylma[10] - ylmcoef[18] * ylma[5]; // l=4,m=1
        ylma[18] = tmp0 * ylma[11] - ylmcoef[18] * ylma[6]; // l=4,m=-1

        tmp0 = ylmcoef[19] * dr[2];
        ylma[19] = tmp0 * ylma[12] - ylmcoef[20] * ylma[7]; // l=4,m=2
        ylma[20] = tmp0 * ylma[13] - ylmcoef[20] * ylma[8]; // l=4,m=-2

        tmp0 = 3.0 * dr[2];
        ylma[21] = tmp0 * ylma[14]; // l=4,m=3
        ylma[22] = tmp0 * ylma[15]; // l=4,m=-3

        tmp0 = ylmcoef[23] * dr[0];
        ylma[23] = ylmcoef[21] * ylma[19] - ylmcoef[22] * ylma[7] - tmp0 * ylma[14]; // l=4,m=4
        ylma[24] = ylmcoef[21] * ylma[20] - ylmcoef[22] * ylma[8] - tmp0 * ylma[15]; // l=4,m=-4
        if (nwl == 4)
            goto YLM_END;

        /***************************
        L = 5
        ***************************/
        ylma[25] = ylmcoef[24] * dr[2] * ylma[16] - ylmcoef[25] * ylma[9]; // l=5,m=0

        tmp0 = ylmcoef[26] * dr[2];
        ylma[26] = tmp0 * ylma[17] - ylmcoef[27] * ylma[10]; // l=5,m=1
        ylma[27] = tmp0 * ylma[18] - ylmcoef[27] * ylma[11]; // l=5,m=-1

        tmp0 = ylmcoef[28] * dr[2];
        ylma[28] = tmp0 * ylma[19] - ylmcoef[29] * ylma[12]; // l=5,m=2
        ylma[29] = tmp0 * ylma[20] - ylmcoef[29] * ylma[13]; // l=5,m=-2

        tmp0 = ylmcoef[30] * dr[2];
        ylma[30] = tmp0 * ylma[21] - ylmcoef[31] * ylma[14]; // l=5,m=3
        ylma[31] = tmp0 * ylma[22] - ylmcoef[31] * ylma[15]; // l=5,m=-3

        tmp0 = ylmcoef[32] * dr[2];
        ylma[32] = tmp0 * ylma[23]; // l=5,m=4
        ylma[33] = tmp0 * ylma[24]; // l=5,m=-4

        tmp0 = ylmcoef[35] * dr[0];
        ylma[34] = ylmcoef[33] * ylma[30] - ylmcoef[34] * ylma[14] - tmp0 * ylma[23]; // l=5,m=5
        ylma[35] = ylmcoef[33] * ylma[31] - ylmcoef[34] * ylma[15] - tmp0 * ylma[24]; // l=5,m=-5
        if (nwl == 5)
            goto YLM_END;
        /*
        // if nwl > 5
        for (int il = 6; il <= nwl; il++)
        {
            int istart = il * il;
            int istart1 = (il - 1) * (il - 1);
            int istart2 = (il - 2) * (il - 2);

            double fac2 = sqrt(4.0 * istart - 1.0);
            double fac4 = sqrt(4.0 * istart1 - 1.0);

            for (int im = 0; im < 2 * il - 1; im++)
            {
                int imm = (im + 1) / 2;
                ylma[istart + im] = fac2 / sqrt((double)istart - imm * imm) * (dr[2] * ylma[istart1 + im] - sqrt((double)istart1 - imm * imm) / fac4 * ylma[istart2 + im]);
            }

            double bl1 = sqrt(2.0 * il / (2.0 * il + 1.0));
            double bl2 = sqrt((2.0 * il - 2.0) / (2.0 * il - 1.0));
            double bl3 = sqrt(2.0) / fac2;

            ylma[istart + 2 * il - 1] = (bl3 * ylma[istart + 2 * il - 5] - bl2 * ylma[istart2 + 2 * il - 5] - 2.0 * dr[0] * ylma[istart1 + 2 * il - 3]) / bl1;
            ylma[istart + 2 * il] = (bl3 * ylma[istart + 2 * il - 4] - bl2 * ylma[istart2 + 2 * il - 4] - 2.0 * dr[0] * ylma[istart1 + 2 * il - 2]) / bl1;
        }*/
    YLM_END:
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
