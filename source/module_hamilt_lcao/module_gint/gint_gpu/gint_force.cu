#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_lcao/module_gint/gint_gpu/gint_force.cuh"
#include "module_hamilt_lcao/module_gint/gint_gpu/gint_force.h"
#include "module_hamilt_lcao/module_gint/gint_gpu/vbatch_matrix_multiple/cuda_tools.cuh"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

// CUDA kernel to calculate psi and force
__global__ void get_psi_force(double *ylmcoef,
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
                        int *atom_iw2_l,
                        int *atom_nw,
                        int nr_max,
                        double *psi_u,
                        double *psir_ylm_right,
                        double *dpsir_ylm_left_x,
                        double *dpsir_ylm_left_y,
                        double *dpsir_ylm_left_z,
                        double *ddpsir_ylm_left_xx,
                        double *ddpsir_ylm_left_xy,
                        double *ddpsir_ylm_left_xz,
                        double *ddpsir_ylm_left_yy,
                        double *ddpsir_ylm_left_yz,
                        double *ddpsir_ylm_left_zz) {
  // Get the size of psi for the current block
  int size = num_psir[blockIdx.x];
  int start_index = psi_size_max * blockIdx.x;
  int end_index = start_index + size;
  start_index += threadIdx.x + blockDim.x * blockIdx.y;
  // Loop over the psi indices for the current block
  for (int index = start_index; index < end_index;
       index += blockDim.x * gridDim.y) {
    // Extract information from input arrays
    double dr[3];
    int index_double = index * 5;
    dr[0] = input_double[index_double];
    dr[1] = input_double[index_double + 1];
    dr[2] = input_double[index_double + 2];
    double distance = input_double[index_double + 3];
    distance = distance * distance;
    double tx = 2.0 * dr[0];
    double ty = 2.0 * dr[1];
    double tz = 2.0 * dr[2];
    double vlbr3_value = input_double[index_double + 4];
    // begin calculation
    double ylma[49]; // Attention!!! At present, we only use L=5 at most. So
                     // (L+1) * (L+1)=36
    double grly[49][3];
    int index_int = index * 2;
    int it = input_int[index_int];
    int dist_tmp = input_int[index_int + 1];

    /***************************
    L = 0
    ***************************/
    ylma[0] = ylmcoef[0]; // l=0, m=0
    grly[0][0] = grly[0][1] = grly[0][2] = 0.0;
    int nwl = ucell_atom_nwl[it];
    double tmp0;
    if (nwl == 0)
      goto YLM_END;

    /***************************
    L = 1
    ***************************/
    ylma[1] = ylmcoef[1] * dr[2]; // l=1, m=0
    grly[1][0] = grly[1][1] = 0.0;
    grly[1][2] = ylmcoef[1];
    ylma[2] = -ylmcoef[1] * dr[0]; // l=1, m=1
    grly[2][1] = grly[2][2] = 0.0;
    grly[2][0] = -ylmcoef[1];
    ylma[3] = -ylmcoef[1] * dr[1]; // l=1, m=-1
    grly[3][0] = grly[3][2] = 0.0;
    grly[3][1] = -ylmcoef[1];
    if (nwl == 1)
      goto YLM_END;

    /***************************
    L = 2
    ***************************/
    ylma[4] = ylmcoef[2] * dr[2] * ylma[1] -
              ylmcoef[3] * ylma[0] * distance; // l=2, m=0
    grly[4][0] =
        ylmcoef[2] * dr[2] * grly[1][0] -
        ylmcoef[3] * (grly[0][0] * distance + ylma[0] * tx); // l=2, m=0
    grly[4][1] =
        ylmcoef[2] * dr[2] * grly[1][1] -
        ylmcoef[3] * (grly[0][1] * distance + ylma[0] * ty); // l=2, m=0
    grly[4][2] =
        ylmcoef[2] * (dr[2] * grly[1][2] + ylma[1]) -
        ylmcoef[3] * (grly[0][2] * distance + ylma[0] * tz); // l=2, m=0

    tmp0 = ylmcoef[4] * dr[2];
    ylma[5] = tmp0 * ylma[2]; // l=2,m=1
    grly[5][0] = tmp0 * grly[2][0];
    grly[5][1] = tmp0 * grly[2][1];
    grly[5][2] = ylmcoef[4] * (ylma[2] + dr[2] * grly[2][2]);

    ylma[6] = tmp0 * ylma[3]; // l=2,m=-1
    grly[6][0] = tmp0 * grly[3][0];
    grly[6][1] = tmp0 * grly[3][1];
    grly[6][2] = ylmcoef[4] * (ylma[3] + dr[2] * grly[3][2]);

    tmp0 = ylmcoef[4] * dr[0];
    ylma[7] = ylmcoef[5] * ylma[4] - ylmcoef[6] * ylma[0] * distance -
              tmp0 * ylma[2]; // l=2,m=2
    grly[7][0] = ylmcoef[5] * grly[4][0] -
                 ylmcoef[6] * (ylma[0] * tx + grly[0][0] * distance) -
                 ylmcoef[4] * (dr[0] * grly[2][0] + ylma[2]);
    grly[7][1] = ylmcoef[5] * grly[4][1] -
                 ylmcoef[6] * (ylma[0] * ty + grly[0][1] * distance) -
                 tmp0 * grly[2][1];
    grly[7][2] = ylmcoef[5] * grly[4][2] -
                 ylmcoef[6] * (ylma[0] * tz + grly[0][2] * distance) -
                 tmp0 * grly[2][2];

    ylma[8] = -tmp0 * ylma[3];
    grly[8][0] = -ylmcoef[4] * (ylma[3] + dr[0] * grly[3][0]);
    grly[8][1] = -tmp0 * grly[3][1];
    grly[8][2] = -tmp0 * grly[3][2];
    if (nwl == 2)
      goto YLM_END;

    /***************************
    L = 3
    ***************************/
    ylma[9] = ylmcoef[7] * dr[2] * ylma[4] -
              ylmcoef[8] * ylma[1] * distance; // l=3, m=0
    grly[9][0] = ylmcoef[7] * dr[2] * grly[4][0] -
                 ylmcoef[8] * (ylma[1] * tx + grly[1][0] * distance);
    grly[9][1] = ylmcoef[7] * dr[2] * grly[4][1] -
                 ylmcoef[8] * (ylma[1] * ty + grly[1][1] * distance);
    grly[9][2] = ylmcoef[7] * (ylma[4] + dr[2] * grly[4][2]) -
                 ylmcoef[8] * (ylma[1] * tz + grly[1][2] * distance);

    tmp0 = ylmcoef[9] * dr[2];
    ylma[10] = tmp0 * ylma[5] - ylmcoef[10] * ylma[2] * distance; // l=3,m=1
    grly[10][0] = tmp0 * grly[5][0] -
                  ylmcoef[10] * (grly[2][0] * distance + ylma[2] * tx);
    grly[10][1] = tmp0 * grly[5][1] -
                  ylmcoef[10] * (grly[2][1] * distance + ylma[2] * ty);
    grly[10][2] = ylmcoef[9] * (dr[2] * grly[5][2] + ylma[5]) -
                  ylmcoef[10] * (grly[2][2] * distance + ylma[2] * tz);

    ylma[11] = tmp0 * ylma[6] - ylmcoef[10] * ylma[3] * distance; // l=3,m=-1
    grly[11][0] = tmp0 * grly[6][0] -
                  ylmcoef[10] * (grly[3][0] * distance + ylma[3] * tx);
    grly[11][1] = tmp0 * grly[6][1] -
                  ylmcoef[10] * (grly[3][1] * distance + ylma[3] * ty);
    grly[11][2] = ylmcoef[9] * (dr[2] * grly[6][2] + ylma[6]) -
                  ylmcoef[10] * (grly[3][2] * distance + ylma[3] * tz);

    tmp0 = ylmcoef[11] * dr[2];
    ylma[12] = tmp0 * ylma[7]; // l=3,m=2
    grly[12][0] = tmp0 * grly[7][0];
    grly[12][1] = tmp0 * grly[7][1];
    grly[12][2] = ylmcoef[11] * (dr[2] * grly[7][2] + ylma[7]);

    ylma[13] = tmp0 * ylma[8]; // l=3,m=-2
    grly[13][0] = tmp0 * grly[8][0];
    grly[13][1] = tmp0 * grly[8][1];
    grly[13][2] = ylmcoef[11] * (dr[2] * grly[8][2] + ylma[8]);

    tmp0 = ylmcoef[14] * dr[0];
    ylma[14] = ylmcoef[12] * ylma[10] - ylmcoef[13] * ylma[2] * distance -
               tmp0 * ylma[7]; // l=3,m=3
    grly[14][0] = ylmcoef[12] * grly[10][0] -
                  ylmcoef[13] * (ylma[2] * tx + grly[2][0] * distance) -
                  ylmcoef[14] * (ylma[7] + dr[0] * grly[7][0]);
    grly[14][1] = ylmcoef[12] * grly[10][1] -
                  ylmcoef[13] * (ylma[2] * ty + grly[2][1] * distance) -
                  tmp0 * grly[7][1];
    grly[14][2] = ylmcoef[12] * grly[10][2] -
                  ylmcoef[13] * (ylma[2] * tz + grly[2][2] * distance) -
                  tmp0 * grly[7][2];

    ylma[15] = ylmcoef[12] * ylma[11] - ylmcoef[13] * ylma[3] * distance -
               tmp0 * ylma[8]; // l=3,m=-3
    grly[15][0] = ylmcoef[12] * grly[11][0] -
                  ylmcoef[13] * (ylma[3] * tx + grly[3][0] * distance) -
                  ylmcoef[14] * (ylma[8] + dr[0] * grly[8][0]);
    grly[15][1] = ylmcoef[12] * grly[11][1] -
                  ylmcoef[13] * (ylma[3] * ty + grly[3][1] * distance) -
                  tmp0 * grly[8][1];
    grly[15][2] = ylmcoef[12] * grly[11][2] -
                  ylmcoef[13] * (ylma[3] * tz + grly[3][2] * distance) -
                  tmp0 * grly[8][2];
    if (nwl == 3)
      goto YLM_END;

    /***************************
    L = 4
    ***************************/
    ylma[16] = ylmcoef[15] * dr[2] * ylma[9] -
               ylmcoef[16] * ylma[4] * distance; // l=4,m=0
    grly[16][0] = ylmcoef[15] * dr[2] * grly[9][0] -
                  ylmcoef[16] * (ylma[4] * tx + grly[4][0] * distance);
    grly[16][1] = ylmcoef[15] * dr[2] * grly[9][1] -
                  ylmcoef[16] * (ylma[4] * ty + grly[4][1] * distance);
    grly[16][2] = ylmcoef[15] * (dr[2] * grly[9][2] + ylma[9]) -
                  ylmcoef[16] * (ylma[4] * tz + grly[4][2] * distance);

    tmp0 = ylmcoef[17] * dr[2];
    ylma[17] = tmp0 * ylma[10] - ylmcoef[18] * ylma[5] * distance; // l=4,m=1
    grly[17][0] = tmp0 * grly[10][0] -
                  ylmcoef[18] * (ylma[5] * tx + grly[5][0] * distance);
    grly[17][1] = tmp0 * grly[10][1] -
                  ylmcoef[18] * (ylma[5] * ty + grly[5][1] * distance);
    grly[17][2] = ylmcoef[17] * (dr[2] * grly[10][2] + ylma[10]) -
                  ylmcoef[18] * (ylma[5] * tz + grly[5][2] * distance);

    ylma[18] = tmp0 * ylma[11] - ylmcoef[18] * ylma[6] * distance; // l=4,m=-1
    grly[18][0] = tmp0 * grly[11][0] -
                  ylmcoef[18] * (ylma[6] * tx + grly[6][0] * distance);
    grly[18][1] = tmp0 * grly[11][1] -
                  ylmcoef[18] * (ylma[6] * ty + grly[6][1] * distance);
    grly[18][2] = ylmcoef[17] * (dr[2] * grly[11][2] + ylma[11]) -
                  ylmcoef[18] * (ylma[6] * tz + grly[6][2] * distance);

    tmp0 = ylmcoef[19] * dr[2];
    ylma[19] = tmp0 * ylma[12] - ylmcoef[20] * ylma[7] * distance; // l=4,m=2
    grly[19][0] = tmp0 * grly[12][0] -
                  ylmcoef[20] * (ylma[7] * tx + grly[7][0] * distance);
    grly[19][1] = tmp0 * grly[12][1] -
                  ylmcoef[20] * (ylma[7] * ty + grly[7][1] * distance);
    grly[19][2] = ylmcoef[19] * (dr[2] * grly[12][2] + ylma[12]) -
                  ylmcoef[20] * (ylma[7] * tz + grly[7][2] * distance);

    ylma[20] = tmp0 * ylma[13] - ylmcoef[20] * ylma[8] * distance; // l=4,m=-2
    grly[20][0] = tmp0 * grly[13][0] -
                  ylmcoef[20] * (ylma[8] * tx + grly[8][0] * distance);
    grly[20][1] = tmp0 * grly[13][1] -
                  ylmcoef[20] * (ylma[8] * ty + grly[8][1] * distance);
    grly[20][2] = ylmcoef[19] * (dr[2] * grly[13][2] + ylma[13]) -
                  ylmcoef[20] * (ylma[8] * tz + grly[8][2] * distance);

    tmp0 = 3.0 * dr[2];
    ylma[21] = tmp0 * ylma[14]; // l=4,m=3
    grly[21][0] = tmp0 * grly[14][0];
    grly[21][1] = tmp0 * grly[14][1];
    grly[21][2] = 3.0 * (dr[2] * grly[14][2] + ylma[14]);

    ylma[22] = tmp0 * ylma[15]; // l=4,m=-3
    grly[22][0] = tmp0 * grly[15][0];
    grly[22][1] = tmp0 * grly[15][1];
    grly[22][2] = 3.0 * (dr[2] * grly[15][2] + ylma[15]);

    tmp0 = ylmcoef[23] * dr[0];
    ylma[23] = ylmcoef[21] * ylma[19] - ylmcoef[22] * ylma[7] * distance -
               tmp0 * ylma[14]; // l=4,m=4
    grly[23][0] = ylmcoef[21] * grly[19][0] -
                  ylmcoef[22] * (ylma[7] * tx + grly[7][0] * distance) -
                  ylmcoef[23] * (dr[0] * grly[14][0] + ylma[14]);
    grly[23][1] = ylmcoef[21] * grly[19][1] -
                  ylmcoef[22] * (ylma[7] * ty + grly[7][1] * distance) -
                  tmp0 * grly[14][1];
    grly[23][2] = ylmcoef[21] * grly[19][2] -
                  ylmcoef[22] * (ylma[7] * tz + grly[7][2] * distance) -
                  tmp0 * grly[14][2];

    ylma[24] = ylmcoef[21] * ylma[20] - ylmcoef[22] * ylma[8] * distance -
               tmp0 * ylma[15]; // l=4,m=-4
    grly[24][0] = ylmcoef[21] * grly[20][0] -
                  ylmcoef[22] * (ylma[8] * tx + grly[8][0] * distance) -
                  ylmcoef[23] * (dr[0] * grly[15][0] + ylma[15]);
    grly[24][1] = ylmcoef[21] * grly[20][1] -
                  ylmcoef[22] * (ylma[8] * ty + grly[8][1] * distance) -
                  tmp0 * grly[15][1];
    grly[24][2] = ylmcoef[21] * grly[20][2] -
                  ylmcoef[22] * (ylma[8] * tz + grly[8][2] * distance) -
                  tmp0 * grly[15][2];
    if (nwl == 4)
      goto YLM_END;

    /***************************
    L = 5
    ***************************/
    ylma[25] = ylmcoef[24] * dr[2] * ylma[16] -
               ylmcoef[25] * ylma[9] * distance; // l=5,m=0
    grly[25][0] = ylmcoef[24] * dr[2] * grly[16][0] -
                  ylmcoef[25] * (ylma[9] * tx + grly[9][0] * distance);
    grly[25][1] = ylmcoef[24] * dr[2] * grly[16][1] -
                  ylmcoef[25] * (ylma[9] * ty + grly[9][1] * distance);
    grly[25][2] = ylmcoef[24] * (dr[2] * grly[16][2] + ylma[16]) -
                  ylmcoef[25] * (ylma[9] * tz + grly[9][2] * distance);

    tmp0 = ylmcoef[26] * dr[2];
    ylma[26] = tmp0 * ylma[17] - ylmcoef[27] * ylma[10] * distance; // l=5,m=1
    grly[26][0] = tmp0 * grly[17][0] -
                  ylmcoef[27] * (ylma[10] * tx + grly[10][0] * distance);
    grly[26][1] = tmp0 * grly[17][1] -
                  ylmcoef[27] * (ylma[10] * ty + grly[10][1] * distance);
    grly[26][2] = ylmcoef[26] * (dr[2] * grly[17][2] + ylma[17]) -
                  ylmcoef[27] * (ylma[10] * tz + grly[10][2] * distance);

    ylma[27] = tmp0 * ylma[18] - ylmcoef[27] * ylma[11] * distance; // l=5,m=-1
    grly[27][0] = tmp0 * grly[18][0] -
                  ylmcoef[27] * (ylma[11] * tx + grly[11][0] * distance);
    grly[27][1] = tmp0 * grly[18][1] -
                  ylmcoef[27] * (ylma[11] * ty + grly[11][1] * distance);
    grly[27][2] = ylmcoef[26] * (dr[2] * grly[18][2] + ylma[18]) -
                  ylmcoef[27] * (ylma[11] * tz + grly[11][2] * distance);

    tmp0 = ylmcoef[28] * dr[2];
    ylma[28] = tmp0 * ylma[19] - ylmcoef[29] * ylma[12] * distance; // l=5,m=2
    grly[28][0] = tmp0 * grly[19][0] -
                  ylmcoef[29] * (ylma[12] * tx + grly[12][0] * distance);
    grly[28][1] = tmp0 * grly[19][1] -
                  ylmcoef[29] * (ylma[12] * ty + grly[12][1] * distance);
    grly[28][2] = ylmcoef[28] * (dr[2] * grly[19][2] + ylma[19]) -
                  ylmcoef[29] * (ylma[12] * tz + grly[12][2] * distance);

    ylma[29] = tmp0 * ylma[20] - ylmcoef[29] * ylma[13] * distance; // l=5,m=-2
    grly[29][0] = tmp0 * grly[20][0] -
                  ylmcoef[29] * (ylma[13] * tx + grly[13][0] * distance);
    grly[29][1] = tmp0 * grly[20][1] -
                  ylmcoef[29] * (ylma[13] * ty + grly[13][1] * distance);
    grly[29][2] = ylmcoef[28] * (dr[2] * grly[20][2] + ylma[20]) -
                  ylmcoef[29] * (ylma[13] * tz + grly[13][2] * distance);

    tmp0 = ylmcoef[30] * dr[2];
    ylma[30] = tmp0 * ylma[21] - ylmcoef[31] * ylma[14] * distance; // l=5,m=3
    grly[30][0] = tmp0 * grly[21][0] -
                  ylmcoef[31] * (grly[14][0] * distance + ylma[14] * tx);
    grly[30][1] = tmp0 * grly[21][1] -
                  ylmcoef[31] * (grly[14][1] * distance + ylma[14] * ty);
    grly[30][2] = ylmcoef[30] * (dr[2] * grly[21][2] + ylma[21]) -
                  ylmcoef[31] * (ylma[14] * tz + grly[14][2] * distance);

    ylma[31] = tmp0 * ylma[22] - ylmcoef[31] * ylma[15] * distance; // l=5,m=-3
    grly[31][0] = tmp0 * grly[22][0] -
                  ylmcoef[31] * (grly[15][0] * distance + ylma[15] * tx);
    grly[31][1] = tmp0 * grly[22][1] -
                  ylmcoef[31] * (grly[15][1] * distance + ylma[15] * ty);
    grly[31][2] = ylmcoef[30] * (dr[2] * grly[22][2] + ylma[22]) -
                  ylmcoef[31] * (ylma[15] * tz + grly[15][2] * distance);

    tmp0 = ylmcoef[32] * dr[2];
    ylma[32] = tmp0 * ylma[23]; // l=5,m=4
    grly[32][0] = tmp0 * grly[23][0];
    grly[32][1] = tmp0 * grly[23][1];
    grly[32][2] = ylmcoef[32] * (ylma[23] + dr[2] * grly[23][2]);

    ylma[33] = tmp0 * ylma[24]; // l=5,m=-4
    grly[33][0] = tmp0 * grly[24][0];
    grly[33][1] = tmp0 * grly[24][1];
    grly[33][2] = ylmcoef[32] * (ylma[24] + dr[2] * grly[24][2]);

    tmp0 = ylmcoef[35] * dr[0];
    ylma[34] = ylmcoef[33] * ylma[30] - ylmcoef[34] * ylma[14] * distance -
               tmp0 * ylma[23]; // l=5,m=5
    grly[34][0] = ylmcoef[33] * grly[30][0] -
                  ylmcoef[34] * (ylma[14] * tx + grly[14][0] * distance) -
                  ylmcoef[35] * (dr[0] * grly[23][0] + ylma[23]);
    grly[34][1] = ylmcoef[33] * grly[30][1] -
                  ylmcoef[34] * (ylma[14] * ty + grly[14][1] * distance) -
                  tmp0 * grly[23][1];
    grly[34][2] = ylmcoef[33] * grly[30][2] -
                  ylmcoef[34] * (ylma[14] * tz + grly[14][2] * distance) -
                  tmp0 * grly[23][2];

    ylma[35] = ylmcoef[33] * ylma[31] - ylmcoef[34] * ylma[15] * distance -
               tmp0 * ylma[24]; // l=5,m=-5
    grly[35][0] = ylmcoef[33] * grly[31][0] -
                  ylmcoef[34] * (ylma[15] * tx + grly[15][0] * distance) -
                  ylmcoef[35] * (dr[0] * grly[24][0] + ylma[24]);
    grly[35][1] = ylmcoef[33] * grly[31][1] -
                  ylmcoef[34] * (ylma[15] * ty + grly[15][1] * distance) -
                  tmp0 * grly[24][1];
    grly[35][2] = ylmcoef[33] * grly[31][2] -
                  ylmcoef[34] * (ylma[15] * tz + grly[15][2] * distance) -
                  tmp0 * grly[24][2];

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
            ylma[istart + im] = fac2 / sqrt((double)istart - imm * imm) * (dr[2]
    * ylma[istart1 + im] - sqrt((double)istart1 - imm * imm) / fac4 *
    ylma[istart2 + im]);
        }

        double bl1 = sqrt(2.0 * il / (2.0 * il + 1.0));
        double bl2 = sqrt((2.0 * il - 2.0) / (2.0 * il - 1.0));
        double bl3 = sqrt(2.0) / fac2;

        ylma[istart + 2 * il - 1] = (bl3 * ylma[istart + 2 * il - 5] - bl2 *
    ylma[istart2 + 2 * il - 5] - 2.0 * dr[0] * ylma[istart1 + 2 * il - 3]) /
    bl1; ylma[istart + 2 * il] = (bl3 * ylma[istart + 2 * il - 4] - bl2 *
    ylma[istart2 + 2 * il - 4] - 2.0 * dr[0] * ylma[istart1 + 2 * il - 2]) /
    bl1;
    }*/
  YLM_END:
    // Calculate normalized position for interpolation
    distance = sqrt(distance);
    const double postion = distance / delta_r_g;
    // Extract integer part and fractional part of the position
    const double ip = static_cast<int>(postion);
    const double x0 = postion - ip;
    const double x1 = 1.0 - x0;
    const double x2 = 2.0 - x0;
    const double x3 = 3.0 - x0;
    const double x12 = x1 * x2 / 6;
    const double x03 = x0 * x3 / 2;
    // Temporary variables for interpolation
    double tmp, dtmp;

    double phi = 0.0;
    // Loop over non-zero elements in atom_nw array
    int it_nw = it * nwmax_g;
    int iw_nr = (it_nw * nr_max + ip) * 2;
    int it_nw_iw = it_nw;
    for (int iw = 0; iw < atom_nw[it]; ++iw) {
      if (atom_iw2_new[it_nw_iw]) {
        // Perform interpolation using cubic B-spline basis functions
        tmp = x12 * (psi_u[iw_nr] * x3 + psi_u[iw_nr + 6] * x0) +
              x03 * (psi_u[iw_nr + 2] * x2 - psi_u[iw_nr + 4] * x1);
        dtmp = x12 * (psi_u[iw_nr + 1] * x3 + psi_u[iw_nr + 7] * x0) +
               x03 * (psi_u[iw_nr + 3] * x2 - psi_u[iw_nr + 5] * x1);
      }
      // Extract information from atom_iw2_* arrays
      const int ll = atom_iw2_l[it_nw_iw];

      const int idx_lm = atom_iw2_ylm[it_nw_iw];

      const double rl = pow(distance, ll);

      // Compute right-hand side of the equation
      psir_ylm_right[dist_tmp] = tmp * ylma[idx_lm] / rl * vlbr3_value;
      // Compute derivatives with respect to spatial coordinates
      const double tmpdphi_rly =
          (dtmp - tmp * ll / distance) / rl * ylma[idx_lm] / distance;
      const double tmprl = tmp / rl;
      dpsir_ylm_left_x[dist_tmp] =
          tmpdphi_rly * dr[0] + tmprl * grly[idx_lm][0];

      dpsir_ylm_left_y[dist_tmp] =
          tmpdphi_rly * dr[1] + tmprl * grly[idx_lm][1];
      dpsir_ylm_left_z[dist_tmp] =
          tmpdphi_rly * dr[2] + tmprl * grly[idx_lm][2];

      ddpsir_ylm_left_xx[dist_tmp] = dpsir_ylm_left_x[dist_tmp] * dr[0];
      ddpsir_ylm_left_xy[dist_tmp] = dpsir_ylm_left_x[dist_tmp] * dr[1];
      ddpsir_ylm_left_xz[dist_tmp] = dpsir_ylm_left_x[dist_tmp] * dr[2];
      ddpsir_ylm_left_yy[dist_tmp] = dpsir_ylm_left_y[dist_tmp] * dr[1];
      ddpsir_ylm_left_yz[dist_tmp] = dpsir_ylm_left_y[dist_tmp] * dr[2];
      ddpsir_ylm_left_zz[dist_tmp] = dpsir_ylm_left_z[dist_tmp] * dr[2];

      // Update loop counters and indices
      dist_tmp += 1;
      iw_nr += nr_max;
      iw_nr += nr_max;
      it_nw_iw++;
    }
  }
}

__global__ void psir_dot_stress(int * n,
                        double **x_array_g,
                        int incx,
                        double **y_array_g,
                        int incy,
                        double **results_g,
                        int batchcount)
{
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for(int i = id; i < batchcount; i += stride){
        double *sum = results_g[i];
        double *x = x_array_g[i];
        double *y = y_array_g[i];

        for(int j = 0; j < n[i]; j++){
            sum[0] += x[j*incx] * y[j*incy];
        }
    }
}

__global__ void dot_product_stress(double *ddpsir_ylm_left_xx,
                                    double *ddpsir_ylm_left_xy,
                                    double *ddpsir_ylm_left_xz,
                                    double *ddpsir_ylm_left_yy,
                                    double *ddpsir_ylm_left_yz,
                                    double *ddpsir_ylm_left_zz,
                                    double *psir_ylm_dm,
                                    double *stress_dot,
                                    int elements_num)
{
  
	__shared__ double cache[256][6]; // == threadsPerBlock
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int cacheIndex = threadIdx.x;
	double tmp[6]= {0.0};
	while (tid < elements_num) {
		tmp[0] += ddpsir_ylm_left_xx[tid] * psir_ylm_dm[tid] * 2;
    tmp[1] += ddpsir_ylm_left_xy[tid] * psir_ylm_dm[tid] * 2;
    tmp[2] += ddpsir_ylm_left_xz[tid] * psir_ylm_dm[tid] * 2;
    tmp[3] += ddpsir_ylm_left_yy[tid] * psir_ylm_dm[tid] * 2;
    tmp[4] += ddpsir_ylm_left_yz[tid] * psir_ylm_dm[tid] * 2;
    tmp[5] += ddpsir_ylm_left_zz[tid] * psir_ylm_dm[tid] * 2;
		tid += blockDim.x * gridDim.x;
	}

  for (int i=0;i<6;i++)
	  cache[cacheIndex][i] = tmp[i];

	__syncthreads();

	int i = blockDim.x / 2;
	while (i != 0) {
		if (cacheIndex < i)
      for (int index=0;index<6;index++)
			  cache[cacheIndex][index] += cache[cacheIndex + i][index];
		__syncthreads();
		i /= 2;
	}
  
	if (cacheIndex == 0)
    for (int index=0;index<6;index++)
		  stress_dot[blockIdx.x+gridDim.x*index] = cache[0][index];
}

__global__ void dot_product_force(double *dpsir_ylm_left_x,
                                  double *dpsir_ylm_left_y,
                                  double *dpsir_ylm_left_z,
                                  double *psir_ylm_dm,
                                  double *force_dot,
                                  int *iat,
                                  int nwmax,
                                  int max_size,
                                  int elements_num)
{
  int tid=threadIdx.x + blockIdx.x * blockDim.x;
  while (tid<elements_num){
  int iat_on_nbz=iat[tid];
  if (iat_on_nbz<=-1) 
  {
    tid += blockDim.x * gridDim.x;
    continue;
  }

  int iat_index =tid*3;
  int dist =tid*nwmax;
  double tmp[3]={0.0};
  
  for (int i=0;i<nwmax;i++)
  {
    tmp[0]+=dpsir_ylm_left_x[dist+i]*psir_ylm_dm[dist+i]*2;
    tmp[1]+=dpsir_ylm_left_y[dist+i]*psir_ylm_dm[dist+i]*2;
    tmp[2]+=dpsir_ylm_left_z[dist+i]*psir_ylm_dm[dist+i]*2;
  }

  for (int i=0;i<3;i++)
     force_dot[iat_index+i]=tmp[i];
  tid += blockDim.x * gridDim.x;
  }
  
}
