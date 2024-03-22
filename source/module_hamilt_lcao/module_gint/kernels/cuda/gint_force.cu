#include "module_hamilt_lcao/module_gint/kernels/cuda/gint_force.cuh"
#include "module_hamilt_lcao/module_gint/kernels/cuda/gint_force.h"
#include "module_hamilt_lcao/module_gint/kernels/cuda/vbatch_matrix_multiple/cuda_tools.cuh"
#include "module_hamilt_lcao/module_gint/kernels/cuda/sph.cuh"
// CUDA kernel to calculate psi and force
namespace lcaoCudaKernel{

/*!
 * \file
 * \brief CUDA kernel to calculate psi and force
 *
 * CUDA kernel that performs calculations on psi and force.
 *
 * \param ylmcoef Pointer to the Ylm coefficients
 * \param delta_r_g Delta r value
 * \param bxyz_g Bxyz value
 * \param nwmax_g Nwmax value
 * \param input_double Array of double input values
 * \param input_int Array of int input values
 * \param num_psir Array containing the number of psi for each block
 * \param psi_size_max Maximum size of psi
 * \param ucell_atom_nwl Array containing Ucell atom nwl values
 * \param atom_iw2_new Array indicating whether atom_iw2 is new
 * \param atom_iw2_ylm Array of atom_iw2 Ylm values
 * \param atom_iw2_l Array of atom_iw2 l values
 * \param atom_nw Array of atom_nw values
 * \param nr_max Maximum nr value
 * \param psi_u Array for psi_u values
 * \param psir_ylm_right Array for psir_ylm_right values
 * \param dpsir_ylm_left_x Array for dpsir_ylm_left_x values
 * \param dpsir_ylm_left_y Array for dpsir_ylm_left_y values
 * \param dpsir_ylm_left_z Array for dpsir_ylm_left_z values
 * \param ddpsir_ylm_left_xx Array for ddpsir_ylm_left_xx values
 * \param ddpsir_ylm_left_xy Array for ddpsir_ylm_left_xy values
 * \param ddpsir_ylm_left_xz Array for ddpsir_ylm_left_xz values
 * \param ddpsir_ylm_left_yy Array for ddpsir_ylm_left_yy values
 * \param ddpsir_ylm_left_yz Array for ddpsir_ylm_left_yz values
 * \param ddpsir_ylm_left_zz Array for ddpsir_ylm_left_zz values
 */

__global__ void
get_psi_force(double *ylmcoef, double delta_r_g, double bxyz_g, double nwmax_g,
              double *input_double, int *input_int, int *num_psir,
              int psi_size_max, int *ucell_atom_nwl, bool *atom_iw2_new,
              int *atom_iw2_ylm, int *atom_iw2_l, int *atom_nw, int nr_max,
              double *psi_u, double *psir_ylm_right, double *dpsir_ylm_left_x,
              double *dpsir_ylm_left_y, double *dpsir_ylm_left_z,
              double *ddpsir_ylm_left_xx, double *ddpsir_ylm_left_xy,
              double *ddpsir_ylm_left_xz, double *ddpsir_ylm_left_yy,
              double *ddpsir_ylm_left_yz, double *ddpsir_ylm_left_zz) {
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
    double vlbr3_value = input_double[index_double + 4];
    // begin calculation
    double ylma[49]; // Attention!!! At present, we only use L=5 at most. So
                     // (L+1) * (L+1)=36
    double grly[49][3];
    int index_int = index * 2;
    int it = input_int[index_int];
    int dist_tmp = input_int[index_int + 1];

    int nwl = ucell_atom_nwl[it];
    spherical_harmonics_d(dr,distance,grly,nwl,ylma,ylmcoef);
  
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

__global__ void psir_dot_stress(int *n, double **x_array_g, int incx,
                                double **y_array_g, int incy,
                                double **results_g, int batchcount) {
  int id = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int i = id; i < batchcount; i += stride) {
    double *sum = results_g[i];
    double *x = x_array_g[i];
    double *y = y_array_g[i];

    for (int j = 0; j < n[i]; j++) {
      sum[0] += x[j * incx] * y[j * incy];
    }
  }
}

/**
 * \brief Compute dot product of stress components and partial derivatives.
 *
 * This CUDA kernel computes the dot product of stress components and partial
 * derivatives based on the input arrays.
 *
 * \param ddpsir_ylm_left_xx Array of ddpsir_ylm_left_xx values.
 * \param ddpsir_ylm_left_xy Array of ddpsir_ylm_left_xy values.
 * \param ddpsir_ylm_left_xz Array of ddpsir_ylm_left_xz values.
 * \param ddpsir_ylm_left_yy Array of ddpsir_ylm_left_yy values.
 * \param ddpsir_ylm_left_yz Array of ddpsir_ylm_left_yz values.
 * \param ddpsir_ylm_left_zz Array of ddpsir_ylm_left_zz values.
 * \param psir_ylm_dm Array of psir_ylm_dm values.
 * \param stress_dot Output array for the dot product of stress components.
 * \param elements_num Number of elements in the input arrays.
 */


__global__ void
dot_product_stress(double *ddpsir_ylm_left_xx, double *ddpsir_ylm_left_xy,
                   double *ddpsir_ylm_left_xz, double *ddpsir_ylm_left_yy,
                   double *ddpsir_ylm_left_yz, double *ddpsir_ylm_left_zz,
                   double *psir_ylm_dm, double *stress_dot, int elements_num) {

  __shared__ double cache[256][6]; // == threadsPerBlock
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int cacheIndex = threadIdx.x;
  double tmp[6] = {0.0};
  while (tid < elements_num) {
    tmp[0] += ddpsir_ylm_left_xx[tid] * psir_ylm_dm[tid] * 2;
    tmp[1] += ddpsir_ylm_left_xy[tid] * psir_ylm_dm[tid] * 2;
    tmp[2] += ddpsir_ylm_left_xz[tid] * psir_ylm_dm[tid] * 2;
    tmp[3] += ddpsir_ylm_left_yy[tid] * psir_ylm_dm[tid] * 2;
    tmp[4] += ddpsir_ylm_left_yz[tid] * psir_ylm_dm[tid] * 2;
    tmp[5] += ddpsir_ylm_left_zz[tid] * psir_ylm_dm[tid] * 2;
    tid += blockDim.x * gridDim.x;
  }

  for (int i = 0; i < 6; i++)
    cache[cacheIndex][i] = tmp[i];

  __syncthreads();

  int i = blockDim.x / 2;
  while (i != 0) {
    if (cacheIndex < i)
      for (int index = 0; index < 6; index++)
        cache[cacheIndex][index] += cache[cacheIndex + i][index];
    __syncthreads();
    i /= 2;
  }

  if (cacheIndex == 0)
    for (int index = 0; index < 6; index++)
      stress_dot[blockIdx.x + gridDim.x * index] = cache[0][index];
}

/**
 * @brief Calculate the dot product force.
 *
 * This function calculates the dot product force based on the provided parameters.
 *
 * @param dpsir_ylm_left_x Pointer to the array of dpsir_ylm_left_x values.
 * @param dpsir_ylm_left_y Pointer to the array of dpsir_ylm_left_y values.
 * @param dpsir_ylm_left_z Pointer to the array of dpsir_ylm_left_z values.
 * @param psir_ylm_dm Pointer to the array of psir_ylm_dm values.
 * @param force_dot Pointer to the array where the calculated force will be stored.
 * @param iat Pointer to the array of iat values.
 * @param nwmax Maximum value for nwmax.
 * @param max_size Maximum size for arrays.
 * @param elements_num Number of elements to process.
 */

__global__ void dot_product_force(double *dpsir_ylm_left_x,
                                  double *dpsir_ylm_left_y,
                                  double *dpsir_ylm_left_z, double *psir_ylm_dm,
                                  double *force_dot, int *iat, int nwmax,
                                  int max_size, int elements_num) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < elements_num) {
    int iat_on_nbz = iat[tid];
    if (iat_on_nbz <= -1) {
      tid += blockDim.x * gridDim.x;
      continue;
    }

    int iat_index = tid * 3;
    int dist = tid * nwmax;
    double tmp[3] = {0.0};

    for (int i = 0; i < nwmax; i++) {
      tmp[0] += dpsir_ylm_left_x[dist + i] * psir_ylm_dm[dist + i] * 2;
      tmp[1] += dpsir_ylm_left_y[dist + i] * psir_ylm_dm[dist + i] * 2;
      tmp[2] += dpsir_ylm_left_z[dist + i] * psir_ylm_dm[dist + i] * 2;
    }

    for (int i = 0; i < 3; i++)
      force_dot[iat_index + i] = tmp[i];
    tid += blockDim.x * gridDim.x;
  }
}

} // namespace lcaoCudaKernel
