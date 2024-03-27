#include <fstream>
#include <sstream>
#include <omp.h>

#include "kernels/cuda/gint_force.cuh"
#include "gint_force.h"
#include "module_base/ylm.h"
#include "module_hamilt_lcao/module_gint/gint_tools.h"
#include "kernels/cuda/vbatch_matrix_multiple/cuda_tools.cuh"

namespace GintKernel{

// Function to calculate forces using GPU-accelerated gamma point Gint
/**
 * @brief Calculate forces and stresses for the `gint_gamma_force_gpu` function.
 *
 * This function calculates forces and stresses based on given parameters.
 *
 * @param dm A pointer to the HContainer<double> object.
 * @param vfactor The scaling factor for some calculation.
 * @param vlocal A pointer to an array of doubles.
 * @param force A pointer to an array to store the calculated forces.
 * @param stress A pointer to an array to store the calculated stresses.
 * @param nczp An integer representing a parameter.
 * @param ylmcoef_now A pointer to an array of doubles representing Ylm coefficients.
 * @param gridt A reference to a Grid_Technique object.
 */
void gint_gamma_force_gpu(hamilt::HContainer<double> *dm, const double vfactor,
                          const double *vlocal, double *force, double *stress,
                          const int nczp, const double *ylmcoef_now,
                          const double dr, const double *rcut,
                          const Grid_Technique &gridt,
                          const UnitCell &ucell) 
{
  const int nbz = gridt.nbzp;
  const int lgd = gridt.lgd;
  const int max_size = gridt.max_atom;
  const int nwmax = ucell.nwmax;
  const int bxyz = gridt.bxyz;
  /* compute the dm matrix from Hcontainer */
  double *dm_matrix_h = new double[lgd * lgd];
  ModuleBase::GlobalFunc::ZEROS(dm_matrix_h, lgd * lgd);
  for (int iat1 = 0; iat1 < ucell.nat; iat1++) {
    for (int iat2 = 0; iat2 < ucell.nat; iat2++) {
      int it1 = ucell.iat2it[iat1];
      int it2 = ucell.iat2it[iat2];
      int lo1 = gridt.trace_lo[ucell.itiaiw2iwt(
          it1, ucell.iat2ia[iat1], 0)];
      int lo2 = gridt.trace_lo[ucell.itiaiw2iwt(
          it2, ucell.iat2ia[iat2], 0)];

      hamilt::AtomPair<double> *tmp_ap = dm->find_pair(iat1, iat2);
      int orb_index = 0;
      if (tmp_ap == NULL) {
        continue;
      }
      for (int orb_i = 0; orb_i < tmp_ap->get_row_size(); orb_i++) {
        for (int orb_j = 0; orb_j < tmp_ap->get_col_size(); orb_j++) {

          dm_matrix_h[(lo1 + orb_i) * lgd + (lo2 + orb_j)] =
              tmp_ap->get_pointer(0)[orb_index];
          orb_index++;
        }
      }
    }
  }
  double *dm_matrix_g;
  checkCuda(cudaMalloc((void **)&dm_matrix_g, lgd * lgd * sizeof(double)));
  checkCuda(cudaMemcpy(dm_matrix_g, dm_matrix_h, lgd * lgd * sizeof(double),
                       cudaMemcpyHostToDevice));

  /*prepare the force and stress dot parameter */
  const int threadsPerBlock = 256;
  const int blocksPerGrid =
      std::min(64, (gridt.psir_size + threadsPerBlock - 1) / threadsPerBlock);
  const int threadsPerBlock_force = 256;
  const int blocksPerGrid_force =
      std::min(64, (nbz * bxyz * max_size + threadsPerBlock_force - 1) /
                       threadsPerBlock_force);
  double *stress_dot = new double[6 * blocksPerGrid];

  for (int i = 0; i < 6 * blocksPerGrid; i++) {
    stress_dot[i] = 0.0;
  }
  double *stress_dot_global_g;
  checkCuda(cudaMalloc((void **)&stress_dot_global_g,
                       6 * blocksPerGrid * gridt.nstreams * sizeof(double)));
  checkCuda(cudaMemset(stress_dot_global_g, 0,
                       6 * blocksPerGrid * gridt.nstreams * sizeof(double)));

  /* cudaMalloc global parameter in device*/
  double *force_dot_global_g;
  checkCuda(
      cudaMalloc((void **)&force_dot_global_g,
                 3 * nbz * bxyz * max_size * gridt.nstreams * sizeof(double)));
  checkCuda(
      cudaMemset(force_dot_global_g, 0,
                 3 * nbz * bxyz * max_size * gridt.nstreams * sizeof(double)));

  int *iat_global_g;
  checkCuda(cudaMalloc((void **)&iat_global_g,
                       nbz * bxyz * max_size * gridt.nstreams * sizeof(int)));
  checkCuda(cudaMemset(iat_global_g, 0,
                       nbz * bxyz * max_size * gridt.nstreams * sizeof(int)));

  /*cuda stream allocate */
  for (int i = 0; i < gridt.nstreams; i++) {
    checkCuda(cudaStreamSynchronize(gridt.streams[i]));
  }

  int iter_num = 0;
  for (int i = 0; i < gridt.nbx; i++) {
    for (int j = 0; j < gridt.nby; j++) {
      /* psi parameter allocate */
      int stream_num = iter_num % gridt.nstreams;
      double *psi_input_double =
          &gridt.psi_input_double_global[gridt.psi_size_max * stream_num * 5];
      int *psi_input_int =
          &gridt.psi_input_int_global[gridt.psi_size_max * stream_num * 2];
      int *num_psir = &gridt.num_psir_global[nbz * stream_num];
      int *atom_pair_A_m =
          &gridt.atom_pair_left_info_global[gridt.atom_pair_size_over_nbz *
                                            stream_num];
      int *atom_pair_B_n =
          &gridt.atom_pair_right_info_global[gridt.atom_pair_size_over_nbz *
                                             stream_num];
      int *atom_pair_k =
          &gridt.atom_pair_k_info_global[gridt.atom_pair_size_over_nbz *
                                         stream_num];
      int *atom_pair_lda =
          &gridt.atom_pair_lda_global[gridt.atom_pair_size_over_nbz *
                                      stream_num];
      int *atom_pair_ldb =
          &gridt.atom_pair_ldb_global[gridt.atom_pair_size_over_nbz *
                                      stream_num];
      int *atom_pair_ldc =
          &gridt.atom_pair_ldc_global[gridt.atom_pair_size_over_nbz *
                                      stream_num];

      double *psi_input_double_g =
          &gridt.psi_input_double_global_g[gridt.psi_size_max * stream_num * 5];
      int *psi_input_int_g =
          &gridt.psi_input_int_global_g[gridt.psi_size_max * stream_num * 2];
      int *num_psir_g = &gridt.num_psir_global_g[nbz * stream_num];
      double *psir_ylm_dm_g =
          &gridt.psir_ylm_dm_global_g[gridt.psir_size * stream_num];
      double *psir_ylm_right_g =
          &gridt.psir_ylm_right_global_g[gridt.psir_size * stream_num];
      double *dpsir_ylm_left_x_g =
          &gridt.dpsir_ylm_left_x_global_g[gridt.psir_size * stream_num];
      double *dpsir_ylm_left_y_g =
          &gridt.dpsir_ylm_left_y_global_g[gridt.psir_size * stream_num];
      double *dpsir_ylm_left_z_g =
          &gridt.dpsir_ylm_left_z_global_g[gridt.psir_size * stream_num];
      double *ddpsir_ylm_left_xx_g =
          &gridt.ddpsir_ylm_left_xx_global_g[gridt.psir_size * stream_num];
      double *ddpsir_ylm_left_xy_g =
          &gridt.ddpsir_ylm_left_xy_global_g[gridt.psir_size * stream_num];
      double *ddpsir_ylm_left_xz_g =
          &gridt.ddpsir_ylm_left_xz_global_g[gridt.psir_size * stream_num];
      double *ddpsir_ylm_left_yy_g =
          &gridt.ddpsir_ylm_left_yy_global_g[gridt.psir_size * stream_num];
      double *ddpsir_ylm_left_yz_g =
          &gridt.ddpsir_ylm_left_yz_global_g[gridt.psir_size * stream_num];
      double *ddpsir_ylm_left_zz_g =
          &gridt.ddpsir_ylm_left_zz_global_g[gridt.psir_size * stream_num];

      int *atom_pair_A_m_g =
          &gridt.atom_pair_left_info_global_g[gridt.atom_pair_size_over_nbz *
                                              stream_num];
      int *atom_pair_B_n_g =
          &gridt.atom_pair_right_info_global_g[gridt.atom_pair_size_over_nbz *
                                               stream_num];
      int *atom_pair_k_g =
          &gridt.atom_pair_k_info_global_g[gridt.atom_pair_size_over_nbz *
                                           stream_num];

      int *atom_pair_lda_g =
          &gridt.atom_pair_lda_global_g[gridt.atom_pair_size_over_nbz *
                                        stream_num];
      int *atom_pair_ldb_g =
          &gridt.atom_pair_ldb_global_g[gridt.atom_pair_size_over_nbz *
                                        stream_num];
      int *atom_pair_ldc_g =
          &gridt.atom_pair_ldc_global_g[gridt.atom_pair_size_over_nbz *
                                        stream_num];

      double **atom_pair_mat_A_array =
          &gridt.atom_pair_left_global[gridt.atom_pair_size_over_nbz *
                                       stream_num];
      double **atom_pair_mat_B_array =
          &gridt.atom_pair_right_global[gridt.atom_pair_size_over_nbz *
                                        stream_num];
      double **atom_pair_mat_C_array =
          &gridt.atom_pair_output_global[gridt.atom_pair_size_over_nbz *
                                         stream_num];

      double **atom_pair_mat_A_array_g =
          &gridt.atom_pair_left_global_g[gridt.atom_pair_size_over_nbz *
                                         stream_num];
      double **atom_pair_mat_B_array_g =
          &gridt.atom_pair_right_global_g[gridt.atom_pair_size_over_nbz *
                                          stream_num];
      double **atom_pair_mat_C_array_g =
          &gridt.atom_pair_output_global_g[gridt.atom_pair_size_over_nbz *
                                           stream_num];

      double *stress_dot_g =
          &stress_dot_global_g[6 * blocksPerGrid * stream_num];
      double *force_dot_g =
          &force_dot_global_g[3 * nbz * bxyz * max_size * stream_num];

      int *iat_g = &iat_global_g[nbz * bxyz * max_size * stream_num];
      int *iat = new int[nbz * bxyz * max_size];
      for (int index = 0; index < nbz * bxyz * max_size; index++) {
        iat[index] = -max_size - 1;
      }
      double *force_h = new double[3 * nbz * bxyz * max_size];
      for (int index = 0; index < 3 * nbz * bxyz * max_size; index++) {
        force_h[index] = 0.0;
      }

      int max_m = 0;
      int max_n = 0;
      int atom_pair_num = 0;
      checkCuda(cudaStreamSynchronize(gridt.streams[stream_num]));
      /*gpu task compute in CPU */
      gtask_force(
          gridt, rcut, ucell, i, j, gridt.psi_size_max_per_z, max_size, nczp, vfactor,
          vlocal, iat, psi_input_double, psi_input_int, num_psir, lgd,
          psir_ylm_right_g, psir_ylm_dm_g, dm_matrix_g, atom_pair_A_m,
          atom_pair_B_n, atom_pair_k, atom_pair_lda, atom_pair_ldb,
          atom_pair_ldc, atom_pair_mat_A_array, atom_pair_mat_B_array,
          atom_pair_mat_C_array, max_m, max_n, atom_pair_num);

      /*variables memcpy to gpu host*/
      checkCuda(cudaMemcpyAsync(psi_input_double_g, psi_input_double,
                                gridt.psi_size_max * 5 * sizeof(double),
                                cudaMemcpyHostToDevice,
                                gridt.streams[stream_num]));
      checkCuda(cudaMemcpyAsync(
          psi_input_int_g, psi_input_int, gridt.psi_size_max * 2 * sizeof(int),
          cudaMemcpyHostToDevice, gridt.streams[stream_num]));
      checkCuda(cudaMemcpyAsync(num_psir_g, num_psir, nbz * sizeof(int),
                                cudaMemcpyHostToDevice,
                                gridt.streams[stream_num]));
      checkCuda(cudaMemcpyAsync(iat_g, iat, nbz * bxyz * max_size * sizeof(int),
                                cudaMemcpyHostToDevice,
                                gridt.streams[stream_num]));

      checkCuda(cudaMemcpyAsync(atom_pair_A_m_g, atom_pair_A_m,
                                gridt.atom_pair_size_over_nbz * sizeof(int),
                                cudaMemcpyHostToDevice,
                                gridt.streams[stream_num]));
      checkCuda(cudaMemcpyAsync(atom_pair_B_n_g, atom_pair_B_n,
                                gridt.atom_pair_size_over_nbz * sizeof(int),
                                cudaMemcpyHostToDevice,
                                gridt.streams[stream_num]));
      checkCuda(cudaMemcpyAsync(atom_pair_k_g, atom_pair_k,
                                gridt.atom_pair_size_over_nbz * sizeof(int),
                                cudaMemcpyHostToDevice,
                                gridt.streams[stream_num]));
      checkCuda(cudaMemcpyAsync(atom_pair_lda_g, atom_pair_lda,
                                gridt.atom_pair_size_over_nbz * sizeof(int),
                                cudaMemcpyHostToDevice,
                                gridt.streams[stream_num]));
      checkCuda(cudaMemcpyAsync(atom_pair_ldb_g, atom_pair_ldb,
                                gridt.atom_pair_size_over_nbz * sizeof(int),
                                cudaMemcpyHostToDevice,
                                gridt.streams[stream_num]));
      checkCuda(cudaMemcpyAsync(atom_pair_ldc_g, atom_pair_ldc,
                                gridt.atom_pair_size_over_nbz * sizeof(int),
                                cudaMemcpyHostToDevice,
                                gridt.streams[stream_num]));

      checkCuda(
          cudaMemcpyAsync(atom_pair_mat_A_array_g, atom_pair_mat_A_array,
                          gridt.atom_pair_size_over_nbz * sizeof(double *),
                          cudaMemcpyHostToDevice, gridt.streams[stream_num]));
      checkCuda(
          cudaMemcpyAsync(atom_pair_mat_B_array_g, atom_pair_mat_B_array,
                          gridt.atom_pair_size_over_nbz * sizeof(double *),
                          cudaMemcpyHostToDevice, gridt.streams[stream_num]));
      checkCuda(
          cudaMemcpyAsync(atom_pair_mat_C_array_g, atom_pair_mat_C_array,
                          gridt.atom_pair_size_over_nbz * sizeof(double *),
                          cudaMemcpyHostToDevice, gridt.streams[stream_num]));

      checkCuda(cudaMemsetAsync(psir_ylm_dm_g, 0,
                                gridt.psir_size * sizeof(double),
                                gridt.streams[stream_num]));
      checkCuda(cudaMemsetAsync(psir_ylm_right_g, 0,
                                gridt.psir_size * sizeof(double),
                                gridt.streams[stream_num]));
      checkCuda(cudaMemsetAsync(dpsir_ylm_left_x_g, 0,
                                gridt.psir_size * sizeof(double),
                                gridt.streams[stream_num]));
      checkCuda(cudaMemsetAsync(dpsir_ylm_left_y_g, 0,
                                gridt.psir_size * sizeof(double),
                                gridt.streams[stream_num]));
      checkCuda(cudaMemsetAsync(dpsir_ylm_left_z_g, 0,
                                gridt.psir_size * sizeof(double),
                                gridt.streams[stream_num]));
      checkCuda(cudaMemsetAsync(ddpsir_ylm_left_xx_g, 0,
                                gridt.psir_size * sizeof(double),
                                gridt.streams[stream_num]));
      checkCuda(cudaMemsetAsync(ddpsir_ylm_left_xy_g, 0,
                                gridt.psir_size * sizeof(double),
                                gridt.streams[stream_num]));
      checkCuda(cudaMemsetAsync(ddpsir_ylm_left_xz_g, 0,
                                gridt.psir_size * sizeof(double),
                                gridt.streams[stream_num]));
      checkCuda(cudaMemsetAsync(ddpsir_ylm_left_yy_g, 0,
                                gridt.psir_size * sizeof(double),
                                gridt.streams[stream_num]));
      checkCuda(cudaMemsetAsync(ddpsir_ylm_left_yz_g, 0,
                                gridt.psir_size * sizeof(double),
                                gridt.streams[stream_num]));
      checkCuda(cudaMemsetAsync(ddpsir_ylm_left_zz_g, 0,
                                gridt.psir_size * sizeof(double),
                                gridt.streams[stream_num]));
      checkCuda(cudaMemsetAsync(stress_dot_g, 0,
                                6 * blocksPerGrid * sizeof(double),
                                gridt.streams[stream_num]));
      checkCuda(cudaMemsetAsync(force_dot_g, 0,
                                3 * nbz * bxyz * max_size * sizeof(double),
                                gridt.streams[stream_num]));
      dim3 grid_psi(nbz, 8);
      dim3 block_psi(64);
      /* cuda stream compute and Multiplication of multinomial matrices */
      get_psi_force<<<grid_psi, block_psi, 0, gridt.streams[stream_num]>>>(
          gridt.ylmcoef_g, dr, gridt.bxyz,
          ucell.nwmax, psi_input_double_g, psi_input_int_g, num_psir_g,
          gridt.psi_size_max_per_z, gridt.ucell_atom_nwl_g,
          gridt.atom_iw2_new_g, gridt.atom_iw2_ylm_g, gridt.atom_iw2_l_g,
          gridt.atom_nw_g, gridt.nr_max, gridt.psi_u_g, psir_ylm_right_g,
          dpsir_ylm_left_x_g, dpsir_ylm_left_y_g, dpsir_ylm_left_z_g,
          ddpsir_ylm_left_xx_g, ddpsir_ylm_left_xy_g, ddpsir_ylm_left_xz_g,
          ddpsir_ylm_left_yy_g, ddpsir_ylm_left_yz_g, ddpsir_ylm_left_zz_g);
      checkCudaLastError();
      gridt.fastest_matrix_mul(
          max_m, max_n, atom_pair_A_m_g, atom_pair_B_n_g, atom_pair_k_g,
          atom_pair_mat_A_array_g, atom_pair_lda_g, atom_pair_mat_B_array_g,
          atom_pair_ldb_g, atom_pair_mat_C_array_g, atom_pair_ldc_g,
          atom_pair_num, gridt.streams[stream_num], nullptr);
      
      /* force compute in GPU */
      dim3 grid_dot_force(blocksPerGrid_force);
      dim3 block_dot_force(threadsPerBlock_force);
      dim3 grid_dot(blocksPerGrid);
      dim3 block_dot(threadsPerBlock);
      dot_product_force<<<grid_dot_force, block_dot_force, 0,
                          gridt.streams[stream_num]>>>(
          dpsir_ylm_left_x_g, dpsir_ylm_left_y_g, dpsir_ylm_left_z_g,
          psir_ylm_dm_g, force_dot_g, iat_g, nwmax, max_size,
          gridt.psir_size / nwmax);

      checkCuda(cudaMemcpy(force_h, force_dot_g,
                           3 * nbz * bxyz * max_size * sizeof(double),
                           cudaMemcpyDeviceToHost));
      for (int index1 = 0; index1 < nbz * bxyz * max_size; index1++) {
        int iat1 = iat[index1];
        if (iat1 >= 0) {
          for (int index2 = 0; index2 < 3; index2++) {
            force[iat1 * 3 + index2] += force_h[index1 * 3 + index2];
          }
        }
      }
      /*stress compute in GPU host*/
      dot_product_stress<<<grid_dot, block_dot, 0, gridt.streams[stream_num]>>>(
          ddpsir_ylm_left_xx_g, ddpsir_ylm_left_xy_g, ddpsir_ylm_left_xz_g,
          ddpsir_ylm_left_yy_g, ddpsir_ylm_left_yz_g, ddpsir_ylm_left_zz_g,
          psir_ylm_dm_g, stress_dot_g, gridt.psir_size);

      checkCuda(cudaMemcpy(stress_dot, stress_dot_g,
                           6 * blocksPerGrid * sizeof(double),
                           cudaMemcpyDeviceToHost));
      for (int i = 0; i < 6; i++) {
        for (int index = 0; index < blocksPerGrid; index++) {
          stress[i] += stress_dot[i * blocksPerGrid + index];
        }
      }
      /*free variables in GPU host*/
      delete[] iat;
      delete[] force_h;
      iter_num++;
    }
  }
  /*free variables in CPU host*/
  delete[] stress_dot;
  delete[] dm_matrix_h;
  checkCuda(cudaFree(dm_matrix_g));
  checkCuda(cudaFree(stress_dot_global_g));
  checkCuda(cudaFree(force_dot_global_g));
  checkCuda(cudaFree(iat_global_g));
  for (int i = 0; i < gridt.nstreams; i++) {
    checkCuda(cudaStreamSynchronize(gridt.streams[i]));
  }
}

} // namespace GintKernel
