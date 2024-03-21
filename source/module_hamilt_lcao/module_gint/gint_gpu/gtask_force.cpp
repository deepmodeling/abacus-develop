#include "gint_force.h"
#include "module_base/ylm.h"
#include "module_hamilt_lcao/module_gint/gint_tools.h"
#include "omp.h"
namespace lcaoCudaKernel{
/**
 * @brief Description of the function.
 *
 * Detailed description of the function.
 *
 * @param gridt The Grid_Technique object.
 * @param i The integer parameter i.
 * @param j The integer parameter j.
 * @param psi_size_max The maximum size of psi.
 * @param max_size The maximum size.
 * @param nczp The nczp parameter.
 * @param vfactor The vfactor parameter.
 * @param vlocal_global_value The array of vlocal_global_value.
 * @param iat_per_nbz The array of iat_per_nbz.
 * @param psi_input_double The double array of psi_input.
 * @param psi_input_int The integer array of psi_input.
 * @param num_psir The array of num_psir.
 * @param lgd The lgd parameter.
 * @param psir_ylm_g The double array of psir_ylm_g.
 * @param psir_zeros_g The double array of psir_zeros_g.
 * @param dm_matrix_g The double array of dm_matrix_g.
 * @param mat_m The array of mat_m.
 * @param mat_n The array of mat_n.
 * @param mat_k The array of mat_k.
 * @param mat_lda The array of mat_lda.
 * @param mat_ldb The array of mat_ldb.
 * @param mat_ldc The array of mat_ldc.
 * @param mat_A The pointer to mat_A.
 * @param mat_B The pointer to mat_B.
 * @param mat_C The pointer to mat_C.
 * @param max_m The reference to max_m.
 * @param max_n The reference to max_n.
 * @param atom_pair_num The reference to atom_pair_num.
 */
void gpu_task_generator_force(
    const Grid_Technique &gridt, const LCAO_Orbitals &ORB,
    const UnitCell &ucell, const int i, const int j,
    const int psi_size_max, const int max_size, const int nczp,
    const double vfactor, const double *vlocal_global_value, int *iat_per_nbz,
    double *psi_input_double, int *psi_input_int, int *num_psir, const int lgd,
    double *psir_ylm_g, double *psir_zeros_g, double *dm_matrix_g, int *mat_m,
    int *mat_n, int *mat_k, int *mat_lda, int *mat_ldb, int *mat_ldc,
    double **mat_A, double **mat_B, double **mat_C, int &max_m, int &max_n,
    int &atom_pair_num) {
  const int grid_index_ij = i * gridt.nby * gridt.nbzp + j * gridt.nbzp;
  const int nwmax = ucell.nwmax;
  bool *gpu_mat_cal_flag = new bool[max_size * gridt.nbzp];

  for (int i = 0; i < max_size * gridt.nbzp; i++) {
    gpu_mat_cal_flag[i] = false;
  }
  // psir generate
  for (int z_index = 0; z_index < gridt.nbzp; z_index++) {
    int num_get_psi = 0;
    int grid_index = grid_index_ij + z_index;
    int num_psi_pos = psi_size_max * z_index;
    int calc_flag_index = max_size * z_index;
    int bcell_start_index = gridt.bcell_start[grid_index];
    int na_grid = gridt.how_many_atoms[grid_index];

    for (int id = 0; id < na_grid; id++) {
      int ib = 0;
      int mcell_index = bcell_start_index + id;
      int imcell = gridt.which_bigcell[mcell_index];
      int iat = gridt.which_atom[mcell_index];
      int it_temp = ucell.iat2it[iat];
      int start_ind_grid = gridt.start_ind[grid_index];

      for (int bx_index = 0; bx_index < gridt.bx; bx_index++) {
        for (int by_index = 0; by_index < gridt.by; by_index++) {
          for (int bz_index = 0; bz_index < gridt.bz; bz_index++) {
            double dr_temp[3];
            dr_temp[0] = gridt.meshcell_pos[ib][0] +
                         gridt.meshball_positions[imcell][0] -
                         gridt.tau_in_bigcell[iat][0];
            dr_temp[1] = gridt.meshcell_pos[ib][1] +
                         gridt.meshball_positions[imcell][1] -
                         gridt.tau_in_bigcell[iat][1];
            dr_temp[2] = gridt.meshcell_pos[ib][2] +
                         gridt.meshball_positions[imcell][2] -
                         gridt.tau_in_bigcell[iat][2];

            double distance =
                sqrt(dr_temp[0] * dr_temp[0] + dr_temp[1] * dr_temp[1] +
                     dr_temp[2] * dr_temp[2]);
            if (distance <= ORB.Phi[it_temp].getRcut()) {
              gpu_mat_cal_flag[calc_flag_index + id] = true;
              int pos_temp_double = num_psi_pos + num_get_psi;
              int pos_temp_int = pos_temp_double * 2;
              pos_temp_double *= 5;
              if (distance < 1.0E-9)
                distance += 1.0E-9;
              psi_input_double[pos_temp_double] = dr_temp[0];
              psi_input_double[pos_temp_double + 1] = dr_temp[1];
              psi_input_double[pos_temp_double + 2] = dr_temp[2];
              psi_input_double[pos_temp_double + 3] = distance;
              int vindex_global = bx_index * gridt.ncy * nczp +
                                  by_index * nczp + bz_index + start_ind_grid;
              psi_input_double[pos_temp_double + 4] =
                  vlocal_global_value[vindex_global] * vfactor;

              psi_input_int[pos_temp_int] = it_temp;
              psi_input_int[pos_temp_int + 1] =
                  (z_index * gridt.bxyz + ib) * max_size * nwmax + id * nwmax;
              iat_per_nbz[z_index * gridt.bxyz * max_size + ib * max_size +
                          id] = iat;
              // printf("the na_grid is %d,the id is %d the ib is %d the iat is
              // %d the nwmax is %d\n",na_grid,id,ib,iat,nwmax);
              num_get_psi++;
            }
            ib++;
          }
        }
      }
    }
    num_psir[z_index] = num_get_psi;
  }

  // TODO:Separate the following code into a single function
  int tid = 0;
  max_m = 0;
  max_n = 0;

  for (int z_index = 0; z_index < gridt.nbzp; z_index++) {
    int grid_index = grid_index_ij + z_index;
    int calc_flag_index = max_size * z_index;
    int bcell_start_index = gridt.bcell_start[grid_index];
    int bcell_start_psir = z_index * gridt.bxyz * max_size * nwmax;

    for (int atom1 = 0; atom1 < gridt.how_many_atoms[grid_index]; atom1++) {
      if (!gpu_mat_cal_flag[calc_flag_index + atom1]) {
        continue;
      }
      int mcell_index1 = bcell_start_index + atom1;
      int iat1 = gridt.which_atom[mcell_index1];
      int it1 = ucell.iat2it[iat1];
      int lo1 = gridt.trace_lo[ucell.itiaiw2iwt(
          it1, ucell.iat2ia[iat1], 0)];
      int nw1 = ucell.atoms[it1].nw;

      for (int atom2 = 0; atom2 < gridt.how_many_atoms[grid_index]; atom2++) {
        if (!gpu_mat_cal_flag[calc_flag_index + atom2]) {
          continue;
        }
        int mcell_index2 = bcell_start_index + atom2;
        int iat2 = gridt.which_atom[mcell_index2];
        int it2 = ucell.iat2it[iat2];
        int lo2 = gridt.trace_lo[ucell.itiaiw2iwt(
            it2, ucell.iat2ia[iat2], 0)];
        int nw2 = ucell.atoms[it2].nw;

        int mat_A_idx = bcell_start_psir + atom2 * nwmax;
        int mat_B_idx = lgd * lo1 + lo2;
        int mat_C_idx = bcell_start_psir + atom1 * nwmax;
        mat_m[tid] = gridt.bxyz;
        mat_n[tid] = nw1;
        mat_k[tid] = nw2;
        mat_lda[tid] = nwmax * max_size;
        mat_ldb[tid] = lgd;
        mat_ldc[tid] = nwmax * max_size;
        mat_A[tid] = psir_ylm_g + mat_A_idx;
        mat_B[tid] = dm_matrix_g + mat_B_idx;
        mat_C[tid] = psir_zeros_g + mat_C_idx;

        if (mat_m[tid] > max_m) {
          max_m = mat_m[tid];
        }

        if (mat_n[tid] > max_n) {
          max_n = mat_n[tid];
        }

        tid++;
      }
    }
  }
  atom_pair_num = tid;

  delete[] gpu_mat_cal_flag;
}
}