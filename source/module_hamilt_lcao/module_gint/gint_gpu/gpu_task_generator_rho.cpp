#include "gint_rho.h"
#include "omp.h"

#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_lcao/module_gint/gint_tools.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

void gpu_task_generator_rho(const Grid_Technique &GridT, 
                            const int i, const int j,
                            const int psi_size_max, const int max_size,
                            const int nczp,
                            double *psi_input_double, int *psi_input_int,
                            int *num_psir,
                            const int lgd,
                            double *psir_ylm_g,
                            double *psir_zeros_g,
                            double *dm_matrix_g,
                            int *mat_m,
                            int *mat_n,
                            int *mat_k,
                            int *mat_lda,
                            int *mat_ldb,
                            int *mat_ldc,
                            double **mat_A,
                            double **mat_B,
                            double **mat_C,
                            int &max_m,
                            int &max_n,
                            int &atom_pair_num,
                            double *rho_g,
                            double **vec_l,
                            double **vec_r,
                            double **dot_product,
                            int *vec_len,
                            int &dot_count 
                            ) 
{ 
  const int grid_index_ij = i * GridT.nby * GridT.nbzp + j * GridT.nbzp;
  const int nwmax = GlobalC::ucell.nwmax;
  bool *gpu_mat_cal_flag = new bool[max_size * GridT.nbzp];
  int * start_idx_psi = new int[max_size * GridT.nbzp];
  for (int i = 0; i < max_size * GridT.nbzp; i++)
  {  
    gpu_mat_cal_flag[i] = false;
  }
  dot_count = 0;

  // psir generate
  for (int z_index = 0; z_index < GridT.nbzp; z_index++) {
    int num_get_psi = 0;
    int grid_index = grid_index_ij + z_index;
    int num_psi_pos = psi_size_max * z_index;
    int calc_flag_index = max_size * z_index;
    int bcell_start_index = GridT.bcell_start[grid_index];
    int atom_offset_psi = 0;

    for (int id = 0; id < GridT.how_many_atoms[grid_index]; id++) {
      int ib = 0;
      int mcell_index = bcell_start_index + id;
      int imcell = GridT.which_bigcell[mcell_index];
      int iat = GridT.which_atom[mcell_index];
      int it_temp = GlobalC::ucell.iat2it[iat];
      int start_ind_grid = GridT.start_ind[grid_index];
      int nw = GlobalC::ucell.atoms[it_temp].nw;
      start_idx_psi[z_index*max_size+id] = z_index * max_size * GridT.bxyz * nwmax+
                                                atom_offset_psi;

      for (int bx_index = 0; bx_index < GridT.bx; bx_index++) {
        for (int by_index = 0; by_index < GridT.by; by_index++) {
          for (int bz_index = 0; bz_index < GridT.bz; bz_index++) {
            double dr_temp[3];
            dr_temp[0] = GridT.meshcell_pos[ib][0] +
                         GridT.meshball_positions[imcell][0] -
                         GridT.tau_in_bigcell[iat][0];
            dr_temp[1] = GridT.meshcell_pos[ib][1] +
                         GridT.meshball_positions[imcell][1] -
                         GridT.tau_in_bigcell[iat][1];
            dr_temp[2] = GridT.meshcell_pos[ib][2] +
                         GridT.meshball_positions[imcell][2] -
                         GridT.tau_in_bigcell[iat][2];

            double distance =
                sqrt(dr_temp[0] * dr_temp[0] + dr_temp[1] * dr_temp[1] +
                     dr_temp[2] * dr_temp[2]);
            if (distance <= GlobalC::ORB.Phi[it_temp].getRcut()) {
              gpu_mat_cal_flag[calc_flag_index + id] = true;
              int pos_temp_double = num_psi_pos + num_get_psi;
              int pos_temp_int = pos_temp_double * 2;
              pos_temp_double *= 5;
              if (distance < 1.0E-9)
                distance += 1.0E-9;
              psi_input_double[pos_temp_double] = dr_temp[0] / distance;
              psi_input_double[pos_temp_double + 1] = dr_temp[1] / distance;
              psi_input_double[pos_temp_double + 2] = dr_temp[2] / distance;
              psi_input_double[pos_temp_double + 3] = distance;

              int vindex_global = bx_index * GridT.ncy * nczp + by_index * nczp +
                                  bz_index + start_ind_grid;
              psi_input_double[pos_temp_double + 4] = 0;  // to be removed

              psi_input_int[pos_temp_int] = it_temp;
              psi_input_int[pos_temp_int + 1] =
                  start_idx_psi[z_index * max_size + id] + ib;
              num_get_psi++;
            }
            ib++;
          }
        }
      }
      atom_offset_psi += nw * GridT.bxyz;
    }
    num_psir[z_index] = num_get_psi;
  }


  //TODO:Separate the following code into a single function
  int tid = 0;
  max_m = 0;
  max_n = 0;

  for (int z_index = 0; z_index < GridT.nbzp; z_index++) {
    int grid_index = grid_index_ij + z_index;
    int calc_flag_index = max_size * z_index;
    int bcell_start_index = GridT.bcell_start[grid_index];
    int nw_total = 0;

    for (int atom1 = 0; atom1 < GridT.how_many_atoms[grid_index]; atom1++) {
      int mcell_index1 = bcell_start_index + atom1;
      int iat1 = GridT.which_atom[mcell_index1];
      int it1 = GlobalC::ucell.iat2it[iat1];
      int lo1=GridT.trace_lo[GlobalC::ucell.itiaiw2iwt(
                it1, GlobalC::ucell.iat2ia[iat1],0)];
      int nw1 = GlobalC::ucell.atoms[it1].nw;
      nw_total += nw1;
      if(!gpu_mat_cal_flag[calc_flag_index + atom1]){
        continue;
      }
      
      for(int atom2 = 0; atom2 < GridT.how_many_atoms[grid_index]; atom2++) {
        if(!gpu_mat_cal_flag[calc_flag_index + atom2]){
        continue;
        }
        int mcell_index2 = bcell_start_index + atom2;
        int iat2 = GridT.which_atom[mcell_index2];
        int it2 = GlobalC::ucell.iat2it[iat2];
        int lo2=GridT.trace_lo[GlobalC::ucell.itiaiw2iwt(
                it2, GlobalC::ucell.iat2ia[iat2],0)];
        int nw2 = GlobalC::ucell.atoms[it2].nw;

        int mat_A_idx = lgd * lo1 + lo2;
        int mat_B_idx = start_idx_psi[z_index*max_size+atom2];
        int mat_C_idx = start_idx_psi[z_index*max_size+atom1];

        mat_m[tid] = nw1;
        mat_n[tid] = GridT.bxyz;
        mat_k[tid] = nw2;
        mat_lda[tid] = lgd;
        mat_ldb[tid] = GridT.bxyz;
        mat_ldc[tid] = GridT.bxyz;
        mat_A[tid] = dm_matrix_g + mat_A_idx;
        mat_B[tid] = psir_ylm_g + mat_B_idx;
        mat_C[tid] = psir_zeros_g + mat_C_idx;

        if(mat_m[tid] > max_m){
          max_m = mat_m[tid];
        }

        if(mat_n[tid] > max_n){
          max_n = mat_n[tid];
        }
        
        tid++;
      }

    }

    //generate dot tasks
    double *vec_l_start = psir_ylm_g + start_idx_psi[z_index*max_size];
    double *vec_r_start = psir_zeros_g + start_idx_psi[z_index*max_size];
    int* vindex = Gint_Tools::get_vindex(GridT.bxyz, GridT.bx, GridT.by, GridT.bz,
                        nczp, GridT.start_ind[grid_index], GridT.ncy*nczp);
    for(int i = 0; i < GridT.bxyz; i++){
      vec_l[dot_count] = vec_l_start + i;
      vec_r[dot_count] = vec_r_start + i;
      dot_product[dot_count] = rho_g + vindex[i];
      vec_len[dot_count] = nw_total;
      dot_count++;
    }
  }
  atom_pair_num = tid;

  delete[] gpu_mat_cal_flag;
  delete[] start_idx_psi;
} 

// rho calculation tasks generate
// void rho_cal_task( const Grid_Technique &GridT,
//                       const int i, const int j,
//                       const int max_size,
//                       const int lgd,
//                       const bool *gpu_mat_cal_flag,
//                       const int *start_idx_psi,
//                       const double *psir_ylm,
//                       const double *psir_zeros,
//                       const double *dm_matrix,
//                       int *mat_m,
//                       int *mat_n,
//                       int *mat_k,
//                       int *mat_lda,
//                       int *mat_ldb,
//                       int *mat_ldc,
//                       double **mat_A,
//                       double **mat_B,
//                       double **mat_C,
//                       int &max_m,
//                       int &max_n,
//                       int &atom_pair_num
//                       )
// {
//   int tid = 0;
//   max_m = 0;
//   max_n = 0;
//   const int grid_index_ij =
//             i * GridT.nby * GridT.nbzp + j * GridT.nbzp;
//   for (int z_index = 0; z_index < GridT.nbzp; z_index++) {
//     int grid_index = grid_index_ij + z_index;
//     int calc_flag_index = max_size * z_index;
//     int bcell_start_index = GridT.bcell_start[grid_index];

//     for (int atom1 = 0; atom1 < GridT.how_many_atoms[grid_index]; atom1++) {
//       if(!gpu_mat_cal_flag[calc_flag_index + atom1]){
//         continue;
//       }
//       int mcell_index1 = bcell_start_index + atom1;
//       int iat1 = GridT.which_atom[mcell_index1];
//       int it1 = GlobalC::ucell.iat2it[iat1];
//       int lo1=GridT.trace_lo[GlobalC::ucell.itiaiw2iwt(
//                 it1, GlobalC::ucell.iat2ia[iat1],0)];
//       int nw1 = GlobalC::ucell.atoms[it1].nw;
      
//       for(int atom2 = 0; atom2 < GridT.how_many_atoms[grid_index]; atom2++) {
//         if(!gpu_mat_cal_flag[calc_flag_index + atom2]){
//         continue;
//         }
//         int mcell_index2 = bcell_start_index + atom2;
//         int iat2 = GridT.which_atom[mcell_index2];
//         int it2 = GlobalC::ucell.iat2it[iat2];
//         int lo2=GridT.trace_lo[GlobalC::ucell.itiaiw2iwt(
//                 it2, GlobalC::ucell.iat2ia[iat2],0)];
//         int nw2 = GlobalC::ucell.atoms[it2].nw;

//         int mat_A_idx = lgd * lo1 + lo2;
//         int mat_B_idx = start_idx_psi[z_index*max_size+atom2];
//         int mat_C_idx = start_idx_psi[z_index*max_size+atom1];

//         mat_m[tid] = nw1;
//         mat_n[tid] = GridT.bxyz;
//         mat_k[tid] = nw2;
//         mat_lda[tid] = lgd;
//         mat_ldb[tid] = GridT.bxyz;
//         mat_ldc[tid] = GridT.bxyz;
//         mat_A[tid] = dm_matrix + mat_A_idx;
//         mat_B[tid] = psir_ylm + mat_B_idx;
//         mat_C[tid] = psir_zeros + mat_C_idx;

//         if(mat_m[tid] > max_m){
//           max_m = mat_m[tid];
//         }

//         if(mat_n[tid] > max_n){
//           max_n = mat_n[tid];
//         }
        
//         tid++;
//       }

//     }
//   }
//   atom_pair_num = tid;
// }


