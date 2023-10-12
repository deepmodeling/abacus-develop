#include "gint_vl.h"
#include "omp.h"

#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_lcao/module_gint/gint_tools.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

void gpu_task_generate_vlocal(const Grid_Technique &GridT, const int i,
                              const int j, const int bx, const int by,
                              const int bz, const int bxyz,
                              const int atom_pair_size_of_meshcell,
                              const int psi_size_max, const int max_size,
                              const int ncx, const int ncy, const int nczp,
                              const double vfactor, const int *start_ind,
                              const double *vlocal_global_value,
                              double *psi_input_double, int *psi_input_int,
                              int *num_psir, int *atom_pair_input_info,
                              int *num_atom_pair) {

  const int grid_index_ij = i * GridT.nby * GridT.nbzp + j * GridT.nbzp;
  const int nwmax = GlobalC::ucell.nwmax;
#pragma omp parallel for
  for (int z_index = 0; z_index < GridT.nbzp; z_index++) {
    int num_get_psi = 0;
    int grid_index = grid_index_ij + z_index;
    int num_psi_pos = psi_size_max * z_index;
    int bcell_start_index = GridT.bcell_start[grid_index];

    for (int id = 0; id < GridT.how_many_atoms[grid_index]; id++) {
      int ib = 0;
      int mcell_index = bcell_start_index + id;
      int imcell = GridT.which_bigcell[mcell_index];
      int iat = GridT.which_atom[mcell_index];
      int it_temp = GlobalC::ucell.iat2it[iat];
      int start_ind_grid = start_ind[grid_index];
      for (int bx_index = 0; bx_index < bx; bx_index++) {
        for (int by_index = 0; by_index < by; by_index++) {
          for (int bz_index = 0; bz_index < bz; bz_index++) {
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
              int pos_temp_double = num_psi_pos + num_get_psi;
              int pos_temp_int = pos_temp_double * 2;
              pos_temp_double *= 5;
              if (distance < 1.0E-9)
                distance += 1.0E-9;
              psi_input_double[pos_temp_double] = dr_temp[0] / distance;
              psi_input_double[pos_temp_double + 1] = dr_temp[1] / distance;
              psi_input_double[pos_temp_double + 2] = dr_temp[2] / distance;
              psi_input_double[pos_temp_double + 3] = distance;

              int vindex_global = bx_index * ncy * nczp + by_index * nczp +
                                  bz_index + start_ind_grid;
              psi_input_double[pos_temp_double + 4] =
                  vlocal_global_value[vindex_global] * vfactor;

              psi_input_int[pos_temp_int] = it_temp;
              psi_input_int[pos_temp_int + 1] =
                  ((z_index * max_size + id) * bxyz) * nwmax + ib;
              num_get_psi++;
            }
            ib++;
          }
        }
      }
    }
    num_psir[z_index] = num_get_psi;

    int atom_pair_index_in_nbz = atom_pair_size_of_meshcell * z_index;
    int atom_pair_index_in_meshcell = 0;
    int atom_num = GridT.how_many_atoms[grid_index];
    for (int atom1 = 0; atom1 < atom_num; atom1++) {

      int iat1 = GridT.which_atom[bcell_start_index + atom1];
      int it1 = GlobalC::ucell.iat2it[iat1];
      int lo1 = GridT.trace_lo[GlobalC::ucell.itiaiw2iwt(
          it1, GlobalC::ucell.iat2ia[iat1], 0)];

      for (int atom2 = 0; atom2 < atom_num; atom2++) {
        int iat2 = GridT.which_atom[bcell_start_index + atom2];
        int it2 = GlobalC::ucell.iat2it[iat2];
        int lo2 = GridT.trace_lo[GlobalC::ucell.itiaiw2iwt(
            it2, GlobalC::ucell.iat2ia[iat2], 0)];
        if (lo1 <= lo2) {
          int atom_pair_index =
              atom_pair_index_in_nbz + atom_pair_index_in_meshcell;
          atom_pair_input_info[atom_pair_index] = atom1;
          atom_pair_input_info[atom_pair_index + 1] = atom2;
          atom_pair_input_info[atom_pair_index + 2] =
              GlobalC::ucell.atoms[it1].nw;
          atom_pair_input_info[atom_pair_index + 3] =
              GlobalC::ucell.atoms[it2].nw;
          atom_pair_input_info[atom_pair_index + 4] = lo1;
          atom_pair_input_info[atom_pair_index + 5] = lo2;
          atom_pair_index_in_meshcell += 6;
        }
      }
    }
    num_atom_pair[z_index] = atom_pair_index_in_meshcell;
  }
}