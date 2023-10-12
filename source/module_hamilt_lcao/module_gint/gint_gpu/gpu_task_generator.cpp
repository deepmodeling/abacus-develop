#include "gint_vl.h"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_lcao/module_gint/gint_tools.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

void gpu_task_generate_vlocal(const Grid_Technique &GridT, int i, int j, int bx,
                              int by, int bz, int bxyz,
                              int atom_pair_size_of_meshcell, int psi_size_max,
                              int max_size, const int ncx, const int ncy,
                              const int nczp, double *dr, int *it,
                              int *psir_ylm_start, int *ib_index, int *num_psir,
                              int *start_ind, int *vindex,
                              int *atom_pair_input_info, int *num_atom_pair) {
  const int nbx = GridT.nbx;
  const int nby = GridT.nby;
  const int nbz = GridT.nbzp;
  const int nwmax = GlobalC::ucell.nwmax;
  const int namax = GlobalC::ucell.namax;
  const int ntype = GlobalC::ucell.ntype;
  for (int z_index = 0; z_index < nbz; z_index++) {
    int num_get_psi = 0;
    int grid_index = i * nby * nbz + j * nbz + z_index;
    int num_psi_pos = psi_size_max * z_index;
    for (int id = 0; id < GridT.how_many_atoms[grid_index]; id++) {
      for (int ib = 0; ib < bxyz; ib++) {
        int mcell_index = GridT.bcell_start[grid_index] + id;
        int imcell = GridT.which_bigcell[mcell_index];
        int iat = GridT.which_atom[mcell_index];
        int it_temp = GlobalC::ucell.iat2it[iat];
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
          int pos_temp = num_psi_pos + num_get_psi;
          if (distance < 1.0E-9)
            distance += 1.0E-9;
          dr[pos_temp * 4] = dr_temp[0] / distance;
          dr[pos_temp * 4 + 1] = dr_temp[1] / distance;
          dr[pos_temp * 4 + 2] = dr_temp[2] / distance;
          dr[pos_temp * 4 + 3] = distance;
          it[pos_temp] = it_temp;
          int dist_tmp = z_index * bxyz * max_size + ib * max_size + id;
          psir_ylm_start[pos_temp] = dist_tmp;
          ib_index[pos_temp] = ib;
          num_get_psi++;
        }
      }
    }
    num_psir[z_index] = num_get_psi;

    int vindex_temp = z_index * bxyz;
    for (int bx_index = 0; bx_index < bx; bx_index++) {
      for (int by_index = 0; by_index < by; by_index++) {
        for (int bz_index = 0; bz_index < bz; bz_index++) {
          int vindex_global = bx_index * ncy * nczp + by_index * nczp +
                              bz_index + start_ind[grid_index];
          vindex[vindex_temp] = vindex_global;
          vindex_temp++;
        }
      }
    }

    int atom_pair_index_in_nbz = atom_pair_size_of_meshcell * z_index;
    int atom_pair_index_in_meshcell = 0;
    for (int atom1 = 0; atom1 < GridT.how_many_atoms[grid_index]; atom1++) {
      for (int atom2 = 0; atom2 < GridT.how_many_atoms[grid_index]; atom2++) {
        int iat1 = GridT.which_atom[GridT.bcell_start[grid_index] + atom1];
        int iat2 = GridT.which_atom[GridT.bcell_start[grid_index] + atom2];
        int it1 = GlobalC::ucell.iat2it[iat1];
        int it2 = GlobalC::ucell.iat2it[iat2];
        int lo1 = GridT.trace_lo[GlobalC::ucell.itiaiw2iwt(
            it1, GlobalC::ucell.iat2ia[iat1], 0)];
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