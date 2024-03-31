#include <omp.h>

#include "gint_vl.h"
#include "module_base/ylm.h"
#include "module_hamilt_lcao/module_gint/gint_tools.h"
namespace GintKernel
{

void gtask_vlocal(const Grid_Technique& gridt,
                  const double* rcut,
                  const UnitCell& ucell,
                  const int i,
                  const int j,
                  const int max_size,
                  const int nczp,
                  const double vfactor,
                  const double* vlocal_global_value,
                  double* psir_ylm_left,
                  double* psir_ylm_right,
                  double* psi_input_double,
                  int* input_int,
                  int* num_psir,
                  int* atom_pair_A_m,
                  int* atom_pair_B_n,
                  int* atom_pair_lda,
                  int* atom_pair_ldb,
                  int* atom_pair_ldc,
                  double** atom_pair_mat_A,
                  double** atom_pair_mat_B,
                  double** atom_pair_mat_C,
                  int& atom_pair_num,
                  int& max_m,
                  int& max_n)
{

    const int grid_index_ij = i * gridt.nby * gridt.nbzp + j * gridt.nbzp;
    const int nwmax = ucell.nwmax;
    bool* gpu_matrix_calc_flag = new bool[max_size * gridt.nbzp];
    for (int i = 0; i < max_size * gridt.nbzp; i++)
        gpu_matrix_calc_flag[i] = false;

    for (int z_index = 0; z_index < gridt.nbzp; z_index++)
    {
        int num_get_psi = 0;
        int grid_index = grid_index_ij + z_index;
        int num_psi_pos = gridt.psi_size_max_z * z_index;
        int calc_flag_index = max_size * z_index;
        int bcell_start_index = gridt.bcell_start[grid_index];

        for (int id = 0; id < gridt.how_many_atoms[grid_index]; id++)
        {
            int ib = 0;
            int mcell_index = bcell_start_index + id;
            int imcell = gridt.which_bigcell[mcell_index];
            int iat = gridt.which_atom[mcell_index];
            int it_temp = ucell.iat2it[iat];
            int start_ind_grid = gridt.start_ind[grid_index];
            for (int bx_index = 0; bx_index < gridt.bx; bx_index++)
            {
                for (int by_index = 0; by_index < gridt.by; by_index++)
                {
                    for (int bz_index = 0; bz_index < gridt.bz; bz_index++)
                    {
                        double dr_temp[3];
                        dr_temp[0] = gridt.meshcell_pos[ib][0]
                                     + gridt.meshball_positions[imcell][0]
                                     - gridt.tau_in_bigcell[iat][0];
                        dr_temp[1] = gridt.meshcell_pos[ib][1]
                                     + gridt.meshball_positions[imcell][1]
                                     - gridt.tau_in_bigcell[iat][1];
                        dr_temp[2] = gridt.meshcell_pos[ib][2]
                                     + gridt.meshball_positions[imcell][2]
                                     - gridt.tau_in_bigcell[iat][2];

                        double distance = sqrt(dr_temp[0] * dr_temp[0]
                                               + dr_temp[1] * dr_temp[1]
                                               + dr_temp[2] * dr_temp[2]);
                        if (distance <= rcut[it_temp])
                        {
                            gpu_matrix_calc_flag[calc_flag_index + id] = true;
                            int pos_temp_double = num_psi_pos + num_get_psi;
                            int pos_temp_int = pos_temp_double * 2;
                            pos_temp_double *= 5;
                            if (distance < 1.0E-9)
                                distance += 1.0E-9;
                            psi_input_double[pos_temp_double]
                                = dr_temp[0] / distance;
                            psi_input_double[pos_temp_double + 1]
                                = dr_temp[1] / distance;
                            psi_input_double[pos_temp_double + 2]
                                = dr_temp[2] / distance;
                            psi_input_double[pos_temp_double + 3] = distance;

                            int vindex_global = bx_index * gridt.ncy * nczp
                                                + by_index * nczp + bz_index
                                                + start_ind_grid;
                            psi_input_double[pos_temp_double + 4]
                                = vlocal_global_value[vindex_global] * vfactor;

                            input_int[pos_temp_int] = it_temp;
                            input_int[pos_temp_int + 1]
                                = ((z_index * max_size + id) * gridt.bxyz)
                                      * nwmax
                                  + ib;
                            num_get_psi++;
                        }
                        ib++;
                    }
                }
            }
        }
        num_psir[z_index] = num_get_psi;
    }

    atom_pair_num = 0;
    max_m = 0;
    max_n = 0;
    for (int z_index = 0; z_index < gridt.nbzp; z_index++)
    {
        int grid_index = grid_index_ij + z_index;
        int atom_num = gridt.how_many_atoms[grid_index];
        int vldr3_index = z_index * max_size * nwmax * gridt.bxyz;
        int bcell_start_index = gridt.bcell_start[grid_index];
        int calc_flag_index = max_size * z_index;
        for (int atom1 = 0; atom1 < atom_num; atom1++)
        {

            int iat1 = gridt.which_atom[bcell_start_index + atom1];
            int it1 = ucell.iat2it[iat1];
            int lo1
                = gridt.trace_lo[ucell.itiaiw2iwt(it1, ucell.iat2ia[iat1], 0)];
            if (gpu_matrix_calc_flag[calc_flag_index + atom1] == false)
                continue;

            for (int atom2 = 0; atom2 < atom_num; atom2++)
            {
                if (gpu_matrix_calc_flag[calc_flag_index + atom2] == false)
                    continue;

                int iat2 = gridt.which_atom[bcell_start_index + atom2];
                int it2 = ucell.iat2it[iat2];
                int lo2 = gridt.trace_lo[ucell.itiaiw2iwt(it2,
                                                          ucell.iat2ia[iat2],
                                                          0)];
                if (lo1 <= lo2)
                {
                    int atom_pair_nw
                        = ucell.atoms[it1].nw * ucell.atoms[it2].nw;
                    if (gridt.grid_vlocal_g[iat1 * ucell.nat + iat2] == nullptr)
                    {
                        // Note that this situation occurs here because the
                        // logic in hcontainer and
                        //  grid integration is different.
                        //  In hcontainer, it is iat1<=iat2, and in grid
                        //  integral, it is lo1<=lo2. This is not entirely
                        //  equivalent in practice. We need to investigate
                        //  what's going on later.
                        //  TODO
                        continue;
                        // std::cout << "Error: GridVlocal did not malloc" <<
                        // std::endl;
                    }

                    int calc_index1 = vldr3_index + atom1 * nwmax * gridt.bxyz;
                    int calc_index2 = vldr3_index + atom2 * nwmax * gridt.bxyz;

                    atom_pair_mat_A[atom_pair_num]
                        = psir_ylm_left + calc_index1;
                    atom_pair_mat_B[atom_pair_num]
                        = psir_ylm_right + calc_index2;
                    atom_pair_mat_C[atom_pair_num]
                        = gridt.grid_vlocal_g[iat1 * ucell.nat + iat2];

                    atom_pair_lda[atom_pair_num] = gridt.bxyz;
                    atom_pair_ldb[atom_pair_num] = gridt.bxyz;
                    atom_pair_ldc[atom_pair_num] = ucell.atoms[it2].nw;

                    atom_pair_A_m[atom_pair_num] = ucell.atoms[it1].nw;
                    atom_pair_B_n[atom_pair_num] = ucell.atoms[it2].nw;
                    if (atom_pair_A_m[atom_pair_num] > max_m)
                    {
                        max_m = atom_pair_A_m[atom_pair_num];
                    }
                    if (atom_pair_B_n[atom_pair_num] > max_n)
                    {
                        max_n = atom_pair_B_n[atom_pair_num];
                    }
                    atom_pair_num++;
                }
            }
        }
    }

    delete[] gpu_matrix_calc_flag;
}

} // namespace GintKernel